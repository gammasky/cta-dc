# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Make sky model.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
from collections import OrderedDict
from astropy import units as u
from astropy.coordinates import SkyCoord

__all__ = [
    'GPSSkyModelSourcesBright',
    'GPSSkyModelSourcesFaint',
    'make_all_sky_models',
]

log = logging.getLogger(__name__)


ENERGY_BANDS = dict(red=[0.1, 1] * u.TeV,
                    green=[1, 10] * u.TeV,
                    blue=[10, 100] * u.TeV)

WCS_SPEC = dict(nxpix=9400,
                nypix=500,
                binsz=0.02,
                xref=341,
                yref=0,
                proj='CAR',
                coordsys='GAL')

class SkyModelMixin:
    """
    Helper class for sky models.

    This way of structuring code is called a "mix-in" class in Python.
    """

    def write_xml(self):
        xml = self.xml
        filename = self.filename

        log.info('Writing {}'.format(filename))
        with open(filename, 'w') as fh:
            fh.write(xml)


class GPSSkyModelSourcesBright(SkyModelMixin):
    """
    GPS sky model - bright sources component.

    https://github.com/gammapy/gamma-cat/
    """
    filename = 'sky_model/ctadc_skymodel_gps_sources_bright.xml'

    def __init__(self):
        pass

    def make(self):
        self.read_gammacat()
        self.make_xml()

    def read_gammacat(self):
        from gammapy.catalog import SourceCatalogGammaCat
        self.gammacat = SourceCatalogGammaCat()

    def make_xml(self):
        """Make XML version of sky model."""
        from .extern.xml import gammacat_to_xml
        self.xml = gammacat_to_xml(self.gammacat)


    def make_rgb_gammacat(self, energy_bands=ENERGY_BANDS, wcs_spec=WCS_SPEC):
        """
        Make RGB sky model image based on on the gamma-cat catalog.

        Parameters
        ----------
        energy_bands : dict
            Specification of the energy bands for the R, G and B channel
            of the output image.
        wcs_spec : dict

        """
        from gammapy.image import SkyImage, SkyImageList
        from gammapy.catalog import select_sky_box, SourceCatalogGammaCat

        # reload gammacat because it's modified later
        gammacat = SourceCatalogGammaCat()

        # set up empty images
        red = SkyImage.empty(name='red', **wcs_spec)
        green = SkyImage.empty(name='green', **wcs_spec)
        blue = SkyImage.empty(name='blue', **wcs_spec)
        flux_rgb = SkyImageList([red, green, blue])

        # filter catalog so that center of the source is within image boundaries
        footprint = red.footprint('corner')

        lon_lim = footprint['lower right'].l.wrap_at('180d'), footprint['lower left'].l
        lat_lim = footprint['lower left'].b, footprint['upper left'].b

        gammacat.table = select_sky_box(gammacat.table, lon_lim, lat_lim, frame='galactic')

        for source in gammacat:
            for channel in ['red', 'green', 'blue']:
                emin, emax = energy_bands[channel]
                try:
                    morph = source.morphology_model(emin, emax)
                except ValueError:
                    continue

                # Evaluate morphology model on cutout
                try:
                    x, y = morph.x_mean, morph.y_mean
                except AttributeError:
                    x, y = morph.x_0, morph.y_0
                pos = SkyCoord(x, y, unit='deg', frame='galactic')

                # TODO: handle bounding box according to source size, now just
                # use default of 4 x 4 deg
                size = (4 * u.deg, 4 * u.deg)
                cutout = flux_rgb[channel].cutout(pos, size=size)

                c = cutout.coordinates()
                l, b = c.galactic.l.wrap_at('180d'), c.galactic.b
                cutout.data = morph(l.deg, b.deg)
                flux_rgb[channel].paste(cutout)

        for channel in ['red', 'green', 'blue']:
            filename = 'sky_model/images/{}.fits'.format(channel)
            log.info('Writing {}'.format(filename))
            flux_rgb[channel].write(filename, clobber=True)


class GPSSkyModelSourcesFaint(SkyModelMixin):
    """
    GPS sky model - faint sources component.

    """
    filename = 'sky_model/ctadc_skymodel_gps_sources_faint.xml'

    def __init__(self):
        pass

    def make(self):
        self.make_xml()

    def make_xml(self):
        self.xml = '<sources></sources>'  # a dummy


def make_all_sky_models():
    """Make all sky models.
    """
    gps_sources_bright = GPSSkyModelSourcesBright()
    gps_sources_bright.make()
    gps_sources_bright.write_xml()
    gps_sources_bright.make_rgb_gammacat()

    gps_sources_faint = GPSSkyModelSourcesFaint()
    gps_sources_faint.make()
    gps_sources_faint.write_xml()
