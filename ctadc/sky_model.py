# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Make sky model.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
from collections import OrderedDict
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.image import SkyImage, SkyImageList
from gammapy.catalog import select_sky_box, SourceCatalogGammaCat
from .extern.xml import gammacat_to_xml

__all__ = [
    'GPSSkyModelSourcesBright',
    'GPSSkyModelSourcesFaint',
    'make_sky_models_xml',
    'make_sky_models_images',
]

log = logging.getLogger(__name__)


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
        self.gammacat = SourceCatalogGammaCat()

    def make_xml(self):
        """Make XML version of sky model."""
        self.xml = gammacat_to_xml(self.gammacat)


class SkyModelImageMaker:
    """
    Make sky model images.
    """

    IMAGE_SPECS_DEFAULT = [
        dict(name='red', energy_band=[0.1, 1] * u.TeV),
        dict(name='green', energy_band=[1, 10] * u.TeV),
        dict(name='blue', energy_band=[10, 100] * u.TeV),
    ]

    WCS_SPEC_DEFAULT = OrderedDict(
        nxpix=9400, nypix=500, binsz=0.02,
        xref=341, yref=0,
        proj='CAR', coordsys='GAL',
    )

    def __init__(self, image_specs=None, wcs_spec=None):
        if image_specs is None:
            image_specs = self.IMAGE_SPECS_DEFAULT.copy()

        if wcs_spec is None:
            wcs_spec = self.WCS_SPEC_DEFAULT.copy()

        self.image_specs = image_specs
        self.wcs_spec = wcs_spec

    def make_all_images(self):
        for image_spec in self.image_specs:
            self.make_flux_image(image_spec)

    def make_flux_image(self, image_spec):
        name = image_spec['name']
        energy_band = image_spec['energy_band']
        log.info('Making image: name={}, energy_band={}'.format(name, energy_band))

        gammacat = SourceCatalogGammaCat()
        image = SkyImage.empty(name=name, dtype='float32', **self.wcs_spec)

        # filter catalog so that center of the source is within image boundaries
        footprint = image.footprint('corner')
        lon_lim = footprint['lower right'].l.wrap_at('180d'), footprint['lower left'].l
        lat_lim = footprint['lower left'].b, footprint['upper left'].b
        gammacat.table['GLON'] = gammacat.table['glon']
        gammacat.table['GLAT'] = gammacat.table['glat']
        gammacat.table = select_sky_box(gammacat.table, lon_lim, lat_lim, frame='galactic')

        for source in gammacat:
            self.add_source_to_image(image, source, energy_band)

        filename = 'sky_model/images/ctadc_skymodel_gps_sources_bright_{}.fits.gz'.format(name)
        log.info('Writing {}'.format(filename))
        image.write(filename, clobber=True)

    @staticmethod
    def add_source_to_image(image, source, energy_band):
        try:
            morph = source.spatial_model(*energy_band)
        except ValueError:
            log.warning('Missing spatial model: source_id={}, common_name={}'
                        ''.format(source.data['source_id'], source.data['common_name']))
            return

        # Evaluate morphology model on cutout
        try:
            x, y = morph.x_mean, morph.y_mean
        except AttributeError:
            # TODO: these different cases should be handles in Gammapy spatial model classes!
            x, y = morph.x_0, morph.y_0

        pos = SkyCoord(x, y, unit='deg', frame='galactic')

        # TODO: handle bounding box according to source size, now just
        # use default of 4 x 4 deg
        size = (4 * u.deg, 4 * u.deg)
        cutout = image.cutout(pos, size=size)

        c = cutout.coordinates()
        l, b = c.galactic.l.wrap_at('180d'), c.galactic.b
        cutout.data = morph(l.deg, b.deg)
        image.paste(cutout)


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


def make_sky_models_xml():
    """Make all sky model XML files.
    """
    gps_sources_bright = GPSSkyModelSourcesBright()
    gps_sources_bright.make()
    gps_sources_bright.write_xml()

    gps_sources_faint = GPSSkyModelSourcesFaint()
    gps_sources_faint.make()
    gps_sources_faint.write_xml()


def make_sky_models_images():
    """Make all sky model image files.
    """
    maker = SkyModelImageMaker()
    maker.make_all_images()

    # TODO: combined RGB image.
