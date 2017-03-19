# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Make sky model.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
import textwrap
from collections import OrderedDict
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from astropy.table import vstack as table_vstack
from gammapy.image import SkyImage
from gammapy.catalog import select_sky_box, SourceCatalogGammaCat
from gammapy.astro.population import make_base_catalog_galactic, add_observed_parameters

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
    filename = 'sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat.xml'

    def __init__(self):
        pass

    def make(self):
        self.read_gammacat()
        self.apply_selections()
        self.make_xml()

    def read_gammacat(self):
        self.gammacat = SourceCatalogGammaCat()

    def apply_selections(self):
        return
        # TODO: Select on source_class = Galactic once available in gamma-cat
        # For now do a quick selection in GLAT
        # import IPython; IPython.embed(); 1/0
        # if np.abs(source.data['glat']) > 5 * u.deg:
        #     continue

        # For debugging, just do a few sources:
        # if source_idx == 3: break

    def make_xml(self):
        """Make XML version of sky model."""
        cat = self.gammacat

        # TODO: temp hack. Remove or move to gamma-cat creation.
        # import IPython; IPython.embed()
        idx = np.where(cat.table['morph_type'] == '')[0]
        cat.table.remove_rows(idx)

        source_library = cat.to_source_library()
        header = textwrap.dedent("""
        <!-- Bright sources for CTA-1DC from gamma-cat -->
        """)
        self.xml = source_library.to_xml(
            title='gammacat',
            header=header,
        )


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

    def make_rgb_image(self, smooth_radius='0.1 deg'):
        from astropy.visualization import make_lupton_rgb

        image_r = SkyImage.read(self._filename('red'))
        image_g = SkyImage.read(self._filename('green'))
        image_b = SkyImage.read(self._filename('blue'))
        images = [image_r, image_g, image_b]

        for image in images:
            # TODO: for now, we need to flip in y direction
            # https://github.com/matplotlib/matplotlib/issues/7656
            # https://github.com/astropy/astropy/pull/5535#issuecomment-267127610
            image.data = np.flipud(image.data)

            image.data = image.smooth(radius=Angle(smooth_radius)).data
            image.data /= image.data.max()
            image.data = np.power(image.data, 0.15)
            log.info('RGB max values: {}'.format(image.data.max()))

        minimum = min([_.data.min() for _ in images])
        log.info('RGB minimum: {}'.format(minimum))

        filename = 'sky_model/images/ctadc_skymodel_gps_sources_bright_rgb.png'
        log.info('Writing {}'.format(filename))
        rgb = make_lupton_rgb(
            image_r=image_r.data,
            image_g=image_g.data,
            image_b=image_b.data,
            minimum=-0,
            stretch=1,
            Q=0,
            filename=filename,
        )
        means = rgb.mean(axis=0).mean(axis=0)
        log.info('RBG image means: {}'.format(means))
        # import IPython; IPython.embed()

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

        filename = self._filename(name)
        log.info('Writing {}'.format(filename))
        image.write(filename, clobber=True)

    @staticmethod
    def _filename(name):
        return 'sky_model/images/ctadc_skymodel_gps_sources_bright_{}.fits.gz'.format(name)

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

        # TODO: handle bounding box according to source size, now just use a large constant bbox
        size = (5 * u.deg, 5 * u.deg)
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
        self.make_pwn()
        self.make_snr()
        self.make_total()
        self.make_xml()

    def make_pwn(self):
        # TODO: the PWN code is at the moment elsewhere. Merge!
        max_age = 1E6 * u.yr
        SN_rate = 3. / (100. * u.yr)
        n_sources = max_age * SN_rate
        table = make_base_catalog_galactic(
            n_sources=n_sources,
            rad_dis='L06',
            vel_dis='F06B',
            max_age=max_age,
            spiralarms=True,
        )
        add_observed_parameters(table)

        # Example how to make a selection using a boolean mask
        # TODO: put something useful!
        mask = table['distance'] < 3 * u.kpc
        table = table[mask]

        # Example how to make a selection using row indices
        idx = [0, 42, 99]
        table = table[idx]

        self.table_pwn = table

    def make_snr(self):
        self.table_snr = Table()

    def make_total(self):
        self.table_total = table_vstack([
            self.table_pwn,
            self.table_snr,
        ])

    def make_xml(self):
        self.xml = '<sources></sources>'  # a dummy

    def write_table(self):
        filename = 'sky_model/ctadc_skymodel_gps_sources_faint.ecsv'
        log.info('Writing {}'.format(filename))
        self.table_total.write(filename, format='ascii.ecsv')


def make_sky_models_xml():
    """Make all sky model XML files.
    """
    gps_sources_bright = GPSSkyModelSourcesBright()
    gps_sources_bright.make()
    gps_sources_bright.write_xml()

    # gps_sources_faint = GPSSkyModelSourcesFaint()
    # gps_sources_faint.make()
    # gps_sources_faint.write_xml()
    # gps_sources_faint.write_table()


def make_sky_models_images():
    """Make all sky model image files.
    """
    maker = SkyModelImageMaker()
    # maker.make_all_images()
    maker.make_rgb_image()
