"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
import logging
from collections import OrderedDict
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
import gammalib
from gammapy.utils.energy import Energy
from gammapy.spectrum import CrabSpectrum

log = logging.getLogger(__name__)


class Histogram:
    """
    Helper class to work with 1-dim histograms.

    Splits histogram computation from plotting,
    and introduces some conveniences.
    """

    def __init__(self, bins, vals):
        self.bins = bins
        self.vals = vals

    @classmethod
    def from_points(cls, bins, points):
        vals = np.histogram(points, bins=bins)[0]
        return cls(bins=bins, vals=vals)

    def copy(self):
        return self.__class__(
            bins=self.bins.copy(),
            vals=self.vals.copy(),
        )

    @property
    def bin_centers(self):
        return 0.5 * (self.bins[:-1] + self.bins[1:])

    def smooth(self, sigma=None):
        if sigma is None:
            vals = self.vals
        else:
            vals = gaussian_filter1d(self.vals.astype(float), sigma=sigma)
        return self.__class__(bins=self.bins, vals=vals)

    def normalise(self, method='max', value=1):
        if method == 'max':
            scale = float(value) / np.max(self.vals)
        elif method == 'int':
            integral = np.sum(self.vals * np.diff(self.bins))
            scale = float(value) / integral
        elif method == 'int log':
            integral = np.sum(self.vals * np.diff(np.log10(self.bins)))
            scale = float(value) / integral
        else:
            raise ValueError('Invalid method: {}'.format(method))

        return self.__class__(bins=self.bins, vals=scale * self.vals)

    def plot_hist(self, ax, **kwargs):
        """Plot a pre-binned histogram"""
        bins = self.bins
        vals = self.vals
        # Plotting a histogram nicely ourselves is nontrivial
        # Example from http://stackoverflow.com/a/18611135/498873
        left, right, width = bins[:-1], bins[1:], np.diff(bins)
        x = np.array([left, right]).T.flatten()
        y = np.array([vals, vals]).T.flatten()
        ax.plot(x, y, **kwargs)
        # ax.bar(
        #     left=left, height=hist, width=width,
        #     label=component['tag'], alpha=0.7,
        #     align='edge', color='none',
        # )

    def plot_smooth(self, ax, **kwargs):
        x = self.bin_centers
        y = self.vals
        ax.plot(x, y, **kwargs)


class HistogramStack:
    def __init__(self, hists=None):
        self.hists = hists if hists else []

    def copy(self):
        return [h.copy() for h in self.hists]

    def add(self, hist):
        self.hists.append(hist)

    @property
    def total(self):
        bins = self.hists[0].bins
        vals = np.vstack([h.vals for h in self.hists]).sum(axis=0)
        return Histogram(bins=bins, vals=vals)


def make_source_tables():
    log.info('Starting make_source_tables')
    parts = [
        dict(tag='gamma-cat', filename='gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'),
        dict(tag='image_sources', filename='image_sources/ctadc_skymodel_gps_sources_images.xml'),
        dict(tag='pwn', filename='pwn/ctadc_skymodel_gps_sources_pwn.xml'),
        dict(tag='snr', filename='snrs/ctadc_skymodel_gps_sources_snr_2.xml'),
        dict(tag='binaries', filename='binaries/ctadc_skymodel_gps_sources_binaries.xml'),
        dict(tag='pulsars', filename='pulsars/ctadc_skymodel_gps_sources_pulsars.xml'),
    ]
    for data in parts:
        filename = '../sky_model/' + data['filename']
        log.info('Reading {}'.format(filename))
        data['models'] = gammalib.GModels(filename)
        table = make_source_table(data)
        filename = filename.replace('.xml', '_summary.ecsv')
        log.info('Writing {}'.format(filename))
        table.write(filename, format='ascii.ecsv', overwrite=True)


def make_source_table(data):
    models = data['models']
    rows = []
    for source_idx, model in enumerate(models):
        # if source_idx == 10: break

        row = OrderedDict()
        row['file'] = data['tag']
        row['name'] = model.name()
        add_source_info_spatial(model.spatial(), row)
        add_source_info_spectral(model.spectral(), row)
        rows.append(row)

    meta = OrderedDict()
    meta['tag'] = data['tag']
    meta['filename'] = data['filename']
    table = Table(rows=rows, meta=meta, names=list(rows[0].keys()))
    # table.info('stats')
    return table


def add_source_info_spatial(spatial, row):
    spatial_type = spatial.type()

    if spatial_type in ['PointSource', 'SkyDirFunction']:
        spatial_type = 'PointSource'
        ra = spatial.ra()
        dec = spatial.dec()
        size = np.nan
    elif spatial_type in ['DiffuseMap', 'SpatialMap']:
        spatial_type = 'DiffuseMap'
        ra, dec, glon, glat = np.nan, np.nan, np.nan, np.nan
        size = np.nan
    elif spatial_type in ['RadialGaussian', 'GaussFunction']:
        spatial_type = 'RadialGaussian'
        ra = spatial.ra()
        dec = spatial.dec()
        size = spatial.sigma()
    elif spatial_type == 'RadialShell':
        ra = spatial.ra()
        dec = spatial.dec()
        size = spatial.radius()
    else:
        raise ValueError('Invalid: {}'.format(spatial_type))

    pos = SkyCoord(ra, dec, unit='deg').galactic
    glon = pos.l.deg
    glat = pos.b.deg

    row['spatial_type'] = spatial_type
    row['ra'] = ra
    row['dec'] = dec
    row['glon'] = glon
    row['glat'] = glat
    row['size'] = size


def add_source_info_spectral(spectral, row):
    row['spectral_type'] = spectral.type()

    emin = gammalib.GEnergy(1, 'TeV')
    emax = gammalib.GEnergy(10, 'TeV')
    flux = spectral.flux(emin, emax)
    row['flux_1_10'] = flux


class GPSSkyModel:
    def __init__(self, data):
        self.data = data
        self.tags = [
            'gamma-cat',
            'image_sources',
            'pwn',
            'snr',
            'binaries',
            'pulsars',
        ]

    @classmethod
    def load_sky_models(cls):
        data = []

        tag = 'gamma-cat'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table))

        tag = 'image_sources'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table))

        tag = 'pwn'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.ecsv'
        table_in = Table.read(filename, format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table, table_in=table_in))

        tag = 'snr'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        filename = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.ecsv'
        table_in = Table.read(filename, format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table, table_in=table_in))

        tag = 'binaries'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table))

        tag = 'pulsars'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pulsars/ctadc_skymodel_gps_sources_pulsars.xml'
        models = gammalib.GModels(filename)
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        data.append(dict(tag=tag, filename=filename, models=models, table=table))

        return cls(data=data)

    def get_component(self, tag):
        for component in self.data:
            if tag == component['tag']:
                return component
        raise ValueError('Invalid component tag: {}'.format(tag))

    def get_components(self, tags=None, skip=None):
        if tags is None:
            tags = self.tags.copy()
        if skip is not None:
            for tag in skip:
                tags.remove(tag)

        for tag in tags:
            yield self.get_component(tag)

    def print_summary(self):
        rows = []
        for component in self.get_components():
            row = OrderedDict()
            row['tag'] = component['tag']
            row['n_sources'] = len(component['models'])
            row['flux_1_10'] = compute_total_flux(component['models'])
            rows.append(row)

        table = Table(rows=rows, names=rows[0].keys())
        filename = 'ctadc_skymodel_gps_sources_summary.csv'
        log.info('Writing {}'.format(filename))
        table.write(filename, format='ascii.csv', overwrite=True)

    def plot_sky_positions(self):
        fig, ax = plt.subplots(figsize=(15, 5))
        for component in self.get_components():
            table = component['table']
            x = Angle(table['glon'], 'deg').wrap_at('180d').deg
            y = table['glat']
            if component['tag'] in ['pwn', 'snr']:
                s, alpha = 10, 0.3
                alpha
            else:
                s, alpha = 30, 0.8
            ax.scatter(x, y, label=component['tag'], s=s, alpha=alpha)

        ax.legend(loc='best')
        ax.set_xlim(180, -180)
        ax.set_ylim(-60, 60)
        ax.set_xlabel('GLON (deg)')
        ax.set_ylabel('GLAT (deg)')
        ax.grid()
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_sky_positions.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)

        ax.set_ylim(-8, 8)
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_sky_positions_gps.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_glon(self):
        fig, ax = plt.subplots()
        bins = np.arange(-180, 181, 1)
        for component in self.get_components(skip=['image_sources']):
            table = component['table']
            points = Angle(table['glon'], 'deg').wrap_at('180d').deg
            hist = Histogram.from_points(bins=bins, points=points)
            hist = hist.normalise(method='int', value=1)
            hist = hist.smooth(sigma=3)
            hist.plot_smooth(
                ax=ax, label=component['tag'], alpha=0.8,
            )

        ax.legend(loc='best')
        ax.set_xlim(180, -180)
        ax.set_xlabel('GLON (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_glon.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_glat(self):
        fig, ax = plt.subplots()
        bins = np.arange(-10, 10.1, 0.1)
        for component in self.get_components(skip=['image_sources']):
            table = component['table']
            points = Angle(table['glat'], 'deg').deg
            hist = Histogram.from_points(bins=bins, points=points)
            hist = hist.smooth(sigma=3)
            hist = hist.normalise(method='int', value=1)
            hist.plot_smooth(
                ax=ax, label=component['tag'], alpha=0.8,
            )

        ax.legend(loc='best')
        ax.set_xlim(-10, 10)
        ax.set_xlabel('GLAT (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_glat.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_size(self):
        fig, ax = plt.subplots()
        bins = np.arange(0, 3, 0.05)
        for component in self.get_components(tags=['gamma-cat', 'pwn', 'snr']):
            table = component['table']
            vals = table['size']
            if component['tag'] == 'gamma-cat':
                vals = np.nan_to_num(vals)  # put gamma-cat point sources at size=0
            # TODO: this is still failing because the size column was changed in the SNR table.
            # Change back in that table, or here in this script.
            ax.hist(
                vals, bins=bins, label=component['tag'], histtype='step',
                alpha=0.8, normed=True,
            )

        ax.legend(loc='best')
        # ax.set_xlim(bins[0], bins-[1])
        ax.set_xlim(0, 2)
        ax.set_xlabel('Source apparent size (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_size.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_size_physical(self):
        fig, ax = plt.subplots()
        bins = 30  # np.arange(0, 3, 0.05)
        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            vals = table['size_physical']
            ax.hist(
                vals, bins=bins, label=component['tag'], histtype='step',
                alpha=0.8, normed=True,
            )

        ax.legend(loc='best')
        # ax.set_xlim(bins[0], bins-[1])
        # ax.set_xlim(0, 2)
        ax.set_xlabel('Source physical size (pc)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_physical_size.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def check_snr_size(self):
        snr = self.get_component('snr')

        size_in = snr['table_in']['size'].quantity.to('deg')
        size_out = snr['table']['size'] * u.deg

        fig, ax = plt.subplots()
        bins = np.arange(0, 5, 0.05)
        ax.hist(size_in.value, bins=bins, label=['snr input'], histtype='step', normed=True)
        ax.hist(size_out.value, bins=bins, label=['snr output'], histtype='step', normed=True)

        ax.legend(loc='best')
        ax.set_xlabel('Source apparent size (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_size_snr_check.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_xy(self):
        fig, ax = plt.subplots(figsize=(7, 7))

        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            x = table['galactocentric_x'].quantity.to('kpc').value
            y = table['galactocentric_y'].quantity.to('kpc').value
            ax.scatter(x, y, label=component['tag'], s=10, alpha=0.5)

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric X (kpc)')
        ax.set_ylabel('Galactocentric Y (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_xy.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_xz(self):
        fig, ax = plt.subplots(figsize=(15, 5))

        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            x = table['galactocentric_x'].quantity.to('kpc').value
            z = table['galactocentric_z'].quantity.to('kpc').value
            ax.scatter(x, z, label=component['tag'], s=10, alpha=0.5)

        ax.set_ylim(-3, 3)
        ax.set_xlabel('Galactocentric X (kpc)')
        ax.set_ylabel('Galactocentric Z (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_xz.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_z(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            bins = np.arange(-2, 2, 0.05)
            vals = table['galactocentric_z'].quantity.to('kpc').value
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric Z (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_z.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_r(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            bins = np.arange(0, 20, 1)
            vals = table['galactocentric_r']
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric radius (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_r.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_distance(self):
        fig, ax = plt.subplots()

        bins = 50  # np.arange(0, 20, 1)

        for component in self.get_components(tags=['pwn', 'snr']):
            table = component['table_in']
            vals = table['distance']
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Distance to Earth (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_distance.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def get_logn_logs(self, quantity, variant):
        crab_1_10 = CrabSpectrum().model.integral(1 * u.TeV, 10 * u.TeV).to('cm-2 s-1').value
        bins = np.logspace(-6.1, 0.7, 100)

        hists = []
        for component in self.get_components():
            table = component['table']
            points = table['flux_1_10'] / crab_1_10
            hist = Histogram.from_points(bins=bins, points=points)

            if quantity == 'n':
                pass
            elif quantity == 'f':
                hist.vals = hist.vals * hist.bin_centers
            else:
                raise ValueError('Invalid quantity: {}'.format(quantity))

            if variant == 'diff':
                pass
            elif variant == 'int':
                hist.vals = np.cumsum(hist.vals[::-1])[::-1]
            else:
                raise ValueError('Invalid variant: {}'.format(variant))

            hists.append(hist)

        return HistogramStack(hists)

    def plot_logn_logs(self, quantity, variant, sigma):
        fig, ax = plt.subplots()

        hists = self.get_logn_logs(quantity, variant)
        kwargs = dict(alpha=0.8, lw=2)
        hists.total.smooth(sigma=sigma).plot_smooth(ax, label='total', **kwargs)
        for tag, hist in zip(self.tags, hists.hists):
            hist.smooth(sigma=sigma).plot_smooth(ax, label=tag, **kwargs)

        ax.set_xlabel('Integral flux 1-10 TeV (% Crab)')

        if quantity == 'n':
            ax.set_ylabel('Number of sources')
        elif quantity == 'f':
            ax.set_ylabel('Flux of sources')
        else:
            raise ValueError('Invalid quantity: {}'.format(quantity))

        ax.legend(loc='best')
        ax.grid('on')

        ax.semilogx()
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_log{}_logs_{}.png'.format(quantity, variant)
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)

        ax.loglog()
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_log{}_logs_{}_logscale.png'.format(quantity, variant)
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_all_spectral_models(self):
        for component in self.get_components():
            plot_spectral_models(component)


def plot_spectral_models(component):
    log.info('Starting plot_spectral_models for: {}'.format(component['tag']))
    fig, ax = plt.subplots()

    models = component['models']
    for idx, model in enumerate(models):
        # if idx == 100: break
        # log.debug('Plotting spectral model for: {}'.format(model.name()))
        plot_kwargs = dict(color='black')
        if component['tag'] in ['pwn', 'snr']:
            plot_kwargs['alpha'] = 0.2
            plot_kwargs['lw'] = 1
        else:
            plot_kwargs['alpha'] = 0.3
            plot_kwargs['lw'] = 2
        plot_spectral_model(model.spectral(), ax, plot_kwargs)

    ax.set_title('{} spectra'.format(component['tag']))
    ax.loglog()
    ax.set_xlabel('Energy (TeV)')
    ax.set_ylabel('e2dnde (erg cm-2 s-1)')
    ax.set_ylim(2e-18, 5e-10)
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_spectra_{}.png'.format(component['tag'])
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)
    plt.close(fig)


def plot_spectral_model(model, ax, plot_kwargs):
    energies = Energy.equal_log_spacing(emin=0.02, emax=100, unit='TeV', nbins=40)
    fluxes = []
    for energy in energies:
        energy_gammalib = gammalib.GEnergy(energy.value, 'TeV')
        fluxes.append(model.eval(energy_gammalib))

    dnde = u.Quantity(fluxes, 'cm-2 s-1 MeV-1')
    e2dnde = (energies ** 2 * dnde).to('erg cm-2 s-1')

    ax.plot(energies.value, e2dnde.value, **plot_kwargs)


def compute_total_flux(models):
    flux_total = 0
    for model in models:
        emin = gammalib.GEnergy(1, 'TeV')
        emax = gammalib.GEnergy(10, 'TeV')
        flux = model.spectral().flux(emin, emax)
        flux_total += flux
    return flux_total


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    make_source_tables()

    gps = GPSSkyModel.load_sky_models()
    gps.print_summary()

    gps.plot_sky_positions()
    gps.plot_glon()
    gps.plot_glat()

    gps.plot_size()
    gps.plot_size_physical()

    gps.check_snr_size()

    gps.plot_galactic_xy()
    gps.plot_galactic_xz()
    gps.plot_galactic_z()
    gps.plot_galactic_r()
    gps.plot_distance()

    gps.plot_logn_logs(quantity='n', variant='diff', sigma=2)
    gps.plot_logn_logs(quantity='n', variant='int', sigma=None)
    gps.plot_logn_logs(quantity='f', variant='diff', sigma=2)
    gps.plot_logn_logs(quantity='f', variant='int', sigma=None)

    gps.plot_all_spectral_models()
