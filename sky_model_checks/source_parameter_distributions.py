"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
import logging
from collections import OrderedDict
from pathlib import Path
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
        # import IPython; IPython.embed(); 1 / 0
        add_source_info_spatial(model.spatial(), row)
        add_source_info_spectral(model.spectral(), row)
        rows.append(row)

    meta = OrderedDict(
        tag=data['tag'],
        filename=data['filename'],
    )
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


def load_sky_models():
    log.info('Starting load_sky_models')
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
    table_in['POS_X'] = table_in['x']
    table_in['POS_Y'] = table_in['y']
    table_in['POS_Z'] = table_in['z']
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

    return data


def compute_total_flux(models):
    flux_total = 0
    for model in models:
        emin = gammalib.GEnergy(1, 'TeV')
        emax = gammalib.GEnergy(10, 'TeV')
        flux = model.spectral().flux(emin, emax)
        flux_total += flux
    return flux_total


def print_skymodel_summary(data):
    log.info('Starting print_skymodel_summary')

    table = Table()
    table['tag'] = [_['tag'] for _ in data]
    table['n_sources'] = [len(_['models']) for _ in data]
    table['flux_1_10'] = [compute_total_flux(_['models']) for _ in data]

    table.pprint()


def plot_sky_positions(data):
    log.info('Starting plot_sky_positions')
    fig, ax = plt.subplots(figsize=(15, 5))
    for component in data:
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


def plot_glon_distribution(data):
    log.info('Starting plot_glon_distribution')
    fig, ax = plt.subplots()
    bins = np.arange(-180, 181, 5)
    for component in data:
        if component['tag'] == 'image_sources':
            continue
        table = component['table']
        vals = Angle(table['glon'], 'deg').wrap_at('180d').deg
        ax.hist(
            vals, bins=bins, label=component['tag'], histtype='step',
            alpha=0.8, normed=True,
        )

    ax.legend(loc='best')
    ax.set_xlim(180, -180)
    ax.set_xlabel('GLON (deg)')
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_glon.png'
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_glat_distribution(data):
    log.info('Starting plot_glat_distribution')
    fig, ax = plt.subplots()
    bins = np.arange(-10, 10.1, 0.5)
    for component in data:
        if component['tag'] == 'image_sources':
            continue
        table = component['table']
        vals = Angle(table['glat'], 'deg').deg
        ax.hist(
            vals, bins=bins, label=component['tag'], histtype='step',
            alpha=0.8, normed=True,
        )

    ax.legend(loc='best')
    ax.set_xlim(-10, 10)
    ax.set_xlabel('GLAT (deg)')
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_glat.png'
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_size_distribution(data):
    log.info('Starting plot_size_distribution')
    fig, ax = plt.subplots()
    bins = np.arange(0, 3, 0.1)
    for component in data:
        if component['tag'] not in ['gamma-cat', 'pwn', 'snr']:
            continue
        table = component['table']
        vals = Angle(table['size'], 'deg').deg
        ax.hist(
            vals, bins=bins, label=component['tag'], histtype='step',
            alpha=0.8, normed=True,
        )

    ax.legend(loc='best')
    # ax.set_xlim(bins[0], bins-[1])
    ax.set_xlim(0, 2)
    ax.set_xlabel('Source size (deg)')
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_size.png'
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_galactic_xy(data):
    log.info('Starting plot_galactic_xy')
    fig, ax = plt.subplots(figsize=(7, 7))

    for component in data:
        if component['tag'] not in ['snr', 'pwn']:
            continue

        table = component['table_in']
        x = table['POS_X'].quantity.to('kpc').value
        y = table['POS_Y'].quantity.to('kpc').value
        ax.scatter(x, y, label=component['tag'], s=10, alpha=0.5)

    # ax.set_xlim(0, 2)
    ax.set_xlabel('Galactocentric X (kpc)')
    ax.set_xlabel('Galactocentric Y (kpc)')
    ax.legend(loc='best')
    fig.tight_layout()

    filename = 'ctadc_skymodel_gps_sources_galactic_xy.png'
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_galactic_z(data):
    log.info('Starting plot_galactic_z')
    fig, ax = plt.subplots()

    for component in data:
        if component['tag'] not in ['snr', 'pwn']:
            continue

        table = component['table_in']
        bins = np.arange(-2, 2, 0.05)
        vals = table['POS_Z'].quantity.to('kpc').value
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


def plot_logn_logs():
    log.info('Starting plot_logn_logs')
    fig, ax = plt.subplots(figsize=(15, 5))

    crab_1_10 = CrabSpectrum().model.integral(1 * u.TeV, 10 * u.TeV).to('cm-2 s-1').value
    bins = np.logspace(-10, 4, 200)
    left, right, width = bins[:-1], bins[1:], np.diff(bins)

    for component in data:
        table = component['table']
        flux = table['flux_1_10'] / crab_1_10
        hist = np.histogram(flux, bins=bins)[0]
        hist = gaussian_filter1d(hist.astype(float), sigma=2)
        hist = hist / np.max(hist)

        # Plotting a histogram nicely ourselves is nontrivial
        # Example from http://stackoverflow.com/a/18611135/498873
        x = np.array([left, right]).T.flatten()
        y = np.array([hist, hist]).T.flatten()

        ax.plot(x, y, label=component['tag'], alpha=0.8, lw=2)

        # ax.bar(
        #     left=left, height=hist, width=width,
        #     label=component['tag'], alpha=0.7,
        #     align='edge', color='none',
        # )

    # import IPython; IPython.embed(); 1/0
    ax.set_xlabel('Integral flux 1-10 TeV (% Crab)')
    ax.semilogx()
    ax.legend(loc='best')
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_logn_logs_diff.png'
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_all_spectral_models(data):
    log.info('Starting plot_all_spectral_models')
    for component in data:
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
    # import IPython; IPython.embed()
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_spectra_{}.png'.format(component['tag'])
    log.info('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_spectral_model(model, ax, plot_kwargs):
    energies = Energy.equal_log_spacing(emin=0.02, emax=100, unit='TeV', nbins=40)
    fluxes = []
    for energy in energies:
        energy_gammalib = gammalib.GEnergy(energy.value, 'TeV')
        fluxes.append(model.eval(energy_gammalib))

    dnde = u.Quantity(fluxes, 'cm-2 s-1 MeV-1')
    e2dnde = (energies ** 2 * dnde).to('erg cm-2 s-1')

    ax.plot(energies.value, e2dnde.value, **plot_kwargs)
    # 1/0


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)

    # make_source_tables()
    data = load_sky_models()
    # print_skymodel_summary(data)

    # plot_sky_positions(data)
    # plot_glon_distribution(data)
    # plot_glat_distribution(data)

    # plot_size_distribution(data)

    # plot_galactic_xy(data)
    # plot_galactic_z(data)

    # plot_logn_logs()

    plot_all_spectral_models(data)
