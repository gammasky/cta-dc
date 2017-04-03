"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
from collections import OrderedDict
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
import gammalib


def make_source_tables():
    parts = [
        dict(tag='gamma-cat', filename='gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'),
        dict(tag='image_sources', filename='image_sources/ctadc_skymodel_gps_sources_images.xml'),
        dict(tag='pwn', filename='pwn/ctadc_skymodel_gps_sources_pwn.xml'),
        dict(tag='snr', filename='snrs/ctadc_skymodel_gps_sources_snr_2.xml'),
        dict(tag='binaries', filename='binaries/ctadc_skymodel_gps_sources_binaries.xml'),
        # TODO: pulsars
    ]
    for data in parts:
        filename = '../sky_model/' + data['filename']
        print('Reading {}'.format(filename))
        data['models'] = gammalib.GModels(filename)
        table = make_source_table(data)
        filename = filename.replace('.xml', '_summary.ecsv')
        print('Writing {}'.format(filename))
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
    table.info('stats')
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
    data = []

    tag = 'gamma-cat'
    filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'
    models = gammalib.GModels(filename)
    table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
    data.append(dict(tag=tag, filename=filename, models=models, table=table))

    tag = 'image_sources'
    filename = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml'
    models = gammalib.GModels(filename)
    table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
    data.append(dict(tag=tag, filename=filename, models=models, table=table))

    tag = 'pwn'
    filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml'
    models = gammalib.GModels(filename)
    table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
    data.append(dict(tag=tag, filename=filename, models=models, table=table))

    tag = 'snr'
    filename = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml'
    models = gammalib.GModels(filename)
    table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
    data.append(dict(tag=tag, filename=filename, models=models, table=table))

    tag = 'binaries'
    filename = '../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml'
    models = gammalib.GModels(filename)
    table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
    data.append(dict(tag=tag, filename=filename, models=models, table=table))

    # TODO: add pulsars

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
    table = Table()
    table['tag'] = [_['tag'] for _ in data]
    table['n_sources'] = [len(_['models']) for _ in data]
    table['flux_1_10'] = [compute_total_flux(_['models']) for _ in data]

    table.pprint()


def plot_sky_positions(data):
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
    print('Writing {}'.format(filename))
    fig.savefig(filename)

    ax.set_ylim(-8, 8)
    fig.tight_layout()
    filename = 'ctadc_skymodel_gps_sources_sky_positions_gps.png'
    print('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_glon_distribution(data):
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
    print('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_glat_distribution(data):
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
    print('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_size_distribution(data):
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
    print('Writing {}'.format(filename))
    fig.savefig(filename)


def plot_logn_logs():
    # import IPython; IPython.embed(); 1/0
    fig, ax = plt.subplots(figsize=(15, 5))
    for component in data:
        table = component['table']
        table.rename_column('flux_1_10', 'S')
        dist = FluxDistribution(table, label=component['tag'])
        dist.plot_integral_count()

    filename = 'ctadc_skymodel_gps_sources_logn_logs.png'
    print('Writing {}'.format(filename))
    fig.savefig(filename)


if __name__ == '__main__':
    # make_source_tables()
    data = load_sky_models()
    # print_skymodel_summary(data)

    # plot_sky_positions(data)
    # plot_glon_distribution(data)
    # plot_glat_distribution(data)

    plot_size_distribution(data)

    # plot_logn_logs()
