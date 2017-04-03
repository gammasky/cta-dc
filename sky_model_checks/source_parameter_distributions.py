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
from astropy.coordinates import SkyCoord
import gammalib


def make_source_tables():
    filenames = [
        'gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml',
        'image_sources/ctadc_skymodel_gps_sources_images.xml',
        'pwn/ctadc_skymodel_gps_sources_pwn.xml',
        'snrs/ctadc_skymodel_gps_sources_snr_2.xml',
        'binaries/ctadc_skymodel_gps_sources_binaries.xml',
    ]
    for filename in filenames:
        filename = '../sky_model/' + filename
        print('Reading {}'.format(filename))
        models = gammalib.GModels(filename)

        table = make_source_table(models)
        filename = filename.replace('.xml', '_summary.ecsv')
        print('Writing {}'.format(filename))
        table.write(filename, format='ascii.ecsv', overwrite=True)


def make_source_table(models):
    rows = []
    for source_idx, model in enumerate(models):
        if source_idx == 10: break

        row = OrderedDict()
        row['name'] = model.name()
        # import IPython; IPython.embed(); 1 / 0
        make_source_table_spatial(model.spatial(), row)
        rows.append(row)

    meta = OrderedDict(
        info='This is TODO',
    )
    table = Table(rows=rows, meta=meta, names=list(rows[0].keys()))
    return table


def make_source_table_spatial(spatial, row):
    if isinstance(spatial, gammalib.GModelSpatialDiffuseMap):
        ra, dec, glon, glat = np.nan, np.nan, np.nan, np.nan
    else:
        ra = spatial.ra()
        dec = spatial.dec()
        pos = SkyCoord(ra, dec, unit='deg').galactic
        glon = pos.l.deg
        glat = pos.b.deg

    row['ra'] = ra
    row['dec'] = dec
    row['glon'] = glon
    row['glat'] = glat


def load_sky_models():
    data = []

    tag = 'gamma-cat'
    filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

    tag = 'image_sources'
    filename = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

    tag = 'pwn'
    filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

    tag = 'snr'
    filename = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

    tag = 'binaries'
    filename = '../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

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
    fig, ax = plt.subplot()

    filename = 'ctadc_skymodel_gps_sources_sky_positions.png'
    print('Writing {}'.format(filename))
    fig.savefig(filename)


if __name__ == '__main__':
    make_source_tables()
    # data = load_sky_models()
    # print_skymodel_summary(data)
    # plot_sky_positions(data)
