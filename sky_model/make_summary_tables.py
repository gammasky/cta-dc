"""Make ECSV format tables summarising parameters for all CTA 1DC sources.

This is useful mostly for checking / plotting the models,
because it's quick and simple to access and analysis the info,
and some info (like position in the Galaxy)
isn't available in the model XML files.
"""
import logging
from collections import OrderedDict
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import gammalib

log = logging.getLogger(__name__)


def make_all_summary_tables():
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
        filename = data['filename']
        log.info('Reading {}'.format(filename))
        data['models'] = gammalib.GModels(filename)
        table = make_summary_table(data)
        filename = filename.replace('.xml', '_summary.ecsv')
        log.info('Writing {}'.format(filename))
        table.write(filename, format='ascii.ecsv', overwrite=True)


def make_summary_table(data):
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
    flux_crab = flux
    row['int_flux_above_1TeV'] = flux


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    make_all_summary_tables()
