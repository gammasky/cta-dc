"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
import numpy as np
from astropy.table import Table
import gammalib


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


if __name__ == '__main__':
    data = load_sky_models()
    print_skymodel_summary(data)
