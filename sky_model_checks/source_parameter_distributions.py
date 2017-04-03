"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
from astropy.table import Table


def load_sky_models():
    import gammalib
    data = []

    tag = 'gamma-cat'
    filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'
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

    tag = 'image_sources'
    filename = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml'
    models = gammalib.GModels(filename)
    data.append(dict(tag=tag, filename=filename, models=models))

    return data


def print_skymodel_summary(data):
    table = Table()
    table['tag'] = [_['tag'] for _ in data]
    table['n_sources'] = [len(_['models']) for _ in data]
    table.pprint()


if __name__ == '__main__':
    data = load_sky_models()
    print_skymodel_summary(data)
