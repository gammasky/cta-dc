"""
Convert Gamma-cat to model XML format.

Use Gammalib for now to make sure the XML format is OK.
"""
import gammalib
from gammapy.catalog import SourceCatalogGammaCat

def gammacat_source_to_gammalib_model(source):
    # import IPython; IPython.embed(); 1/0

    gammapy_spectral = source.spectral_model
    gammalib_spectral = gammacat_source_to_gammalib_model_spectral(gammapy_spectral)

    gammapy_spatial = source.spatial_model()
    gammalib_spatial = gammacat_source_to_gammalib_model_spatial(gammapy_spatial)

    gammalib_model = gammalib.GModelSky(gammalib_spatial, gammalib_spectral)
    gammalib_model.name(source.name)
    return gammalib_model


def gammacat_source_to_gammalib_model_spatial(gammapy_spatial):
    pos = gammalib.GSkyDir()

    kind = gammapy_spatial.__class__.__name__
    if kind == 'Gaussian2D':
        pos.lb_deg(gammapy_spatial.x_mean.value, gammapy_spatial.y_mean.value)
        sigma = gammapy_spatial.x_stddev.value
        # TODO: use gammalib.GModelSpatialEllipticalGauss
        gammalib_spatial = gammalib.GModelSpatialRadialGauss(pos, sigma)
    elif kind == 'Shell2D':
        # import IPython; IPython.embed(); 1/0
        # gammalib_spatial = gammalib.GModelSpatialPointSource(pos)
        pos.lb_deg(gammapy_spatial.x_0.value, gammapy_spatial.y_0.value)
        radius = gammapy_spatial.r_in.value # TODO: is this inner or outer radius?
        width= gammapy_spatial.width.value
        gammalib_spatial = gammalib.GModelSpatialRadialShell(pos, radius, width)
    else:
        raise NotImplementedError('Not supported: {}'.format(kind))

    return gammalib_spatial


def gammacat_source_to_gammalib_model_spectral(gammapy_spectral):
    kind = gammapy_spectral.__class__.__name__
    if kind == 'PowerLaw':
        prefactor = gammapy_spectral.parameters['amplitude'].quantity.to('cm-2 s-1 MeV-1').value
        index = gammapy_spectral.parameters['index'].value
        pivot = gammalib.GEnergy(gammapy_spectral.parameters['reference'].quantity.to('MeV').value, 'MeV')
        gammalib_spectral = gammalib.GModelSpectralPlaw(prefactor, index, pivot)
    elif kind == 'PowerLaw2':
        # import IPython; IPython.embed(); 1/0
        flux = gammapy_spectral.parameters['amplitude'].quantity.to('cm-2 s-1').value
        index = gammapy_spectral.parameters['index'].value
        emin = gammalib.GEnergy(gammapy_spectral.parameters['emin'].quantity.to('MeV').value, 'MeV')
        emax = gammalib.GEnergy(gammapy_spectral.parameters['emax'].quantity.to('MeV').value, 'MeV')
        gammalib_spectral = gammalib.GModelSpectralPlawPhotonFlux(flux, index, emin, emax)
    elif kind == 'ExponentialCutoffPowerLaw':
        prefactor = gammapy_spectral.parameters['amplitude'].quantity.to('cm-2 s-1 MeV-1').value
        index = gammapy_spectral.parameters['index'].value
        pivot = gammalib.GEnergy(gammapy_spectral.parameters['reference'].quantity.to('MeV').value, 'MeV')
        cutoff = gammalib.GEnergy(1. / gammapy_spectral.parameters['lambda_'].quantity.to('MeV-1').value, 'MeV')
        gammalib_spectral = gammalib.GModelSpectralExpPlaw(prefactor, index, pivot, cutoff)
    else:
        raise NotImplementedError('Not supported: {}'.format(kind))
    return gammalib_spectral


def gammacat_to_xml_version1():
    cat = SourceCatalogGammaCat()

    models = gammalib.GModels()
    for source_idx in range(len(cat.table)):
        # For debugging
        # if source_idx > 5: break

        source = cat[source_idx]
        # print('\n\n*** Source: {} ***\n'.format(source.name))
        try:
            model = gammacat_source_to_gammalib_model(source)
            models.append(model)
        except ValueError as exc:
            # import IPython; IPython.embed(); 1/0
            if source.data['where'] != 'egal':
                print('Skipped non-EGAL source: {} (gammacat source_id={})'.format(source.name, source.data['source_id']))
            # print(exc)
            # raise

    filename = 'ctadc_skymodel_gps_sources_gamma-cat2.xml'
    print('Writing {}'.format(filename))
    models.save(filename)


def gammacat_to_xml_version2():
    cat = SourceCatalogGammaCat()
    xml = ''

    for source_idx in range(len(cat.table)):
        source = cat[source_idx]
        xml += source_to_xml(source)

    return xml


if __name__ == '__main__':
    gammacat_to_xml_version1()
