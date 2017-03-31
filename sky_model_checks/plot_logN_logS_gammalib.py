

def compute_fluxes():
    import gammalib
    filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml'
    models = gammalib.GModels(filename)
    fluxes = []
    for model in models:
        # import IPython;IPython.embed(); 1/0
        emin = gammalib.GEnergy(1, 'TeV')
        emax = gammalib.GEnergy(10, 'TeV')
        flux = model.spectral().flux(emin, emax)
        fluxes.append(flux)

    return fluxes


def make_plot(fluxes):
    import matplotlib.pyplot as plt
    plt.hist(fluxes)
    filename = 'logN_logS_gammalib.png'
    print('Writing {}'.format(filename))
    plt.savefig(filename)


if __name__ == '__main__':
    fluxes = compute_fluxes()
    make_plot(fluxes)
