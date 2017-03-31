




def compute_fluxes(filename):
    import gammalib
    crab_integral_flux_1to10TeV = 2.5520883250253037e-11

    models = gammalib.GModels(filename)
    fluxes = []
    fluxes_cu = []
    names = []
    for model in models:
        print(model)
        #import IPython;IPython.embed();
        name = model.name()
        names.append(name)
        emin = gammalib.GEnergy(1, 'TeV')
        emax = gammalib.GEnergy(10, 'TeV')
        flux = model.spectral().flux(emin, emax)
        flux_cu = (flux/crab_integral_flux_1to10TeV)*100
        fluxes.append(flux)
        fluxes_cu.append(flux_cu)
    return names, fluxes, fluxes_cu


def make_plot(num_plot, color_scale, fluxes_cu,fluxes_cu_2,fluxes_cu_3,fluxes_cu_4,fluxes_cu_5):
    import matplotlib.pyplot as plt
    import numpy as np
    num_bins = 100
    flux_min = 0.07
    flux_max = 110
    plt.figure()
    bins = np.logspace(np.log10(flux_min), np.log10(flux_max), num_bins)
    hist, bins_edges = np.histogram(fluxes_cu, bins=bins)
    y = np.insert(hist[::-1].astype('float64').cumsum(), 0, 0.01)
    plt.step(bins[::-1], y, color=color_scale[0], lw=2)
    if num_plot > 1:
        hist_2, bins_edges_2 = np.histogram(fluxes_cu_2, bins=bins)
        y_2 = np.insert(hist_2[::-1].astype('float64').cumsum(), 0, 0.01)
        plt.step(bins[::-1], y_2, color=color_scale[1], lw=2)
    if num_plot > 2:
        hist_3, bins_edges_3 = np.histogram(fluxes_cu_3, bins=bins)
        y_3 = np.insert(hist_3[::-1].astype('float64').cumsum(), 0, 0.01)
        plt.step(bins[::-1], y_3, color=color_scale[2], lw=2)
    if num_plot > 3:
        hist_4, bins_edges_4 = np.histogram(fluxes_cu_4, bins=bins)
        y_4 = np.insert(hist_4[::-1].astype('float64').cumsum(), 0, 0.01)
        plt.step(bins[::-1], y_4, color=color_scale[3], lw=2)
    if num_plot > 4:
        hist_5, bins_edges_5 = np.histogram(fluxes_cu_5, bins=bins)
        y_5 = np.insert(hist_5[::-1].astype('float64').cumsum(), 0, 0.01)
        plt.step(bins[::-1], y_5, color=color_scale[4], lw=2)

    plt.loglog()
    plt.ylim(0.8, 5000)
    plt.xlim(0.05, 100)
    plt.ylabel('#')
    plt.xlabel('F (1-10 TeV) [cu]')

    plt.plot([10, 10], [flux_min, 100], 'k-', lw=1, color='black')
    plt.plot([0.2, 0.2], [flux_min, 5000], 'k--', lw=1, color='black')


    #plt.hist(fluxes)
    filename = 'logN_logS_gammalib.png'
    print('Writing {}'.format(filename))
    plt.savefig(filename)

 #hist_all, bin_edges_sim_all = np.histogram(merged_table['int_flux_above_1TeV_cu'], bins=bins)

#    y = np.insert(hist[::-1].astype('float64').cumsum(), 0, 0.01)



if __name__ == '__main__':
    num_plot = 5
    color_scale = ['red','blue','green','magenta','black']

    filename_pwn = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml'
    names_pwn, fluxes_pwn, fluxes_cu_pwn = compute_fluxes(filename_pwn)

    filename_binaries = '../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml'
    names_binaries, fluxes_binaries, fluxes_cu_binaries = compute_fluxes(filename_binaries)

    filename_snr = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml'
    names_snr, fluxes_snr, fluxes_cu_snr = compute_fluxes(filename_snr)

    filename_gammacat = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat.xml'
    #names_gammacat, fluxes_gammacat, fluxes_cu_gammacat = compute_fluxes(filename_gammacat)

    filename_images = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml'
    names_images, fluxes_images, fluxes_cu_images = compute_fluxes(filename_images)

    make_plot(num_plot, color_scale, fluxes_cu_pwn, fluxes_cu_snr, fluxes_cu_binaries, fluxes_cu_images, fluxes_cu_pwn)
