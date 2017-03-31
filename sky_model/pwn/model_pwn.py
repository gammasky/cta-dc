"""
Make synthetic PNW population for CTA 1DC
"""
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from gammapy.astro.population import make_base_catalog_galactic, add_observed_parameters
from gammapy.utils.random import get_random_state
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum


def define_flux_crab_above_energy(emin=1 * u.TeV, emax=10 * u.TeV):
    crab = CrabSpectrum('meyer').model
    crabMAGIC = LogParabola(amplitude=3.23e-11 * u.Unit('cm-2 s-1 TeV-1'), reference=1 * u.TeV, alpha=2.47, beta=0.24)
    crab_flux_above_1TeV = crabMAGIC.integral(emin=emin, emax=emax)
    crab_flux_above_1TeV_model = crab.integral(emin=emin, emax=emax)
    return crab_flux_above_1TeV, crab_flux_above_1TeV_model


def flux_amplitude_from_energy_flux(alpha, beta, energy_flux):
    spec = LogParabola(
        amplitude=1 * u.Unit('cm-2 s-1 TeV-1'),
        reference=1 * u.TeV,
        alpha=alpha,
        beta=-beta,
    )

    # Assumption before for `energy_flux` was energy band 1 to 10 TeV
    emin_int = 1 * u.TeV
    emax_int = 10 * u.TeV
    pivot_energy = 1 * u.TeV
    energy_flux_standard_candle = spec.energy_flux(emin=emin_int, emax=emax_int)
    flux_at_1TeV = energy_flux / energy_flux_standard_candle * u.Unit(
        'cm-2 s-1 TeV-1')  # * spec.parameters['amplitude'].quantity
    flux_at_1TeV = flux_at_1TeV.to('cm-2 s-1 MeV-1')

    flux_above_1TeV = LogParabola(amplitude=flux_at_1TeV, reference=pivot_energy, alpha=alpha, beta=beta)
    flux_above_1TeV = flux_above_1TeV.integral(emin=emin_int, emax=10 * u.TeV)
    flux_above_1TeV = flux_above_1TeV.to('cm-2 s-1')
    # evaluating Crab flux at and above 1 TeV by using MAGIC Crab spectrum from JHEA 2015
    crab_flux_above_1TeV, crab_flux_above_1TeV_model = define_flux_crab_above_energy(emin=1 * u.TeV, emax=10 * u.TeV)
    crabMAGIC = LogParabola(amplitude=3.23e-11 * u.Unit('cm-2 s-1 TeV-1'), reference=1 * u.TeV, alpha=2.47, beta=0.24)
    # plt.figure()
    # min, max = [50*u.GeV,50*u.TeV]
    # crabMAGIC.plot(min,max)
    # crab = CrabSpectrum('meyer').model
    crab_flux_above_1TeV = crab_flux_above_1TeV.to('cm-2 s-1')
    crab_flux_above_1TeV_model = crab_flux_above_1TeV_model.to('cm-2 s-1')
    crab_flux_at_1TeV = crabMAGIC(pivot_energy).to('MeV-1 cm-2 s-1')

    # computing flux at and above 1 TeV in crab units
    flux_at_1TeV_cu = (flux_at_1TeV / crab_flux_at_1TeV).to('%')
    flux_above_1TeV_cu = (flux_above_1TeV / crab_flux_above_1TeV).to('%')
    flux_above_1TeV_cu_model = (flux_above_1TeV / crab_flux_above_1TeV_model).to('%')

    # print(flux_at_1TeV, crab_flux_at_1TeV, flux_at_1TeV_cu)
    # print('  ', flux_above_1TeV, crab_flux_above_1TeV, flux_above_1TeV_cu, alpha, beta)
    # print('  ', flux_above_1TeV, crab_flux_above_1TeV, flux_above_1TeV_cu, flux_above_1TeV_cu_model)

    return flux_at_1TeV, flux_at_1TeV_cu, flux_above_1TeV, flux_above_1TeV_cu


def make_pwn_table(
        n_sources=650, random_state=0,
        mean_extension=0.13, sigma_extension=0.1, mean_intrinsic=20, sigma_intrinsic=5,  # intrinsic_extension=50,
        mean_index_alpha=2.1, sigma_index_alpha=0.3, max_index_beta=0.5,
        mean_logluminosity=33.5, sigma_logluminosity=1.0, ):
    """Make a catalog of PWN.

    to be defined
    1. how many from LogN-LogS
    2. spatial distribution: spiral arms
    3. flux: Log10(luminosity) gaussian
    4. morphlogy: Gaussian distribution with the mean and the sigma derived from
    the the Gaus distribution of the size of the known PWNe
    5. spectral: logparabola (how to choose the parameters)

    Parameters
    -----------
    n_sources: int
        number of sources defined from the logN-logS of the PWN population
    rad_dis : callable
    rad_dis : callable
        Radial surface density distribution of sources.
        YusifovKucuk2004 is radial distribution of the surface density of pulsars in the galaxy
    mean_extension: float deg
        mean of the Gaussian distribution of the angular extension of the known PWN
    sigma_extension: float deg
        sigma of the Gaussian distribution of the angular extension of the known PWN
    intrinsic_extension: int in pc
        this is needed to compute the distance
    mean_index_alpha: float
        mean of the Gaus distribution of the alpha index under the assumption of a logparabola spectral model
    sigma_index_alpha: float
        mean of the Gaus distribution of the alpha index under the assumption of a logparabola spectral model
    max_index_beta: float
        maximum value of the beta index under the assumption of a logparabola spectral model.
        Assuming an uniform distribution for beta
    """
    random_state = get_random_state(random_state)

    table = make_base_catalog_galactic(
        n_sources=n_sources,
        max_age=0 * u.year,
        random_state=random_state,
    )
    table = add_observed_parameters(table)


    # Draw angular extension, then compute physical extension
    # intrinsic_extension = random_state.normal(mean_intrinsic,sigma_intrinsic,n_sources)
    intrinsic_extension = random_state.uniform(3, 60, n_sources)

    intrinsic_extension = intrinsic_extension * u.pc
    constant = 0.000291 / u.arcmin
    angular_extension = intrinsic_extension / (table['distance'] * constant)
    angular_extension = angular_extension.to('deg') / 2.0
    for idx in range(len(table)):
        if (angular_extension[idx] < 0.01 * u.Unit('deg')):
            # print('PROBLEM: ', angular_extension[idx])
            angular_extension[idx] = 0.01 * u.Unit('deg')

    for idx in range(len(table)):
        if (angular_extension[idx] < 0):
            print('ROBLEM: ', angular_extension[idx], intrinsic_extension[idx])
    # print('size ', intrinsic_extension, angular_extension, table['distance'])

    # Draw log parabola spectral model parameters
    alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
    # random_state.uniform(1.5, 2.6, n_sources)
    beta = random_state.uniform(0, max_index_beta, n_sources)

    # Define the luminosity

    logluminosity = random_state.normal(mean_logluminosity, sigma_logluminosity, n_sources)
    for idx in range(len(table)):
        if (logluminosity[idx] > 35):
            logluminosity[idx] = random_state.normal(mean_logluminosity, sigma_logluminosity, 1)
    luminosity = (10 ** logluminosity) * u.erg / u.second

    distance = table['distance'].quantity

    # integral sed between 1 and 10 TeV
    energy_flux = luminosity / (4 * np.pi * distance ** 2)
    # print(energy_flux)
    energy_flux = energy_flux.to('TeV cm-2 s-1')

    vals = []
    vals_cu = []
    int_vals = []
    int_vals_cu = []
    name = []
    for idx in range(len(table)):
        val, val_cu, int_val, int_val_cu = flux_amplitude_from_energy_flux(
            alpha[idx], beta[idx], energy_flux[idx],
        )
        vals.append(val)
        vals_cu.append(val_cu)
        int_vals.append(int_val)
        int_vals_cu.append(int_val_cu)
        source_name = 'pwn_{}'.format(idx)
        name.append(source_name)
    norm = u.Quantity(vals)
    norm_cu = u.Quantity(vals_cu)
    int_flux = u.Quantity(int_vals)
    int_flux_cu = u.Quantity(int_vals_cu)

    # import IPython; IPython.embed()
    table = table[['x', 'y', 'z', 'spiralarm', 'GLON', 'GLAT', 'RA', 'DEC', 'distance']]
    table['sigma'] = Column(angular_extension, description='Angular extension (deg)', unit='deg')
    table['spec_alpha'] = Column(alpha, description='Spectral model parameter (log parabola)')
    table['spec_beta'] = Column(beta, description='Spectral model parameter (log parabola)')
    table['spec_norm'] = Column(norm, description='Spectral model norm parameter (log parabola)', unit='MeV-1 s-1 cm-2')
    table['spec_norm_cu'] = Column(norm_cu, description='Spectral model norm parameter (log parabola) in crab units')
    table['int_flux_above_1TeV'] = Column(int_flux, description='Integral flux above 1 TeV ', unit='s-1 cm-2')
    table['int_flux_above_1TeV_cu'] = Column(int_flux_cu, description='Integral flux above 1 TeV in crab units')
    table['source_name'] = Column(name, description='source name')
    return table


if __name__ == '__main__':
    table = make_pwn_table()

    filename = 'ctadc_skymodel_gps_sources_pwn.ecsv'
    print('Writing {}'.format(filename))
    table.write(filename, format='ascii.ecsv', overwrite=True)
