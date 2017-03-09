import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.table import Table, Column
from astropy.modeling.models import Gaussian2D
from gammapy.astro.population import make_base_catalog_galactic, add_observed_parameters
from gammapy.utils.random import get_random_state
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum

def flux_amplitude_from_energy_flux(alpha, beta, energy_flux):
    spec = LogParabola(
        amplitude=1 * u.Unit('cm-2 s-1 TeV-1'),
        reference=1 * u.TeV,
        alpha=alpha,
        beta=beta,
    )

    # Assumption before for `energy_flux` was energy band 1 to 10 TeV
    energy_flux_standard_candle = spec.energy_flux(emin=1 * u.TeV, emax=10 * u.TeV)
    amplitude = energy_flux / energy_flux_standard_candle * u.Unit('cm-2 s-1 TeV-1')#* spec.parameters['amplitude'].quantity
    amplitude = amplitude.to('cm-2 s-1 MeV-1')
    crab = CrabSpectrum('meyer').model
    crab = CrabSpectrum('meyer').model
    crab_energy = 1*u.TeV
    crab_flux_at_1TeV =  crab(crab_energy).to('MeV-1 cm-2 s-1')
    amplitude_crab = amplitude/crab_flux_at_1TeV * 100

    print(energy_flux_standard_candle, energy_flux, amplitude, crab_flux_at_1TeV,amplitude_crab)

    return amplitude, amplitude_crab


def make_pwn_table(
        n_sources=450, random_state=0,
        mean_extension=0.13, sigma_extension=0.1, mean_intrinsic = 20, sigma_intrinsic=10, #intrinsic_extension=50,
        mean_index_alpha=1.8, sigma_index_alpha=0.27, max_index_beta=0.5,
        mean_logluminosity=33.9, sigma_logluminosity=0.5,
):
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
    table = make_base_catalog_galactic(
        n_sources=n_sources,
        max_age=0 * u.year,
    )
    table = add_observed_parameters(table)

    random_state = get_random_state(random_state)

    # Draw angular extension, then compute physical extension
    intrinsic_extension = random_state.normal(mean_intrinsic,sigma_intrinsic,n_sources)
    intrinsic_extension = intrinsic_extension * u.pc
    constant = 0.000291 / u.arcmin
    angular_extension = intrinsic_extension/(table['distance'] * constant)
    angular_extension = angular_extension.to('deg')/2.0

    print('size ', intrinsic_extension, angular_extension, table['distance'])

    # Draw log parabola spectral model parameters
    alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
    beta = random_state.uniform(0, max_index_beta, n_sources)

    # Define the luminosity

    logluminosity = random_state.normal(mean_logluminosity, sigma_logluminosity, n_sources)
    luminosity = (10 ** logluminosity) * u.erg / u.second

    distance = table['distance'].quantity

    # integral sed between 1 and 10 TeV
    energy_flux = luminosity / (4 * np.pi * distance ** 2)
    #print(energy_flux)
    energy_flux = energy_flux.to('TeV cm-2 s-1')

    vals = []
    vals_crab = []
    name = []
    for idx in range(len(table)):
        val, val_crab = flux_amplitude_from_energy_flux(
            alpha[idx], beta[idx], energy_flux[idx],
        )
        vals.append(val)
        vals_crab.append(val_crab)
        source_name = 'pwn_{}'.format(idx)
        name.append(source_name)
    norm = u.Quantity(vals)
    norm_crab = u.Quantity(vals_crab)



    # import IPython; IPython.embed()
    table = table[['x', 'y', 'z', 'spiralarm', 'GLON', 'GLAT', 'RA', 'DEC', 'distance']]
    table['sigma'] = Column(angular_extension, description='Angular extension (deg)', unit='deg')
    table['spec_alpha'] = Column(alpha, description='Spectral model parameter (log parabola)')
    table['spec_beta'] = Column(beta, description='Spectral model parameter (log parabola)')
    table['spec_norm'] = Column(norm, description='Spectral model norm parameter (log parabola)', unit='MeV-1 s-1 cm-2')
    table['spec_norm_crab'] = Column(norm_crab, description='Spectral model norm parameter (log parabola) in crab units')

    table['source_name'] = Column(name, description='source name')
    return table

if __name__ == '__main__':
    table = make_pwn_table()

    filename = 'ctadc_skymodel_gps_sources_pwn.ecsv'
    print('Writing {}'.format(filename))
    table.write(filename, format='ascii.ecsv', overwrite=True)
