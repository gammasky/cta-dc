import astropy.units as u
from astropy.table import Table, Column
from astropy.units import Quantity
from astropy.coordinates import SkyCoord, spherical_to_cartesian

from gammapy.astro.population.simulate import astrometry
from gammapy.astro.population import Exponential, FaucherSpiral, RMIN, RMAX, ZMIN, ZMAX, radial_distributions
from astropy.modeling.models import Gaussian2D
from gammapy.utils.random import get_random_state
from gammapy.spectrum import SpectrumObservation, SpectrumFit, models


def make_base_pwn_catalog(n_sources = 375, spiral_arms=True, rad_dis='YK04',
                          mean_extension=0.13, sigma_extension=0.2, intrinsic_extension = 50,
                          mean_index_alpha=1.8, sigma_index_alpha=0.27, max_index_beta):

    """
        Make a catalog of PWN

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



    ### simulate position. Returns x, y in cartesian coordinates
    if isinstance(rad_dis, str):
        rad_dis = radial_distributions[rad_dis]

    # Draw r from the radial distribution
    r = Quantity(draw(RMIN.value, RMAX.value, n_sources, pdf(rad_dis())), 'kpc')

    # Apply spiralarm modelling, if not something else should be added
    if spiral_arms:
        r, theta, spiralarm = FaucherSpiral()(r)

    #Compute cartesian coordinates
    x, y = astrometry.cartesian(r, theta)

    random_state = get_random_state(random_state)

    # Define spectral model. I want a logparabola.

    alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
    beta = random_state.uniform(0, max_index_beta, n_sources)


    #Define the luminosity

    logluminosity = random_state.normal(mean_logluminosity,sigma_logluminosity,n_sources)

    luminosity = power(10,logluminosity)

    luminosity = luminosity * u.erg / u.second


    #Define a 2DGauss morphology
    angular_extension = random_state.normal(mean_extension, sigma_extension, n_sources)
    #how do I define the unit?
    angular_extension = angular_extension[angular_extension < 0] = 0.005
    #### SIGMA NEEDS TO BE POSITIVE

    intrinsic_extension = intrinsic_extension * u.pc
    angular_extension_arcmin = angular_extension.to('arcmin')
    # like this it is not an array!
    constant = 0.000291 / u.arcmin
    distance = intrinsic_extension / (angular_extension_arcmin * constant)

    #integral sed between 1 and 10 TeV
    sed = luminosity / (4 * np.pi * (distance ** 2))
    sed =  sed.to('TeV cm**-2 s**-1')


    #ampl = 1. / (2 * np.pi * (sigma) ** 2)
    ampl = draw(flux_min, flux_max, n_sources, power_law,
                index=flux_index)

    spatial_model = Gaussian2D(ampl, x, y, x_stddev=sigma, y_stddev=sigma)






    table = Table()
    table['x'] = Column(x, unit='kpc', description='Galactocentric x coordinate')
    table['y'] = Column(y, unit='kpc', description='Galactocentric y coordinate')

