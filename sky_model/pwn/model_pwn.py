"""
Make synthetic PNW population for CTA 1DC
"""
import numpy as np
import astropy.units as u
from astropy.table import Table, Column, vstack
from gammapy.astro.population import make_base_catalog_galactic, add_observed_parameters
from gammapy.utils.random import get_random_state
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum


def compute_glat_glon_distance(table):
    from astropy.coordinates import SkyCoord
    from gammapy.utils import coordinates as astrometry
    x, y, z = table['x'].quantity, table['y'].quantity, table['z'].quantity
    distance, glon, glat = astrometry.galactic(x, y, z)
    phys_size = table['physical_size'].quantity
    coordinate = SkyCoord(glon, glat, unit='deg', frame='galactic').transform_to('icrs')
    ra, dec = coordinate.ra.deg, coordinate.dec.deg
    r = np.sqrt(table['x'] ** 2 + table['y'] ** 2)


    constant = 0.000291 / u.arcmin

    size=phys_size/((distance.to('pc'))*constant)
    size = size.to('deg')

    for idx in range(len(size)):
        if size[idx] < 0.01 * u.Unit('deg'):
            print('PROBLEM: ', size[idx])
            size[idx] = 0.01 * u.Unit('deg')




    # Add columns to table

    table['distance'] = Column(distance, unit='kpc',
                               description='Distance observer to source center')
    table['GLON'] = Column(glon, unit='deg', description='Galactic longitude')
    table['GLAT'] = Column(glat, unit='deg', description='Galactic latitude')
    table['RA'] = Column(ra, unit='deg')
    table['DEC'] = Column(dec, unit='deg')
    table['size'] = Column(size, unit='deg')
    table['galactocentric_r'] = Column(r, unit='kpc', description='Galactocentric radius in the xy plan')

    return table


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
    flux_at_1TeV = energy_flux / energy_flux_standard_candle * u.Unit('cm-2 s-1 TeV-1')
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


# def make_pwn_table(
#         n_sources=650, random_state=0,
#         mean_extension=0.13, sigma_extension=0.1, mean_intrinsic=20, sigma_intrinsic=5,  # intrinsic_extension=50,
#         mean_index_alpha=2.1, sigma_index_alpha=0.3, max_index_beta=0.5,
#         mean_logluminosity=33.5, sigma_logluminosity=1.0, ):
#     """Make a catalog of PWN.
#
#     to be defined
#     1. how many from LogN-LogS
#     2. spatial distribution: spiral arms
#     3. flux: Log10(luminosity) gaussian
#     4. morphlogy: Gaussian distribution with the mean and the sigma derived from
#     the the Gaus distribution of the size of the known PWNe
#     5. spectral: logparabola (how to choose the parameters)
#
#     Parameters
#     -----------
#     n_sources: int
#         number of sources defined from the logN-logS of the PWN population
#     rad_dis : callable
#     rad_dis : callable
#         Radial surface density distribution of sources.
#         YusifovKucuk2004 is radial distribution of the surface density of pulsars in the galaxy
#     mean_extension: float deg
#         mean of the Gaussian distribution of the angular extension of the known PWN
#     sigma_extension: float deg
#         sigma of the Gaussian distribution of the angular extension of the known PWN
#     intrinsic_extension: int in pc
#         this is needed to compute the distance
#     mean_index_alpha: float
#         mean of the Gaus distribution of the alpha index under the assumption of a logparabola spectral model
#     sigma_index_alpha: float
#         mean of the Gaus distribution of the alpha index under the assumption of a logparabola spectral model
#     max_index_beta: float
#         maximum value of the beta index under the assumption of a logparabola spectral model.
#         Assuming an uniform distribution for beta
#     """
#     random_state = get_random_state(random_state)
#
#     table = make_base_catalog_galactic(
#         n_sources=n_sources,
#         max_age=0 * u.year,
#         random_state=random_state,
#     )
#     table = add_observed_parameters(table)
#
#     table['distance'] = table['distance'].quantity.to('kpc')
#     table['distance'].unit = 'kpc'
#
#     # Draw angular extension, then compute physical extension
#     # intrinsic_extension = random_state.normal(mean_intrinsic,sigma_intrinsic,n_sources)
#     intrinsic_extension = random_state.uniform(3, 60, n_sources)
#
#     intrinsic_extension = intrinsic_extension * u.pc
#     constant = 0.000291 / u.arcmin
#     angular_extension = intrinsic_extension / (table['distance'].quantity.to('pc') * constant)
#     angular_extension = angular_extension.to('deg') / 2.0
#     for idx in range(len(table)):
#         if angular_extension[idx] < 0.01 * u.Unit('deg'):
#             # print('PROBLEM: ', angular_extension[idx])
#             angular_extension[idx] = 0.01 * u.Unit('deg')
#
#     for idx in range(len(table)):
#         if angular_extension[idx] < 0:
#             print('ROBLEM: ', angular_extension[idx], intrinsic_extension[idx])
#     # print('size ', intrinsic_extension, angular_extension, table['distance'])
#
#     # Draw log parabola spectral model parameters
#     alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
#     # random_state.uniform(1.5, 2.6, n_sources)
#     beta = random_state.uniform(0, max_index_beta, n_sources)
#
#     # Define the luminosity
#
#     logluminosity = random_state.normal(mean_logluminosity, sigma_logluminosity, n_sources)
#     for idx in range(len(table)):
#         if logluminosity[idx] > 35:
#             logluminosity[idx] = random_state.normal(mean_logluminosity, sigma_logluminosity, 1)
#     luminosity = (10 ** logluminosity) * u.erg / u.second
#
#     distance = table['distance'].quantity
#     # We compute intrinsic physical extension from the angular extension here,
#     # to make sure the two quantities are consistent
#     size_physical = distance * np.tan(angular_extension.to('rad').value)
#
#     # integral sed between 1 and 10 TeV
#     energy_flux = luminosity / (4 * np.pi * distance ** 2)
#     energy_flux = energy_flux.to('TeV cm-2 s-1')
#
#     vals = []
#     vals_cu = []
#     int_vals = []
#     int_vals_cu = []
#     name = []
#     for idx in range(len(table)):
#         val, val_cu, int_val, int_val_cu = flux_amplitude_from_energy_flux(
#             alpha[idx], beta[idx], energy_flux[idx],
#         )
#         vals.append(val)
#         vals_cu.append(val_cu)
#         int_vals.append(int_val)
#         int_vals_cu.append(int_val_cu)
#         source_name = 'pwn_{}'.format(idx)
#         name.append(source_name)
#     norm = u.Quantity(vals)
#     norm_cu = u.Quantity(vals_cu)
#     int_flux = u.Quantity(int_vals)
#     int_flux_cu = u.Quantity(int_vals_cu)
#
#     table = table[['x', 'y', 'z', 'spiralarm', 'GLON', 'GLAT', 'RA', 'DEC', 'distance']]
#     table['sigma'] = Column(angular_extension, description='Angular size', unit='deg')
#     table['size_physical'] = Column(size_physical.to('pc').value, description='Physical size', unit='pc')
#     table['spec_alpha'] = Column(alpha, description='Spectral model parameter (log parabola)')
#     table['spec_beta'] = Column(beta, description='Spectral model parameter (log parabola)')
#     table['spec_norm'] = Column(norm, description='Spectral model norm parameter (log parabola)', unit='MeV-1 s-1 cm-2')
#     table['spec_norm_cu'] = Column(norm_cu, description='Spectral model norm parameter (log parabola) in crab units')
#     table['int_flux_above_1TeV'] = Column(int_flux, description='Integral flux above 1 TeV ', unit='s-1 cm-2')
#     table['int_flux_above_1TeV_cu'] = Column(int_flux_cu, description='Integral flux above 1 TeV in crab units')
#     table['source_name'] = Column(name, description='source name')
#
#     return table


def polish_pwn_table(table):
    table.rename_column('x', 'galactocentric_x')
    table.rename_column('y', 'galactocentric_y')
    table.rename_column('z', 'galactocentric_z')
    table.rename_column('size', 'sigma')
    table['galactocentric_x'].format = '%.5f'
    table['galactocentric_y'].format = '%.5f'
    table['galactocentric_z'].format = '%.5f'
    table['GLON'].format = '%.5f'
    table['GLAT'].format = '%.5f'
    table['RA'].format = '%.5f'
    table['DEC'].format = '%.5f'
    table['distance'].format = '%.5f'
    table['sigma'].format = '%.5f'
    table['physical_size'].format = '%.5f'
    table['spec_alpha'].format = '%5g'
    table['spec_beta'].format = '%5g'
    table['spec_norm'].format = '%5g'
    table['spec_norm_cu'].format = '%5g'
    table['int_flux_above_1TeV'].format = '%5g'
    table['int_flux_above_1TeV_cu'].format = '%5g'
    table['galactocentric_r'].format = '%0.4f'

    return table

def make_composites(random_state, min_frac_radius=0.1, max_frac_radius=0.7, type=1):
    """
    Exchange some PWN to be composites,
    i.e. coupled to the SNRs
    """
    filename = '../snrs/ctadc_skymodel_gps_sources_snr_2.ecsv'
    table_snrs = Table.read(filename, format='ascii.ecsv')

   # print(table_snrs)
   # table_snrs.info()

    #Do we want to select only one type of SNR for the composite
    # mask_type = table_snrs['type'] == 1
    # table_snrs_type = table_snrs[mask_type]

    n_snrs = len(table_snrs)

    n_composites = int(n_snrs/3)
    frac_radius_composite = random_state.uniform(min_frac_radius, max_frac_radius, n_composites)


    x_composites = []
    y_composites = []
    z_composites = []
    size_composites = []
    size_snrs = []
    frac = []
    type = []
    for id in range(1, int(n_composites)):
        x_snr = table_snrs[id]['galactocentric_x']
        y_snr = table_snrs[id]['galactocentric_y']
        z_snr = table_snrs[id]['galactocentric_z']
        size_phys_snr = table_snrs[id]['size_physical']
        x_composites.append(x_snr)
        y_composites.append(y_snr)
        z_composites.append(z_snr)
        size_phys_composite = frac_radius_composite[id]*size_phys_snr
        size_composites.append(size_phys_composite)
        size_snrs.append(size_phys_snr)
        frac.append(frac_radius_composite[id])
        type.append('composite')

    table_composite = Table()
    table_composite['x'] = Column(x_composites, description='Galactocentric x coordinate', unit='kpc')
    table_composite['y'] = Column(y_composites, description='Galactocentric y coordinate', unit='kpc')
    table_composite['z'] = Column(z_composites, description='Galactocentric z coordinate', unit='kpc')
    table_composite['physical_size'] = Column(size_composites, description='physical size', unit='pc')
    table_composite['type'] = Column(type, description='type of PWN')
    #table_composite['size_snr'] = Column(size_snrs, description='physical size', unit='pc')
    #table_composite['frac'] = Column(frac, description='physical size', unit='pc')

    #print('------------------------------------------')
    #print(table_composite)

    return table_composite

def make_pwn_pos(random_state,
            n_sources, min_intrinsic_extension=3, max_intrinsic_extension=30):

    print('n_sources ', n_sources)

    table = make_base_catalog_galactic(
        n_sources=n_sources,
        max_age=0 * u.year,
        random_state=random_state,
    )

    physical_size = random_state.uniform(min_intrinsic_extension, max_intrinsic_extension, n_sources)
    physical_size = u.Quantity(physical_size, 'pc')
    type = []
    for iii in range(1, len(physical_size)):
        type.append('isolated')

    table.remove_column('morph_type')
    table.remove_column('age')
    table.remove_column('n_ISM')
    table.remove_column('v_abs')
    table.remove_column('vx')
    table.remove_column('vy')
    table.remove_column('vz')
    table.remove_column('x_birth')
    table.remove_column('y_birth')
    table.remove_column('z_birth')
    table['physical_size'] = Column(physical_size,
                                          description='physical size', unit='pc')
    table['type'] = Column(type,
                                    description='type of PWN', unit='pc')

    return table

def add_spectra(table, random_state,
                mean_index_alpha=2.1, sigma_index_alpha=0.3, max_index_beta=0.5,
                mean_logluminosity=33.5, sigma_logluminosity=1.0):

    alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
    beta = random_state.uniform(0, max_index_beta, n_sources)

    # Define the luminosity

    logluminosity = random_state.normal(mean_logluminosity, sigma_logluminosity, n_sources)
    for idx in range(len(table)):
        if logluminosity[idx] > 35:
            logluminosity[idx] = random_state.normal(mean_logluminosity, sigma_logluminosity, 1)
    luminosity = (10 ** logluminosity) * u.erg / u.second

    distance = table['distance'].to('cm')
    # integral sed between 1 and 10 TeV
    energy_flux = luminosity / (4 * np.pi * distance * distance)
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


    table['spec_alpha'] = Column(alpha, description='Spectral model parameter (log parabola)')
    table['spec_beta'] = Column(beta, description='Spectral model parameter (log parabola)')
    table['spec_norm'] = Column(norm, description='Spectral model norm parameter (log parabola)', unit='MeV-1 s-1 cm-2')
    table['spec_norm_cu'] = Column(norm_cu, description='Spectral model norm parameter (log parabola) in crab units')
    table['int_flux_above_1TeV'] = Column(int_flux, description='Integral flux above 1 TeV ', unit='s-1 cm-2')
    table['int_flux_above_1TeV_cu'] = Column(int_flux_cu, description='Integral flux above 1 TeV in crab units')

    return table

if __name__ == '__main__':
    # table = make_pwn_table()
    # table = polish_pwn_table(table)

    n_sources = 650
    random_seed=0
    random_state = get_random_state(random_seed)

    table_composites = make_composites(random_state=random_state)
    n_isolated_pwn = n_sources - len(table_composites)
    table_isolated = make_pwn_pos(n_sources=n_isolated_pwn, random_state=random_state)

    table = vstack([table_composites, table_isolated])
    compute_glat_glon_distance(table)

    table = add_spectra(table, random_state=random_state)

    polish_pwn_table(table)

    table.pprint()

    filename = 'ctadc_skymodel_gps_sources_pwn.ecsv'
    print('Writing {}'.format(filename))
    table.write(filename, format='ascii.ecsv', overwrite=True)
