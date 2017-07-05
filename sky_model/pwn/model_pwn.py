"""
Make synthetic PNW population for CTA 1DC
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
from gammapy.astro.population import make_base_catalog_galactic
from gammapy.utils.random import get_random_state
from gammapy.utils import coordinates as astrometry
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum
from gammapy.utils.random import sample_powerlaw

def compute_glat_glon_distance(table):
    x, y, z = table['x'].quantity, table['y'].quantity, table['z'].quantity
    distance, glon, glat = astrometry.galactic(x, y, z)
    phys_size = table['size_physical'].quantity
    coordinate = SkyCoord(glon, glat, unit='deg', frame='galactic').transform_to('icrs')
    ra, dec = coordinate.ra.deg, coordinate.dec.deg
    r = np.sqrt(table['x'] ** 2 + table['y'] ** 2)

    constant = 0.000291 / u.arcmin

    size = phys_size / ((distance.to('pc')) * constant)
    # Clip size because we have a min / max for allowed size in XML / science tools
    size = np.clip(size.to('deg').value, 0.01, 100) * u.deg

    # Add columns to table

    table['distance'] = Column(distance, unit='kpc', description='Distance observer to source center')
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
    #print('crab_flux_above_1TeV ', crab_flux_above_1TeV )
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
    #print('std: ', energy_flux_standard_candle, 'energy_flux: ', energy_flux, flux_at_1TeV)
    flux_at_1TeV = flux_at_1TeV.to('cm-2 s-1 MeV-1')

    spec2 = LogParabola(amplitude=flux_at_1TeV, reference=pivot_energy, alpha=alpha, beta=beta)
    energy_flux_above_1TeV = spec2.energy_flux(emin=emin_int, emax=10 * u.TeV)
   # print('c',energy_flux_above_1TeV.to('TeV cm-2 s-1'),energy_flux_above_1TeV.to('TeV cm-2 s-1')/(flux_at_1TeV.to('cm-2 s-1 TeV-1')/)) )
    flux_above_1TeV = spec2.integral(emin=emin_int, emax=10 * u.TeV)
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
    #print(crab_flux_above_1TeV, flux_above_1TeV, flux_above_1TeV / crab_flux_above_1TeV, flux_above_1TeV_cu  )
    flux_above_1TeV_cu_model = (flux_above_1TeV / crab_flux_above_1TeV_model).to('%')

    return flux_at_1TeV, flux_at_1TeV_cu, flux_above_1TeV, flux_above_1TeV_cu


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
    table['size_physical'].format = '%.5f'
    table['spec_alpha'].format = '%5g'
    table['spec_beta'].format = '%5g'
    table['spec_norm'].format = '%5g'
    table['spec_norm_cu'].format = '%5g'
    table['int_flux_above_1TeV'].format = '%5g'
    table['int_flux_above_1TeV_cu'].format = '%5g'
    table['galactocentric_r'].format = '%0.4f'

    return table


def select_those_to_removed(table, tag='pwn'):
    #[20, 3, 3, 1, 6, 9]
    skip_array_pwn = [20,3,3,1,6,9]
    skip_array_composite = [20,0,0,0,0,0]

    if (tag == 'pwn'):
        skip_array = skip_array_pwn
    if (tag == 'composite'):
        skip_array = skip_array_composite

    name = []
    skip = []
    more_than_10cu = 0
    between_10_and_8cu = 0
    between_8_and_6cu = 0
    between_6_and_4cu = 0
    between_4_and_2cu = 0
    between_2_and_1cu = 0
    for row in table:
        if (tag == 'pwn'):
            source_name = 'pwn_{}'.format(row.index)
        if (tag == 'composite'):
            source_name = 'composite_{}'.format(row.index)
        name.append(source_name)
        if (row['int_flux_above_1TeV_cu'] > 10):
            print('crab: ', row['int_flux_above_1TeV_cu'], row['int_flux_above_1TeV'])
            more_than_10cu += 1
            if (more_than_10cu <= skip_array[0]):
                skip.append(1)
            else:
                skip.append(0)
        if (row['int_flux_above_1TeV_cu'] > 8 and row['int_flux_above_1TeV_cu'] < 10):
            between_10_and_8cu += 1
            #print('8-10    ', row.index, row['int_flux_above_1TeV_cu'])
            if (between_10_and_8cu <= skip_array[1]):
                skip.append(1)
            else:
                skip.append(0)
        if (row['int_flux_above_1TeV_cu'] < 8 and row['int_flux_above_1TeV_cu'] > 6):
            #print('4-6: ', remove_or_not_3, row.index, row['int_flux_above_1TeV_cu'])
            between_8_and_6cu += 1
            if (between_8_and_6cu <= skip_array[2]):
                skip.append(1)
            else:
                skip.append(0)
        if (row['int_flux_above_1TeV_cu'] < 6 and row['int_flux_above_1TeV_cu'] > 4):
            #print('4-6: ', remove_or_not_3, row.index, row['int_flux_above_1TeV_cu'])
            between_6_and_4cu += 1
            if (between_6_and_4cu <= skip_array[3]):
                skip.append(1)
            else:
                skip.append(0)
        if (row['int_flux_above_1TeV_cu'] < 4 and row['int_flux_above_1TeV_cu'] > 2):
            #print('2-4: ', remove_or_not_4, row.index, row['int_flux_above_1TeV_cu'])
            between_4_and_2cu += 1
            if (between_4_and_2cu <= skip_array[4]):
                skip.append(1)
            else:
                skip.append(0)
        if (row['int_flux_above_1TeV_cu'] < 2 and row['int_flux_above_1TeV_cu'] > 1):
            #print('1-2: ', remove_or_not_5, row.index, row['int_flux_above_1TeV_cu'])
            between_2_and_1cu += 1
            if (between_2_and_1cu <= skip_array[5]):
                skip.append(1)
            else:
                skip.append(0)
        elif (row['int_flux_above_1TeV_cu'] < 1):
            skip.append(0)


    print('more than 10cu: ', more_than_10cu)
    print('between 10 and 8 cu: ', between_10_and_8cu)
    print('between 8 and 6 cu: ', between_8_and_6cu)
    print('between 6 and 4 cu: ', between_6_and_4cu)
    print('between 4 and 2 cu: ', between_10_and_8cu)
    print('between 2 and 1 cu: ', between_2_and_1cu)
    table['skip'] = Column(skip, description='skip boolean: 0 take, 1 skip')
    table['source_name'] = Column(name, description='source name')
    return table


def make_composites(random_state, min_frac_radius=0.1, max_frac_radius=0.7, type=1):
    """
    Exchange some PWN to be composites,
    i.e. coupled to the SNRs
    """
    filename = '../snrs/ctadc_skymodel_gps_sources_snr_2.ecsv'
    table_snrs = Table.read(filename, format='ascii.ecsv')

    # Do we want to select only one type of SNR for the composite
    # mask_type = table_snrs['type'] == 1
    # table_snrs_type = table_snrs[mask_type]

    n_snrs = len(table_snrs)

    n_composites = int(n_snrs / 3)
    frac_radius_composite = random_state.uniform(min_frac_radius, max_frac_radius, n_composites)
    #offset = random_state.norm(0.1, 0.1, n_composites )

    x_composites = []
    y_composites = []
    z_composites = []
    size_composites = []
    size_snrs = []
    frac = []
    type = []
    for id in range(1, int(n_composites)):
        x_snr = table_snrs[id+500]['galactocentric_x']
        y_snr = table_snrs[id+500]['galactocentric_y']
        z_snr = table_snrs[id+500]['galactocentric_z']
        size_phys_snr = table_snrs[id]['size_physical']
        x_composites.append(x_snr)
        y_composites.append(y_snr)
        z_composites.append(z_snr)
        size_phys_composite = frac_radius_composite[id] * size_phys_snr
        size_composites.append(size_phys_composite)
        size_snrs.append(size_phys_snr)
        frac.append(frac_radius_composite[id])
        type.append('composite')

    table_composite = Table()
    table_composite['x'] = Column(x_composites, description='Galactocentric x coordinate', unit='kpc')
    table_composite['y'] = Column(y_composites, description='Galactocentric y coordinate', unit='kpc')
    table_composite['z'] = Column(z_composites, description='Galactocentric z coordinate', unit='kpc')

    table_composite['size_physical'] = Column(size_composites, description='physical size', unit='pc')
    #table_composite['size_snr'] = Column(size_snrs, description='physical size', unit='pc')
    #table_composite['frac'] = Column(frac, description='physical size', unit='pc')

    #print('------------------------------------------')
    #print(table_composite)
    table_composite['type'] = Column(type, description='type of PWN')

    return table_composite


def make_pwn_pos(random_state,
                 n_sources, min_intrinsic_extension=3, max_intrinsic_extension=50,
                 mean_intrinsic_extension=15, sigma_intrinsic_extension=10, ):
   # print('n_sources ', n_sources)

    table = make_base_catalog_galactic(
        n_sources=n_sources,
        max_age=0 * u.year,
        random_state=random_state,
    )


    size_physical = random_state.normal(mean_intrinsic_extension,sigma_intrinsic_extension, n_sources)


    for iii in range(0, len(size_physical)):
        if (size_physical[iii] < 0):
            print(size_physical[iii])
            #size_physical[iii] = random_state.normal(mean_intrinsic_extension, sigma_intrinsic_extension, n_sources)
            #if (size_physical[iii] < 0):
            size_physical[iii] = 0.1


    #physical_size = random_state.uniform(min_intrinsic_extension, max_intrinsic_extension, n_sources)
    size_physical = u.Quantity(size_physical, 'pc')
    type = []
    for iii in range(0, len(size_physical)):
        type.append('isolated')

    #size_physical = random_state.uniform(min_intrinsic_extension, max_intrinsic_extension, n_sources)
    #size_physical = u.Quantity(size_physical, 'pc')


    table.remove_column('age')
    table.remove_column('n_ISM')
    table.remove_column('v_abs')
    table.remove_column('vx')
    table.remove_column('vy')
    table.remove_column('vz')
    table.remove_column('x_birth')
    table.remove_column('y_birth')
    table.remove_column('z_birth')

    table['type'] = Column(type, description='type of PWN')

    table['size_physical'] = Column(size_physical, description='Physical size', unit='pc')


    return table


def add_spectra(table, random_state,
                mean_index_alpha=2.1, sigma_index_alpha=0.2, max_index_beta=0.4,
                mean_logluminosity=33.5, sigma_logluminosity=1.0):

    n_sources = len(table)
    alpha = random_state.normal(mean_index_alpha, sigma_index_alpha, n_sources)
    beta = random_state.uniform(0.1, max_index_beta, n_sources)

    #Define the luminosity
    luminosity = sample_powerlaw(
        x_min=2.5e33,
        x_max=1e37,
        gamma=1.3,
        size=n_sources,
    )
    luminosity = luminosity * u.erg / u.second

   # logluminosity = random_state.normal(mean_logluminosity, sigma_logluminosity, n_sources)
   # for idx in range(len(table)):
   #     if logluminosity[idx] > 35:
   #         logluminosity[idx] = random_state.normal(mean_logluminosity, sigma_logluminosity, 1)
   # luminosity = (10 ** logluminosity) * u.erg / u.second

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
        energy_flux_check = (LogParabola(amplitude=val,
                                         alpha=alpha[idx],
                                         beta=beta[idx],
                                         reference=1 * u.TeV,
                                         ).energy_flux(emin=1 * u.TeV, emax=10 * u.TeV)).to('TeV cm-2 s-1')

        #lum = energy_flux_check
        #print(energy_flux[idx], energy_flux_check)
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
    table['luminosity']= Column(luminosity, description='Intrinsic source luminosity', unit='erg s-1')

    # luminosity_check = []
    # for row in table:
    #      energy_flux_check = (LogParabola(amplitude=row['spec_norm'] * u.Unit('MeV-1 s-1 cm-2'),
    #                                       alpha=row['spec_alpha'],
    #                                       beta=row['spec_beta'],
    #                                       reference=1 * u.TeV,
    #                          ).energy_flux(emin=1 * u.TeV, emax=10 * u.TeV)).to('erg cm-2 s-1')
    #      print(energy_flux_check)
    #      dist = u.Quantity(row['distance'], 'kpc')
    #      lum_check = energy_flux_check * 4 * np.pi * (dist.to('cm')) ** 2
    #      luminosity_check.append(lum_check.value)
    #
    # table['luminosity_check'] = Column(luminosity_check, description='Intrinsic source luminosity check', unit='erg s-1')
    #
    # print(table)

    return table


def main():
    n_sources = 850
    random_seed = 0
    random_state = get_random_state(random_seed)

    table_composites = make_composites(random_state=random_state)
    #table_composites.pprint()

    n_isolated_pwn = n_sources - len(table_composites)
    table_isolated = make_pwn_pos(n_sources=n_isolated_pwn, random_state=random_state)
    #table_isolated.pprint()


    #table = vstack([table_composites, table_isolated])

    compute_glat_glon_distance(table_composites)

    table_composites = add_spectra(table_composites, random_state=random_state)
    polish_pwn_table(table_composites)

    select_those_to_removed(table_composites, tag='composite')

    filename_composite = 'ctadc_skymodel_gps_sources_composite.ecsv'
    print('Writing {}'.format(filename_composite))
    table_composites.write(filename_composite, format='ascii.ecsv', overwrite=True)


    compute_glat_glon_distance(table_isolated)
    table_isolated = add_spectra(table_isolated, random_state=random_state)
    polish_pwn_table(table_isolated)

    for row in table_isolated:
        if (row['int_flux_above_1TeV_cu']>10):
            print(row['int_flux_above_1TeV'], row['int_flux_above_1TeV_cu'])
    select_those_to_removed(table_isolated, tag='pwn')

    filename = 'ctadc_skymodel_gps_sources_pwn.ecsv'
    print('Writing {}'.format(filename))
    table_isolated.write(filename, format='ascii.ecsv', overwrite=True)


if __name__ == '__main__':
    main()
