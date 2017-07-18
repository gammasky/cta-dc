"""
Reformat data for simulated SNRs by Pierre Cristofari from TXT to ECSV
"""
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from gammapy.utils.coordinates import galactic as compute_galactic_coordinates
from astropy.units import Quantity
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum

def define_flux_crab_above_energy(emin=1 * u.TeV, emax=10 * u.TeV):
    crab = CrabSpectrum('meyer').model
    crabMAGIC = LogParabola(amplitude=3.23e-11 * u.Unit('cm-2 s-1 TeV-1'), reference=1 * u.TeV, alpha=2.47, beta=0.24)
    crab_flux_above_1TeV = crabMAGIC.integral(emin=emin, emax=emax)
    crab_flux_above_1TeV_model = crab.integral(emin=emin, emax=emax)

    return crab_flux_above_1TeV, crab_flux_above_1TeV_model

def read_txt_files(version):
    filename = 'ctadc_skymodel_gps_sources_snr_{}.txt'.format(version)
    print('Reading {}'.format(filename))
    data = np.genfromtxt(filename, defaultfmt='%0.6f', names=True)

    # Trick to read the Energy grid
    data2 = np.genfromtxt(filename, defaultfmt='%0.6f')

    EGRID = data2[0]
    EGRID = EGRID[~np.isnan(EGRID)]
    size = len(data.dtype.names)
    number_of_SNRs = len(data)  # number of SNRs

    # Reformat spectral data into arrays
    test = [None] * number_of_SNRs
    for j in range(0, number_of_SNRs):
        test[j] = [data[j][0]]
        for i in range(1, size):
            test[j] = np.append(test[j], data[j][i])

    temp = test[0]
    for j in range(1, number_of_SNRs):
        temp = np.vstack([temp, test[j]])

    NAMES = data.dtype.names[0]
    for x in range(1, 9):
        NAMES = np.append(NAMES, data.dtype.names[x])
    for x in range(9, size):
        NAMES = np.append(NAMES, EGRID[x - 9])

    t = Table(temp, names=NAMES)

    # Parameters used for this simulation
    t.meta['simulation_parameters'] = [
        dict(name='alpha', value=4.1, description='Spectral index'),
        dict(name='Kep', value=1e-2, description='Electron to proton ratio'),
        dict(name='xi', value=0.1, description='CR efficiency'),
    ]

    t.meta['energy_array'] = EGRID.tolist()

    t[NAMES[0]].unit = ''
    t[NAMES[0]].format = '%i'
    t[NAMES[0]].description = 'Simulation number'

    t[NAMES[1]].unit = 'kpc'
    t[NAMES[1]].format = '%0.4f'
    t[NAMES[1]].description = 'X Position'

    t[NAMES[2]].unit = 'kpc'
    t[NAMES[2]].format = '%0.4f'
    t[NAMES[2]].description = 'Y Position'

    t[NAMES[3]].unit = 'kpc'
    t[NAMES[3]].format = '%0.4f'
    t[NAMES[3]].description = 'Z Position'

    t[NAMES[4]].unit = u.Unit('cm-3')
    t[NAMES[4]].format = '%0.4f'
    t[NAMES[4]].description = 'density'

    t[NAMES[5]].unit = ''
    t[NAMES[5]].format = '%i'
    t[NAMES[5]].description = 'type of progenitor'

    t[NAMES[6]].unit = 'kyear'
    t[NAMES[6]].format = '%0.4f'
    t[NAMES[6]].description = 'age of the SNR'

    t[NAMES[7]].unit = 'arcmin'
    t[NAMES[7]].format = '%0.4f'
    t[NAMES[7]].description = 'apparent size of the SNR'

    t[NAMES[8]].unit = 'pc'
    t[NAMES[8]].format = '%0.4f'
    t[NAMES[8]].description = 'SNR Radius'

    for x in range(9, size):
        t[NAMES[x]].unit = 'cm-2 s-1 TeV-1'
        t[NAMES[x]].description = 'Differential spectrum'

    return t


def add_extra_info(table):
    # Change Pierre's XYZ to the one used in Gammapy at the moment
    # This was checked to be correct in https://github.com/gammasky/cta-dc/issues/17
    table['galactocentric_x'] = Column(table['POS_Y'].data, unit='kpc', description='Galactocentric X', format='%0.5f')
    table['galactocentric_y'] = Column(-table['POS_X'].data, unit='kpc', description='Galactocentric Y', format='%0.5f')
    table['galactocentric_z'] = Column(table['POS_Z'].data, unit='kpc', description='Galactocentric Y', format='%0.5f')
    table.remove_columns(['POS_X', 'POS_Y', 'POS_Z'])

    table.rename_column('Radius', 'size_physical')
    table.rename_column('size', 'sigma')

    r = np.sqrt(table['galactocentric_x'] ** 2 + table['galactocentric_y'] ** 2)
    table['galactocentric_r'] = Column(r, unit='kpc', description='Galactocentric radius in the xy plan')

    distance, glon, glat = compute_galactic_coordinates(
        x=table['galactocentric_x'].quantity,
        y=table['galactocentric_y'].quantity,
        z=table['galactocentric_z'].quantity,
    )

    table['distance'] = Column(distance, unit='kpc', description='Distance from Earth')
    table['distance'].format = '%.5f'

    table['GLON'] = Column(glon, unit='deg', description='Galactic longitude')
    table['GLON'].format = '%.5f'
    table['GLAT'] = Column(glat, unit='deg', description='Galactic latitude')
    table['GLAT'].format = '%.5f'
    table['skip'] = Column(0, description='Skip boolean, 1 skip 0 keep')
    return table


def add_sed_columns(table):
    energy_array = np.array(table.meta['energy_array'])
    sed_energy = np.tile(energy_array, reps=(len(table), 1))

    #table.info()

    # Copy over fluxes into array column
    sed_dnde = np.empty_like(sed_energy)
    for col_idx in range(50):
        sed_dnde[:, col_idx] = table.columns[6 + col_idx]

    table['sed_energy'] = u.Quantity(sed_energy, 'TeV').to('MeV')
    table['sed_dnde'] = u.Quantity(sed_dnde, 'cm-2 s-1 TeV-1').to('cm-2 s-1 MeV-1')

    return table

def make_spectral_point_selection(row):
    # Jurgen requested that we remove nodes with zero or very low flux
    # so that it works for ctools.
    # So here we remove the ones below an arbirtrary low threshold

    # In addition we noticed that some SNRs have all fluxes very low
    # We remove these super faint SNRs completely.

    mask = row['sed_dnde'] > 1e-20
    sed_energy = row['sed_energy'][mask]
    sed_dnde = row['sed_dnde'][mask]

    sum = 0
    for idx in range(len(sed_dnde) - 1):
        if (sed_energy[idx]> 1000000 and sed_energy[idx]<10000000):
            # print (idx, sed_energy[idx])
            bin = sed_energy[idx + 1] - sed_energy[idx]
            dnde = bin * sed_dnde[idx]
            sum = + dnde

    #print(sum)
    keep = (mask.sum() > 3)

    return dict(
        sed_energy=sed_energy,
        sed_dnde=sed_dnde,
        keep=keep,
        sum=sum,
    )

def select_those_to_keep(table):
    keep = []
    flux_1_10 = []
    flux_1_10_cu = []
    idx = 0
    crab_flux_above_1TeV, crab_flux_above_1TeV_model = define_flux_crab_above_energy();
    crab_flux = crab_flux_above_1TeV.value
    idx = 0

    print('---------------------------------------')
    for row in table:

        spec = make_spectral_point_selection(row)
        keep.append(spec['keep'])
        flux_1_10.append(spec['sum'])
        flux_1_10_cu.append(spec['sum']/crab_flux*100)
        if (flux_1_10_cu[idx]>1):
            print(flux_1_10_cu[idx])
        idx += 1
        if not spec['keep']:
            continue

    table['ctools_compatible'] = Column(keep, description='boolean indicating ctools compatible')
    table['int_flux_above_1TeV'] = Column(flux_1_10, description='integral flux between 1 and 10 TeV', unit= 'cm-2 s-1')
    table['int_flux_above_1TeV_cu'] = Column(flux_1_10_cu, description='integral flux between 1 and 10 TeV', unit= 'cm-2 s-1')
    return table

if __name__ == '__main__':
    for version in [1, 2]:
        table = read_txt_files(version=version)
        #print(table.info())
        table = add_sed_columns(table)
        table = select_those_to_keep(table)

        print(table)
        table = add_extra_info(table)


        #filename = 'ctadc_skymodel_gps_sources_snr_{}_robi.ecsv'.format(version)
        #print('Writing {}'.format(filename))
        #table.write(filename, format='ascii.ecsv', overwrite=True)
