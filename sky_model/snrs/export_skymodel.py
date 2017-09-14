"""
Reformat data for simulated SNRs by Pierre Cristofari from TXT to ECSV
"""
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from gammapy.utils.coordinates import galactic as compute_galactic_coordinates
import numpy as np

def read_txt_files(version):
    filename = 'test_roberta_september_{}.txt'.format(version)
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
    for x in range(1, 10):
        NAMES = np.append(NAMES, data.dtype.names[x])
    for x in range(10, size):
        NAMES = np.append(NAMES, EGRID[x - 10])

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

    t[NAMES[9]].unit = 'cm-2 s-1'
    t[NAMES[9]].format = '%0.4e'
    t[NAMES[9]].description = 'Integral flux 1-10TeV'

    for x in range(10, size):
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
    table['flux_1_10'] = Column(table['flux_1_10'].data, unit='cm-2 s-1', description='Integral flux from 1 to 10 TeV',
                                format='%0.4e')

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

    table['glon'] = Column(glon, unit='deg', description='Galactic longitude')
    table['glon'].format = '%.5f'
    table['glat'] = Column(glat, unit='deg', description='Galactic latitude')
    table['glat'].format = '%.5f'
    zarray = np.zeros(len(table['distance']))
    table['skip'] = Column(zarray, description='Skip boolean, 1 skip 0 keep')
    table['skip'].format = '%.d'

    #print(table['flux_1_10'])
    #print(table.info())
    return table


if __name__ == '__main__':
    for version in [1]:
        table = read_txt_files(version=version)
        table = add_extra_info(table)

        filename = 'ctadc_skymodel_gps_sources_snr_{}.ecsv'.format(version)
        print('Writing {}'.format(filename))
        table.write(filename, format='ascii.ecsv', overwrite=True)
