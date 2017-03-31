"""
Reformat data for simulated SNRs by Pierre Cristofari from TXT to ECSV
"""
import numpy as np
import astropy.units as u
from astropy.table import Table


def read_txt_files():
    filename = 'ctadc_skymodel_gps_sources_snr_2.txt'
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
    t[NAMES[0]].description = 'Simulation_number'

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
        t[NAMES[x]].unit = 'TeV^{-1}cm^{-2}s^{-1}'
        t[NAMES[x]].description = 'Differential spectrum'

    return t


if __name__ == '__main__':
    table = read_txt_files()
    filename = 'ctadc_skymodel_gps_sources_snr_2.ecsv'
    print('Writing {}'.format(filename))
    table.write(filename, format='ascii.ecsv', overwrite=True)
