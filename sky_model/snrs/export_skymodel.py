##############################################################################################
############# Jan 2017 - export simulated SNRs for Roberta and Christoph        ##############
##############################################################################################
## Export '(Results.txt)'  from .txt to .ecsv 



import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
# from matplotlib.ticker import ScalarFormatter
# import matplotlib.ticker as mticker
# import matplotlib.ticker as ticker
import astropy.units as u

######################################################################
# READ ORIGINAL .TXT FILES

data = np.genfromtxt('test_skymodel.txt', defaultfmt='%0.6f', names=True)
# Trick to read the Energy grid
data2 = np.genfromtxt('test_skymodel.txt', defaultfmt='%0.6f')
EGRID = data2[0]
EGRID = EGRID[~np.isnan(EGRID)]
size = len(data.dtype.names)
number_of_SNRs = len(data)  # number of SNRs

######################################################################
# BUILD ARRAYS BEFORE WRITING IN .ECSV
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

######################################################################
# WRITE IN .ECSV

from astropy.table import Table

t = Table(temp, names=NAMES)

t.meta[
    'hello'] = 'Parameters used for this simulation: spectral index alpha=4.1; Electron to proton ratio Kep=10-2; CR efficiency xi=0.1'
t.meta['energy_array'] = EGRID

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

print(t)
t.write('SNRs_SIMULATED_SKYMODEL_TEST.ecsv', format='ascii.ecsv')
