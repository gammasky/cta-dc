from astropy.table import Table, Column
import matplotlib.pyplot as plt
import numpy as np

infile = "/Users/ltibaldo/GitHub/cta-dc/sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.ecsv"

snrmodel = Table.read(infile,format='ascii.ecsv')
names = ['Ia','IIP','Ib/c','IIb']

z = [np.array([]),
     np.array([]),
     np.array([]),
     np.array([])]

for obj in snrmodel:
    type = int(obj['type']) - 1
    z[type] = np.append(z[type],obj['galactocentric_z'])

N = len(snrmodel)
for s, zval in enumerate(z):
    print(s+1,float(len(zval))/N)

ax = plt.subplot()
ax.set_xlabel("z (kpc)")

for s, zval in enumerate(z):
    plt.hist(zval,histtype = 'step',label = names[s])

x = np.linspace(-1.5,1.5,100)
ax.plot(x,140*np.exp(-x**2/(2*(0.146/2.355)**2)),linestyle='--',label="Gauss FWHM 146 pc")
ax.plot(x,140*np.exp(-x**2/(2*(1./2.355)**2)),linestyle='--',label="Gauss FWHM 1 kpc")

ax.legend(loc='best')
plt.savefig("SNR_z.png")