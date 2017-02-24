import numpy as np
from astropy.table import Table, vstack
from gammapy.catalog import SourceCatalogGammaCat, SourceCatalogHGPS
import matplotlib.pyplot as plt
import os
from gammapy.spectrum import CrabSpectrum
from astropy import units as u
from astropy.table import Column
from scipy import optimize
from scipy.stats import norm

flux_min = 0.5
flux_max = 100
bins = np.logspace(np.log10(flux_min),np.log10(flux_max),30)

gammacat = SourceCatalogGammaCat()

cat_table = gammacat.table
#print(cat_table.info())

mask_pwn = cat_table['classes'] == 'pwn'
cat_pwn = cat_table[mask_pwn]
#cat_pwn.pprint()
mask_nan= np.isfinite(cat_pwn['spec_flux_above_1TeV'])
#np.isnan
cat_pwn= cat_pwn[mask_nan]

#print(cat_pwn['common_name','spec_flux_above_1TeV_crab'])



#hist, bin_edges = np.histogram(cat_pwn['spec_flux_above_1TeV_crab'],bins=bins)

#print(hist.sum())

#y = np.insert(hist[::-1].astype('float64').cumsum(),0, 0.01)
##reverse array bins[::-1]
#plt.step(bins[::-1],y)
#plt.loglog()
#plt.savefig('gammacat.png')

hgpscat = SourceCatalogHGPS(os.environ['HGPS'])
hgpscat_table = hgpscat.table
emin = 1
emax = 100
crab = CrabSpectrum('meyer').model
flux_above_1TeV_crab = crab.integral(emin * u.TeV, emax * u.TeV)

#print(hgpscat_table.info())
mask_pwn_ = hgpscat_table['Source_Name'] == 'HESS J1834-087'
hgpscat_ = hgpscat_table[mask_pwn_]
print(hgpscat_)

mask_pwn_hgps = hgpscat_table['Source_Class'] == 'PWN'
hgpscat_pwn = hgpscat_table[mask_pwn_hgps]
#hgpscat_pwn.pprint()
print(hgpscat_pwn['Source_Name','Flux_Spec_Int_1TeV'])

num_pwn = len(hgpscat_pwn)
print('num_pwn = ',num_pwn)


mask_composite_hgps = hgpscat_table['Source_Class'] == 'Composite'
hgpscat_composite = hgpscat_table[mask_composite_hgps]
#print(hgpscat_composite['Source_Name','Flux_Spec_Int_1TeV'])
idxx1 = np.where(hgpscat_composite['Source_Name'] == 'HESS J1714-385')[0]
print(len(idxx1),idxx1[0])
hgpscat_composite.remove_row(int(idxx1[0]))
print(hgpscat_composite['Source_Name','Flux_Spec_Int_1TeV'])
num_composite = len(hgpscat_composite)
print('num composite = ',num_composite)

mask_unid_hgps = hgpscat_table['Source_Class'] == 'Unid'
hgpscat_unid = hgpscat_table[mask_unid_hgps]
#hgpscat_unid.pprint()
num_unid = len(hgpscat_unid)
#print(hgpscat_unid['Source_Name','Flux_Spec_Int_1TeV'])

idx1 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1018-589B')[0]
hgpscat_unid.remove_row(int(idx1[0]))
idx2 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1023-575')[0]
hgpscat_unid.remove_row(int(idx2[0]))
idx3 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1457-593')[0]
hgpscat_unid.remove_row(int(idx3[0]))
idx4 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1503-582')[0]
hgpscat_unid.remove_row(int(idx4[0]))
idx5 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1614-518')[0]
hgpscat_unid.remove_row(int(idx5[0]))
idx6 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1634-472')[0]
hgpscat_unid.remove_row(int(idx6[0]))
idx7 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1646-458')[0]
hgpscat_unid.remove_row(int(idx7[0]))
idx8 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1641-463')[0]
hgpscat_unid.remove_row(int(idx8[0]))
idx9 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1729-345')[0]
hgpscat_unid.remove_row(int(idx9[0]))
idx10 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1747-248')[0]
hgpscat_unid.remove_row(int(idx10[0]))
idx11 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1800-240')[0]
hgpscat_unid.remove_row(int(idx11[0]))
idx12 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1808-204')[0]
hgpscat_unid.remove_row(int(idx12[0]))
idx13 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1832-085')[0]
hgpscat_unid.remove_row(int(idx13[0]))
idx14 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1832-093')[0]
hgpscat_unid.remove_row(int(idx14[0]))
idx15 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1843-033')[0]
hgpscat_unid.remove_row(int(idx15[0]))
idx16 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1848-018')[0]
hgpscat_unid.remove_row(int(idx16[0]))
idx17 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1923+141')[0]
hgpscat_unid.remove_row(int(idx17[0]))
idx18 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1943+213')[0]
hgpscat_unid.remove_row(int(idx18[0]))
idx19 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1702-420')[0]
hgpscat_unid.remove_row(int(idx19[0]))
idx20 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1713-381')[0]
hgpscat_unid.remove_row(int(idx20[0]))
idx21 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1507-622')[0]
hgpscat_unid.remove_row(int(idx21[0]))
idx22 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1844-030')[0]
hgpscat_unid.remove_row(int(idx22[0]))
idx23 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1746-285')[0]
hgpscat_unid.remove_row(int(idx23[0]))
idx24 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1745-303')[0]
hgpscat_unid.remove_row(int(idx24[0]))
idx25 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1741-302')[0]
hgpscat_unid.remove_row(int(idx25[0]))
idx26 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1745-290')[0]
hgpscat_unid.remove_row(int(idx26[0]))
idx27 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1746-308')[0]
hgpscat_unid.remove_row(int(idx27[0]))


print(hgpscat_unid['Source_Name','Flux_Spec_Int_1TeV'])
print('num_unid = ',len(hgpscat_unid))

hgpscat_pwn_extended = vstack([hgpscat_pwn, hgpscat_unid, hgpscat_composite])
#hgpscat_pwn_extended.pprint()
#hgpscat_pwn_extended.info()

#print(hgpscat_pwn_extended['Source_Name','Flux_Spec_Int_1TeV'])
num_pwn_extended = len(hgpscat_pwn_extended)
print('num_tot = ',num_pwn_extended)

col_crab = np.ndarray(shape=(num_pwn_extended, 1), dtype=float, order='F')

counter = 0
for row in hgpscat_pwn_extended['Flux_Spec_Int_1TeV']:
    flux_cu = ((hgpscat_pwn_extended[counter]['Flux_Spec_Int_1TeV'])/flux_above_1TeV_crab.value)*100
    #print(flux_cu)
    col_crab[counter]= flux_cu
    counter += 1

print('flux_crab_above_1TEV: ',flux_above_1TeV_crab.value)


new_column = Column(col_crab, name='spec_flux_above_1TeV_crab')
hgpscat_pwn_extended.add_column(new_column, index=10)
print(hgpscat_pwn_extended['Source_Name','spec_flux_above_1TeV_crab','Spatial_Model','Size'])


################ LogN-LogS ####################################################################

hist_hgps, bin_edges_hgps = np.histogram(hgpscat_pwn_extended['spec_flux_above_1TeV_crab'],bins=bins)


y_hgps = np.insert(hist_hgps[::-1].astype('float64').cumsum()+0.001,0, 0.01)
y_hgps_err = np.sqrt(y_hgps)
#print(' leng array: ',len(hist_hgps), hist_hgps.sum())

num_bins = len(y_hgps)

for i in range(1, num_bins):
    #print(bin_edges_hgps[num_bins - i], '  ',(y_hgps[i]))
    print(np.log10(bin_edges_hgps[num_bins - i]), '  ', np.log10(y_hgps[i]))#, ' ', bin_edges_hgps[num_bins - i], '  ',(y_hgps[i]))

logx = np.log10(bin_edges_hgps[::-1])
logy = np.log10(y_hgps)
logyerr = y_hgps_err / y_hgps


rangelogx = []
rangelogy = []
rangelogyerr = []

for i in range(1, 15):
    #print(i, logx[i-1],logy[i-1])

    rangelogx.append(logx[i-1])
    rangelogy.append(logy[i-1])
    rangelogyerr.append(logyerr[i-1])
    print(type(rangelogx[i-1]))

#print(len(rangelogx))
for i in range(1, len(rangelogx)+1):
    print(i, '  ', rangelogx[i-1], ' ', rangelogy[i-1])


fitfunc = lambda p, x: p[0] - p[1] * x
errfunc = lambda p, x, y: (y - fitfunc(p, x))

pinit = [2.2, 1.0]
out = optimize.leastsq(errfunc, pinit, args=(logx, logy))

#p0 = np.power(10,out[0][0])
#p1 = out[0][1]
#print(p0,' ',p1)

#

import powerlaw
fit = powerlaw.Fit(hgpscat_pwn_extended['spec_flux_above_1TeV_crab'].data, xmin=10, xmax=100)

print('fit :', fit.power_law.alpha)
N = (hgpscat_pwn_extended['spec_flux_above_1TeV_crab'].data > 10).sum()
print('N: ', N)


p0 = 2.28
p1 = 1.12

x0 = np.log10(53)
x1 = np.log10(10)
x2 = np.log10(0.2)
y1 = p0 -x1*p1
y0 = p0-x0*p1

y2= p0-x2*p1

#print(x0,x1,y0,y1)
#reverse array bins[::-1]
plt.figure()
plt.step(bins[::-1],y_hgps,color='r')
plt.loglog()
plt.ylim(flux_min,5000)
plt.ylabel('N (F>1 TeV')
plt.xlabel('F>1 TeV [cu]')


plt.plot([8, 8], [flux_min, 100], 'k--', lw=1, color='b')
plt.plot([0.2, 0.2], [flux_min, 5000], 'k--', lw=1, color='b')
plt.plot([0.2, 53], [np.power(10,y2), np.power(10,y0)], 'k--', lw=1, color='g')
plt.savefig('hgps.png')

##############################################
print('#########  SIZE DISTRIB ##################')
bins_size = np.linspace(0,1,num=10)
mask_nan_size= np.isfinite(hgpscat_pwn_extended['Size'])
#np.isnan
hgpscat_pwn_extended_size = hgpscat_pwn_extended[mask_nan_size]

print(hgpscat_pwn_extended_size['Spatial_Model','Size'])

print(type(hgpscat_pwn_extended_size['Size']))

hist_hgps_size, bin_edges_hgps_size = np.histogram(hgpscat_pwn_extended_size['Size'],bins=bins_size,normed=True)

bin_center = (bins_size[:-1] + bins_size[1:]) / 2


for i in range(1, (len(hist_hgps_size)+1)):
    print(i-1, '  ', hist_hgps_size[i-1], ' ', bin_center[i-1], bins_size[i-1])

plt.figure()

#print(type(bins_size),len(bin_s),len(hist_hgps_size))


#plt.plot(bin_center,hist_hgps_size, 'ro') #where='mid'
#plt.bar(bin_center, hist_hgps_size, align='center')
plt.step(bin_center, hist_hgps_size, where='post')
rv = norm()
mean, sigma = norm.fit(hgpscat_pwn_extended_size['Size'],loc=0.2)
x = np.linspace(0,1,100)
y = norm.pdf(x,mean, sigma)
plt.plot(x,y)
print(mean, sigma)
plt.savefig('size.png')

import IPython; IPython.embed()