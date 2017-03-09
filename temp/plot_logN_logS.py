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


def define_flux_crab_above_energy(emin=1* u.TeV, emax=100* u.TeV):
    crab = CrabSpectrum('meyer').model
    flux_crab = crab.integral(emin * u.TeV, emax * u.TeV)
    return flux_crab


def prepare_the_pwn_table(flux_min=0.5, flux_max=100):

    gammacat = SourceCatalogGammaCat()
    cat_table = gammacat.table
    mask_pwn = cat_table['classes'] == 'pwn'
    cat_pwn = cat_table[mask_pwn]
    cat_pwn.pprint()
    mask_nan= np.isfinite(cat_pwn['spec_flux_above_1TeV'])
    cat_pwn= cat_pwn[mask_nan]



    hgpscat = SourceCatalogHGPS(os.environ['HGPS'])
    hgpscat_table = hgpscat.table
    mask_pwn_hgps = hgpscat_table['Source_Class'] == 'PWN'
    hgpscat_pwn = hgpscat_table[mask_pwn_hgps]
    #hgpscat_pwn.pprint()
    print('----------------- PWN -----------------------')
    print(hgpscat_pwn['Source_Name','Flux_Spec_Int_1TeV'])
    num_pwn = len(hgpscat_pwn)

    print('num_pwn = ', num_pwn)



    mask_composite_hgps = hgpscat_table['Source_Class'] == 'Composite'
    hgpscat_composite = hgpscat_table[mask_composite_hgps]
    # print(hgpscat_composite['Source_Name','Flux_Spec_Int_1TeV'])
    idxx1 = np.where(hgpscat_composite['Source_Name'] == 'HESS J1714-385')[0]
    #print(len(idxx1), idxx1[0])
    hgpscat_composite.remove_row(int(idxx1[0]))
    print('----------------- composite -----------------------')
    print(hgpscat_composite['Source_Name', 'Flux_Spec_Int_1TeV'])
    num_composite = len(hgpscat_composite)

    print('num composite = ', num_composite)

    mask_unid_hgps = hgpscat_table['Source_Class'] == 'Unid'
    hgpscat_unid = hgpscat_table[mask_unid_hgps]
    # hgpscat_unid.pprint()
    #print(hgpscat_unid['Source_Name','Flux_Spec_Int_1TeV'])

    idx1 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1018-589B')[0]
    hgpscat_unid.remove_row(int(idx1[0]))
    idx2 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1708-410')[0]
    hgpscat_unid.remove_row(int(idx2[0]))
    idx2 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1852-000')[0]
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
    idx27 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1841-055')[0]
    hgpscat_unid.remove_row(int(idx27[0]))
    idx28 = np.where(hgpscat_unid['Source_Name'] == 'HESS J1626-490')[0]
    hgpscat_unid.remove_row(int(idx28[0]))


    #print(hgpscat_unid['Source_Name','Flux_Spec_Int_1TeV'])
    num_unid = len(hgpscat_unid)
    print('num_unid = ',num_unid)

    hgpscat_pwn_extended = vstack([hgpscat_pwn, hgpscat_unid, hgpscat_composite])
    print('----------------- extended -----------------------')
    print(hgpscat_pwn_extended['Source_Name','Flux_Spec_Int_1TeV'])
    num_pwn_extended = len(hgpscat_pwn_extended)
    print('num_tot = ',num_pwn_extended)

    flux_crab_above_1TeV = define_flux_crab_above_energy(emin=1,emax=100)

    flux_crab = []
    sigma = []
    for id in range(len(hgpscat_pwn_extended)):
        #print((hgpscat_pwn_extended[id]['Source_Name']),(hgpscat_pwn_extended[id]['Flux_Spec_Int_1TeV']))

        flux_cu = ((hgpscat_pwn_extended[id]['Flux_Spec_Int_1TeV']) / flux_crab_above_1TeV.value) * 100
        flux_crab.append(flux_cu)
        sigma.append(hgpscat_pwn_extended[id]['Size']/2.0)

    flux_crab_perc = u.Quantity(flux_crab)
    sigma_d = u.Quantity(sigma)
    hgpscat_pwn_extended['flux_crab'] = Column(flux_crab_perc, description='Flux above 1 TeV in crab units')
    hgpscat_pwn_extended['sigma'] = Column(sigma_d, description='radius of angular extension')
    #print(hgpscat_pwn_extended['flux_crab'])

    return hgpscat_pwn_extended


def plot_logN_logS(table_pwn, new_table_sim_pwn, merged_table, flux_min, flux_max):

    print('--------------------------------------- table ready now plotting ----------------------------')
    #print(table_pwn['Source_Name', 'Flux_Spec_Int_1TeV'])

    num_bins = 100
    bins = np.logspace(np.log10(flux_min), np.log10(flux_max), num_bins)

    hist, bin_edges = np.histogram(table_pwn['flux_crab'], bins=bins)
    #  print(table_sim_pwn['spec_norm_crab'])
    hist_sim, bin_edges_sim = np.histogram(new_table_sim_pwn['spec_norm_crab'], bins=bins)

    hist_merged, bin_edges_sim_merged = np.histogram(merged_table['spec_norm_crab'], bins=bins)

    y = np.insert(hist[::-1].astype('float64').cumsum(), 0, 0.01)
    y_err = np.sqrt(y)
    print(' lenght hgps array: ', len(hist), hist.sum())

    #print(y)

    y_sim = np.insert(hist_sim[::-1].astype('float64').cumsum(), 0, 0.01)
    y_err_sim = np.sqrt(y)
    print(' lenght sim array: ', len(hist_sim), hist_sim.sum())

    #print(y_sim)

    y_merged = np.insert(hist_merged[::-1].astype('float64').cumsum(), 0, 0.01)
    y_err_merged = np.sqrt(y)
    print(' lenght merged array: ', len(hist_merged), hist_merged.sum())

    #for i in range(1, num_bins):
    #        print(bin_edges_sim[num_bins-i],'-',bin_edges_sim[num_bins-i-1],y[i],y_sim[i], y_merged[i])



 #   for i in range(1, num_bins):
  #      # print(bin_edges_hgps[num_bins - i], '  ',(y_hgps[i]))
  #      print(np.log10(bin_edges_hgps[num_bins - i]), '  ',
  #            np.log10(y_hgps[i]))  # , ' ', bin_edges_hgps[num_bins - i], '  ',(y_hgps[i]))

   # logx = np.log10(bin_edges_hgps[::-1])
   ## logy = np.log10(y_hgps)
    #logyerr = y_hgps_err / y_hgps


    p0 = 2.28
    p1 = 1.20
    #p1 = 1.24

    x0 = np.log10(53)
    x1 = np.log10(10)
    x2 = np.log10(0.2)
    y1 = p0 - x1 * p1
    y0 = p0 - x0 * p1

    y2 = p0 - x2 * p1

    #plotting
    plt.figure()
    plt.step(bins[::-1], y, color='r',lw=2)
    #plt.step(bins[::-1], y_sim, color='black',lw=1)
    plt.step(bins[::-1], y_merged, color='b',lw=2)
    plt.loglog()
    plt.ylim(0.8, 5000)
    plt.xlim(0.05, 100)
    plt.ylabel('N (F>1 TeV)')
    plt.xlabel('F>1 TeV [cu]')

    plt.plot([10, 10], [flux_min, 100], 'k-', lw=1, color='black')
    plt.plot([0.2, 0.2], [flux_min, 5000], 'k--', lw=1, color='black')
    plt.plot([0.2, 53], [np.power(10, y2), np.power(10, y0)], 'k--', lw=2, color='b')
    plt.savefig('logNlogS.png')



def remove_bright_pwn(table_sim_pwn):

    remove_or_not, remove_or_not_2,remove_or_not_3, remove_or_not_4, remove_or_not_5 = 0, 0, 0, 0, 0

    lenght_table = len(table_sim_pwn)
    #print('len table: ', len(table_sim_pwn))
    remove_idx = []
    for row in range(1, lenght_table):
        if (table_sim_pwn[row]['spec_norm_crab']>10):
            print('crab: ',row,table_sim_pwn[row]['spec_norm_crab'],table_sim_pwn[row]['source_name'])
            remove_idx.append(row)
            continue


        if (table_sim_pwn[row]['spec_norm_crab']>8):
            if (remove_or_not < 8):
                print('8-10: ', row, table_sim_pwn[row]['spec_norm_crab'], table_sim_pwn[row]['source_name'])
                remove_idx.append(row)
                remove_or_not += 1

        if (table_sim_pwn[row]['spec_norm_crab'] > 6 and table_sim_pwn[row]['spec_norm_crab'] < 8):
            if (remove_or_not_3 < 5):
                print('6-8: ', row, table_sim_pwn[row]['spec_norm_crab'],table_sim_pwn[row]['source_name'])
                remove_idx.append(row)
                remove_or_not_3 += 1

        if (table_sim_pwn[row]['spec_norm_crab'] > 4 and table_sim_pwn[row]['spec_norm_crab'] < 6):
            if (remove_or_not_4 < 5):
                print('4-6: ', row, table_sim_pwn[row]['spec_norm_crab'],table_sim_pwn[row]['source_name'])
                remove_idx.append(row)
                remove_or_not_4 += 1

        if (table_sim_pwn[row]['spec_norm_crab'] > 2 and table_sim_pwn[row]['spec_norm_crab'] < 4):
            if (remove_or_not_5 < 4):
                print('2-4: ', row, table_sim_pwn[row]['spec_norm_crab'],table_sim_pwn[row]['source_name'])
                remove_idx.append(row)
                remove_or_not_5 += 1

       # if (table_sim_pwn[row]['spec_norm_crab'] > 1 and table_sim_pwn[row]['spec_norm_crab'] < 2):
       #     if (remove_or_not_2 < 1):
       #         print('1-2: ', row, table_sim_pwn[row]['spec_norm_crab'],table_sim_pwn[row]['source_name'])
       #         remove_idx.append(row)
       #         remove_or_not_2 += 1


    print('quanti: ',len(remove_idx))
    print(remove_idx)
    for idx in range(0, len(remove_idx)):
        source_name = 'pwn_{}'.format(remove_idx[idx])
        id = np.where(table_sim_pwn['source_name']== source_name)[0]
        #print(remove_idx[idx],source_name, id)
        table_sim_pwn.remove_row(int(id[0]))


    table_sim_pwn_reduced = table_sim_pwn['source_name','spec_norm_crab','sigma']

    #hgpscat_unid.remove_row(int(idx7[0]))

    #print(table_sim_pwn['spec_norm_crab'])

    return table_sim_pwn, table_sim_pwn_reduced


def plot_size_distrib(table_pwn, table_sim, merged_table):

    bins_size = np.linspace(0, 1.0, num=100)
    bin_center = (bins_size[:-1] + bins_size[1:]) / 2

    print('here')
    table_pwn.pprint()

    #for row in range(1, len(table_pwn) ):
     #   if (table_pwn[row]['sigma'] > 0.1):
     #       print('sigma:', row, table_pwn[row]['sigma'], table_pwn[row]['source_name'])


    mask_nan_hgps = np.isfinite(table_pwn['sigma'])
    # np.isnan
    table_pwn_size = table_pwn[mask_nan_hgps]

    table_pwn_size.pprint()
    mask_nan_merged = np.isfinite(merged_table['sigma'])
    # np.isnan
    merged_table_size = merged_table[mask_nan_merged]


    hist_size_hgps, bin_edges_hgps_size = np.histogram(table_pwn_size['sigma'], bins=bins_size, normed=True)
    hist_size_sim, bin_edges_sim_size = np.histogram(table_sim['sigma'], bins=bins_size, normed=True)
    hist_size_merged, bin_edges_merged_size = np.histogram(merged_table_size['sigma'], bins=bins_size, normed=True)

    print(hist_size_hgps)

    #for i in range(0, len(hist_size_sim) ):
    #    print(i , ' tt ', hist_size_sim[i], ' ', bin_center[i], bins_size[i])

    plt.figure()


    # plt.plot(bin_center,hist_hgps_size, 'ro') #where='mid'
    # plt.bar(bin_center, hist_hgps_size, align='center')
    plt.step(bin_center, hist_size_hgps, where='post')
    plt.step(bin_center, hist_size_sim, where='post', color='b')
    plt.step(bin_center, hist_size_merged, where='post', color='r')
    rv = norm()
    mean, sigma = norm.fit(table_pwn_reduced['sigma'], loc=0.2)
    x = np.linspace(0, 1, 100)
    y = norm.pdf(x, mean, sigma)
    plt.plot(x, y)
    print(mean, sigma)
    plt.savefig('size.png')



if __name__ == '__main__':
    table_pwn = prepare_the_pwn_table(flux_min=0.07,flux_max=100)

    print('-----------------------------------------------------')
    table_pwn_reduced = table_pwn['Source_Name','flux_crab','sigma']
    table_pwn_reduced.rename_column('flux_crab', 'spec_norm_crab')
    table_pwn_reduced.rename_column('Source_Name', 'source_name')
    #table_pwn_reduced.rename_column('Size', 'sigma')

    print('-----------------------------------------------------')
    table_sim_pwn = Table.read('ctadc_skymodel_gps_sources_pwn.ecsv', format='ascii.ecsv')
    new_table_sim_pwn, new_table_sim_pwn_reduced = remove_bright_pwn(table_sim_pwn)
    merged_table = vstack([table_pwn_reduced, new_table_sim_pwn_reduced])
    merged_table.pprint()
    plot_logN_logS(table_pwn, new_table_sim_pwn, merged_table, flux_min=0.07, flux_max=100)

    plot_size_distrib(table_pwn=table_pwn_reduced, table_sim=new_table_sim_pwn_reduced, merged_table=merged_table)