from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
import astropy.units as u
from astropy.units import Quantity
from gammapy.spectrum.models import LogParabola
from gammapy.spectrum import CrabSpectrum

def define_flux_crab_above_energy(emin=1 * u.TeV, emax=10 * u.TeV):
    crab = CrabSpectrum('meyer').model
    crabMAGIC = LogParabola(amplitude=3.23e-11 * u.Unit('cm-2 s-1 TeV-1'), reference=1 * u.TeV, alpha=2.47, beta=0.24)
    crab_flux_above_1TeV = crabMAGIC.integral(emin=emin, emax=emax)
    crab_flux_above_1TeV_model = crab.integral(emin=emin, emax=emax)

    return crab_flux_above_1TeV, crab_flux_above_1TeV_model


def load_tables():
    #filename_summary = 'ctadc_skymodel_gps_sources_snr_2_summary.ecsv'
    #table_summary = Table.read(filename_summary, format='ascii.ecsv')
    filename_orig = 'ctadc_skymodel_gps_sources_snr_1.ecsv'
    table_orig = Table.read(filename_orig, format='ascii.ecsv')

    return table_orig



def select_needed_parameters(table_o):
    crab_flux_above_1TeV, crab_flux_above_1TeV_model = define_flux_crab_above_energy();
    crab_flux = crab_flux_above_1TeV.value
    q = Quantity(0.000291, unit='arcmin**-1')
    flux_1_10_cu = []
    distance = []
    age = []
    ism_density = []
    skip_known = []
    galactocentric_x = []
    galactocentric_y = []
    galactocentric_z = []
    galactocentric_r = []
    size_physical = []
    size_degree = []
    glat = []
    glon = []
    flux_q = []
    #for row in table_o:
    #    fl_cu = ()
    #   #print(fl_cu)
    #    flux_1_10_cu.append(fl_cu)


    print(table_o.info())
    for idx in range(len(table_o)):
        if table_o[idx]['skip'] == 1:
            #print(row['keep'])
            continue;

        skip_known.append(table_o[idx]['skip'])
        dis = Quantity(table_o[idx]['distance'], 'kpc')
        distance.append(table_o[idx]['distance'])
        age.append(table_o[idx]['age'])
        ism_density.append(table_o[idx]['n0'])
        ang_size = Quantity(table_o[idx]['sigma'], 'arcmin')
        size_degree.append(ang_size.to('degree').value)
        flux_q.append(table_o[idx]['flux_1_10'])
        print(table_o[idx]['flux_1_10'])
        # if (idx <10):
        #     print(ang_size, ang_size.to('degree'), dis, dis.to('pc'), size_ph, q)
        galactocentric_x.append(table_o[idx]['galactocentric_x'])
        galactocentric_y.append(table_o[idx]['galactocentric_y'])
        galactocentric_z.append(table_o[idx]['galactocentric_z'])
        galactocentric_r.append(table_o[idx]['galactocentric_r'])
        glon.append(table_o[idx]['glon'])
        glat.append(table_o[idx]['glat'])
        size_physical.append(table_o[idx]['size_physical'])
        print(table_o[idx]['size_physical'], table_o[idx]['flux_1_10'])

    tab = Table([distance])
    tab['distance']=Column(distance, description='distance', unit='kpc')
    tab['GLAT'] = Column(glat, description='Latitude')
    tab['GLON'] = Column(glon, description='Latitude')
    tab['galactocentric_x'] = Column(galactocentric_x, description='galactocentric_x', unit='kpc')
    tab['galactocentric_y'] = Column(galactocentric_y, description='galactocentric_x', unit='kpc')
    tab['galactocentric_z'] = Column(galactocentric_z, description='galactocentric_x', unit='kpc')
    tab['galactocentric_r'] = Column(galactocentric_r, description='galactocentric_x', unit='kpc')
    tab['age'] = Column(age, description='age', unit='kyr')
    tab['n0'] = Column(ism_density, description='n0 column density', unit='cm-3')
    tab['int_flux_above_1TeV_cu'] = Column(flux_1_10_cu, description='integral flux in crab unit', unit='%')
    tab['skip'] = Column(skip_known, description='skip because already known')
    tab['size_physical'] = Column(size_physical, description='intrinsic physical size', unit='pc')
    tab['sigma'] = Column(size_degree, description='angular size', unit='deg')
    #tab['int_flux_above_1TeV_cu'] = Column(flux_1_10_cu, description='Integral flux between 1 and 10 TeV in crab units')
    tab['flux_1_10'] = Column(flux_q, description='Integral flux between 1 and 10 TeV', unit='cm-2 s-1')

    tab.remove_column('col0')
    print('----------------------------')
    print(tab)
    return tab

if __name__ == '__main__':
    table_o = load_tables()
    print(table_o.info())
    #print(table_s.info())

    print(table_o['flux_1_10'])

    #tab = select_needed_parameters(table_o)
    #print(tab.info())
  #  tab.rename_column('glat','GLAT')
  #  tab.rename_column('glon', 'GLON')
   # tab.rename_column('size_ang','sigma')
    #print(table_s['sigma'], table_s['size'])

    #print(table_s.info())


    #filename = 'ctadc_skymodel_gps_sources_snr_1_summary_all.ecsv'
    #print('Writing {}'.format(filename))
    #tab.write(filename, format='ascii.ecsv', overwrite=True)






