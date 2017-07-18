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
    filename_summary = 'ctadc_skymodel_gps_sources_snr_2_summary.ecsv'
    table_summary = Table.read(filename_summary, format='ascii.ecsv')
    filename_orig = 'ctadc_skymodel_gps_sources_snr_2_keep.ecsv'
    table_orig = Table.read(filename_orig, format='ascii.ecsv')

    return table_summary, table_orig



def select_needed_parameters(table_s, table_o):
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
    for row in table_s:
        fl_cu = (row['flux_1_10']/crab_flux)
        #print(fl_cu)
        flux_1_10_cu.append(fl_cu)

    print(table_o.info())
    for idx in range(len(table_o)):
        if table_o[idx]['keep'] == False:
            #print(row['keep'])
            continue;
        skip_known.append(table_o[idx]['skip'])
        dis = Quantity(table_o[idx]['distance'], 'kpc')
        distance.append(table_o[idx]['distance'])
        age.append(table_o[idx]['age'])
        ism_density.append(table_o[idx]['n0'])
        ang_size = Quantity(table_o[idx]['sigma'], 'arcmin')
        size_degree.append(ang_size.to('degree').value)

        # if (idx <10):
        #     print(ang_size, ang_size.to('degree'), dis, dis.to('pc'), size_ph, q)
        galactocentric_x.append(table_o[idx]['galactocentric_x'])
        galactocentric_y.append(table_o[idx]['galactocentric_y'])
        galactocentric_z.append(table_o[idx]['galactocentric_z'])
        galactocentric_r.append(table_o[idx]['galactocentric_r'])
        size_physical.append(table_o[idx]['size_physical'])

    #print(' distanec', len(distance))


    table_s['distance'] = Column(distance, description='distance', unit='kpc')
    table_s['galactocentric_x'] = Column(galactocentric_x, description='galactocentric_x', unit='kpc')
    table_s['galactocentric_y'] = Column(galactocentric_y, description='galactocentric_x', unit='kpc')
    table_s['galactocentric_z'] = Column(galactocentric_z, description='galactocentric_x', unit='kpc')
    table_s['galactocentric_r'] = Column(galactocentric_r, description='galactocentric_x', unit='kpc')
    table_s['age'] = Column(age, description='age', unit='kyr')
    table_s['n0'] = Column(ism_density, description='n0 column density', unit='cm-3')
    table_s['int_flux_above_1TeV_cu'] = Column(flux_1_10_cu, description='integral flux in crab unit', unit='%')
    table_s['skip'] = Column(skip_known, description='skip because already known')
    table_s['size_physical'] = Column(size_physical, description='intrinsic physical size', unit='pc')
    table_s['size_ang'] = Column(size_degree, description='angular size', unit='deg')
    return table_s

if __name__ == '__main__':
    table_s, table_o = load_tables()
    #print(table_o.info())
    #print(table_s.info())



    select_needed_parameters(table_s, table_o)
    #print(table_s.info())
    table_s.rename_column('glat','GLAT')
    table_s.rename_column('glon', 'GLON')
    table_s.rename_column('size_ang','sigma')
    print(table_s['sigma'], table_s['size'])

    #print(table_s.info())


    filename = 'ctadc_skymodel_gps_sources_snr_2_summary_all.ecsv'
    print('Writing {}'.format(filename))
    table_s.write(filename, format='ascii.ecsv', overwrite=True)






