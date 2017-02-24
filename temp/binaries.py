from astropy.table import Table


lightcurve = Table.read('../data/N5.csv', format='ascii.no_header', delimiter=';', comment='#')

lightcurve.rename_column('col1', 'PHASE')
lightcurve.rename_column('col2', 'day')
lightcurve.rename_column('col3', 'index')
lightcurve.rename_column('col4', 'flux_at_1_TeV')
lightcurve.rename_column('col5', 'cutoff')
lightcurve.rename_column('col7', 'flux_above_350GeV')
lightcurve.rename_column('col8', 'phase1')
lightcurve.rename_column('col9', 'flux_above_1TeV')
lightcurve.rename_column('col10', 'flux_above_10GeV')
lightcurve.rename_column('col11', 'NORM_above1TeV')
lightcurve.rename_column('col12', 'NORM')
lightcurve.pprint()


lightcurve = lightcurve[['PHASE','NORM']]

print(lightcurve)

lightcurve.write('../data/Binary_N5.fits', overwrite=True)


