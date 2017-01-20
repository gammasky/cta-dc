"""
Convert CTA IRF files from ROOT to FITS
"""
from glob import glob
import sys
sys.path.append('/Users/deil/code/gammapy-extra/datasets/cta/')
from export_cta_prod2_root_to_fits import cta_perf_root_to_fits

# cta_perf_root_to_fits('CTA-Performance-South-0.5h_20170119.root', 'CTA-Performance-South-0.5h_20170119.fits.gz')

for root_filename in glob('*.root'):
    fits_filename = root_filename.replace('.root', '.fits.gz')
    cta_perf_root_to_fits(root_filename, fits_filename)


