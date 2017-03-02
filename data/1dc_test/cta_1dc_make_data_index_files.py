"""
Make an HDU and observation index file for the CTA 1DC dataset.

Format is described here:
http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/index.html

Written by Christoph Deil, March 2, 2017
"""
import logging
from glob import glob
from astropy.table import Table

logging.basicConfig(level='DEBUG')
log = logging.getLogger()


"""
  <observation name="GPS" id="000362" instrument="CTA">
    <parameter name="EventList" file="../data/gps_baseline_000001.fits" />
    <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
  </observation>
"""


def make_observation_table():
    """
    Make a table of observations.

    TODO: the observation table should really be provided as an input to 1DC.
    The only reason we make one here is because none was provided so far.
    """
    rows = []
    filenames = glob('data/*.fits')
    for filename in filenames[:2]:
        log.debug('Reading {}'.format(filename))
        table = Table.read(filename)
        print(table.meta)
        row = dict(
            OBS_ID=42
        )
        rows.append(row)

    obs_table = Table(rows=rows)
    obs_table.meta['dataset'] = 'CTA 1DC test data'
    filename = 'obs-index.fits.gz'
    print('Writing {}'.format(filename))
    obs_table.write(filename, overwrite=True)


def make_index_files():
    filename = "data/obs_gps_baseline.xml"


if __name__ == '__main__':
    make_observation_table()
