"""
Make an HDU and observation index tables for the CTA 1DC dataset.

Format is described here:
http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/index.html
"""
from collections import OrderedDict
import logging
from glob import glob
from pathlib import Path
import subprocess
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack as table_vstack
from astropy.coordinates import SkyCoord

log = logging.getLogger()

BASE_PATH = Path('1dc/1dc')


def get_events_file_info(filename):
    log.debug('Reading {}'.format(filename))
    header = fits.open(filename)['EVENTS'].header

    info = OrderedDict()

    info['OBS_ID'] = header['OBS_ID']
    info['RA_PNT'] = header['RA_PNT']
    info['DEC_PNT'] = header['DEC_PNT']

    pos = SkyCoord(info['RA_PNT'], info['DEC_PNT'], unit='deg').galactic
    info['GLON_PNT'] = pos.l.deg
    info['GLAT_PNT'] = pos.b.deg

    info['ZEN_PNT'] = 90 - float(header['ALT_PNT'])
    info['ALT_PNT'] = header['ALT_PNT']
    info['AZ_PNT'] = header['AZ_PNT']
    info['ONTIME'] = header['ONTIME']
    info['LIVETIME'] = header['LIVETIME']
    info['DEADC'] = header['DEADC']
    info['TSTART'] = header['TSTART']
    info['TSTOP'] = header['TSTOP']
    info['DATE_OBS'] = header['DATE_OBS']
    info['TIME_OBS'] = header['TIME_OBS']
    info['DATE_END'] = header['DATE_END']
    info['TIME_END'] = header['TIME_END']

    info['OBJECT'] = header['OBJECT']
    info['CALDB'] = header['CALDB']
    info['IRF'] = header['IRF']

    # Not part of the spec, but good to know from which file the info comes
    info['EVENTS_FILENAME'] = filename

    # info['EVENT_COUNT'] = len(events)
    # info['EVENT_TIME_MIN'] = events['TIME'].min()
    # info['EVENT_TIME_MAX'] = events['TIME'].max()
    # info['EVENT_ENERGY_MIN'] = events['ENERGY'].min()
    # info['EVENT_ENERGY_MAX'] = events['ENERGY'].max()

    # gti = Table.read(filename, hdu='GTI')
    # info['GTI_START'] = gti['START'][0]
    # info['GTI_STOP'] = gti['STOP'][0]

    return info


class ObservationDefinition:
    """Helper class to create a row for the HDU index table."""

    def __init__(self, data):
        expected_keys = {
            'obs_id',
            'dataset',
        }
        if set(data.keys()) != expected_keys:
            log.error(data)
            raise ValueError('No no no...')

        self.data = data

    def make_hdu_index_rows(self):
        yield self.make_hdu_index_entry_events()
        yield self.make_hdu_index_entry_gti()
        yield self.make_hdu_index_entry_aeff()
        yield self.make_hdu_index_entry_edisp()
        yield self.make_hdu_index_entry_psf()
        yield self.make_hdu_index_entry_bkg()

    def make_hdu_index_entry_events(self):
        return dict(
            OBS_ID=self.data['obs_id'],
            HDU_TYPE='events',
            HDU_CLASS='events',
            FILE_DIR=self.events_dir,
            FILE_NAME=self.events_filename,
            HDU_NAME='EVENTS',
        )

    def make_hdu_index_entry_gti(self):
        row = self.make_hdu_index_entry_events()
        row['HDU_TYPE'] = 'gti'
        row['HDU_CLASS'] = 'gti'
        return row

    def make_hdu_index_entry_aeff(self):
        return dict(
            OBS_ID=self.data['obs_id'],
            HDU_TYPE='aeff',
            HDU_CLASS='aeff_2d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='EFFECTIVE AREA',
        )

    def make_hdu_index_entry_edisp(self):
        return dict(
            OBS_ID=self.data['obs_id'],
            HDU_TYPE='edisp',
            HDU_CLASS='edisp_2d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='ENERGY DISPERSION',
        )

    def make_hdu_index_entry_psf(self):
        return dict(
            OBS_ID=self.data['obs_id'],
            HDU_TYPE='psf',
            HDU_CLASS='psf_3gauss',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='POINT SPREAD FUNCTION',
        )

    def make_hdu_index_entry_bkg(self):
        return dict(
            OBS_ID=self.data['obs_id'],
            HDU_TYPE='bkg',
            HDU_CLASS='bkg_3d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='BACKGROUND',
        )

    @property
    def base_dir(self):
        """Base dir for all files, what Jürgen calls CTADATA.

        We use relative paths instead of relying on that env variable.
        TODO: good choice or change?
        """
        return '../..'

    @property
    def events_dir(self):
        return self.base_dir + '/data/baseline/' + self.data['dataset']

    @property
    def events_filename(self):
        return '{}_baseline_{:06d}.fits'.format(
            self.data['dataset'],
            self.data['obs_id'],
        )

    @property
    def irf_dir(self):
        return self.base_dir + '/caldb/data/cta/1dc/bcf/South_z20_50h'

    @property
    def irf_filename(self):
        return 'irf_file.fits'


def add_provenance(meta):
    meta['observer'] = 'CTA first data challenge (1DC)'


def make_observation_index_table(dataset, out_dir, max_rows=-1, progress_bar=True):
    """
    Make observation index table.

    Format: http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/obs_index/index.html
    """
    log.info('Gathering observation index info from events for: {}'.format(dataset))

    glob_pattern = str(BASE_PATH / 'data/baseline/{}/*.fits'.format(dataset))
    log.debug('glob pattern: {}'.format(glob_pattern))
    filenames = list(glob(glob_pattern))
    log.debug('Number of files matching: {}'.format(len(filenames)))

    if max_rows > 0:
        log.warning('Selecting subset of observations: first {}'.format(max_rows))
        filenames = filenames[:max_rows]

    if progress_bar:
        from tqdm import tqdm
        filenames = tqdm(filenames)

    rows = []
    for filename in filenames:
        row = get_events_file_info(filename)
        rows.append(row)

    names = list(rows[0].keys())
    obs_table = Table(rows=rows, names=names)
    add_provenance(obs_table.meta)
    obs_table.meta['dataset'] = dataset

    filename = out_dir / 'obs-index.fits.gz'
    log.info('Writing {}'.format(filename))
    obs_table.write(filename, overwrite=True)


def make_hdu_index_table(dataset, out_dir, max_rows=-1):
    """Make HDU index table.

    This is not a general solution, it has some hard-coded stuff
    and needs the observation index file to be there already.
    """
    filename = out_dir / 'obs-index.fits.gz'
    log.info('Reading {}'.format(filename))
    obs_table = Table.read(filename)

    if max_rows > 0:
        log.warning('Selecting subset of observations: first {}'.format(max_rows))
        obs_table = obs_table[:max_rows]

    rows = []
    for obs_table_row in obs_table:
        obs_def = ObservationDefinition(data=dict(
            dataset=dataset,
            obs_id=obs_table_row['OBS_ID']
        ))

        rows.extend(obs_def.make_hdu_index_rows())

    names = ['OBS_ID', 'HDU_TYPE', 'HDU_CLASS', 'FILE_DIR', 'FILE_NAME', 'HDU_NAME']

    hdu_table = Table(rows=rows, names=names)
    add_provenance(hdu_table.meta, dataset)
    hdu_table.meta['dataset'] = dataset

    filename = out_dir / 'hdu-index.fits.gz'
    log.info('Writing {}'.format(filename))
    hdu_table.write(filename, overwrite=True)


def make_concatenated_index_files():
    """Make index files for all observations combined."""
    datasets = ['agn', 'egal', 'gc', 'gps']

    (BASE_PATH / 'index/all').mkdir(exist_ok=True)

    table = table_vstack([
        Table.read(BASE_PATH / 'index' / dataset / 'obs-index.fits.gz')
        for dataset in datasets
    ], metadata_conflicts='silent')
    table.meta = OrderedDict()
    add_provenance(table.meta)
    table.meta['dataset'] = 'all'

    filename = BASE_PATH / 'index/all/obs-index.fits.gz'
    log.info('Writing {}'.format(filename))
    table.write(filename, overwrite=True)

    table = table_vstack([
        Table.read(BASE_PATH / 'index' / dataset / 'hdu-index.fits.gz')
        for dataset in datasets
    ], metadata_conflicts='silent')
    table.meta = OrderedDict()
    add_provenance(table.meta)
    table.meta['dataset'] = 'all'

    filename = BASE_PATH / 'index/all/hdu-index.fits.gz'
    log.info('Writing {}'.format(filename))
    table.write(filename, overwrite=True)


def make_tarball():
    """Make index file tarball, ready for upload to CTA server."""
    cmd = 'cd 1dc; tar zcvf index.tar.gz 1dc/index'
    log.info('Executing: {}'.format(cmd))
    subprocess.call(cmd, shell=True)


def main():
    # Config options
    out_base_dir = BASE_PATH / 'index'
    datasets = ['agn', 'egal', 'gc', 'gps']
    max_rows = -1
    progress_bar = True
    loglevel = 'INFO'
    tarball = True

    # For debugging the script, use these options:
    # out_base_dir = BASE_PATH / 'index-test'
    # datasets = ['agn', 'gps']
    # max_rows = 10
    # tarball = False

    # Execute steps
    logging.basicConfig(level=loglevel)

    # for dataset in datasets:
    #     out_dir = out_base_dir / dataset
    #     if not out_dir.is_dir():
    #         log.info('Making directory: {}'.format(out_dir))
    #         out_dir.mkdir(parents=True)
    #
    #     make_observation_index_table(dataset, out_dir, max_rows=max_rows, progress_bar=progress_bar)
    #     make_hdu_index_table(dataset, out_dir, max_rows=max_rows)

    make_concatenated_index_files()

    if tarball:
        make_tarball()


if __name__ == '__main__':
    main()
