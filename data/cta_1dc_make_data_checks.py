"""
TODO:
- check for consistency against XML observation lists
- check if AGN observations are wobble observations

"""
import logging
from pathlib import Path
import tempfile
import hashlib
import subprocess
import ruamel.yaml
from astropy.io import fits
from astropy.table import Table
from gammapy.extern import xmltodict
from gammapy.data import DataStore

BASE_PATH = Path('1dc/1dc')

log = logging.getLogger()

# {'@name': 'AGN', '@id': '510000', '@instrument': 'CTA',
#  'parameter': [OrderedDict([('@name', 'Pointing'), ('@ra', '35.665'), ('@dec', '43.0355')]),
#                OrderedDict([('@name', 'EnergyBoundaries'), ('@emin', '30000'), ('@emax', '50000000')]),
#                OrderedDict([('@name', 'GoodTimeIntervals'), ('@tmin', '662774400'), ('@tmax', '662776200')]),
#                OrderedDict([('@name', 'TimeReference'), ('@mjdrefi', '51544'), ('@mjdreff', '0.5'), ('@timeunit', 's'),
#                             ('@timesys', 'TT'), ('@timeref', 'LOCAL')]),
#                OrderedDict([('@name', 'RegionOfInterest'), ('@ra', '35.665'), ('@dec', '43.0355'), ('@rad', '5')]),
#                OrderedDict([('@name', 'Deadtime'), ('@deadc', '0.98')]),
#                OrderedDict([('@name', 'Calibration'), ('@database', '1dc'), ('@response', 'North_z20_50h')]),
#                OrderedDict([('@name', 'EventList'), ('@file', '$CTADATA/data/baselin


# From https://forge.in2p3.fr/projects/data-challenge-1-dc-1/wiki#Observation-pattern
docs_obs_infos = dict(
    gps=dict(n_obs=3270),
    gc=dict(n_obs=1671),
    egal=dict(n_obs=1271),
    agn=dict(n_obs=1920),
)


def write_yaml(data, path):
    """Helper function to write data to a YAML file."""
    path = Path(path)
    log.info('Writing {}'.format(path))
    with path.open('w') as fh:
        ruamel.yaml.round_trip_dump(data, fh)


class IndexFileChecker:
    def __init__(self, dataset):
        self.dataset = dataset

    def run(self):
        self.check_against_docs()
        self.check_against_xml()

    def check_against_docs(self):
        obs_table = self.obs_table
        n_obs_actual = len(obs_table)
        n_obs_expected = docs_obs_infos[self.dataset]['n_obs']
        if n_obs_actual != n_obs_expected:
            log.error(f'dataset={self.dataset!r}, n_obs_actual={n_obs_actual}, n_obs_expected={n_obs_expected}')

    def check_against_xml(self):
        xml_list = self.xml_list
        obs_table = self.obs_table
        assert len(xml_list) == len(obs_table)

        for xml_data, obs_data in zip(xml_list, obs_table):
            self.check_observation(xml_data, obs_data)

    @staticmethod
    def check_observation(xml_data, obs_data):
        if int(xml_data['@id']) != int(obs_data['OBS_ID']):
            log.error(f'Mismatch: xml obs_id: {xml_data["@id"]}, obs_data: {obs_data["OBS_ID"]}')
        # print(xml_data['parameter'][6])
        # print(obs_data['IRF'])
        assert xml_data['parameter'][6]['@response'] == obs_data['IRF']

    @property
    def xml_list(self):
        filename = BASE_PATH / f'obs/obs_{self.dataset}_baseline.xml'
        log.debug(f'Reading {filename}')
        data = xmltodict.parse(filename.read_text())
        return data['observation_list']['observation']

    @property
    def obs_table(self):
        filename = BASE_PATH / 'index' / self.dataset / 'obs-index.fits.gz'
        log.debug(f'Reading {filename}')
        return Table.read(str(filename))


def check_index_files():
    datasets = ['agn', 'egal', 'gc', 'gps']
    for dataset in datasets:
        IndexFileChecker(dataset).run()


def check_composite_index_files():
    n_obs = 8132

    filename = BASE_PATH / 'index/all/obs-index.fits.gz'
    log.debug(f'Reading {filename}')
    table = Table.read(filename)
    if len(table) != n_obs:
        log.error(f'Incorrect number of rows: actual={len(table)}, expected={n_obs}')

    filename = BASE_PATH / 'index/all/hdu-index.fits.gz'
    log.debug(f'Reading {filename}')
    table = Table.read(filename)
    if len(table) != 6 * n_obs:
        log.error(f'Incorrect number of rows: actual={len(table)}, expected={6 * n_obs}')


def check_index_files_checksums(dirname):
    log.info(f'Checking checksums for folder: {dirname}')
    for filename, md5_expected in [
        # dict(zip(['filename', 'md5'], _.split())) for _ in
        _.split() for _ in
        Path('checks/refs/checksums.txt').read_text().splitlines()
    ]:
        path = Path(dirname) / filename
        md5_actual = hashlib.md5(path.read_bytes()).hexdigest()
        if md5_actual == md5_expected:
            log.debug(f'Checksum OK: {filename} {md5_expected}')
        else:
            log.error(f'Checksum mismatch: {filename} expected={md5_expected} actual={md5_actual}')


def check_checksums():
    """Check that local files and tarball and checksums under version control are consistent."""
    check_index_files_checksums('1dc')

    with tempfile.TemporaryDirectory() as tmpdirname:
        cmd = f'tar zxf 1dc/index.tar.gz -C {tmpdirname}'
        subprocess.call(cmd, shell=True)
        check_index_files_checksums(tmpdirname)


class IndexToTextDumper:
    """Create a text representation of the index files.

    We commit this to the git repo and use git diff
    to review fixes / changes. This avoids regressions.
    """
    datasets = ['agn', 'all', 'egal', 'gc', 'gps']

    def run(self):
        for dataset in self.datasets:
            for which in ['hdu', 'obs']:
                OUT_PATH = Path('checks/refs')
                path_in = BASE_PATH / 'index' / dataset / f'{which}-index.fits.gz'
                path_out = OUT_PATH / dataset / f'{which}-index.txt'
                path_out.parent.mkdir(exist_ok=True, parents=True)

                self.dump_header(path_in, path_out)

                path_out = OUT_PATH / dataset / f'{which}-index.ecsv'
                self.dump_data(path_in, path_out)

    @staticmethod
    def dump_header(path_in, path_out):
        txt = ''
        hdu_list = fits.open(str(path_in))
        for hdu in hdu_list:
            txt += hdu.header.tostring(sep='\n', padding=False)
            txt += '\n\n**********\n\n'
        log.info(f'Writing {path_out}')
        path_out.write_text(txt)

    @staticmethod
    def dump_data(path_in, path_out):
        table = Table.read(str(path_in))
        log.info(f'Writing {path_out}')
        table.write(str(path_out), format='ascii.ecsv', overwrite=True)


class DataStoreChecks:
    datasets = ['agn', 'all', 'egal', 'gc', 'gps']

    def run(self):
        for dataset in self.datasets:
            self.check(dataset)

    @staticmethod
    def check(dataset):
        path = BASE_PATH / 'index' / dataset
        ds = DataStore.from_dir(path)

        # Errors are always the same
        # We only check the first and last run
        ds.obs_table = ds.obs_table[[0, -1]]
        ds.hdu_table = ds.hdu_table[[0, 1, 2, 3, 4, 5, -6, -5, -4, -3, -2, -1]]

        results = ds.check()
        results = (_ for _ in results if _['level'] not in {'debug'})
        results = list(results)

        filename = Path('checks/refs') / dataset / 'check.yaml'
        write_yaml(results, filename)


if __name__ == '__main__':
    logging.basicConfig(level='INFO')
    # IndexToTextDumper().run()
    # check_index_files()
    # check_composite_index_files()
    # check_checksums()
    DataStoreChecks().run()
