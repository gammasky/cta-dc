"""
TODO:
- check for consistency against XML observation lists
- check if AGN observations are wobble observations

"""
import logging
from pathlib import Path
from astropy.table import Table
import xmltodict

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


class IndexFileChecker:
    def __init__(self, dataset):
        self.dataset = dataset

    def run(self):
        self.check_against_xml()

    def check_against_xml(self):
        xml_list = self.xml_list
        # print(dict(xml_data[0]))
        obs_table = self.obs_table
        assert len(xml_list) == len(obs_table)

        for xml_data, obs_data in zip(xml_list, obs_table):
            self.check_observation(xml_data, obs_data)

    @staticmethod
    def check_observation(xml_data, obs_data):
        assert int(xml_data['@id']) == int(obs_data['OBS_ID'])
        # print(xml_data['parameter'][6])
        # print(obs_data['IRF'])
        assert xml_data['parameter'][6]['@response'] == obs_data['IRF']

    @property
    def xml_list(self):
        filename = BASE_PATH / 'obs/obs_{}_baseline.xml'.format(self.dataset)
        log.info('Reading {}'.format(filename))
        data = xmltodict.parse(filename.read_text())
        return data['observation_list']['observation']

    @property
    def obs_table(self):
        filename = BASE_PATH / 'index' / self.dataset / 'obs-index.fits.gz'
        return Table.read(str(filename))


def check_tarball():
    # TODO: implement checks if file is OK
    # filename = BASE_PATH / '1dc/index.tar.gz'
    pass


def main():
    # datasets = ['agn', 'egal', 'gc', 'gps']
    datasets = ['agn']

    for dataset in datasets:
        IndexFileChecker(dataset).run()


if __name__ == '__main__':
    main()
