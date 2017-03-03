"""
Make an HDU and observation index file for the CTA 1DC dataset.

Format is described here:
http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/index.html

Written by Christoph Deil, March 2, 2017
"""
from collections import OrderedDict
import logging
from glob import glob
from pprint import pprint
import numpy as np
from astropy.table import Table

logging.basicConfig(level='DEBUG')
log = logging.getLogger()


def make_observation_index_table():
    """
    Make observation index table.

    Format: http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/obs_index/index.html

    TODO: the observation table should really be provided as an input to 1DC.
    The only reason we make one here is because none was provided so far.
    """
    rows = []
    filenames = glob('data/*.fits')
    for filename in filenames[:]:
        log.debug('Reading {}'.format(filename))

        events = Table.read(filename, hdu='EVENTS')
        # print('\n\nEVENTS\n\n')
        # print(events.meta)

        # gti = Table.read(filename, hdu='GTI')
        # print('\n\nGTI\n\n')
        # pprint(gti.meta)
        # print(gti)

        row = OrderedDict(
            OBS_ID=events.meta['OBS_ID'],
            RA_PNT=events.meta['RA_PNT'],
            DEC_PNT=events.meta['DEC_PNT'],
            ZEN_PNT=90 - float(events.meta['ALT_PNT']),
            ALT_PNT=events.meta['ALT_PNT'],
            AZ_PNT=events.meta['AZ_PNT'],
            ONTIME=events.meta['ONTIME'],
            LIVETIME=events.meta['LIVETIME'],
            DEADC=events.meta['DEADC'],
            TSTART=events.meta['TSTART'],
            TSTOP=events.meta['TSTOP'],
            DATE_OBS=events.meta['DATE_OBS'],
            TIME_OBS=events.meta['TIME_OBS'],
            DATE_END=events.meta['DATE_END'],
            TIME_END=events.meta['TIME_END'],
        )

        # Not part of the spec, but good to know from which file the info comes
        row['EVENTS_FILENME'] = filename

        rows.append(row)

    names = list(rows[0].keys())
    obs_table = Table(rows=rows, names=names)
    obs_table.meta['dataset'] = 'CTA 1DC test data'
    obs_table.info()

    # TODO: remove this temp hack as soon as OBS_ID is filled by Jürgen in a good way.
    obs_table['OBS_ID'] = 362 + np.arange(len(obs_table))

    filename = 'obs-index.fits.gz'
    log.info('Writing {}'.format(filename))
    obs_table.write(filename, overwrite=True)


class ObservationDefinition:
    """
    Helper class to make HDU index table entries
    from ctools observation definitions.
    """

    def __init__(self, data):
        self.data = data

    @classmethod
    def from_xml_dict(cls, xml_dict):
        # TODO: reformat xml_dict here
        data = xml_dict
        return cls(data=data)

    def make_hdu_index_rows(self):
        rows = []
        rows.append(self.make_hdu_index_entry_events())
        rows.append(self.make_hdu_index_entry_gti())
        rows.append(self.make_hdu_index_entry_aeff())
        rows.append(self.make_hdu_index_entry_edisp())
        rows.append(self.make_hdu_index_entry_psf())
        rows.append(self.make_hdu_index_entry_bkg())
        return rows

    def make_hdu_index_entry_events(self):
        return OrderedDict(
            OBS_ID=self.obs_id,
            HDU_TYPE='events',
            HDU_CLASS='events',
            FILE_DIR=self.events_dir,
            FILE_NAME=self.events_filename,
            HDU_NAME='EVENTS',
        )

    def make_hdu_index_entry_gti(self):
        row = self.make_hdu_index_entry_events()
        row['HDU_TYPE'] = 'GTI'
        row['HDU_CLASS'] = 'GTI'
        return row

    def make_hdu_index_entry_aeff(self):
        return OrderedDict(
            OBS_ID=self.obs_id,
            HDU_TYPE='aeff',
            HDU_CLASS='aeff_2d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='EFFECTIVE AREA',
        )

    def make_hdu_index_entry_edisp(self):
        return OrderedDict(
            OBS_ID=self.obs_id,
            HDU_TYPE='edisp',
            HDU_CLASS='edisp_2d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='ENERGY DISPERSION',
        )

    def make_hdu_index_entry_psf(self):
        return OrderedDict(
            OBS_ID=self.obs_id,
            HDU_TYPE='psf',
            HDU_CLASS='psf_3gauss',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='POINT SPREAD FUNCTION',
        )

    def make_hdu_index_entry_bkg(self):
        return OrderedDict(
            OBS_ID=self.obs_id,
            HDU_TYPE='bkg',
            HDU_CLASS='bkg_3d',
            FILE_DIR=self.irf_dir,
            FILE_NAME=self.irf_filename,
            HDU_NAME='BACKGROUND',
        )

    @property
    def obs_id(self):
        return int(self.data['@id'])

    @property
    def events_dir(self):
        return 'data'

    @property
    def events_filename(self):
        return self.data['parameter'][0]['@file'].split('/')[-1]

    @property
    def irf_dir(self):
        return 'caldb/data/cta/prod3b/bcf/South_z20_50h'

    @property
    def irf_filename(self):
        return 'irf_file.fits'


# caldb/data/cta/prod3b/bcf/South_z20_50h

"""
  <observation name="GPS" id="000362" instrument="CTA">
    <parameter name="EventList" file="../data/gps_baseline_000001.fits" />
    <parameter name="Calibration" database="prod3b" response="South_z20_50h" />
  </observation>
"""


# tag = 'gps_baseline'
# database = 'prod3b'
# response = 'South_z20_50h'


class ObservationDefinitionList:
    """
    Helper class to reformat ctools observation lists to HDU index tables.

    There's no specification, but some info and examples are
    scattered throughout the ctools docs:
    http://cta.irap.omp.eu/ctools-devel/

    The goal is to have an HDU index table in this standard format:
    http://gamma-astro-data-formats.readthedocs.io/en/latest/data_storage/hdu_index/index.html

    Parameters
    ----------
    data : list
        List of `ObservationDefinition` objects
    """

    def __init__(self, data):
        self.data = data

    @classmethod
    def read(cls, filename):
        """Read XML format."""
        import xmltodict
        with open(filename) as fh:
            text = fh.read()
        content = xmltodict.parse(text)
        xml_list = content['observation_list']['observation']
        obs_defs = []
        for xml_dict in xml_list:
            obs_def = ObservationDefinition.from_xml_dict(xml_dict)
            obs_defs.append(obs_def)

        return cls(data=obs_defs)

    def make_hdu_index_table(self):
        rows = []
        for obs_def in self.data:
            rows.extend(obs_def.make_hdu_index_rows())

        names = list(rows[0].keys())
        table = Table(rows=rows, names=names)
        table.meta['dataset'] = 'CTA 1DC test data'
        # table.info()
        # table[:10].pprint()
        return table


def make_hdu_index_table():
    """
    Make HDU index table.

    TODO: the observation table should really be provided as an input to 1DC.
    The only reason we make one here is because none was provided so far.
    """
    # TODO: should we go via the obs index table or the ctools XML obsdef list?
    filename = 'obs-index.fits.gz'
    log.info('Reading {}'.format(filename))
    obs_table = Table.read(filename)

    filename = 'data/obs_gps_baseline.xml'
    log.info('Reading {}'.format(filename))
    obs_def_list = ObservationDefinitionList.read(filename)
    hdu_table = obs_def_list.make_hdu_index_table()

    filename = 'hdu-index.fits.gz'
    log.info('Writing {}'.format(filename))
    hdu_table.write(filename, overwrite=True)


if __name__ == '__main__':
    # make_observation_index_table()
    make_hdu_index_table()
