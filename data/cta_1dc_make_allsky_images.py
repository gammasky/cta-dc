"""Make all-sky images for 1DC datasets"""
import logging
from gammapy.image import SkyImage
from gammapy.data import DataStore

log = logging.getLogger(__name__)


def make_counts_image(dataset, max_runs):
    log.info(f'Making all-sky image for dataset: {dataset}')
    data_store = DataStore.from_dir('1dc/1dc/index/{}'.format(dataset))

    image = SkyImage.empty(
        nxpix=3600, nypix=1800, binsz=0.1,
        xref=0, yref=0,
        proj='AIT', coordsys='GAL',
    )

    # Define energy band
    # energy_band = Energy([1, 10], 'TeV')

    obs_ids = list(data_store.obs_table['OBS_ID'])
    obs_ids = obs_ids[:max_runs] if max_runs > 0 else obs_ids

    # Show a progress bar
    from tqdm import tqdm
    obs_ids = tqdm(obs_ids)

    for obs_id in obs_ids:
        events = data_store.obs(obs_id).events
        image.fill_events(events)

    image.data = image.data.astype('float32')

    filename = 'checks/images/allsky_counts_{}.fits.gz'.format(dataset)
    log.info(f'Writing {filename}')
    image.write(filename, overwrite=True)


if __name__ == '__main__':
    logging.basicConfig(level='INFO')

    datasets = ['agn', 'egal', 'gc', 'gps']
    # datasets = ['agn', 'egal']
    max_runs = -1

    logging.basicConfig(level='INFO')

    for dataset in datasets:
        make_counts_image(dataset, max_runs)
