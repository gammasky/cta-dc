"""
First attempt to analyse 1DC data with Gammapy
"""
from gammapy.data import DataStore, DataStoreObservation, HDULocation


class CTA1DCDatastore(DataStore):
    def __init__(self, dataset, database, response):
        self.dataset = dataset
        self.database = database
        self.response = response


class CTA1DCObservation(DataStoreObservation):
    def __init__(self, obs_id, data_store):
        self.obs_id = obs_id
        self.data_store = data_store

    def load(self, hdu_type=None, hdu_class=None):
        location = self.location(hdu_type=hdu_type, hdu_class=hdu_class)
        return location.load()

    def location(self, hdu_type=None, hdu_class=None):
        if hdu_type == 'events':
            filename = self.data_store.get_filename(obs_id=obs_id, hdu_type=hdu_type)
        elif hdu_type == 'aeff':
            raise NotImplementedError
        elif hdu_type == 'psf':
            raise NotImplementedError

        return filename


class CTA1DCLocation(HDULocation):



def make_sky_image(config):
    pass


if __name__ == '__main__':
    datastore = CTA1DCDatastore(
        dataset='gps_baseline',
        database="prod3b",
        response="South_z20_50h",
    )

    for obs_id in [1, 2]:
        obs = datastore.obs(obs_id=obs_id)
        print(obs.events)

    # config = dict(
    #     dataset=dataset,
    #     obs_id=[1, 10, 20],
    # )
    # make_sky_image(config)
