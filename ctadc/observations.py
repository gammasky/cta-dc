"""
Make observation lists.
"""
import logging
import numpy as np
from astropy.table import Table, Column

__all__ = [
    'make_obslist_gps',
    'make_all_obslists',
]

log = logging.getLogger(__name__)


def make_obslist_gps(t_obs, n_obs_per_grid_position, glon, glat):
    """Make Galactic plane survey (GPS) observation lists.

    Parameters
    ----------
    t_obs : float
        Duration of one observation (hours)
    n_obs_per_grid_position : int
        Number of observations per grid position
    glon : `~numpy.array`
        GPS Galactic longitude grid positions (deg)
    glat : list of float
        GPS Galactic latitude grid positions (deg)

    Returns
    -------
    table : `~astropy.table.Table`
        Table with parameters for each observation.
    """
    n_pointings = len(glon) * len(glat)
    t_obs_total = t_obs * n_obs_per_grid_position * n_pointings

    table = Table()
    m = table.meta
    m['description'] = 'Simulated Galactic plane survey (GPS) observation list table.'
    m['origin'] = 'Generated with a script from https://github.com/gammasky/cta-dc'
    m['t_obs'] = '{} hour'.format(t_obs)
    m['t_obs_total'] = '{} hour'.format(t_obs_total)
    m['n_obs_per_grid_position'] = n_obs_per_grid_position
    m['n_pointings'] = n_pointings
    m['glon_range'] = '{}'.format([glon[0], glon[-1]])
    m['glat_range'] = '{}'.format([glat[0], glat[-1]])

    # Make a grid of lon / lat pointings
    glon, glat = np.meshgrid(glon, glat)
    glon = np.repeat(glon.flatten(), n_obs_per_grid_position)
    glat = np.repeat(glat.flatten(), n_obs_per_grid_position)

    t_obs = np.ones_like(glon) * t_obs
    zenith = np.random.uniform(low=0, high=50, size=glon.size)

    table['glon'] = Column(glon, unit='deg', description='Galactic longitude pointing position')
    table['glat'] = Column(glat, unit='deg', description='Galactic latitude pointing position')
    table['t_obs'] = Column(t_obs, unit='hour', description='Observation time')
    table['zenith'] = Column(zenith, unit='deg', description='Mean zenith angle of observation')

    return table


def make_all_obslists():
    """Make all obsevation lists.
    """
    glon = np.arange(start=-100, stop=100.1, step=2)
    glat = np.arange(start=-2, stop=+2.1, step=2)
    table = make_obslist_gps(t_obs=0.5, n_obs_per_grid_position=5, glon=glon, glat=glat)

    # TODO: add other obs patterns, or wobble obs on sources, or more obs on Galactic center region

    filename = 'observations/ctadc_observations_gps.ecsv'
    log.info('Writing {}'.format(filename))
    table.write(filename, format='ascii.ecsv')
