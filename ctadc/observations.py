"""
Make observation lists.
"""
import numpy as np
from astropy.table import Table, Column

__all__ = [
    'make_obslists_gps',
]


def make_obslists_gps(
        n_obs=2000,
        t_obs=0.5,
        glon_range=(-100, 100),
        glat_vals=[-2, 0, +2]
):
    """Make Galactic plane survey (GPS) observation lists.

    Returns
    -------
    table : `~astropy.table.Table`
        Table with parameters for each observation.
    """
    n_glat_rows = len(glat_vals)
    n_obs_per_glat_row = n_obs

    # Simulate parameters
    glon = np.ones(n_obs) *


    # Store results in table
    table = Table()
    table['glon'] = 42
    table['glat'] = Column(glon)
    table['t_obs'] = Column(t_obs, unit='hour', description='Observation time')
    m = table.meta
    m['description'] = 'Simulated Galactic plane survey (GPS) observation list table.'
    m['origin'] = 'Generated with a script from https://github.com/gammasky/cta-dc'
    m['n_obs'] = n_obs
    m['t_obs'] = t_obs
    m['glon_range'] = glon_range
    m['glat_vals'] = glat_vals

    return table