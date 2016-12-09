# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Make sky model.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging

__all__ = [
    'GPSSkyModel',
    'make_all_sky_models',
]

log = logging.getLogger(__name__)


class GPSSkyModel:
    """
    GPS (Galactic plane survey) sky model.
    """

    def __init__(self):
        pass

    def make(self):
        pass

    def to_xml(self):
        xml = '<sources></sources>'  # a dummy
        return xml


def make_all_sky_models():
    """Make all sky models.
    """
    gps = GPSSkyModel()
    gps.make()
    xml = gps.to_xml()

    filename = 'sky_model/ctadc_skymodel_gps.xml'
    log.info('Writing {}'.format(filename))
    with open(filename, 'w') as fh:
        fh.write(xml)
