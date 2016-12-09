# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Make sky model.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
from gammapy.catalog import SourceCatalogGammaCat

__all__ = [
    'GPSSkyModelSourcesBright',
    'GPSSkyModelSourcesFaint',
    'make_all_sky_models',
]

log = logging.getLogger(__name__)


class SkyModelMixin:
    """
    Helper class for sky models.

    This way of structuring code is called a "mix-in" class in Python.
    """
    def write_xml(self):
        xml = self.xml
        filename = self.filename

        log.info('Writing {}'.format(filename))
        with open(filename, 'w') as fh:
            fh.write(xml)


class GPSSkyModelSourcesBright(SkyModelMixin):
    """
    GPS sky model - bright sources component.

    https://github.com/gammapy/gamma-cat/
    """
    filename = 'sky_model/ctadc_skymodel_gps_sources_bright.xml'

    def __init__(self):
        pass

    def make(self):
        self.make_xml()

    def make_xml(self):
        self.xml = '<sources></sources>'  # a dummy


class GPSSkyModelSourcesFaint(SkyModelMixin):
    """
    GPS sky model - faint sources component.

    """
    filename = 'sky_model/ctadc_skymodel_gps_sources_faint.xml'

    def __init__(self):
        pass

    def make(self):
        self.make_xml()

    def make_xml(self):
        self.xml = '<sources></sources>'  # a dummy



def make_all_sky_models():
    """Make all sky models.
    """
    gps_sources_bright = GPSSkyModelSourcesBright()
    gps_sources_bright.make()
    gps_sources_bright.write_xml()

    gps_sources_faint = GPSSkyModelSourcesFaint()
    gps_sources_faint.make()
    gps_sources_faint.write_xml()
