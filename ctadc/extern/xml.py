"""
Temp file with helper functions for model XML serialisation.

Will be implemented properly in Gammapy in the coming weeks.

I'm using string templating, not `xmltodict` for now
https://github.com/martinblech/xmltodict#roundtripping
because I couldn't figure out how to make list of dicts work:

    >>> import xmltodict
    >>> xmltodict.unparse(dict(a=[dict(x=1), dict(y=2)]))
"""
import logging
import textwrap
import numpy as np
import astropy.units as u
from .models import SourceModel, MissingGammaCatDataError

__all__ = [
    'gammacat_to_xml'
]

log = logging.getLogger(__name__)


SOURCE_LIBRARY_TEMPLATE = """\
<?xml version="1.0" ?>
{header}
<source_library title="{title}">

{sources_xml}
</source_library>
"""


def gammacat_to_xml(gammacat):
    sources_xml = []
    for source_idx in range(len(gammacat.table)):
        # For debugging, just do a few sources:
        # if source_idx == 3: break

        source = gammacat[source_idx]

        # TODO: Select on source_class = Galactic once available in gamma-cat
        # For now do a quick selection in GLAT
        # import IPython; IPython.embed(); 1/0
        if np.abs(source.data['glat']) > 5 * u.deg:
            continue

        try:
            source_model = SourceModel.from_gammacat(source)
            source_xml = source_model.to_xml()
        except MissingGammaCatDataError:
            log.warning('Skipping source {} (missing data in gamma-cat)'.format(source.name))
            continue

        sources_xml.append(source_xml)

    header = textwrap.dedent("""
    <!-- Bright sources for CTA-1DC from gamma-cat -->
    """)

    data = dict(
        title='gammacat',
        header=header,
        sources_xml='\n'.join(sources_xml),
    )
    return SOURCE_LIBRARY_TEMPLATE.format(**data)
