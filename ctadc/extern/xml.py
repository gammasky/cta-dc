"""
Temp file with helper functions for model XML serialisation.

Will be implemented properly in Gammapy in the coming weeks.

I'm using string templating, not `xmltodict` for now
https://github.com/martinblech/xmltodict#roundtripping
because I couldn't figure out how to make list of dicts work:

    >>> import xmltodict
    >>> xmltodict.unparse(dict(a=[dict(x=1), dict(y=2)]))
"""
import textwrap

__all__ = [
    'gammacat_to_xml'
]

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
        if source_idx == 3: break

        source = gammacat[source_idx]
        source_xml = gammacat_source_to_xml(source)
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


SOURCE_TEMPLATE = """\
<source name="{name}" type="{source_type}">
    <spectrum type="{spectrum_type}">
        <parameter name="Prefactor" value="{prefactor}"/>
        <parameter name="Index" value="{index}"/>
        <parameter name="Scale" value="{scale}"/>
    </spectrum>
    <spatialModel type="{spatial_type}">
        <parameter name="GLON" value="{glon}"/>
        <parameter name="GLAT" value="{glat}"/>
    </spatialModel>
</source>
"""


def gammacat_source_to_xml(source):
    # import IPython; IPython.embed(); 1/0
    d = source.data
    data = dict(
        name=d['common_name'],
        source_type='PointSource',

        spectrum_type='PowerLaw',
        prefactor=42,
        index=99,
        scale=1,

        spatial_type='SkyDirFunction',
        glon=d['glon'],
        glat=d['glat'],
    )
    return SOURCE_TEMPLATE.format(**data)
