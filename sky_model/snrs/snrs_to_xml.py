"""
Convert SNRs to XML format.
"""
from astropy.table import Table

SOURCE_LIBRARY_TEMPLATE = """
<?xml version="1.0" standalone="no"?>

<source_library title="CTA 1DC simulated supernova remnants">

{xml_sources}

</source_library>
"""

SOURCE_TEMPLATE = """
 <source name="{source_name}" type="PointSource">
    <spectrum type="PowerLaw">
        <parameter name="Prefactor" value="6.500000190867716e-13" unit="1 / (cm2 s TeV)"/>
        <parameter name="Index" value="-2.6500000953674316" unit=""/>
        <parameter name="Scale" value="1.0" unit="TeV"/>
    </spectrum>
    <spatialModel type="Gauss">
        <parameter name="GLON" value="{glon}" unit="deg"/>
        <parameter name="GLAT" value="{glat}" unit="deg"/>
        <parameter name="Sigma" value="{sigma}" unit="deg"/>
    </spatialModel>
 </source>
"""


def make_snr_xml():
    filename = 'SNRs_SIMULATED_SKYMODEL_TEST.ecsv'
    table = Table.read(filename, format='ascii.ecsv')

    xml_sources = ''
    for idx in [1, 2, 3]:
        pars = dict(
            source_name='snr_{}'.format(idx),
            glon=42,
            glat=43,
            sigma=99,
        )
        xml_sources += SOURCE_TEMPLATE.format(**pars)

    xml = SOURCE_LIBRARY_TEMPLATE.format(xml_sources=xml_sources)

    filename = 'ctadc_skymodel_gps_sources_snr.xml'
    print('Writing {}'.format(filename))
    with open(filename, 'w') as fh:
        fh.write(xml)


if __name__ == '__main__':
    make_snr_xml()
