"""
Convert SNRs to XML format.
"""
import numpy as np
import astropy.units as u
from astropy.table import Table

SOURCE_LIBRARY_TEMPLATE = """\
<?xml version="1.0" standalone="no"?>

<source_library title="CTA 1DC simulated supernova remnants">

{xml_sources}

</source_library>
"""

SOURCE_TEMPLATE = """
 <source name="{source_name}" type="ExtendedSource">
{xml_spectral}
{xml_spatial}
 </source>
"""

SPATIAL_TEMPLATE = """\
    <spatialModel type="RadialShell">
        <parameter name="GLON" value="{glon:.5f}" unit="deg"/>
        <parameter name="GLAT" value="{glat:.5f}" unit="deg"/>
        <parameter name="Radius" value="{radius:.5f}" unit="deg"/>
        <parameter name="Width" value="{width:.5f}" unit="deg"/>
    </spatialModel>"""

SPECTRUM_TEMPLATE = """\
    <spectrum type="NodeFunction">
{xml_spectrum_nodes}\
    </spectrum>"""

SPECTRUM_NODE_TEMPLATE = """\
      <node>
        <parameter name="Energy" value="{energy:.5f}" unit="TeV"/>
        <parameter name="Intensity" value="{dnde:.5g}" unit="cm-2 s-1 TeV-1"/>
      </node>
"""


def make_table_spectrum_xml(energy_list, dnde_list):
    xml_spectrum_nodes = ''
    for energy, dnde in zip(energy_list, dnde_list):
        xml_spectrum_nodes += SPECTRUM_NODE_TEMPLATE.format(
            energy=energy,
            dnde=dnde,
        )

    return SPECTRUM_TEMPLATE.format(xml_spectrum_nodes=xml_spectrum_nodes)


def make_snr_xml(data):
    table = data['table']
    # print(table.colnames)
    # table.info()
    # table.info('stats')

    xml_sources = ''
    for row in table[:3]:

        xml_spectral = make_table_spectrum_xml(
            energy_list=data['energy'],
            dnde_list=data['dnde'],
        )

        # Arbitrary assumption on width of the SNR shell
        # TODO: ask Pierre what we should put!
        # Also check if radius is inner or outer!
        width_fraction = 0
        radius = u.Quantity(row['size'], 'arcmin').to('deg')
        width = width_fraction * radius
        xml_spatial = SPATIAL_TEMPLATE.format(
            glon=42,
            glat=43,
            radius=radius.value,
            width=width.value,
        )

        source_name = 'snr_{}'.format(row.index)
        xml_source = SOURCE_TEMPLATE.format(
            source_name=source_name,
            xml_spectral=xml_spectral,
            xml_spatial=xml_spatial,
        )

        xml_sources += xml_source

    xml = SOURCE_LIBRARY_TEMPLATE.format(xml_sources=xml_sources)

    filename = 'ctadc_skymodel_gps_sources_snr.xml'
    print('Writing {}'.format(filename))
    with open(filename, 'w') as fh:
        fh.write(xml)


def read_snr_data():
    filename = 'ctadc_skymodel_gps_sources_snr.ecsv'
    table = Table.read(filename, format='ascii.ecsv')
    data = dict(
        table=table,
        energy=np.linspace(1, 3, 5),
        dnde=np.logspace(-10, -8, 5),
    )
    return data


if __name__ == '__main__':
    data = read_snr_data()
    make_snr_xml(data)
