"""
Convert SNRs to XML format.
"""
import numpy as np
import astropy.units as u
from astropy.table import Table
from gammapy.utils.coordinates import galactic as compute_galactic_coordinates

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


def make_table_spectrum_xml(sed_energy, sed_dnde):
    print(sed_energy.shape)
    print(sed_dnde.shape)
    xml_spectrum_nodes = ''
    for energy, dnde in zip(sed_energy, sed_dnde):
        xml_spectrum_nodes += SPECTRUM_NODE_TEMPLATE.format(
            energy=energy,
            dnde=dnde,
        )

    return SPECTRUM_TEMPLATE.format(xml_spectrum_nodes=xml_spectrum_nodes)


def make_snr_xml(table):
    # print(table.colnames)
    # table.info()
    # table.info('stats')

    xml_sources = ''
    for row in table[:2]:

        xml_spectral = make_table_spectrum_xml(
            sed_energy=row['sed_energy'],
            sed_dnde=row['sed_dnde'],
        )

        # Arbitrary assumption on width of the SNR shell
        # TODO: ask Pierre what we should put!
        # Also check if radius is inner or outer!
        width_fraction = 0
        radius = u.Quantity(row['size'], 'arcmin').to('deg')
        width = width_fraction * radius

        xml_spatial = SPATIAL_TEMPLATE.format(
            glon=row['glon'],
            glat=row['glat'],
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

    distance, glon, glat = compute_galactic_coordinates(
        x=table['POS_X'].quantity,
        y=table['POS_Y'].quantity,
        z=table['POS_Z'].quantity,
    )
    table['distance'] = distance
    table['glat'] = glat
    table['glon'] = glon

    sed_energy = np.tile([1, 2, 3], reps=(len(table), 1))
    table['sed_energy'] = sed_energy

    sed_dnde = np.tile([1, 2, 3], reps=(len(table), 1))
    table['sed_dnde'] = sed_dnde

    return table


if __name__ == '__main__':
    table = read_snr_data()
    make_snr_xml(table)
