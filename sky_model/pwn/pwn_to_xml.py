"""
Convert PWN to XML format.
"""
from pathlib import Path
from astropy.table import Table

SOURCE_LIBRARY_TEMPLATE = """\
<?xml version="1.0" standalone="no"?>

<source_library title="CTA 1DC simulated pulsar wind nebulae">

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
    <spatialModel type="RadialGaussian">
        <parameter name="GLON" scale="1.0" value="{glon:.5f}" min="-360" max="360" free="0"/>
        <parameter name="GLAT" scale="1.0" value="{glat:.5f}" min="-90" max="90" free="0"/>
        <parameter name="Sigma" scale="1.0" value="{sigma:.5f}" min="0.01" max="10"  free="1"/>
    </spatialModel>"""

SPECTRUM_TEMPLATE = """\
    <spectrum type="LogParabola">
        <parameter name="Prefactor" scale="1e-20" value="{norm:5g}"  min="1e-07" max="10000.0" free="1"/>
        <parameter name="Index"     scale="-1"    value="{index:.5f}"  min="0.0"   max="+5.0"   free="1"/>
        <parameter name="Curvature" scale="-1"    value="{curvature:.5f}"  min="-5.0"   max="+5.0"   free="1"/>
        <parameter name="PivotEnergy" scale="1e6"   value="{energy:.1f}" min="0.01"  max="1000.0" free="0"/>
      </spectrum>

"""



def make_pwn_xml(table):
    xml_sources = ''

    for row in table:
        if (row['skip'] == 1):
            continue;
        else:
            xml_spectral = SPECTRUM_TEMPLATE.format(
                norm=row['spec_norm'] / 1e-20,
                index=row['spec_alpha'],
                curvature=row['spec_beta'],
                energy=1.0
            )

            xml_spatial = SPATIAL_TEMPLATE.format(
                glon=row['GLON'],
                glat=row['GLAT'],
                sigma=row['sigma']
            )

            source_name = row['source_name']

            xml_source = SOURCE_TEMPLATE.format(
                source_name=source_name,
                xml_spectral=xml_spectral,
                xml_spatial=xml_spatial
            )

            xml_sources += xml_source

    return SOURCE_LIBRARY_TEMPLATE.format(xml_sources=xml_sources)


if __name__ == '__main__':



    filename = 'ctadc_skymodel_gps_sources_pwn.ecsv'
    print('Reading {}'.format(filename))
    table = Table.read(filename, format='ascii.ecsv')
    table.pprint()
    for row in table:
        if (row['skip'] == 1):
            print(row['source_name'], row['int_flux_above_1TeV_cu'])

    xml = make_pwn_xml(table)

    filename = 'ctadc_skymodel_gps_sources_pwn.xml'
    print('Writing {}'.format(filename))
    print(Path(filename))

    Path(filename).write_text(xml)


    filename_composite = 'ctadc_skymodel_gps_sources_composite.ecsv'
    print('Reading {}'.format(filename_composite))
    table_composite = Table.read(filename_composite, format='ascii.ecsv')
    table_composite.pprint()

    xml_composite = make_pwn_xml(table_composite)
    filename_composite = 'ctadc_skymodel_gps_sources_composite.xml'
    print('Writing {}'.format(filename_composite))
    Path(filename_composite).write_text(xml_composite)