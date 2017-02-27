"""
Convert PWN to XML format.
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
    <spatialModel type="RadialGaussian">
        <parameter name="GLON" value="{glon:.5f}" unit="deg"/>
        <parameter name="GLAT" value="{glat:.5f}" unit="deg"/>
        <parameter name="Sigma" value="{sigma:.5f}" unit="deg"/>
    </spatialModel>"""


SPECTRUM_TEMPLATE = """\
    <spectrum type="LogParabola">
        <parameter name="Prefactor" value="{norm:.3e}" unit="MeV-1 cm-2 s-1"/>
        <parameter name="Index" value="{index:.3f}" />
        <parameter name="Curvature" value="{curvature:.3f}"/>
        <parameter name="PivotEnergy" value="{energy:.5f}" unit="MeV"/>
      </node>
"""


# def make_table_spectrum_xml(sed_energy, sed_dnde):
#     xml_spectrum_nodes = ''
#     for energy, dnde in zip(sed_energy, sed_dnde):
#         xml_spectrum_nodes += SPECTRUM_NODE_TEMPLATE.format(
#             energy=energy,
#             dnde=dnde,
#         )
#
#     return SPECTRUM_TEMPLATE.format(xml_spectrum_nodes=xml_spectrum_nodes)


def make_pwn_xml(table):

    xml_sources = ''
    #for row in table[:2]:
    for row in table:
        if (row['sigma']>0.5):
            print(row['sigma'])


        xml_spectral = SPECTRUM_TEMPLATE.format(
            norm=row['spec_norm'],
            index=row['spec_alpha'],
            curvature = row['spec_beta'],
            energy = 1e06
            )

        # Arbitrary assumption on width of the  shell
        # TODO: ask Pierre what we should put!
        # Also check if radius is inner or outer!
        #width_fraction = 0
        #radius = u.Quantity(row['size'], 'arcmin').to('deg')
        #width = width_fraction * radius

        xml_spatial = SPATIAL_TEMPLATE.format(
            glon=row['GLON'],
            glat=row['GLAT'],
            sigma=row['sigma']
            )

       # import IPython; IPython.embed();
        source_name='pwn_{}'.format(row.index)
        #print(source_name)

        xml_source = SOURCE_TEMPLATE.format(
             source_name=source_name,
             xml_spectral=xml_spectral,
             xml_spatial=xml_spatial
             )

        xml_sources += xml_source
        #print(xml_sources)

    xml = SOURCE_LIBRARY_TEMPLATE.format(xml_sources=xml_sources)

    #print(xml)
    filename = 'ctadc_skymodel_gps_sources_pwn.xml'
    print('Writing {}'.format(filename))
    with open(filename, 'w') as fh:
         fh.write(xml)


# def read__data():
#     filename = 'ctadc_skymodel_gps_sources_.esv'
#     print('Reading {}'.format(filename))
#     table = Table.read(filename, format='ascii.ecsv')
#
#     distance, glon, glat = compute_galactic_coordinates(
#         x=table['POS_X'].quantity,
#         y=table['POS_Y'].quantity,
#         z=table['POS_Z'].quantity,
#     )
#     table['distance'] = distance
#     table['glat'] = glat
#     table['glon'] = glon
#
#     energy_array = np.array(table.meta['energy_array'])
#     sed_energy = np.tile(energy_array, reps=(len(table), 1))
#     table['sed_energy'] = sed_energy
#
#     # Copy over fluxes into array column
#     sed_dnde = np.empty_like(sed_energy)
#     for col_idx in range(50):
#         sed_dnde[:, col_idx] = table.columns[9 + col_idx]
#
#     table['sed_dnde'] = sed_dnde
#
#     return table


if __name__ == '__main__':
    table = Table.read('ctadc_skymodel_gps_sources_pwn.ecsv', format='ascii.ecsv')
    table.pprint()
    make_pwn_xml(table)
