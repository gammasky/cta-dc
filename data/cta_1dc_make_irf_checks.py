"""Make some quick check plots for the IRFs."""
from pathlib import Path
import logging
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from gammapy.irf import EffectiveAreaTable2D, EnergyDispersion2D, EnergyDependentMultiGaussPSF

log = logging.getLogger(__name__)


def check_aeff(label):
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    aeff = EffectiveAreaTable2D.read(irf_file, hdu='EFFECTIVE AREA')
    aeff.peek()
    filename = 'checks/irfs/' + label + '_aeff.png'
    log.info(f'Writing {filename}')
    plt.savefig(filename)


def check_edisp(label):
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    edisp = EnergyDispersion2D.read(irf_file, hdu='ENERGY DISPERSION')
    edisp.peek()
    filename = 'checks/irfs/' + label + '_edisp.png'
    log.info(f'Writing {filename}')
    plt.savefig(filename)


def check_psf(label):
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    edisp = EnergyDependentMultiGaussPSF.read(irf_file, hdu='POINT SPREAD FUNCTION')
    edisp.peek()
    filename = 'checks/irfs/' + label + '_psf.png'
    log.info(f'Writing {filename}')
    plt.savefig(filename)


def check_all():
    # table = Table.read('1dc/1dc/caldb/data/cta/1dc/caldb.indx')

    Path('checks/irfs').mkdir(exist_ok=True)

    labels = [
        'North_z20_50h', 'North_z40_50h',
        'South_z20_50h', 'South_z40_50h',
    ]
    for label in labels:
        check_aeff(label)
        check_edisp(label)
        check_psf(label)


if __name__ == '__main__':
    logging.basicConfig(level='INFO')

    check_all()
