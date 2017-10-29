"""Make some quick check plots for the IRFs."""
from pathlib import Path
import logging
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.table import Table
from gammapy.irf import EffectiveAreaTable2D
from gammapy.irf import EnergyDispersion2D
from gammapy.irf import EnergyDependentMultiGaussPSF
from gammapy.irf import Background3D

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


def check_bkg(label):
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    plt.clf()

    bkg = Background3D.read(irf_file, hdu='BACKGROUND')
    table = bkg.data.data

    raise NotImplementedError
    # Columns:
    # BGD float32 (21, 36, 36) 1/s/MeV/sr
    # ENERG_LO float32        (21,)        TeV
    # DETX_LO float32        (36,)        deg
    # DETY_LO float32        (36,)        deg
    for idx_detx in [18, 21, 22, 23, 24, 25, 26, 27]:
        detx = table['DETX_LO'].data.squeeze()[idx_detx]
        dety = table['DETY_LO'].data.squeeze()[18]
        energy = table['ENERG_LO'].data.squeeze()
        bkg = table['BGD'].data.squeeze()[:, idx_detx, 18]
        txt = f'detx={detx:.1f}, dety={dety:.1f}'
        plt.plot(energy, bkg, label=txt)

    plt.legend(loc='lower left')
    plt.xlabel('Energy (TeV)')
    plt.xlim(0.01, 1)
    plt.ylim(1e-6, 1)
    plt.ylabel('Background rate (s-1 MeV-1 sr-1)')
    plt.loglog()

    filename = 'checks/irfs/' + label + '_bkg.png'
    log.info(f'Writing {filename}')
    plt.savefig(filename)


def check_all():
    Path('checks/irfs').mkdir(exist_ok=True)

    labels = [
        'North_z20_50h', 'North_z40_50h',
        'South_z20_50h', 'South_z40_50h',
    ]
    for label in labels:
        check_aeff(label)
        check_edisp(label)
        check_psf(label)
        # check_bkg(label)


if __name__ == '__main__':
    logging.basicConfig(level='INFO')

    check_all()
