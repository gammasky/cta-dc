"""Make some quick check plots for the IRFs."""
from pathlib import Path
import logging
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.table import Table
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


def check_bkg_rate(label):
    # Gammapy background3d class currently not working (well)
    # So we'll just work with the table data directly to make some plots
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    plt.clf()

    table = Table.read(irf_file, hdu='BACKGROUND')
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

    filename = 'checks/irfs/' + label + '_bkg_rate.png'
    log.info(f'Writing {filename}')
    plt.savefig(filename)


def check_bkg_count(label):
    """Expected normalised background counts histograms.

    This is to check the following observed distribution:
    https://forge.in2p3.fr/boards/236/topics/1824?r=2057#message-2057
    """
    # Gammapy background3d class currently not working (well)
    # So we'll just work with the table data directly to make some plots
    irf_file = '1dc/1dc/caldb/data/cta/1dc/bcf/' + label + '/irf_file.fits'
    log.info(f'Reading {irf_file}')

    plt.clf()

    table = Table.read(irf_file, hdu='BACKGROUND')
    # Columns:
    # BGD float32 (21, 36, 36) 1/s/MeV/sr
    # ENERG_LO float32        (21,)        TeV
    # DETX_LO float32        (36,)        deg
    # DETY_LO float32        (36,)        deg

    # dety = table['DETY_LO'].data.squeeze()[18]
    # print(dety)  # this shows dety == 0.0

    for idx_detx in [18, 21, 22, 23, 24, 25, 26, 27]:
        detx = table['DETX_LO'].data.squeeze()[idx_detx]
        energy = table['ENERG_LO'].data.squeeze()
        bkg = table['BGD'].data.squeeze()[:, idx_detx, 18]
        val = bkg * energy  # this is to account for equal-width log-energy bins.
        val /= val.sum()
        txt = f'offset={detx:.1f}'
        plt.plot(energy, val, label=txt)

    plt.legend(loc='upper right')
    plt.xlabel('Energy (TeV)')
    plt.xlim(0.01, 1)
    plt.ylim(0, 0.6)
    plt.ylabel('Background distribution for equal-width logE binning')
    plt.semilogx()

    filename = 'checks/irfs/' + label + '_bkg_counts.png'
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
        # check_aeff(label)
        # check_edisp(label)
        # check_psf(label)
        check_bkg_rate(label)
        check_bkg_count(label)


if __name__ == '__main__':
    logging.basicConfig(level='INFO')

    check_all()
