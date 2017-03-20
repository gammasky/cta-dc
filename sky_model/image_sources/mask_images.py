"""
Mask out diffuse emission and noise from images.
"""
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion, write_ds9
from gammapy.image import SkyImage


def mask_rxj1713():
    filename = 'RXJ1713_above2TeV_2016.fits.gz'
    print('Reading {}'.format(filename))
    image = SkyImage.read(filename)

    # The position from SIMBAD for "RX J1713.7-3946" and "SNR G347.3-00.5"
    # are offset by ~ 0.2 deg from the TeV shell center position
    # http://simbad.u-strasbg.fr/simbad/sim-id?Ident=SNR+G347.3-00.5
    # SkyCoord.from_name('SNR G347.3-00.5')

    # Here we put the center position from TeVCat, which was taken from a HESS paper
    # http://tevcat.uchicago.edu/?mode=1&showsrc=84
    region = CircleSkyRegion(
        center=SkyCoord(ra='17h13m33.6s', dec='-39d45m36s'),
        radius=Angle(0.7, 'deg'),
    )

    # filename = 'RXJ1713_mask.reg'
    # print('Writing {}'.format(filename))
    # write_ds9([region], filename=filename)
    #
    # fig, ax, cbar = image.plot()
    # filename = 'RXJ1713.png'
    # print('Writing {}'.format(filename))
    # fig.savefig(filename)

    mask = image.region_mask(region)
    image.data *= mask.data

    filename = 'RXJ1713_above2TeV_2016_masked.fits.gz'
    print('Writing {}'.format(filename))
    image.write(filename, overwrite=True)


if __name__ == '__main__':
    mask_rxj1713()
