"""
Plot some simple source parameter distributions
to illustrate / check the CTA 1DC GPS sky model.
"""
import logging
from collections import OrderedDict
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
#import gammalib
from gammapy.utils.energy import Energy
from gammapy.spectrum import CrabSpectrum
from gammapy.spectrum.models import LogParabola

log = logging.getLogger(__name__)


class Histogram:
    """
    Helper class to work with 1-dim histograms.

    Splits histogram computation from plotting,
    and introduces some conveniences.
    """

    def __init__(self, bins, vals):
        self.bins = bins
        self.vals = vals

    @classmethod
    def from_points(cls, bins, points):
        vals = np.histogram(points, bins=bins)[0]
        return cls(bins=bins, vals=vals)

    def copy(self):
        return self.__class__(
            bins=self.bins.copy(),
            vals=self.vals.copy(),
        )

    @property
    def bin_centers(self):
        return 0.5 * (self.bins[:-1] + self.bins[1:])

    def smooth(self, sigma=None):
        if sigma is None:
            vals = self.vals
        else:
            vals = gaussian_filter1d(self.vals.astype(float), sigma=sigma)
        return self.__class__(bins=self.bins, vals=vals)

    def normalise(self, method='max', value=1):
        if method == 'max':
            scale = float(value) / np.max(self.vals)
        elif method == 'int':
            integral = np.sum(self.vals * np.diff(self.bins))
            scale = float(value) / integral
        elif method == 'int log':
            integral = np.sum(self.vals * np.diff(np.log10(self.bins)))
            scale = float(value) / integral
        else:
            raise ValueError('Invalid method: {}'.format(method))

        return self.__class__(bins=self.bins, vals=scale * self.vals)

    def plot_hist(self, ax, **kwargs):
        """Plot a pre-binned histogram"""
        bins = self.bins
        vals = self.vals
        # Plotting a histogram nicely ourselves is nontrivial
        # Example from http://stackoverflow.com/a/18611135/498873
        left, right, width = bins[:-1], bins[1:], np.diff(bins)
        x = np.array([left, right]).T.flatten()
        y = np.array([vals, vals]).T.flatten()
        ax.plot(x, y, **kwargs)
        # ax.bar(
        #     left=left, height=hist, width=width,
        #     label=component['tag'], alpha=0.7,
        #     align='edge', color='none',
        # )

    def plot_smooth(self, ax, **kwargs):
        x = self.bin_centers
        y = self.vals
        ax.plot(x, y, **kwargs)


class HistogramStack:
    def __init__(self, hists=None):
        self.hists = hists if hists else []

    def copy(self):
        return [h.copy() for h in self.hists]

    def add(self, hist):
        self.hists.append(hist)

    @property
    def total(self):
        bins = self.hists[0].bins
        vals = np.vstack([h.vals for h in self.hists]).sum(axis=0)
        return Histogram(bins=bins, vals=vals)

class GPSSkyModel:
    ALL_TAGS = [
            'gamma-cat',
            'image_sources',
            'pwn',
            'snr',
            'composite',
            'binaries',
            'pulsars',
        ]

    def __init__(self, data, tags=None):
        self.data = data
        self.tags = tags or self.ALL_TAGS

    @classmethod
    def load_tables(cls, tags=None):
        tags = tags or cls.ALL_TAGS
        data = []
        for tag in tags:
            data.append(getattr(cls, '_load_tables_' + tag)())
        return cls(data=data)

    @classmethod
    def _load_tables_gammacat(cls):
        tag = 'gamma-cat'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml'
        table = Table.read(filename.replace('.xml', '_summary.ecsv'), format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table)

    @classmethod
    def _load_tables_pwn(cls):
        tag = 'pwn'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table)

    @classmethod
    def _load_tables_composite(cls):
        tag = 'composite'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_composite.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table)

    @classmethod
    def _load_tables_snr(cls):
        tag = 'snr'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/snr/ctadc_skymodel_gps_sources_snr_2.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table)


    def get_component(self, tag):
        for component in self.data:
            if tag == component['tag']:
                return component
        raise ValueError('Invalid component tag: {}'.format(tag))

    def get_components(self, tags=None, skip=None):
        if tags is None:
            tags = self.tags.copy()
        if skip is not None:
            for tag in skip:
                tags.remove(tag)

        for tag in tags:
            yield self.get_component(tag)

    def plot_luminosity(self):
        fig, ax = plt.subplots()
        fig_skip, ax_skip = plt.subplots()
        bins = 50  # np.arange(0, 20, 1)
        fig_sim, ax_sim = plt.subplots()
        fig_skip_sim, ax_skip_sim = plt.subplots()

        for component in self.get_components(tags=['pwn', 'composite']):
            table = component['table']
            print(table)
            vals = []
            vals_sim = []
            vals_skip_sim = []
            vals_skip = []
            for row in table:
                dist_kpc = u.Quantity(row['distance'], 'kpc')

                dist = 4 * np.pi * (dist_kpc.to('cm')) ** 2

                energy_flux = (LogParabola(amplitude=row['spec_norm'] * u.Unit('MeV-1 cm-2 s-1'),
                                           alpha=row['spec_alpha'],
                                           beta=row['spec_beta'],
                                           reference=1 * u.TeV,
                                           ).energy_flux(
                    emin=1 * u.TeV, emax=10 * u.TeV)).to('TeV cm-2 s-1')
                luminosity_check = energy_flux * dist

                vals_sim.append(np.log10(row['luminosity']))
                vals.append(np.log10(luminosity_check.value))
                if (row['skip'] == 0):
                    vals_skip_sim.append(np.log10(row['luminosity']))
                    vals_skip.append(np.log10(luminosity_check.value))

            ax_sim.hist(
                vals_sim, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
             )

            ax_skip_sim.hist(
                vals_skip_sim, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

            ax_skip.hist(
                vals_skip, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

        ax_sim.set_xlabel('Luminosity [erg/s]')
        ax_sim.legend(loc='best')
        fig_sim.tight_layout()

        ax_skip_sim.set_xlabel('Luminosity [erg/s]')
        ax_skip_sim.legend(loc='best')
        fig_skip_sim.tight_layout()

        ax_skip.set_xlabel('Luminosity [erg/s]')
        ax_skip.legend(loc='best')
        fig_skip.tight_layout()

        ax.set_xlabel('Luminosity [erg/s]')
        ax.legend(loc='best')
        fig.tight_layout()

        filename_sim = 'ctadc_skymodel_gps_sources_luminosity_simulated.png'
        log.info('Writing {}'.format(filename_sim))
        fig_sim.savefig(filename_sim)
        plt.close(fig)

        filename_skip_sim = 'ctadc_skymodel_gps_sources_luminosity_simulated_selection.png'
        log.info('Writing {}'.format(filename_skip_sim))
        fig_skip_sim.savefig(filename_skip_sim)
        plt.close(fig_skip_sim)

        filename = 'ctadc_skymodel_gps_sources_luminosity.png'
        #log.info('Writing {}'.format(filename_check)
        fig.savefig(filename)
        plt.close(fig)

        filename_skip = 'ctadc_skymodel_gps_sources_luminosity_selection.png'
        # log.info('Writing {}'.format(filename_check)
        fig_skip.savefig(filename_skip)
        plt.close(fig_skip)


    def plot_galactic_z(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['pwn', 'composite']):
            table = component['table']
            print(table.info())
            bins = np.arange(-2, 2, 0.05)
            vals = []
            many = 0
            for row in table:
                if (row['skip'] == 1):
                    many += 1
                    continue;
                else:
                    vals.append(row['galactocentric_z'])
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'],
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric Z (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_z.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)


    gps = GPSSkyModel.load_tables(tags=['pwn', 'composite'])

    gps.plot_luminosity()
    gps.plot_galactic_z()
   #gps.print_summary()
