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
            'gammacat',
            #'image_sources',
            'pwn',
            'snr',
            'composite',
            #'binaries',
            #'pulsars',
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
        tag = 'gammacat'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2_summary.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table, color='black')
        print(table)

    @classmethod
    def _load_tables_images(cls):
        tag = 'image_sources'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/image_sources/ctadc_skymodel_gps_sources_images.ecvs'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table, color='magenta')

    @classmethod
    def _load_tables_pwn(cls):
        tag = 'pwn'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table, color='green')

    @classmethod
    def _load_tables_composite(cls):
        tag = 'composite'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/pwn/ctadc_skymodel_gps_sources_composite.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table, color='black')

    @classmethod
    def _load_tables_snr(cls):
        tag = 'snr'
        log.debug('Reading {}'.format(tag))
        filename = '../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.ecsv'
        table = Table.read(filename, format='ascii.ecsv')
        return dict(tag=tag, filename=filename, table=table, color='blue')

   # @classmethod
   # def _load_tables_snr(cls):
   #     tag = 'binaries'
   #     log.debug('Reading {}'.format(tag))
   #     filename = '../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.ecsv'
   #     table = Table.read(filename, format='ascii.ecsv')
   #     return dict(tag=tag, filename=filename, table=table, color='orange')


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
            vals = []
            vals_sim = []
            vals_skip_sim = []
            vals_skip = []
            for row in table:
                # dist_kpc = u.Quantity(row['distance'], 'kpc')
                #
                #
                #
                # energy_flux = (LogParabola(amplitude=row['spec_norm'] * u.Unit('MeV-1 cm-2 s-1'),
                #                            alpha=row['spec_alpha'],
                #                            beta=row['spec_beta'],
                #                            reference=1 * u.TeV,
                #                            ).energy_flux(
                #     emin=1 * u.TeV, emax=10 * u.TeV)).to('erg cm-2 s-1')
                # luminosity_check = energy_flux * 4 * np.pi * (dist_kpc.to('cm')) ** 2

                vals_sim.append(np.log10(row['luminosity']))
                #vals.append(np.log10(row['luminosity_check']))
                if (row['skip'] == 0):
                    vals_skip_sim.append(np.log10(row['luminosity']))
                    #vals_skip.append(np.log10(luminosity_check.value))

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


        filename_sim = 'ctadc_skymodel_gps_sources_luminosity_simulated.png'
        log.info('Writing {}'.format(filename_sim))
        fig_sim.savefig(filename_sim)
        plt.close(fig)

        filename_skip_sim = 'ctadc_skymodel_gps_sources_luminosity_simulated_selection.png'
        log.info('Writing {}'.format(filename_skip_sim))
        fig_skip_sim.savefig(filename_skip_sim)
        plt.close(fig_skip_sim)

        filename = 'ctadc_skymodel_gps_sources_luminosity.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

        filename_skip = 'ctadc_skymodel_gps_sources_luminosity_selection.png'
        log.info('Writing {}'.format(filename_skip))
        fig_skip.savefig(filename_skip)
        plt.close(fig_skip)

    def get_logn_logs(self, quantity, variant):

        bins = np.logspace(-3.1, 0.7, 100)

        hists = []
        for component in self.get_components():
            table = component['table']
            points = table['int_flux_above_1TeV_cu']
            hist = Histogram.from_points(bins=bins, points=points)

            if quantity == 'n':
                pass
            elif quantity == 'f':
                hist.vals = hist.vals * hist.bin_centers
            else:
                raise ValueError('Invalid quantity: {}'.format(quantity))

            if variant == 'diff':
                pass
            elif variant == 'int':
                hist.vals = np.cumsum(hist.vals[::-1])[::-1]
            else:
                raise ValueError('Invalid variant: {}'.format(variant))

            hists.append(hist)
                            # print('size: ', hists.size())

            return HistogramStack(hists)

    def plot_galactic_z(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['pwn', 'composite', 'snr']):
            table = component['table']

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
                alpha=0.8, normed=True, label=component['tag'], color=component['color'],
                )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric Z (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_z.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)


    def plot_galactic_r(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['pwn', 'composite', 'snr']):
            table = component['table']
            bins = np.arange(0, 20, 1)
            vals = []
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    vals.append(row['galactocentric_r'])
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'], color=component['color'],
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Galactocentric radius (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_r.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_xz(self):
        fig, ax = plt.subplots(figsize=(15, 5))
        x = []
        z = []
        for component in self.get_components(tags=['composite','snr', 'pwn']):
            table = component['table']
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    xx = u.Quantity(row['galactocentric_x'],'kpc')
                    zz = u.Quantity(row['galactocentric_z'], 'kpc')
                    x.append(xx.value)
                    z.append(zz.value)
            ax.scatter(x, z, label=component['tag'], s=10, alpha=0.5, color=component['color'])

        ax.set_ylim(-3, 3)
        ax.set_xlabel('Galactocentric X (kpc)')
        ax.set_ylabel('Galactocentric Z (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_xz.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_galactic_xy(self):
        fig, ax = plt.subplots()
        x = []
        y = []
        for component in self.get_components(tags=['snr','composite','pwn']):
            table = component['table']
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    xx = u.Quantity(row['galactocentric_x'], 'kpc')
                    yy = u.Quantity(row['galactocentric_y'], 'kpc')
                    x.append(xx.value)
                    y.append(yy.value)
            ax.scatter(x, y, label=component['tag'], s=10, alpha=0.5, color=component['color'])

        ax.set_ylim(-20, 20)
        ax.set_xlabel('Galactocentric X (kpc)')
        ax.set_ylabel('Galactocentric Y (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_xy.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_distance(self):
        fig, ax = plt.subplots()
        bins = 50  # np.arange(0, 20, 1)

        for component in self.get_components(tags=['snr', 'composite','pwn']):
            table = component['table']
            vals = []
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    vals.append(row['distance'])
            ax.hist(
                vals, bins=bins, histtype='step',
                alpha=0.8, normed=True, label=component['tag'], color=component['color']
            )

        # ax.set_xlim(0, 2)
        ax.set_xlabel('Distance to Earth (kpc)')
        ax.legend(loc='best')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_distance.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)



    def plot_glon(self):
        fig, ax = plt.subplots()
        bins = np.arange(-180, 181, 1)

        for component in self.get_components():#skip=['image_sources']):

            table = component['table']
            vals = []
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    points = Angle(row['GLON'], 'deg').wrap_at('180d').deg
                    vals.append(points)
            #print(vals)
            hist = Histogram.from_points(bins=bins, points=vals)
            hist = hist.normalise(method='int', value=1)
            hist = hist.smooth(sigma=3)
            hist.plot_smooth(
                ax=ax, label=component['tag'], alpha=0.8,
            )

        ax.legend(loc='best')
        ax.set_xlim(180, -180)
        ax.set_xlabel('GLON (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_glon.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)


    def plot_glat(self):
        fig, ax = plt.subplots()
        bins = np.arange(-180, 181, 1)

        for component in self.get_components():#skip=['image_sources']):

            table = component['table']
            vals = []
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    points = Angle(table['GLAT'], 'deg').deg
                    vals.append(points)
            hist = Histogram.from_points(bins=bins, points=vals)
            hist = hist.normalise(method='int', value=1)
            hist = hist.smooth(sigma=3)
            hist.plot_smooth(
                ax=ax, label=component['tag'], alpha=0.8,
            )

        ax.legend(loc='best')
        ax.set_xlim(180, -180)
        ax.set_xlabel('GLAT (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_glat.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_sky_positions(self):
        fig, ax = plt.subplots(figsize=(15, 5))
        for component in self.get_components():
            table = component['table']
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    x = Angle(table['GLON'], 'deg').wrap_at('180d').deg
                    y = table['GLAT']
            if component['tag'] in ['pwn', 'snr', 'composite']:
                s = 10
                alpha = 0.3
            else:
                s = 30
                alpha = 0.8
            ax.scatter(x, y, label=component['tag'], color=component['color'], s=s, alpha=alpha)

        ax.legend(loc='best')
        ax.set_xlim(180, -180)
        ax.set_ylim(-60, 60)
        ax.set_xlabel('GLON (deg)')
        ax.set_ylabel('GLAT (deg)')
        ax.grid()
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_sky_positions.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)

        ax.set_ylim(-8, 8)
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_sky_positions_gps.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_size(self):
        fig, ax = plt.subplots()
        bins = np.arange(0, 5, 0.01)
        for component in self.get_components(tags=['composite', 'pwn']):
            table = component['table']
            vals = []
            for row in table:
                if (row['skip'] == 1):
                    continue;
                else:
                    vals.append(row['sigma'])
            ax.hist(
                vals, bins=bins, label=component['tag'], histtype='step',
                alpha=0.8, normed=True,
            )

        ax.legend(loc='best')
        # ax.set_xlim(bins[0], bins-[1])
        ax.set_xlim(0, 2)
        ax.set_xlabel('Source apparent size (deg)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_size.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_size_physical(self):
        fig, ax = plt.subplots()
        bins = 30  # np.arange(0, 3, 0.05)
        for component in self.get_components(tags=['pwn', 'snr', 'composite']):
            table = component['table']
            vals = []
            for row in table:
                #if (row['skip'] == 1):
                #    continue;
                #else:
                vals.append(row['size_physical'])
            ax.hist(
                vals, bins=bins, label=component['tag'], histtype='step',
                alpha=0.8, normed=True,
            )

        ax.legend(loc='best')
        # ax.set_xlim(bins[0], bins-[1])
        # ax.set_xlim(0, 2)
        ax.set_xlabel('Source physical size (pc)')
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_physical_size.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)


    def plot_flux_vs_size(self):
        fig, ax = plt.subplots()

        for component in self.get_components(tags=['composite', 'pwn']):
            table = component['table']
            print('here: ',component['tag'], len(table), component['color'])
            int_flux = []
            size = []
            for row in table:
               # if (row['skip'] == 1):
               #     continue;
                #else:
                sigma = u.Quantity(row['sigma'], 'deg')
                flux = u.Quantity((row['int_flux_above_1TeV_cu']), '%')
                size.append(sigma.value)
                int_flux.append(flux.value)
            ax.scatter(size, int_flux, label=component['tag'], s=10, alpha=0.8, color=component['color'])

        ax.grid()
        ax.set_ylim(0.001, 120)
        ax.set_xlim(0.01,0.8)
        ax.semilogy()
        ax.set_xlabel('size (deg)')
        ax.set_ylabel('Integral flux %cu')
        ax.legend(loc='lower right')
        fig.tight_layout()

        filename = 'ctadc_skymodel_gps_sources_galactic_flux_size.png'
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

    def plot_spectral_models(self):


        for component in self.get_components(tags=['pwn', 'composite']):
            fig, ax = plt.subplots()
            table = component['table']
            tag = component['tag']
            vals = []
            idx = 0

            energies = Energy.equal_log_spacing(emin=0.02, emax=100, unit='TeV', nbins=40)
            for row in table:
                idx += 1
                if (idx < 100):
                    spec = LogParabola(amplitude=row['spec_norm'] * u.Unit('MeV-1 s-1 cm-2'),
                                       alpha=row['spec_alpha'],
                                       reference=1 * u.TeV,
                                       beta=row['spec_beta'],
                                       )
                    fluxes = []
                    e2dnde = []
                    for energy in energies:
                        dnde = spec.evaluate(energy=energy,
                                             amplitude=row['spec_norm'] * u.Unit('MeV-1 s-1 cm-2'),
                                             alpha=row['spec_alpha'],
                                             reference=1 * u.TeV,
                                             beta=row['spec_beta'],
                                             )
                        fluxes.append(dnde.value)
                        e2dnde.append(((energy ** 2 * dnde).to('erg cm-2 s-1')).value)
                    ax.plot(energies.value, e2dnde, color='black', alpha=0.2,lw=2)
                else:
                    break

            ax.set_title('{} spectra'.format(component['tag']))
            ax.loglog()
            ax.set_xlabel('Energy (TeV)')
            ax.set_ylabel('e2dnde (erg cm-2 s-1)')
            ax.set_ylim(2e-18, 5e-10)

            fig.tight_layout()
            filename = 'ctadc_skymodel_gps_sources_spectra_{}.png'.format(component['tag'])
            log.info('Writing {}'.format(filename))
            fig.savefig(filename)
            plt.close(fig)

    def get_logn_logs(self, quantity, variant, tags):
        print('ere')
        #bins = np.logspace(-3.1, 0.7, 1000)
        bins = np.logspace(-0.5, 2.1, 1000)
        crab_1_10 = CrabSpectrum().model.integral(1 * u.TeV, 10 * u.TeV).to('cm-2 s-1').value
        hists = []
        hists_skip = []

        for component in self.get_components(tags=tags):
            table = component['table']
            tag = component['tag']
            print('starting from ', tag)
            if tag in {'pwn', 'composite'}:
                fluxes = []
                fluxes_skip = []
                for row in table:
                    flux = row['int_flux_above_1TeV_cu']
                    fluxes.append(flux)
                    if (row['skip'] == 0):
                        fluxes_skip.append(flux)

            if ( tag == 'gammacat'):
                fluxes = []
                fluxes_skip = []
                for row in table:
                    flux = (row['flux_1_10']/crab_1_10)*100
                    fluxes.append(flux)
                    fluxes_skip.append(flux)


            hist = Histogram.from_points(bins=bins, points=fluxes)
            hist_skip = Histogram.from_points(bins=bins, points=fluxes_skip)

            if quantity == 'n':
                pass
            elif quantity == 'f':
                hist.vals = hist.vals * hist.bin_centers
                hist_skip.vals = hist_skip.vals * hist_skip.bin_centers
            else:
                raise ValueError('Invalid quantity: {}'.format(quantity))

            if variant == 'diff':
                pass
            elif variant == 'int':
                hist.vals = np.cumsum(hist.vals[::-1])[::-1]
                hist_skip.vals = np.cumsum(hist_skip.vals[::-1])[::-1]
            else:
                raise ValueError('Invalid variant: {}'.format(variant))

            hists.append(hist)
            hists_skip.append(hist_skip)


        return  HistogramStack(hists_skip), HistogramStack(hists)

    def plot_logn_logs(self, quantity, variant, sigma, tags=['gammacat', 'pwn', 'composite']):

        colors=['blue','red','green','black']
        fig, ax = plt.subplots()

        hists_skip, hists = self.get_logn_logs(quantity, variant, tags)
        hists_skip_pwn, hists_pwn = self.get_logn_logs(quantity, variant, tags=['composite', 'pwn'])
        #kwargs = dict(alpha=0.8, lw=2)
        hists_skip_pwn.total.smooth(sigma=sigma).plot_smooth(
            ax, alpha=0.8, lw=2, color='gray', linestyle='dashed',)
        hists_pwn.total.smooth(sigma=sigma).plot_smooth(
            ax, label='pwn+composite', alpha=0.8, lw=2, color='gray',  )
        #hists_skip.total.smooth(sigma=sigma).plot_smooth(ax, label='total', alpha=0.8, lw=1, color='orange')

        idx_skip = 0
        for tag, hist in zip(tags, hists_skip.hists):
            idx_skip += 1
            hist.smooth(sigma=sigma).plot_smooth(ax, alpha=0.8, lw=2, linestyle='dashed', color=colors[idx_skip],)
        idx = 0
        for tag, hist in zip(tags, hists.hists):
            idx +=1
            hist.smooth(sigma=sigma).plot_smooth(ax, alpha=0.8, lw=2, color=colors[idx], label=tag,)

        ax.set_xlabel('Integral flux 1-10 TeV (% Crab)')
        if quantity == 'n':
            ax.set_ylabel('Number of sources')
        elif quantity == 'f':
            ax.set_ylabel('Flux of sources')
        else:
            raise ValueError('Invalid quantity: {}'.format(quantity))

        ax.legend(loc='best')
        ax.grid('on')

        ax.semilogx()
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_log{}_logs_{}.png'.format(quantity, variant)
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)

        ax.loglog()
        ax.set_ylim(0,500)
       # ax.set_xlim(1, 200)
        fig.tight_layout()
        filename = 'ctadc_skymodel_gps_sources_log{}_logs_{}_logscale.png'.format(quantity, variant)
        log.info('Writing {}'.format(filename))
        fig.savefig(filename)
        plt.close(fig)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)


    gps = GPSSkyModel.load_tables(tags=['gammacat', 'pwn', 'composite', 'snr'])

    gps.plot_luminosity()
    gps.plot_galactic_z()
    gps.plot_galactic_r()
    gps.plot_galactic_xz()
    gps.plot_galactic_xy()
    gps.plot_distance()
    gps.plot_glon()
    gps.plot_glat()
    gps.plot_sky_positions()
    gps.plot_size()
    gps.plot_size_physical()
    gps.plot_flux_vs_size()
    gps.plot_logn_logs(quantity = 'n', variant = 'diff', sigma = 2)
    gps.plot_logn_logs(quantity='n', variant='int', sigma=None)
    gps.plot_logn_logs(quantity='f', variant='diff', sigma=2)
    gps.plot_logn_logs(quantity='f', variant='int', sigma=None)
    #
   #gps.print_summary()
