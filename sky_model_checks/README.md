# The CTA-1DC GPS sky model

## Summary

```
     tag      n_sources     flux_1_10    
------------- --------- -----------------
    gamma-cat        56 1.49638858908e-10
image_sources        11 6.27654769846e-11
          pwn       609 7.49925383064e-11
          snr       826 4.20904896842e-11
     binaries        11 9.62148259309e-12
      pulsars        15 1.57700148891e-13
```

## Figures

### Position

* As seen from Earth:
  * [ctadc_skymodel_gps_sources_sky_positions.png](ctadc_skymodel_gps_sources_sky_positions.png)
  * [ctadc_skymodel_gps_sources_sky_positions_gps.png](ctadc_skymodel_gps_sources_sky_positions_gps.png)
  * GLON: [ctadc_skymodel_gps_sources_glon.png](ctadc_skymodel_gps_sources_glon.png)
  * GLAT: [ctadc_skymodel_gps_sources_glat.png](ctadc_skymodel_gps_sources_glat.png)
  * Distance: TODO

* Positions in the Galaxy (only for simulated sources: PWN, SNR):
  * [ctadc_skymodel_gps_sources_galactic_xy.png](ctadc_skymodel_gps_sources_galactic_xy.png)
  * [ctadc_skymodel_gps_sources_galactic_xz.png](ctadc_skymodel_gps_sources_galactic_xz.png)
  * [ctadc_skymodel_gps_sources_galactic_z.png](ctadc_skymodel_gps_sources_galactic_z.png)
  * [ctadc_skymodel_gps_sources_galactic_r.png](ctadc_skymodel_gps_sources_galactic_r.png)

### Size and integral flux

* Source sizes:
  * It's radius for SNRs (TODO: outer?) and Gaussian sigma for Gaussian shapes
  * Apparent size (deg): [ctadc_skymodel_gps_sources_size.png](ctadc_skymodel_gps_sources_size.png)
  * Physical size (parsec): [ctadc_skymodel_gps_sources_physical_size.png](ctadc_skymodel_gps_sources_physical_size.png)
  * I found it weird that there are only few small-size SNRs in the apparent size.
    So I checked that that the angular size distribution as given in Pierre's table
    matches with the angular size distribution we have in the XML sky model,
    i.e. that there's no bug in our scripts:
    [ctadc_skymodel_gps_sources_size_snr_check.png](ctadc_skymodel_gps_sources_size_snr_check.png)
  * There is no bug. The fact that there are only few SNRs with small angular size is a result
    of Pierre's model, which contains very few SNRs with small physical size, i.e. below 10 pc.
* logN-logS: [logN_logS_gammalib.png](logN_logS_gammalib.png)
  * Differential (max at 1, smoothed a bit): [ctadc_skymodel_gps_sources_logn_logs_diff.png](ctadc_skymodel_gps_sources_logn_logs_diff.png)
  * Integral: todo
* Source size vs flux

### Spectra

* Spectra for each of the components:
  * [gamma-cat](ctadc_skymodel_gps_sources_spectra_gamma-cat.png)
  * [image_sources](ctadc_skymodel_gps_sources_spectra_image_sources.png)
  * [pwn](ctadc_skymodel_gps_sources_spectra_pwn.png)
  * [snr](ctadc_skymodel_gps_sources_spectra_snr.png)
  * [binaries](ctadc_skymodel_gps_sources_spectra_binaries.png)
  * [pulsars](ctadc_skymodel_gps_sources_spectra_pulsars.png)

* Compute fraction of emission from different components at different energies.
* Which sources are the PeVatrons?
* What about faint sources from gamma-cat that are hard spectrum and have no measurement at high energy.
  Do they become very bright or even brightest in the 10 - 100 TeV energy range?

## Models

These are the model files, with all information on every source:

* gamma-cat: [xml](../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml)
* image_sources: [xml](../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml)
* pwn: [xml](../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml)
* snr: [xml](../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml)
* binaries: [xml](../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml)
* pulsars: [xml](../sky_model/pulsars/ctadc_skymodel_gps_sources_pulsars.xml)
