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

* Sky positions:
  * [ctadc_skymodel_gps_sources_sky_positions.png](ctadc_skymodel_gps_sources_sky_positions.png)
  * [ctadc_skymodel_gps_sources_sky_positions_gps.png](ctadc_skymodel_gps_sources_sky_positions_gps.png)
* GLON: [ctadc_skymodel_gps_sources_glon.png](ctadc_skymodel_gps_sources_glon.png)
* GLAT: [ctadc_skymodel_gps_sources_glat.png](ctadc_skymodel_gps_sources_glat.png)
* Positions in the Galaxy:
  * [ctadc_skymodel_gps_sources_galactic_xy.png](ctadc_skymodel_gps_sources_galactic_xy.png)
  * [ctadc_skymodel_gps_sources_galactic_z.png](ctadc_skymodel_gps_sources_galactic_z.png)
* Distance: TODO

### Size and integral flux

* Source sizes: [ctadc_skymodel_gps_sources_size.png](ctadc_skymodel_gps_sources_size.png)
  * It's radius for SNRs (TODO: outer?) and Gaussian sigma for Gaussian shapes
* logN-logS: [logN_logS_gammalib.png](logN_logS_gammalib.png)
  * Differential (max at 1, smoothed a bit): [ctadc_skymodel_gps_sources_logn_logs_diff.png](ctadc_skymodel_gps_sources_logn_logs_diff.png)
  * Integral: todo
* Source size vs flux

### Spectra

* Plot all or many SEDs in one plot
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
