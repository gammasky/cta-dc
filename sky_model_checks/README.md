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
```

## Figures

### Position

* Sky positions:
  * [ctadc_skymodel_gps_sources_sky_positions.png](ctadc_skymodel_gps_sources_sky_positions.png)
  * [ctadc_skymodel_gps_sources_sky_positions_gps.png](ctadc_skymodel_gps_sources_sky_positions_gps.png)
* GLON: [ctadc_skymodel_gps_sources_glon.png](ctadc_skymodel_gps_sources_glon.png)
* GLAT: [ctadc_skymodel_gps_sources_glat.png](ctadc_skymodel_gps_sources_glat.png)
* Distance: TODO 
* Galactic top-down view: TODO

### Spatial

* Source sizes: [ctadc_skymodel_gps_sources_size.png](ctadc_skymodel_gps_sources_size.png)
  * It's radius for SNRs (TODO: outer?) and Gaussian sigma for Gaussian shapes

### Spectral

* TODO: Differential logN logS
* TODO: Integral logN logS
* TODO: Spectra

## Models

These are the model files, with all information on every source:

* gamma-cat: [xml](../sky_model/gamma-cat/ctadc_skymodel_gps_sources_gamma-cat2.xml)
* image_sources: [xml](../sky_model/image_sources/ctadc_skymodel_gps_sources_images.xml)
* pwn: [xml](../sky_model/pwn/ctadc_skymodel_gps_sources_pwn.xml)
* snr: [xml](../sky_model/snrs/ctadc_skymodel_gps_sources_snr_2.xml)
* binaries: [xml](../sky_model/binaries/ctadc_skymodel_gps_sources_binaries.xml)
