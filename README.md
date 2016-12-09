# cta-dc

CTA data challenge 2017

## What is this?

This is a repository containing code and files for the CTA data challenge.

It is very preliminary. Contributions welcome!

URL: https://github.com/gammasky/cta-dc

## Repository content

Here's an overview of the content of this repository.

For further details, see the `README.md` files in each folder,
e.g. [sky_model/README.md](sky_model/README.md]) contains a
description of the sky model we generated.

- [sky_model](sky_model) - The sky model (sources and diffuse emission components)
- [observations](observations) - Simulated list of observations for the GPS and EGS
- [irf](irf) - CTA instrument response functions
- [data](data) - Generated data (event lists and more)
- [ctadc](ctadc) - Python package with helper functions / classes
- [make.py](make.py) - Main script, that can run all steps (e.g. generate sky model and data)

## Data formats

The following data formats are used:

* `FITS` - "Flexible Image Transport System",
  see [here](https://en.wikipedia.org/wiki/FITS)
* `ECSV` - "Enhanced Character Separated Values table format",
  see [here](https://github.com/astropy/astropy-APEs/blob/master/APE6.rst)
* `YAML` - a human and machine readable data format
  see [here](https://en.wikipedia.org/wiki/YAML)
* `XML` - "Extensible Markup Language",
  see [here](https://en.wikipedia.org/wiki/XML)

Note that these are only "container" formats that can hold almost any data and
information. The specification of data formats for gamma-ray astronomy is work
in progress at http://gamma-astro-data-formats.readthedocs.io/ .
Some of the formats we use (e.g. event lists) are already described there.


## How to run the scripts?

We use Python scripts to produce the sky model and data.

You can access the files with whatever software you like.

E.g. the sky model and data generation is fully scripted and reproducible.
Anyone can use the [make.py](make.py) script to run them.
To see what's available type:

    ./make.py --help

There's no `setup.py` file, i.e. you can't install `ctadc` into your `site-packages`.
Just execute scripts from the top-level repo folder (where `make.py` is located),
and all scripts should work.

### Installation

To run the scripts, you have to install Python 2.7 or 3.4 or later
and the packages listed in `requirements.txt` via `pip`

    pip install -r requirements.txt

The following packages are also available via `conda` (or other package managers,
like `brew` or `macports` on macOS, or `apt-get` or `yum` on Linux):

    conda install click pyyaml pytest numpy scipy matplotlib astropy

And then you can install the rest with `pip`:

    pip install git+https://github.com/regions/regions.git#egg=regions
    pip install git+https://github.com/gammapy/gammapy.git#egg=gammapy
    pip install git+https://github.com/naima/naima.git#egg=naima


Here's some links to the documentation pages for these packages:

- pytest (http://docs.pytest.org/), optional, for automated tests
- Numpy (https://docs.scipy.org/doc/numpy/reference/)
- Scipy (https://docs.scipy.org/doc/scipy/reference/)
- Matplotlib (http://matplotlib.org/), optional, for plotting
- Astropy (http://www.astropy.org/)
- Gammapy (http://docs.gammapy.org/)
- Maybe Naima (http://naima.readthedocs.io/) for SED models
- Maybe GAMERA (http://joachimhahn.github.io/GAMERA/) for GPS population models

## Glossary

Commonly used abbreviations:

- `CTA` - Cherenkov telescope array
- `GPS` - Galactic plane survey
- `EGS` - Extragalactic survey
- `ST` - Science tools
- `IRF` - Instrument response function