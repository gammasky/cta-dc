# cta-dc

CTA data challenge 2017

## What is this?

This is a repository containing code and files for the CTA data challenge.

It is very preliminary. Contributions welcome!

URL: https://github.com/gammasky/cta-dc

## How to use this?

We use Python scripts to produce the sky model and data.

You can access the files with whatever software you like.

To run the scripts, you have to install Python 2.7 or 3.4 or later
and the packages listed in `requirements.txt` via `conda` or `pip`

    pip install -r requirements.txt

## Content

Here's an overview of the content of this repository.

For further details, see the `README.md` files in each folder,
e.g. `sky_model/README.md` contains a description of the sky model
we generated.

- `sky_model` - The sky model (sources and diffuse emission components)
- `observations` - Simulated list of observations for the GPS and EGS
- `irf` - CTA instrument response functions
- `data` - Generated data (event lists and more)
- `ctadc` - Python package with helper functions / classes
- `make.py` - Main script, that can run all steps (e.g. generate sky model and data)

## Glossary

Commonly used abbreviations:

- `CTA` - Cherenkov telescope array
- `GPS` - Galactic plane survey
- `EGS` - Extragalactic survey
- `ST` - Science tools