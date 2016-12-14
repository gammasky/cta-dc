#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
import logging
import click
import ctadc

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


@click.group()
def cli():
    """CTA-DC command line interface script.

    For further information, see `README.md`.
    """


@cli.group('skymodels')
def cli_skymodels():
    """Make files in `sky_model` folder.
    """


@cli_skymodels.command('xml')
def make_skymodels_xml():
    """Make sky model XML files.
    """
    ctadc.sky_model.make_sky_models_xml()


@cli_skymodels.command('images')
def make_skymodels_images():
    """Make sky model images.
    """
    ctadc.sky_model.make_sky_models_images()


@cli.command('observations')
def make_observations():
    """Make files in `observations` folder.
    """
    ctadc.observations.make_all_obslists()


@cli.command('data')
def make_data():
    """Make files in `data` folder.

    (Event lists and other files.)
    """
    ctadc.data.make_all_data()


if __name__ == '__main__':
    cli()
