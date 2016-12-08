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


@cli.command('make-obslists')
def make_obslists():
    """Make observation lists.
    """
    ctadc.make_obslists()


if __name__ == '__main__':
    cli()
