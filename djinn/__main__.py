#! /usr/bin/env python3

import rich_click as click

from djinn.downsample import downsample
from djinn.extract import extract
from djinn.fastq import fastq
from djinn.ncbi import ncbi
from djinn.standardize import standardize

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.COMMAND_GROUPS = {
    "djinn":
        [
            {
                "name": "Commands",
                "commands": ["fastq", "ncbi", "standardize"],
                "panel_styles": {"border_style": "magenta"}
            },
            {
                "name": "Other Tools",
                "commands": ["downsample", "extract"],
                "panel_styles": {"border_style": "blue"}
            }
        ]
}

@click.group(options_metavar='', context_settings={"help_option_names" : []})
@click.version_option("0.0.0", prog_name="djinn", hidden = True)
def cli():
    """
    Convert between linked-read formats and barcode styles
    """

cli.add_command(ncbi)
cli.add_command(extract)
cli.add_command(downsample)
cli.add_command(fastq)
cli.add_command(standardize)
