#! /usr/bin/env python3

import rich_click as click

from djinn.ncbi import ncbi
from djinn.fastq import fastq
from djinn.standardize_bam import standardize_bam
from djinn.standardize_fastq import standardize_fastq

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
                "commands": ["fastq", "ncbi", "standardize-bam", "standardize-fastq"],
                "panel_styles": {"border_style": "magenta"}
            }
        ]
}

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
@click.version_option("0.0.0", prog_name="djinn")
def cli():
    """
    Convert between linked-read formats and barcode styles
    """

cli.add_command(ncbi)
cli.add_command(fastq)
cli.add_command(standardize_bam)
cli.add_command(standardize_fastq)
