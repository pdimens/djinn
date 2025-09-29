#! /usr/bin/env python3

import rich_click as click

from djinn.downsample import downsample
from djinn.fastq import fastq
from djinn.extract import extract
from djinn.hicspoof import hic_spoof
from djinn.isolate import isolate
from djinn.ncbi import ncbi
from djinn.sort import sort
from djinn.standardize import standardize

config = click.RichHelpConfiguration(
    max_width=80,
    theme = "magenta2-slim",
    use_markdown=True,
    show_arguments=False,
    style_options_panel_border = "magenta",
    style_commands_panel_border = "blue",
    style_option_default= "dim",
    style_deprecated="dim red",
    options_table_column_types = ["opt_long", "opt_short", "help"],
    options_table_help_sections = ["required", "help", "default"]
)

@click.group(options_metavar='', context_settings={"help_option_names" : []})
@click.command_panel("File Conversions")
@click.command_panel("Other Tools")
@click.rich_config(config)
@click.version_option("0.0.0", prog_name="djinn", hidden = True)
def cli():
    """
    Convert between linked-read formats and barcode styles
    """

cli.add_command(downsample)
cli.add_command(extract)
cli.add_command(fastq)
cli.add_command(standardize)
cli.add_command(hic_spoof)
cli.add_command(isolate)
cli.add_command(ncbi)
cli.add_command(sort)
