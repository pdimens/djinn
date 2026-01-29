#! /usr/bin/env python3

import rich_click as click
from djinn.fastq import fastq
from djinn.sam import sam

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

@click.group(options_metavar='')
@click.command_panel("Subcommands", commands = ["fastq", "sam"])
@click.rich_config(config)
@click.help_option('--help', hidden = True)
@click.version_option("0.0.0", prog_name="djinn", hidden = True)
def cli():
    """
    Convert between linked-read formats and barcode styles

    Use the subcommands (e.g. `djinn sam ...`) to explore the options
    for FASTQ and SAM files.
    """

cli.add_command(fastq)
cli.add_command(sam)