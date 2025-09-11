#! /usr/bin/env python3

import rich_click as click

import fastq
import ncbi
import standardize_fastq
import standardize_bam

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
@click.version_option("0.0.0", prog_name="djinn")
def cli():
    """
    Convert between linked-read formats and barcode styles
    """



module_docstring = {
    "djinn": [
        {
            "name": "Commands",
            "commands": ["fastq", "ncbi", "standardize-bam", "standardize-fastq"],
            "panel_styles": {"border_style" : "magenta"}
        }
    ]
}

cli.add_command(fastq)
cli.add_command(ncbi)
cli.add_command(standardize_bam)
cli.add_command(standardize_fastq)
