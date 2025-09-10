#! /usr/bin/env python3

import rich_click as click

import fastq
import ncbi
import standardize_fastq
import standardize_bam

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def djinn():
    """
    Convert between linked-read formats and barcode styles
    """

module_docstring = {
    "dhinn": [
        {
            "name": "Commands",
            "commands": ["fastq", "ncbi", "standardize-bam", "standardize-fastq"],
            "panel_styles": {"border_style" : "magenta"}
        }
    ]
}

djinn.add_command(fastq)
djinn.add_command(ncbi)
djinn.add_command(standardize_bam)
djinn.add_command(standardize_fastq)
