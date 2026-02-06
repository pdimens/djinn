from djinn.sam.sample import sample
from djinn.sam.extract import extract
from djinn.sam.ncbi import ncbi
from djinn.sam.singletons import filter_singletons
from djinn.sam.invalid import filter_invalid
from djinn.sam.standardize import standardize
from djinn.sam.sort import sort
from djinn.sam.concat import concat
from djinn.sam.assign_mi import assign_mi
import rich_click as click

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def sam():
    """
    SAM/BAM file conversions and modifications

    Always writes to stdout. Use the subcommands below with the `sam` prefix, e.g.:
    ```bash
    djinn sam sample options... args...
    ```
    """

sam.add_command(assign_mi)
sam.add_command(sample)
sam.add_command(concat)
sam.add_command(extract)
sam.add_command(ncbi)
sam.add_command(filter_invalid)
sam.add_command(filter_singletons)
sam.add_command(sort)
sam.add_command(standardize)