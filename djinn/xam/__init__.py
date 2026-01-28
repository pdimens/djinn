from djinn.xam.sample import sample
from djinn.xam.extract import extract
from djinn.xam.filter_singletons import filter_singletons
from djinn.xam.invalid import filter_invalid
from djinn.xam.standardize import standardize
from djinn.xam.sort import sort
import rich_click as click

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def sam():
    """
    SAM/BAM file conversions and modifications

    Always writes to stdout. Use the subcommands below with the `sam` prefix, e.g.:
    ```bash
    djinn sam downsample options... args...
    ```
    """

sam.add_command(sample)
sam.add_command(extract)
sam.add_command(standardize)
sam.add_command(filter_invalid)
sam.add_command(filter_singletons)
sam.add_command(sort)