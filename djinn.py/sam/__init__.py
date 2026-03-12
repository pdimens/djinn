from .assign_mi import assign_mi
from .concat import concat
from .count import count
from .extract import extract
from .invalid import filter_invalid
from .ncbi import ncbi
from .sample import sample
from .singletons import filter_singletons
from .sort import sort
from .standardize import standardize
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
sam.add_command(concat)
sam.add_command(count)
sam.add_command(extract)
sam.add_command(filter_invalid)
sam.add_command(filter_singletons)
sam.add_command(ncbi)
sam.add_command(sample)
sam.add_command(sort)
sam.add_command(standardize)