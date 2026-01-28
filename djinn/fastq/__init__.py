from djinn.fastq.sample import sample
from djinn.fastq.convert import convert
from djinn.fastq.extract import extract
from djinn.fastq.filter_singletons import filter_singletons
from djinn.fastq.invalid import filter_invalid
from djinn.fastq.standardize import standardize
from djinn.fastq.sort import sort
import rich_click as click

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def fastq():
    """
    FASTQ file conversions and modifications

    Use the subcommands below with the `fastq` prefi, e.g.:
    
    ```bash
    djinn fastq downsample options... args...
    ```
    """

fastq.add_command(sample)
fastq.add_command(extract)
fastq.add_command(convert)
fastq.add_command(standardize)
fastq.add_command(filter_invalid)
fastq.add_command(filter_singletons)
fastq.add_command(sort)