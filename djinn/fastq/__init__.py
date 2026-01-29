from djinn.fastq.sample import sample
from djinn.fastq.convert import convert
from djinn.fastq.extract import extract
from djinn.fastq.singletons import filter_singletons
from djinn.fastq.invalid import filter_invalid
from djinn.fastq.ncbi import ncbi
from djinn.fastq.spoof_hic import spoof_hic
from djinn.fastq.standardize import standardize
from djinn.fastq.sort import sort
import rich_click as click

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def fastq():
    """
    FASTQ file conversions and modifications

    In most cases, you can specify `--threads` if `pigz` is available in your PATH. If `pigz` isn't available in your path,
    djinn will fall back to `gzip` and this value will be ignored. Otherwise, the number of threads will be divided evenly
    between the number of FASTQ output files (therefore use even numbers for paired-end FASTQ files).
    """

fastq.add_command(convert)
fastq.add_command(extract)
fastq.add_command(filter_invalid)
fastq.add_command(filter_singletons)
fastq.add_command(ncbi)
fastq.add_command(sample)
fastq.add_command(sort)
fastq.add_command(spoof_hic)
fastq.add_command(standardize)