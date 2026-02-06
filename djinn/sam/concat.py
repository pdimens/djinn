#! /usr/bin/env python
import os
import sys
import pysam
import rich_click as click
from djinn.utils.file_ops import print_error, validate_sam, which_linkedread_sam
from djinn.utils.barcodes import haplotagging, tellseq, stlfr, tenx

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/concat")
@click.option('--bx', is_flag = True, default = True, show_default = True, help="Rewrite BX tags for uniqueness")
@click.option('-S', '--sam', is_flag = True, default = False, help = 'Output as SAM instead of BAM')
@click.argument('input', nargs = -1, required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def concat(input, bx, sam):
    """
    Molecule-aware file concatenation
    
    Concatenate records from linked-read SAM/BAM files while making sure molecule identification tags (`MI`)
    remain unique for every sample. This is a means of accomplishing the same as \'samtools cat\', except all MI/BX tags are updated
    so individuals don't have overlapping tags (which would mess up all the linked-read info). Use the \'--bx\'
    option to also rewrite BX tags such that they are unique between samples too. Use `djinn sam assign-mi` if your alignments
    do not have `MI` tags.
    """
    # Get the max number of unique haplotagging barcodes
    lrtype = which_linkedread_sam(input[0])
    if lrtype == "none":
        print_error(
            "undetermined chemistry",
            "Unable to determine linked-read chemistry based on barcodes in the first 100 records. Barcodes must conform to one of `haplotagging`, `stlfr`, or `tellseq`"
        )     

    with pysam.AlignmentFile(input[0], require_index=False) as xam_in:
        header = xam_in.header.to_dict()
    # Remove all @PG lines
    if 'PG' in header:
        del header['PG']
    # Add a new @PG line
    sys.argv[0] = os.path.basename(sys.argv[0])
    new_pg_line = {'ID': 'concatenate', 'PN': 'harpy', 'VN': '1.x', 'CL': " ".join(sys.argv)}
    if 'PG' not in header:
        header['PG'] = []
    header['PG'].append(new_pg_line)

    # update RG lines to match output filename name
    header['RG'][0]['ID'] = "concat"
    header['RG'][0]['SM'] = "concat"

    # set up a generator for the bx tags if --bx was invoked
    if bx:
        if lrtype == "tellseq":
            bc_generator = tellseq()
        if lrtype == "10x":
            bc_generator = tenx()
        elif lrtype == "stlfr":
            bc_generator = stlfr()
        else:
            bc_generator = haplotagging()
        
    with pysam.AlignmentFile(sys.stdout.buffer, "w" if sam else "wb", header = header) as bam_out:
        # current available unused MI tag
        MI_NEW = 1
        if bx:
            BX_NEW = bc_generator.next()
        else:
            BX_NEW = None
        # iterate through the bam files
        for xam in input:
            # create MI dict for this sample
            MI_LOCAL = {}
            with pysam.AlignmentFile(xam, require_index=False) as xamfile:
                for record in xamfile.fetch(until_eof=True):
                    try:
                        mi = record.get_tag("MI")
                        # if previously converted for this sample, use that
                        if mi in MI_LOCAL:
                            record.set_tag("MI", MI_LOCAL[mi][0])
                            if bx:
                                record.set_tag("bx", MI_LOCAL[mi][1])
                        else:
                            record.set_tag("MI", MI_NEW)
                            if bx:
                                record.set_tag("BX", BX_NEW)
                            # add to sample conversion dict
                            MI_LOCAL[mi] = [MI_NEW, BX_NEW]
                            # increment to next unique MI
                            MI_NEW += 1
                            if bx:
                                BX_NEW = bc_generator.next()
                    except StopIteration:
                        sys.exit(1)
                    except KeyError:
                        pass
                    bam_out.write(record)
