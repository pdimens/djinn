#! /usr/bin/env python
"""assign molecular identifier (MI) tags to alignments based on distance and barcode"""
import sys
import pysam
import rich_click as click
from djinn.utils.file_ops import validate_sam, which_linkedread_sam, print_error, is_standard

def write_validbx(bam, alnrecord, mol_id):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    mol_id: the "mol_id" entry of a barcode dictionary
    Formats an alignment record to include the MI tag
    and the BX at the end and writes it to the output
    bam file. Replaces existing MI tag, if exists.
    '''
    # get all the tags except MI b/c it's being replaced (if exists)
    # will manually parse BX, so omit that too
    # also remove DI because it's not necessary
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    tags.append(("MI", mol_id))
    _bx = alnrecord.get_tag("BX")
    if "-" in _bx:
        # it's been deconvolved, set it to a DX tag
        tags.append(("DX", _bx))
    bx_clean = _bx.split("-")[0]
    tags.append(("BX", bx_clean))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_invalidbx(bam, alnrecord):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include the BX 
    at the end and writes it to the output
    bam file. Removes MI tag if present.
    '''
    # will not keep an existing MI tag if present
    # also remove DI because it's not necessary
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    _bx = alnrecord.get_tag("BX")
    # if hyphen is present, it's been deconvolved and shouldn't have been
    # and rm the hyphen part
    bx_clean = _bx.split("-")[0]
    tags.append(("BX", bx_clean))
    # update the record's tags
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_missingbx(bam, alnrecord, missing):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include invalid BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    '''
    # removes MI and DI tags, writes new BX tag
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    tags.append(("BX", missing))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/djinn/assign_mi")
@click.option('-c', '--cutoff', type = int, default = 0, show_default = True, help="Distance in base pairs at which alignments with the same barcode should be considered different molecules. If 0, then alignment distance is ignored.")
@click.option('-u', '--keep-unmapped', is_flag = True, default = False, show_default = True, help="Keep unmapped records")
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), callback = validate_sam)
@click.help_option('--help', hidden = True)
def assign_mi(input, cutoff, keep_unmapped):
    """
    Assign an MI:i (Molecular Identifier) tags
    
    Using a distance `cutoff`, assign barcoded alignments to unique molecules (default = 0).
    Unmapped records are discarded in the output unless `--keep-unmapped` is used. Records without a `BX:Z` tag or
    with an invalid barcode (based on `VX` tag) are preserved
    but are not assigned an `MI:i` tag. Input file must be in standard format (BX and VX tags) and
    must be coordinate sorted (such as with `samtools sort`).
    """
    if not is_standard(input):
        print_error(
            "incorrect file format",
            "Input alignment file must be in standard format, meaning it has a barcode in the [green bold]BX[/] tag and a 0/1 validation in the [green bold]VX[/] tag. Use [blue]djinn sam standardize[/] to convert the input into the standard format."
        )
    lr_type = which_linkedread_sam(input)
    if lr_type == "haplotagging":
        MISSING = "A00C00B00D00"
    elif lr_type == "stlfr":
        MISSING = "0_0_0"
    elif lr_type == "tellseq":
        MISSING = "N" * 18
    else:
        print_error(
            "undetermined chemistry",
            "Unable to determine linked-read chemistry based on barcodes in the first 100 records. Barcodes must conform to one of `haplotagging`, `stlfr`, or `tellseq`"
        )
    
    # initialize the dict
    d = {}
    # LAST_CONTIG keeps track of the last contig so we can
    # clear the dict when it's a new contig/chromosome
    LAST_CONTIG = False
    # MI is the name of the current molecule, starting at 1 (0+1)
    MI = 0

    # iniitalize input/output files
    with (
        pysam.AlignmentFile(input, check_sq=False, require_index=False) as alnfile,
        pysam.AlignmentFile(sys.stdout.buffer, "wb", template = alnfile) as outfile
    ):
        for record in alnfile.fetch(until_eof = True):
            chrm = record.reference_name
            # check if the current chromosome is different from the previous one
            # if so, empty the dict (a consideration for RAM usage)
            if LAST_CONTIG and chrm != LAST_CONTIG:
                d = {}
            if not keep_unmapped and record.is_unmapped:
                # skip, don't output
                LAST_CONTIG = chrm
                continue

            try:
                bx = record.get_tag("BX")
                if record.get_tag("VX") == 0:
                    # VX:i:0 is invalid
                    write_invalidbx(outfile, record)
                    LAST_CONTIG = chrm
                    continue
            except KeyError:
                # There is no bx or maybe no vx tag
                write_missingbx(outfile, record, MISSING)
                LAST_CONTIG = chrm
                continue

            aln = record.get_blocks()
            if not aln:
                if keep_unmapped:
                    write_invalidbx(outfile, record)
                LAST_CONTIG = chrm
                continue

            # logic to accommodate split records
            # start position of first alignment
            pos_start = aln[0][0]
            # end position of last alignment
            pos_end   = aln[-1][1]

            # create bx entry if it's not present
            if bx not in d:
                # increment MI b/c it's a new molecule
                MI += 1
                d[bx] = {
                    "lastpos" : pos_end,
                    "current_suffix": 0,
                    "mol_id": MI
                }
                # write and move on
                write_validbx(outfile, record, MI)
                LAST_CONTIG = chrm
                continue

            # store the original barcode as `orig` b/c we might need to suffix it
            orig = bx
            # if there is a suffix, append it to the barcode name
            if d[orig]["current_suffix"] > 0:
                bx = orig + "." + str(d[orig]["current_suffix"])

            # distance from last alignment = current aln start - previous aln end
            dist = pos_start - d[bx]["lastpos"]
            # handle overlapping or out-of-order alignments
            if dist < 0:
                dist = 0
            # if the distance between alignments is > cutoff, it's a different molecule
            # so we'll +1 the suffix of the original barcode and relabel this one as
            # BX + suffix. Since it's a new entry, we initialize it and move on
            if cutoff > 0 and dist > cutoff:
                # increment MI b/c it's a new molecule
                MI += 1
                # increment original barcode's suffix
                d[orig]["current_suffix"] += 1
                bx = orig + "." + str(d[orig]["current_suffix"])
                # add new entry for new suffixed barcode with unique MI
                d[bx] = {
                    "lastpos" : pos_end,
                    "current_suffix": 0,
                    "mol_id": MI
                }
                # write and move on
                write_validbx(outfile, record, MI)
                LAST_CONTIG = chrm
                continue

            if record.is_reverse or (record.is_forward and not record.is_paired):
                # set the last position to be the end of current alignment
                d[bx]["lastpos"] = pos_end

            # if it hasn't moved on by now, it's a record for an
            # existing barcode/molecule. Write the record.
            write_validbx(outfile, record, d[bx]["mol_id"])

            # update the chromosome tracker
            LAST_CONTIG = chrm
