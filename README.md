![djinn logo](https://raw.githubusercontent.com/pdimens/djinn/refs/heads/docs/static/djinn.png)

Convert linked-read FASTQ and BAM files

## Convert files
Djinn converts between linked-read data formats. It supports:
- 10X
- haplotagging
- stLFR
- TELLseq
- standard

You can convert between these formats in terms of FASTQ type or barcode style.

## NCBI submission
NCBI strips out sequence headers from FASTQ submissions, so it would be best to convert your linked-read
FASTQ data into an unaligned BAM file, with the linked-read barcode stored in the `BX` or `BC` tag.
Djinn provides a convenience function to convert to (or from) this format, although it's really just
a basic `samtools` command.
