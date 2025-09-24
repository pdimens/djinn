![djinn logo](https://raw.githubusercontent.com/pdimens/djinn/refs/heads/docs/static/djinn.png)

[![GitHub Release](https://img.shields.io/github/v/release/pdimens/djinn?style=for-the-badge&logo=anaconda&logoColor=ffffff)](https://github.com/pdimens/djinn/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-docs-fbab3a?style=for-the-badge&logo=quicklook&logoColor=ffffff)](https://pdimens.github.io/djinn)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pdimens/djinn/tests.yaml?style=for-the-badge&logo=cachet&logoColor=ffffff)](https://www.youtube.com/watch?v=F1qdBPlK9M4)

Modify and convert linked-read FASTQ and BAM files

## Extract barcodes
Get all the barcodes from the input file(s)

## Sort by barcode
Sort records by barcode rather than name or position

## Convert files
Djinn converts between linked-read data formats. You can convert between formats in terms of FASTQ type or barcode style. It supports:
- 10X
- haplotagging
- stLFR
- TELLseq
- [standard](https://pdimens.github.io/harpy/getting_started/linked_read_data/#linked-read-data-types)

## NCBI submission
NCBI strips out sequence headers from FASTQ submissions, so it would be best to convert your linked-read
FASTQ data into an unaligned BAM file, with the linked-read barcode stored in the `BX` or `BC` tag.
Djinn provides a convenience function to convert to (or from) this format, although it's really just
a basic `samtools` command.

## Hi-C spoofing
This is **extremely experimental** and it converts a paired-end linked-read fastq file pair into one that conforms to Hi-C expectations by mix-matching the R1s and R2s of reads that share a barcode.