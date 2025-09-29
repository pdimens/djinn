---
label: Home
description: Using Djinn to convert linked-read data formats
icon: home
hidden: true
---

<img src="static/djinn.png" width="600">

You're here because you have linked-read data and might want to convert it between different linked-read formats. Many of the conversions 
available in Djinn work for either FASTQ or BAM files. Djinn also includes other tools for _e.g._ extracting barcodes from linked-read
data, sorting by barcodes, etc. Nifty convenience things. 

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
Djinn provides a convenience function to convert to (or from) this format, although we make no effort to hide
the fact it's just one-liner `samtools` commands.

:::info Useless trivia
The original version of these general functions was written while waiting for repairs at a mechanic shop and it was called `lr-switcheroo`.
:::