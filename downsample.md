---
label: Downsample
description: Downsample data by barcode
category: [fastq,bam]
tags: [fastq,bam]
icon: fold-down
order: 10
---

# :icon-fold-down: Downsample data by barcode

While downsampling (subsampling) FASTQ and BAM files is relatively simple with tools such as `awk`, `samtools`, `seqtk`, `seqkit`, etc.,
[!badge corners="pill" text="downsample"] allows you to downsample a BAM file (or paired-end FASTQ) _by barcodes_. That means you can
keep all the reads associated with `d` number of barcodes or `-d` fraction of barcodes (e.g. `-d 0.5` will downsample to 50% of all barcodes).

```bash usage
djinn downsample -d <value> [-i] PREFIX INPUTS
```

```bash example | downsample bam file by 50% of the barcodes, ignoring invalid barcodes 
djinn downsample -t 8 -d 0.5 -b BC sample1.perc50 sample1.bam
```

```bash example | downsample fastq pair, including invalid barcodes
djinn downsample -d 1000 -i sample1.sub1000 sample1.F.fq.gz sample1.R.fq.gz
```

## :icon-terminal: Running Options
{.compact}
| argument                | description                                                                                          |
|:---------------------|:-----------------------------------------------------------------------------------------------------|
| `INPUTS`             | [!badge variant="info" text="required"] One BAM file or both read files from a paired-end FASTQ pair |
| `--downsample` `-d`  | [!badge variant="info" text="required"] Number/fraction of barcodes to downsample to                          |
| `--invalid` `-i`     | Include invalid barcodes in downsampling                                                 |
| `--prefix` `-p`      | Prefix for output files                                                                              |
| `--random-seed`      | Random seed for sampling [!badge variant="secondary" text="optional"]                                |
| `--threads` `-t`     | Number of threads to use (BAM only, default: 10)                                                                              |

