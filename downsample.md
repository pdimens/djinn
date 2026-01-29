---
label: sample
description: Downsample data by barcode
category: [fastq,bam]
tags: [fastq,bam]
icon: fold-down
order: 10
---

# :icon-fold-down: Downsample data by barcode

While downsampling (subsampling) FASTQ and BAM files is relatively simple with tools such as `awk`, `samtools`, `seqtk`, `seqkit`, etc.,
[!badge corners="pill" text="downsample"] allows you to downsample a BAM file (or FASTQ(s)) **by barcodes**. That means you can
keep all the reads associated with `-d` number of barcodes or fraction of barcodes (e.g. `-d 0.5` will downsample to 50% of all barcodes).

### FASTQ
```bash usage
djinn fastq sample -d <value> [-i] PREFIX INPUT
```

```bash example | downsample fastq pair, including invalid barcodes
djinn fastq sample -d 1000 -i sample1.sub1000 sample1.F.fq.gz sample1.R.fq.gz
```

### SAM
```bash example | downsample bam file by 50% of the barcodes, ignoring invalid barcodes 
djinn sam sample -t 8 -d 0.5 -b BC sample1.bam > sample1.half.bam
```

```bash usage
djinn fastq sample -d <value> [-i] PREFIX INPUT
```

## :icon-terminal: Running Options
{.compact}
| argument                | description                                                                                          |
|:---------------------|:-----------------------------------------------------------------------------------------------------|
| `prefix`      | Prefix for output files                                                                              |
| `INPUTS`             | [!badge variant="info" text="required"] One BAM file or both read files from a paired-end FASTQ pair |
| `-c` `--cache-size` | [!badge text="FASTQ only"] [!badge variant="ghost" text="hidden"] number of reads to store before writing (bigger is faster, default: `5000`) |
| `-d` `--downsample`  | [!badge variant="info" text="required"] Number/fraction of barcodes to downsample to                          |
| `-i` `--invalid`    | Include invalid barcodes in downsampling                                                 |
| `-r` `--random-seed`      | Random seed for sampling [!badge variant="secondary" text="optional"]                                |
| `-t` `--threads`    | [!badge text="SAM/BAM only"] Number of threads to use (default: 10)                                                                              |

