---
label: Filter
description: Remove singletons or invalid barcodes
category: [fastq, bam]
tags: [fastq, bam]
icon: filter
order: 10
---

# :icon-filter: Remove singletons or invalid barcodes

You may be interested in retaining only "useful" linked reads, that is, reads with a valid barcode or reads
that are not singletons (i.e. there are more than 2 reads with the same barcode). Like their names suggests, [!badge corners="pill" text="filter-invalid"]() will remove invalid barcodes and [!badge corners="pill" text="filter-singletons"]() will remove singletons.
Each of these commands can optionally retain the removed reads in separate files.


## filter-invalid

```bash usage
djinn filter-invalid [--invalid] PREFIX INPUTS...
```

```bash example | remove invalid barcodes and retain them in separate files
djinn filter-invalid --invalid echidna_valid echidna.R1.fq.gz echidna.R2.fq.gz
```

### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `INPUTS`          | [!badge variant="info" text="required"] FASTQ pair or SAM/BAM file           |
| `-c` `--cache-size` | [!badge text="FASTQ only"] [!badge variant="ghost" text="hidden"] number of reads to store before writing (bigger is faster, default: `5000`) |
| `-i` `--invalid` | write reads with invalid barcodes to `prefix.invalid.bam` or `prefix.invalid.R[12].fq.gz` |

## filter-singletons

```bash usage
djinn filter-singletons [--singletons] PREFIX INPUTS...
```

```bash example | remove singleton barcodes and retain them in separate files
djinn filter-singletons --singletons hyena_filt hyena.R1.fq.gz hyena.R2.fq.gz
```

### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `INPUTS`          | [!badge variant="info" text="required"] FASTQ pair or SAM/BAM file           |
| `-c` `--cache-size` | [!badge text="FASTQ only"] [!badge variant="ghost" text="hidden"] number of reads to store before writing (bigger is faster, default: `5000`) |
| `-s` `--singletons` | write reads with singleton barcodes to `prefix.singletons.bam` or `prefix.singletons.R[12].fq.gz` |
