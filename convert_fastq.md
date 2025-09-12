---
label: Convert FASTQ
description: Convert between linked-read FASTQ data formats
category: [fastq]
tags: [fastq]
icon: arrow-switch
order: 10
---

In the event you need your linked-read data converted into a different linked-read format, we got you covered.
We might disagree on the fragmented format landscape, but that doesn't mean you
shouldn't be able to use your data how and where you want to. This command converts a paired-end read set of FASTQ
files between the common linked-read FASTQ types.

```bash usage
djinn fastq PREFIX TARGET FQ1 FQ2
```

```bash example | tellseq → stlfr
djinn data/orcs_stlfr stlfr data/orcs.R1.fq.gz data/orcs.R2.fq.gz
```

Auto-detects the input format as one of haplotagging, TELLseq, or stLFR (10X is not supported as input),
and converts it to the format provided as the `TARGET` positional argument.
In all cases, a file will be created with the barcode conversion map.

### :icon-move-to-end: Conversion targets

{.compact}
| `TARGET`       | barcode format                                     | example                     |
|:---------------|:---------------------------------------------------|:----------------------------|
| `10x`          | the first N base pairs of R1, given `--barcodes`   |                             |
| `haplotagging` | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
| `stlfr`        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
| `tellseq`      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |


### :icon-terminal: Running Options
{.compact}
| argument          | description                                                                   |
|:------------------|:------------------------------------------------------------------------------|
| `PREFIX`          | [!badge variant="info" text="required"] output filename prefix                |
| `TARGET`          | [!badge variant="info" text="required"] target format for output FASTQ files  |
| `FQ1`             | [!badge variant="info" text="required"] forward reads of FASTQ pair           |
| `FQ2`             | [!badge variant="info" text="required"] reverse reads of FASTQ pair           |

