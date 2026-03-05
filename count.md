---
label: Extract
description: Extract linked-read barcodes
category: [fastq,bam]
tags: [fastq,bam]
icon: fold-down
order: 10
---

# :icon-fold-down: Extract linked-read barcodes

You might want a list of all the unique barcodes in your data and [!badge corners="pill" text="extract"] will do just that.

### FASTQ
```bash usage
djinn fastq extract INPUT...
```

```bash example 
djinn fastq extract sample1.F.fq.gz sample1.R.fq.gz > sample1.bc
```


### SAM
```bash usage
djinn sam extract INPUT
```

```bash example
djinn sam extract sample1.bam > sample1.bc
```
