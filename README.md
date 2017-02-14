Alignment of Sclerotinia mitochondrial genomes
==============================================

This will contain the analysis of the mitochondrial genomes for our 
*Sclerotinia* isolates. This is controlled via makefile.

This is the current structure of the repository as of 2017-02-14 12:00 CST

```
.
|-- BAMS
|-- bt2-index
|-- mitochondria_genome
|-- reads
|-- runfiles
|-- runs
|   |-- BOWTIE2-BUILD
|   |-- FIXMATE
|   |-- MAP-READS
|   |-- SAM-TO-BAM
|   |-- VALIDATE-BAM
|   `-- VALIDATE-SAM
|-- SAMS
`-- scripts
```

Required directories
--------------------

 - **mitochondria_genome/**:  ftp://ftp.broadinstitute.org/pub/annotation/fungi/sclerotinia_sclerotiorum/broad/genomes/sclerotinia_sclerotiorum/
 - **reads/**: Paired-end genomic data of 55 samples of *S. sclerotiorum* from BGI, in \*.fq.gz format
 
 
 Generated directories
 ---------------------
 
 - runfiles/: shell scripts for pre-processing submission scripts (kept in this directory for posterity)
 - *bt2-index/*: genome index files generated via `make index` 
 - *SAMS/*: mapped sam files generated via `make map`
 - *BAMS/*: filtered bam files
 - *runs/*: std out and std err of runs
