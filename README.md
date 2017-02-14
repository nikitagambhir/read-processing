Alignment of Sclerotinia mitochondrial genomes
==============================================

This will contain the analysis of the mitochondrial genomes for our 
*Sclerotinia* isolates. This is controlled via makefile.

To build your analysis, add your data in the Required Directories below
and run

```
make
```

This will generate the genome index, sam files, bam files, g.vcf files
and a final `GATK/res.vcf.gz`.

You can look at the Makefile for details.


The processing of the analyses are based on the workflow presented here:
https://github.com/knausb/bam_processing

This is the current structure of the repository as of 2017-02-14 16:49 CST

```
.
|-- BAMS
|-- bt2-index
|-- GVCF
|-- mitochondria_genome
|-- reads
|-- REF
|-- runfiles
|-- runs
|   |-- BOWTIE2-BUILD
|   |-- FIXMATE
|   |-- MAKE-GATK-REF
|   |-- MAKE-GVCF
|   |-- MAKE-VCF
|   |-- MAP-READS
|   |-- MARK-DUPS
|   |-- SAM-TO-BAM
|   |-- VALIDATE-BAM
|   |-- VALIDATE-DUPS
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
