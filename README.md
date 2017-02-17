Read Processing for paired-end Illumina reads
==============================================

This repository contains a Makefile to process paired-end illumina data 
with bowtie2, samtools, picard, and gatk. Much of the analyses are based
on @knausb's [bam processing workflow][brianflow], but tweaked for a
haploid plant pathogen with an available genome (here we use it for 
*Sclerotinia sclerotiorum*). **This is designed to run on the [HCC SLURM
Cluster][HCC] using the [SLURM_Array][sarray] submission script in the
users $PATH.** There are no guarantees that this will work anywhere else.

Running the workflow
--------------------

To build your analysis, add your data to directories called **reads/** and
**mitochondria_genome/**, edit variables in make like `ROOT_DIR`, and ensure
that you have an environment variable called `$EMAIL` so you can be spammed
every time the processes start and finish. In your shell, type:

```
make
```

This will generate the genome index, sam files, bam files, g.vcf files
and a final `GATK/res.vcf.gz`. You can find a manifest of the files
generated from each sample in [manifest.txt] (Which will be updated randomly
as I work out the bugs and test things).

Adding steps to the workflow
----------------------------

If you want to add steps/rules to the workflow, you should first be familiar
or comfortable with Makefiles. Here are some helpful guides:

 - [Vince Buffalo's Makefiles in Bioinformatics][buffalo-make]
 - [Makefile Style Guide][make-style]
 
One of the things that gets kind of weird about this Makefile as compared to
traditional makefiles is the fact that I wrote it in a really wonky sort of
way where there are many dependencies for a single rule (This may change in
the future). 

Many of the rules in the makefile take the form of:

```make
.PHONY: all out

all : out

SAMPLES     := $(shell ls -1d samples-dir/*)
OUT_SAMPLES := $(patsubst %.in,%.out,$(SAMPLES))

out : $(OUT_SAMPLES)

runs/ARRAY-JOB-NAME/ARRAY-JOB-NAME.sh : $(SAMPLES)
	echo $^ | sed -r 's/([^ ]+?).in */script-to-run.sh \1.in -o \1.out\n/' > $(RUNFILES)/run-script-to-run.txt
	SLURM_Array -c $(RUNFILES)/run-script-to-run.txt \
	            --mail $(EMAIL) \
				-r runs/ARRAY-JOB-NAME \
				-l $(MODULE) \
				--hold \
				-w $(ROOT_DIR)

$(OUT_SAMPLES) : $(SAMPLES) runs/ARRAY-JOB-NAME/ARRAY-JOB-NAME.sh

```

Where the target is a phony target that generates one `OUT_SAMPLE` for every 
`SAMPLE` via `runs/ARRAY-JOB-NAME/ARRAY-JOB-NAME.sh`. This shell script as a
target is generated via [`SLURM_Array`][sarray] as it submits a 
[SLURM array job][arrayjob]. Is this an elegant solution? No. It's a bit
hamfisted, and I will probably change it in the future.

 
[make-style]: http://clarkgrubb.com/makefile-style-guide
[buffalo-make]: https://github.com/vsbuffalo/makefiles-in-bioinfo
[brianflow]: https://github.com/knausb/bam_processing
[HCC]: http://hcc.unl.edu/
[sarray]: https://github.com/zkamvar/SLURM_Array
[arrayjob]: https://slurm.schedmd.com/job_array.html

Required directories
--------------------

 - **mitochondria_genome/**: A gzipped fasta file such as one from here:  ftp://ftp.broadinstitute.org/pub/annotation/fungi/sclerotinia_sclerotiorum/broad/genomes/sclerotinia_sclerotiorum/
 - **reads/**: Paired-end genomic data, in `\*_[12].fq.gz` format
 
 
 Generated directories
 ---------------------
 
 - runfiles/: shell scripts for pre-processing submission scripts (kept in this directory for posterity)
 - *bt2-index/*: genome index files generated via `make index` 
 - *SAMS/*: mapped sam files generated via `make map`
 - *BAMS/*: filtered bam files
 - *GVCF/*: `*.g.vcf` and `*.vcf` files generated via GATK
 - *runs/*: std out and std err of runs
