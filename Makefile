# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Makefile for analyzing mitochondrial genomes
#
# Author: Zhian N. Kamvar
# Licesne: MIT
#
# This makefile contains rules and recipes for mapping, filtering, and
# analyzing mitochondrial genomes of *Sclerotinia sclerotiorum* treated with
# four different fungicides to assess the impact of fungicide stress on genomic
# architecture.
#
# Since this makefile runs on a SLURM cluster, this is tailored specifically for
# the HCC cluster in UNL. For this to work, the script SLURM_Array must be in 
# your path. You can download it here: https://github.com/zkamvar/SLURM_Array
#
# The general pattern of this makefile is that each target takes two steps:
#
# 1. Run a bash script to collect the dependencies into separate lines of a
#    text file. Each line will be an identical command to run in parallel
#    across the cluster.
# 2. The text file is submitted to the cluster with SLURM_Array
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.PHONY: all index help map clean

all: index map

ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) 
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := runfiles # directory to place the files for running
FASTA    := mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_supercontigs.fasta.gz
PREFIX   := Ssc_mito # prefix for the bowtie2 index
READS    := $(shell ls -d reads/*_1.fq.gz | sed 's/_1.fq.gz//g')
RFILES   := $(addsuffix _1.fq.gz, $(READS))
IDX      := $(addprefix $(strip index/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)
MAPPED   := mapped
SAM      := $(patsubst reads/%.sam, $(MAPPED)/%.sam, $(addsuffix .sam, $(READS)))
SAM_VAL  := $(patsubst %.sam, %_stats.txt.gz, $(SAM))
BAM      := $(patsubst mapped/%.sam, bams/%_nsort, $(SAM))
FIXED    := $(patsubst %_nsort, %_fixed.bam, $(BAM))
BAM_VAL  := $(patsubst %_fixed.bam, %_fixed_stats.txt.gz, $(FIXED))

runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh : scripts/make-index.sh $(FASTA)
ifeq ($(wildcard index/.),)
	mkdir index
endif
	$^ $(addprefix index/, $(PREFIX)) 
	SLURM_Array -c $(RUNFILES)/make-index.txt \
		-r runs/BOWTIE2-BUILD \
		-l bowtie/2.2 \
		-w $(ROOT_DIR)

$(IDX) : scripts/make-index.sh $(FASTA) runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh
index : $(FASTA) $(IDX) 

runs/MAP-READS/MAP-READS.sh: scripts/make-alignment.sh $(RFILES) 
ifeq ($(wildcard $(MAPPED)/.),)
	mkdir $(MAPPED)
	mkdir bams
endif
	$< $(addprefix index/, $(PREFIX)) $(MAPPED) $(READS)
	SLURM_Array -c runfiles/make-alignment.txt \
		--mail  $$EMAIL\
		-r runs/MAP-READS \
		-l bowtie/2.2 \
		--hold \
		-w $(ROOT_DIR)

$(SAM) : scripts/make-alignment.sh $(RFILES) runs/MAP-READS/MAP-READS.sh
$(SAM_VAL): $(SAM) runs/VALIDATE-SAM/VALIDATE-SAM.sh
$(BAM) : $(SAM) runs/SAM-TO-BAM/SAM-TO-BAM.sh
$(FIXED) : $(BAM) runs/FIXMATE/FIXMATE.sh
$(BAM_VAL) : $(FIXED) runs/VALIDATE-BAM/VALIDATE-BAM.sh

map : $(IDX) $(SAM) $(SAM_VAL) $(BAM) $(FIXED) $(BAM_VAL)

runs/VALIDATE-SAM/VALIDATE-SAM.sh: $(SAM)
	echo $^ | \
	sed -r 's/'\
	'([^ ]+?).sam *'\
	'/'\
	'samtools stats \1.sam | gzip -c > \1_stats.txt.gz\n'\
	'/g' > runfiles/validate-sam.txt # end
	SLURM_Array -c runfiles/validate-sam.txt \
		--mail  $$EMAIL\
		-r runs/VALIDATE-SAM \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)

runs/SAM-TO-BAM/SAM-TO-BAM.sh: $(SAM)
	echo $^ | \
	sed -r 's/'\
	'$(MAPPED)([^ ]+?).sam *'\
	'/'\
	'samtools view -bSu $(MAPPED)\1.sam | '\
	'samtools sort -n -O bam -o bams\1_nsort -T bams\1_nsort_tmp\n'\
	'/g' > runfiles/sam-to-bam.txt # end
	SLURM_Array -c runfiles/sam-to-bam.txt \
		--mail  $$EMAIL\
		-r runs/SAM-TO-BAM \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)

runs/FIXMATE/FIXMATE.sh: $(BAM)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_nsort *'\
	'@'\
	'zcat $(FASTA) > $$TMPDIR/r.fa; '\
	'samtools fixmate -O bam \1_nsort /dev/stdout | '\
	'samtools sort -O bam -o - -T \1_csort_tmp | '\
	'samtools calmd -b - $$TMPDIR/r.fa > \1_fixed.bam\n'\
	'@g' > runfiles/fixmate.txt # end
	SLURM_Array -c runfiles/fixmate.txt \
		--mail  $$EMAIL\
		-r runs/FIXMATE \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)

runs/VALIDATE-BAM/VALIDATE-BAM.sh: $(FIXED)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_fixed.bam *'\
	'@'\
	'samtools stats \1_fixed.bam | '\
	'gzip -c > \1_fixed_stats.txt.gz\n'\
	'@g' > runfiles/validate-bam.txt # end
	SLURM_Array -c runfiles/validate-bam.txt \
		--mail  $$EMAIL\
		-r runs/VALIDATE-BAM \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)

help :
	@echo $(IDX)
	@echo $(SAM)
	@echo $(SAM_VAL)

runclean.%:
	$(RM) -rfv runs/$*

