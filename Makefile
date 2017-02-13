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

.PHONY: all index help map

all: index map

ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) # cannot do $(shell pwd) since it executes on a different machine
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := runfiles
FASTA    := mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_supercontigs.fasta.gz
PREFIX   := Ssc_mito 
READS    := $(shell ls -d reads/*_1.fq.gz | sed 's/_1.fq.gz//g')
RFILES   := $(addsuffix _1.fq.gz, $(READS))
IDX      := $(addprefix $(strip index/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)
MAPPED   := mapped
SAM      := $(patsubst reads/%.sam, $(MAPPED)/%.sam, $(addsuffix .sam, $(READS)))
SAM_VAL  := $(patsubst %.sam, %_stats.txt.gz, $(SAM))
BAM      := $(patsubst mapped/%.sam, bams/%_nsort, $(SAM))
FIXED    := $(patsubst %_nsort, %_fixed.bam, $(BAM))

# INDEX CREATION
# ==============
# 
# First, the command has to be built and stored in a file called 
# make-index.txt. This file is then sent to the cluster with SLURM_Array
#
# The next rule declares that the bowtie2 index files $(IDX) should are
# dependent on the shell script produced from SLURM_Array, the generating
# script, and the genome FASTA files. If any of these files changes, the index
# files must be rebuilt. 
# -----------------------------------------------------------------------------
#
# Create run script
runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh : scripts/make-index.sh $(FASTA)
ifeq ($(wildcard index/.),)
	mkdir index
endif
	# All of the dependencies are in the correct order here
	$^ $(addprefix index/, $(PREFIX)) 
	SLURM_Array -c $(RUNFILES)/make-index.txt \
		-r runs/BOWTIE2-BUILD \
		-l bowtie/2.2 \
		-w $(ROOT_DIR)
#
# The important files should depend on the files that built it
# (Index Files) -> (submit script) -> (generator script) + (FASTA files)
$(IDX) : scripts/make-index.sh $(FASTA) runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh
#
# .PHONY target
index : $(FASTA) $(IDX) 
#
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# READ MAPPING
# ============
#
# This goes throught he script make-alignment.sh which creates the runscript
# for SLURM_Array
# -----------------------------------------------------------------------------
#
# Create run script
runs/MAP-READS/MAP-READS.sh: scripts/make-alignment.sh $(RFILES) 
ifeq ($(wildcard $(MAPPED)/.),)
	mkdir $(MAPPED)
	mkdir bams
endif
	# The script input here is the only one that's in order
	$< $(addprefix index/, $(PREFIX)) $(MAPPED) $(READS)
	SLURM_Array -c runfiles/make-alignment.txt \
		--mail  $$EMAIL\
		-r runs/MAP-READS \
		-l bowtie/2.2 \
		--hold \
		-w $(ROOT_DIR)
#
# SAM files need to depend on the index files
# (SAM files) -> (Index files)
# (SAM files) -> (submit script) -> (generator script) + (read files)
$(SAM) : scripts/make-alignment.sh $(RFILES) runs/MAP-READS/MAP-READS.sh
$(SAM_VAL): $(SAM) runs/VALIDATE-SAM/VALIDATE-SAM.sh
$(BAM) : $(SAM) runs/SAM-TO-BAM/SAM-TO-BAM.sh
$(FIXED) : $(BAM) runs/FIXMATE/FIXMATE.sh
# 
# .PHONY target map
map : $(IDX) $(SAM) $(SAM_VAL) $(BAM) $(FIXED)
#
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

runs/VALIDATE-SAM/VALIDATE-SAM.sh: $(SAM)
	echo $^ | sed -r 's/([^ ]+?).sam /samtools stats \1.sam | gzip -c > \1_stats.txt.gz\n/g' > runfiles/validate-sam.txt
	SLURM_Array -c runfiles/validate-sam.txt \
		--mail  $$EMAIL\
		-r runs/VALIDATE-SAM \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)

runs/SAM-TO-BAM/SAM-TO-BAM.sh: $(SAM)
	echo $^ | sed -r 's/mapped([^ ]+?).sam */samtools view -bSu mapped\1.sam | samtools sort -n -O bam -o bams\1_nsort -T bams\1_nsort_tmp\n/g' > runfiles/sam-to-bam.txt
	SLURM_Array -c runfiles/sam-to-bam.txt \
		--mail  $$EMAIL\
		-r runs/SAM-TO-BAM \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)
# $SAMT fixmate -O bam bams/${arr[0]}_nsort /dev/stdout | $SAMT sort -O bam -o - -T bams/${arr[0]}_csort_tmp | $SAMT calmd -b - $REF > bams/${arr[0]}_fixed.bam

runs/FIXMATE/FIXMATE.sh: $(BAM)
	echo $^ | sed -r 's@([^ ]+?)_nsort *@samtools fixmate -O bam \1_nsort \/dev\/stdout | samtools sort -O bam -o - -T \1_csort_tmp | samtools calmd -b - $(FASTA) > \1_fixed.bam\n@g' > runfiles/fixmate.txt
	SLURM_Array -c runfiles/fixmate.txt \
		--mail  $$EMAIL\
		-r runs/FIXMATE \
		-l samtools/1.3 \
		--hold \
		-w $(ROOT_DIR)
help :
	@echo $(IDX)
	@echo $(SAM)
	@echo $(SAM_VAL)

