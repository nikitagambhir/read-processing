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
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := $(addprefix $(ROOT_DIR)/, mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_supercontigs.fasta.gz)
PREFIX   := Ssc_mito 
READS    := $(shell ls -d reads/*_1.fq.gz | sed 's/_1.fq.gz//g')
RFILES   := $(addsuffix _1.fq.gz, $(READS))
IDX      := $(addprefix $(strip index/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)
MAPPED   := $(ROOT_DIR)/mapped
SAM      := $(addprefix $(MAPPED)/, $(READS))


# INDEX CREATION
# ==============
# 
# First, the command has to be built and stored in a file called 
# make-index.txt.
# 
# This file is then sent to the cluster with SLURM_Array
#
# The next rule is more of a formality, declaring that the bz2 files should be 
# created.
# -----------------------------------------------------------------------------
#
# Create run script
runfiles/make-index.txt : scripts/make-index.sh $(FASTA)
ifeq ($(wildcard index/.),)
	mkdir index
endif
	$^ $(addprefix $(ROOT_DIR)/index/, $(PREFIX)) 
#
# submit job
runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh : runfiles/make-index.txt
	SLURM_Array -c $(RUNFILES)/make-index.txt \
		--mail  $$EMAIL\
		-r runs/BOWTIE2-BUILD \
		-l bowtie/2.2 \
		-w $(ROOT_DIR)
#
# Formatlity: submit job should create the proper output files.
$(IDX) : runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh
#
# .PHONY target
index : $(IDX)
#
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# READ MAPPING
# ============
#
# This goes throught he script make-alignment.sh which creates the runscript
# for SLURM_Array
#
# SLURM_Array runs it, emails me, and we're all happy.
#
# Like above, the next rule is a formality accounting for the files produced.
# -----------------------------------------------------------------------------
#
# Create run script
runfiles/make-alignment.txt : scripts/make-alignment.sh $(RFILES) $(ROOT_DIR)
	$< $(addprefix $(ROOT_DIR)/index/, $(PREFIX)) $(MAPPED) $(addprefix $(ROOT_DIR)/, $(READS))
#
# submit job
runs/MAP-READS/MAP-READS.sh : runfiles/make-alignment.txt index
ifeq ($(wildcard $(MAPPED)/.),)
	mkdir $(MAPPED)
endif
	SLURM_Array -c $< \
	--mail  $$EMAIL\
	-r runs/MAP-READS \
	-l bowtie/2.2 \
	-w $(ROOT_DIR)
#
# Formality: the submit job should create these output files
$(SAM): runs/MAP-READS/MAP-READS.sh
# 
# .PHONY target map
map: $(SAM)
#
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
help:
	@echo $(IDX)

