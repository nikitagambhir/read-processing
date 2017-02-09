.PHONY: all index help map

all: index map

ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) # cannot do $(shell pwd) since it executes on a different machine
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := $(addprefix $(ROOT_DIR)/, mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz)
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
	bash $^ $(addprefix $(ROOT_DIR)/index/, $(PREFIX)) 
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
# -----------------------------------------------------------------------------
runfiles/make-alignment.txt : scripts/make-alignment.sh $(RFILES) $(ROOT_DIR)
	$< $(addprefix $(ROOT_DIR)/index/, $(PREFIX)) $(MAPPED) $(addprefix $(ROOT_DIR)/, $(READS))

runs/MAP-READS/MAP-READS.sh : runfiles/make-alignment.txt index
ifeq ($(wildcard $(MAPPED)/.),)
	mkdir $(MAPPED)
endif
	SLURM_Array -c $< \
	--mail  $$EMAIL\
	-r runs/MAP-READS \
	-l bowtie/2.2 \
	-w $(ROOT_DIR)

$(SAM): runs/MAP-READS/MAP-READS.sh

map: $(SAM)

help:
	@echo $(IDX)
