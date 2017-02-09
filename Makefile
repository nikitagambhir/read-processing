.PHONY: all align index help 


ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) # cannot do $(shell pwd) since it executes on a different machine
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := $(addprefix $(ROOT_DIR)/, mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz)
PREFIX   := Ssc_mito 
READS    := $(shell ls -d reads/*_1.fq.gz)
IDX      := $(addprefix $(strip index/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)

index : $(IDX)
align: runfiles/make-alignment.txt

# Step 1: create the index for the mitochondrial genome
# 
# First, the command has to be built and stored in a file called make-index.txt.
# This file is then sent to the cluster with SLURM_Array
#
# The next rule is more of a formality, declaring that the bz2 files should be created.
runfiles/make-index.txt : scripts/make-index.sh $(FASTA)
	# Make output directory
ifeq (wildcard index/.),)
	mkdir index
endif
	# Create command file
	bash $^ $(addprefix $(ROOT_DIR)/index/, $(PREFIX)) 
	# Run command with SLURM_Array
	SLURM_Array -c $(RUNFILES)/make-index.txt --mail  $$EMAIL-r runs/BOWTIE2-BUILD -l bowtie/2.2 -w $(ROOT_DIR)

$(IDX) : runfiles/make-index.txt


runfiles/make-alignment.txt: $(READS) $(ROOT_DIR)
	printf " $(addsuffix \\n, $(addprefix $(ROOT_DIR)/, $(READS)))" > runfiles/make-alignment.txt




help:
	@echo $(IDX)
