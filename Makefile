.PHONY: all align index help fasta


ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) # cannot do $(shell pwd) since it executes on a different machine
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz
PREFIX   := Ssc_mito 
READS    := $(shell ls -d reads/*_1.fq.gz)
IDX      := $(addprefix $(strip index/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)

# Step 1: create the index for the mitochondrial genome

runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh : $(FASTA)
	# Make output directory
ifeq (wildcard index/.),)
	mkdir index
endif
	# Create command file
	echo "zcat $< > \$$TMPDIR/tmp.fna; \
	bowtie2-build \
	--seed 20170207 \
	--threads 4 \
	\$$TMPDIR/tmp.fna \
	$(ROOT_DIR)/index/$(PREFIX)" > $(RUNFILES)/make-index.txt
	# Run command with SLURM_Array
	SLURM_Array -c $(RUNFILES)/make-index.txt --mail  $$EMAIL-r runs/BOWTIE2-BUILD -l bowtie/2.2 -w $(ROOT_DIR)

runfiles/make-alignment.txt: $(READS) $(ROOT_DIR)
	printf " $(addsuffix \\n, $(addprefix $(ROOT_DIR)/, $(READS)))" > runfiles/make-alignment.txt

$(IDX) : runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh

align: runfiles/make-alignment.txt
index: $(IDX)	
help:
	@echo $(IDX)
