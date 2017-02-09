.PHONY: all help


TMP      := \$$TMPDIR # $TMPDIR is only valid on running jobs
ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria) # cannot do $(shell pwd) since it executes on a different machine
ROOT_DIR := $(strip $(ROOT_DIR))
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := mitochondira_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz
PREFIX   := sclerotinia_sclerotiorum_mitochondria_2_contigs
READS    := $(shell ls -d reads/*_1.fq.gz)

# Step 1: create the index for the mitochondrial genome

index: mitochondria_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz
	# Make output directory
	mkdir index
	# Create command file
	echo "zcat $< > $(TMP)/tmp.fna; \
	bowtie2-build \
	--seed 20170207 \
	--threads 4 \
	$(TMP)/tmp.fna \
	$(ROOT_DIR)/index/$(<F)" > $(RUNFILES)/make-index.txt
	# Run command with SLURM_Array
	SLURM_Array -c $(RUNFILES)/make-index.txt --mail  $$EMAIL-r runs/BOWTIE2-BUILD -l bowtie/2.2 -w $(ROOT_DIR)

help: $(READS) $(ROOT_DIR)
	printf "" > runfiles/make-alignment.txt
	printf " $(addsuffix \\n, $(addprefix $(ROOT_DIR)/, $(READS)))" >> runfiles/make-alignment.txt
	head runfiles/make-alignment.txt
