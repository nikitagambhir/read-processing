.PHONY: all 


TMP      := \$$TMPDIR
ROOT_DIR := $(shell echo $$WORK/Sclerotinia_mitochondria)
RUNFILES := $(ROOT_DIR)/runfiles
FASTA    := mitochondira_genome/sclerotinia_sclerotiorum_mitochondria_2_contigs.fasta.gz
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
