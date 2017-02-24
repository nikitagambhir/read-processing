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


all: $(FASTA) $(RFILES) vcf 


EMAIL    := $$EMAIL # Set this environmental variable or change it here
ROOT_DIR := $(shell echo $$WORK/$$(basename $$(pwd))) 
ROOT_DIR := $(strip $(ROOT_DIR))
TMP      := \$$TMPDIR
IDX_DIR  := bt2-index
BOWTIE   := bowtie/2.2
SAMTOOLS := samtools/1.3
PICARD   := picard/1.1
GATK     := gatk/3.4
gatk     := \$$GATK
PIC      := \$$PICARD
TRIM_DIR := TRIM
SAM_DIR  := SAMS
BAM_DIR  := BAMS
GVCF_DIR := GVCF
REF_DIR  := REF
RUNFILES := runfiles
FAST_DIR := mitochondria_genome
FASTA    := $(FAST_DIR)/sclerotinia_sclerotiorum_mitochondria_2_supercontigs.fasta.gz
REF_FNA  := $(patsubst $(FAST_DIR)/%.fasta.gz,$(REF_DIR)/%.fasta, $(FASTA))
REF_IDX  := $(patsubst %.fasta,%.dict, $(REF_FNA))
INTERVALS:= $(patsubst %.fasta,%.intervals.txt, $(REF_FNA))
PREFIX   := Ssc_mito # prefix for the bowtie2 index
READS    := $(shell ls -d reads/*_1.fq.gz | sed 's/_1.fq.gz//g')
RFILES   := $(addsuffix _1.fq.gz, $(READS))
IDX      := $(addprefix $(strip $(IDX_DIR)/$(PREFIX)), .1.bz2 .2.bz2 .3.bz2 .4.bz2 .rev.1.bz2 .rev.2.bz2)
SAM      := $(patsubst reads/%.sam, $(SAM_DIR)/%.sam, $(addsuffix .sam, $(READS)))
SAM_VAL  := $(patsubst %.sam, %_stats.txt.gz, $(SAM))
BAM      := $(patsubst $(SAM_DIR)/%.sam, $(BAM_DIR)/%_nsort, $(SAM))
FIXED    := $(patsubst %_nsort, %_fixed.bam, $(BAM))
DUPMRK   := $(patsubst %_nsort, %_dupmrk.bam, $(BAM))
GVCF     := $(patsubst reads/%,$(GVCF_DIR)/%.g.vcf.gz, $(READS))
DUP_VAL  := $(patsubst %_nsort, %_dupmrk_stats.txt.gz, $(BAM))
BAM_VAL  := $(patsubst %_fixed.bam, %_fixed_stats.txt.gz, $(FIXED))
VCF      := $(GVCF_DIR)/res.vcf.gz

joiner = reads/$(1)_1.fq.gz,\
	reads/$(1)_2.fq.gz,\
	$(SAM_DIR)/$(1).sam,\
	$(SAM_DIR)/$(1)_stats.txt.gz,\
	$(BAM_DIR)/$(1)_nsort,\
	$(BAM_DIR)/$(1)_fixed.bam,\
	$(BAM_DIR)/$(1)_fixed_stats.txt.gz,\
	$(BAM_DIR)/$(1)_dupmrk.bam,\
	$(BAM_DIR)/$(1)_dupmrk_stats.txt.gz,\
	$(GVCF_DIR)/$(1).g.vcf.gz\\n

MANIFEST := $(foreach x,$(patsubst reads/%,%, $(READS)),$(call joiner,$(x)))


$(RUNFILES) $(IDX_DIR) $(SAM_DIR) $(BAM_DIR) $(REF_DIR) $(GVCF_DIR) $(TRIM_DIR):
	-mkdir $@
index : $(FASTA) $(REF_FNA) $(INTERVALS) $(IDX) 
map : index $(SAM) $(SAM_VAL) 
bam : map $(BAM) $(FIXED) $(BAM_VAL) 
dup : bam $(DUPMRK) $(DUP_VAL)
vcf : dup $(REF_IDX) $(GVCF) $(VCF)
concat : runs/CONCAT-VCF/CONCAT-VCF.sh

# Unzips the reference genome
$(REF_DIR)/%.fasta : $(FAST_DIR)/%.fasta.gz | $(REF_DIR) $(RUNFILES)
	zcat $^ | sed -r 's/[ ,]+/_/g' > $@
	
# Creates intervals for the final step. Edit the -w parameter to change.
$(REF_DIR)/%.intervals.txt : $(REF_DIR)/%.fasta
	./scripts/make-GATK-intervals.py -f $< -w 10000 -o $@

runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh : scripts/make-index.sh $(REF_FNA) | $(IDX_DIR) $(RUNFILES)
	$^ $(addprefix $(IDX_DIR)/, $(PREFIX)) 
	SLURM_Array -c $(RUNFILES)/make-index.txt \
		-r runs/BOWTIE2-BUILD \
		-l $(BOWTIE) \
		-w $(ROOT_DIR)

$(IDX) : scripts/make-index.sh $(FASTA) runs/BOWTIE2-BUILD/BOWTIE2-BUILD.sh

runs/TRIM-READS/TRIM-READS.sh: $(RFILES) | $(TRIM_DIR)
	echo $(READS) | \
	sed -r 's@'\
	'reads/([^ ]+?) *'\
	'@'\
	'trimmomatic PE -phred33 reads/\1_1.fq.gz reads/\1_2.fq.gz '\
	'-baseout $(TRIM_DIR)/\1.fq.gz '\
	'ILLUMINACLIP:/util/opt/anaconda/2.0/envs/trimmomatic-0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36\n'\
	'@g' > $(RUNFILES)/trimmomatic-run.txt #end
	SLURM_Array -c $(RUNFILES)/trimmomatic-run.txt \
		-r runs/TRIM-READS \
		-l trimmomatic/0.36 \
		--hold \
		-w $(ROOT_DIR)

runs/MAP-READS/MAP-READS.sh: scripts/make-alignment.sh $(RFILES) | $(SAM_DIR) 
	$< $(addprefix $(IDX_DIR)/, $(PREFIX)) $(SAM_DIR) $(READS)
	SLURM_Array -c $(RUNFILES)/make-alignment.txt \
		--mail $(EMAIL) \
		-r runs/MAP-READS \
		-l $(BOWTIE) \
		--hold \
		-w $(ROOT_DIR)

$(SAM) : scripts/make-alignment.sh $(RFILES) runs/MAP-READS/MAP-READS.sh


runs/VALIDATE-SAM/VALIDATE-SAM.sh: $(SAM) | $(SAM_DIR) 
	echo $^ | \
	sed -r 's/'\
	'($(SAM_DIR)[^ ]+?).sam *'\
	'/'\
	'samtools stats \1.sam | gzip -c > \1_stats.txt.gz\n'\
	'/g' > $(RUNFILES)/validate-sam.txt # end
	SLURM_Array -c $(RUNFILES)/validate-sam.txt \
		--mail $(EMAIL) \
		-r runs/VALIDATE-SAM \
		-l $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)

$(SAM_VAL): $(SAM) runs/VALIDATE-SAM/VALIDATE-SAM.sh

# SAMTOOLS SPECIFICATIONS
#
# http://www.htslib.org/doc/
# view
# # -b       output BAM
# # -S       ignored (input format is auto-detected)
# # -u       uncompressed BAM output (implies -b)
#
# sort
# # -n         Sort by read name
# # -o FILE    output file name [stdout]
# # -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or
# # -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)
#
# calmd
# # -u         uncompressed BAM output (for piping)
runs/SAM-TO-BAM/SAM-TO-BAM.sh: $(SAM) | $(BAM_DIR)
	echo $^ | \
	sed -r 's/'\
	'$(SAM_DIR)([^ ]+?).sam *'\
	'/'\
	'samtools view -bSu $(SAM_DIR)\1.sam | '\
	'samtools sort -n -O bam -o $(BAM_DIR)\1_nsort -T $(BAM_DIR)\1_nsort_tmp\n'\
	'/g' > $(RUNFILES)/sam-to-bam.txt # end
	SLURM_Array -c $(RUNFILES)/sam-to-bam.txt \
		--mail $(EMAIL) \
		-r runs/SAM-TO-BAM \
		-l $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)

$(BAM) : $(SAM) runs/SAM-TO-BAM/SAM-TO-BAM.sh
# Fix mate information and add the MD tag.
# # http://samtools.github.io/hts-specs/
# # MD = String for mismatching positions
# # NM = Edit distance to the reference
runs/FIXMATE/FIXMATE.sh: $(BAM) | $(BAM_DIR)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_nsort *'\
	'@'\
	'samtools fixmate -O bam \1_nsort /dev/stdout | '\
	'samtools sort -O bam -o - -T \1_csort_tmp | '\
	'samtools calmd -b - $(REF_FNA) > \1_fixed.bam\n'\
	'@g' > $(RUNFILES)/fixmate.txt # end
	SLURM_Array -c $(RUNFILES)/fixmate.txt \
		--mail $(EMAIL) \
		-r runs/FIXMATE \
		-l $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)

$(FIXED) : $(BAM) runs/FIXMATE/FIXMATE.sh

runs/VALIDATE-BAM/VALIDATE-BAM.sh: $(FIXED) | $(BAM_DIR)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_fixed.bam *'\
	'@'\
	'samtools stats \1_fixed.bam | '\
	'gzip -c > \1_fixed_stats.txt.gz\n'\
	'@g' > $(RUNFILES)/validate-bam.txt # end
	SLURM_Array -c $(RUNFILES)/validate-bam.txt \
		--mail $(EMAIL) \
		-r runs/VALIDATE-BAM \
		-l $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)

$(BAM_VAL) : $(FIXED) runs/VALIDATE-BAM/VALIDATE-BAM.sh

runs/MARK-DUPS/MARK-DUPS.sh: $(FIXED)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_fixed.bam *'\
	'@'\
	'java -Djava.io.tmpdir=$(TMP) '\
	'-jar $(PIC) MarkDuplicates '\
	'I=\1_fixed.bam '\
	'O=\1_dupmrk.bam '\
	'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '\
	'ASSUME_SORTED=true '\
	'M=\1_marked_dup_metrics.txt; '\
	'samtools index \1_dupmrk.bam\n'\
	'@g' > $(RUNFILES)/mark-dups.txt # end
	SLURM_Array -c $(RUNFILES)/mark-dups.txt \
		--mail $(EMAIL) \
		-r runs/MARK-DUPS \
		-l $(PICARD) $(SAMTOOLS) \
		--hold \
		-m 25g \
		-w $(ROOT_DIR)

$(DUPMRK) : $(FIXED) runs/MARK-DUPS/MARK-DUPS.sh

runs/VALIDATE-DUPS/VALIDATE-DUPS.sh: $(DUPMRK)
	echo $^ | \
	sed -r 's@'\
	'([^ ]+?)_dupmrk.bam *'\
	'@'\
	'samtools stats \1_dupmrk.bam | '\
	'gzip -c \1_dupmrk_stats.txt.gz\n'\
	'@g' > $(RUNFILES)/validate-dups.txt # end
	SLURM_Array -c $(RUNFILES)/validate-dups.txt \
		--mail $(EMAIL) \
		-r runs/VALIDATE-DUPS \
		-l $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)	

$(DUP_VAL): $(DUPMRK) runs/VALIDATE-DUPS/VALIDATE-DUPS.sh


# https://www.broadinstitute.org/gatk/documentation/article?id=3893
# # https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
# # https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php
#
# # Note that Haplotypecaller requires an indexed bam.
# # If yours is not, use SAMtools.
#
# # If you're dealing with legacy data you may encounter legacy quality encodings.
# # If you encounter this use:
# #
# # --fix_misencoded_quality_scores
# #
# # In your GATK call. (But only on the offending libraries.)
# # https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_engine_CommandLineGATK.php#--fix_misencoded_quality_scores
# # https://en.wikipedia.org/wiki/FASTQ_format
#
#
# CMD="$JAVA -Djava.io.tmpdir=/data/ -jar $GATK \
#   -T HaplotypeCaller \
#   -R $REF \
#   --emitRefConfidence GVCF \
#   -ploidy 2 \
#   -I bams/${arr[0]}_dupmrk.bam \
#   -o gvcf/${arr[0]}_2n.g.vcf.gz"
#
#
# CMD="$JAVA -Djava.io.tmpdir=/data/ -jar $GATK \
#    -T HaplotypeCaller \
#    -R $REF \
#    --emitRefConfidence GVCF \
#    -ploidy 3 \
#    -I bams/${arr[0]}_dupmrk.bam \
#    -o gvcf/${arr[0]}_3n.g.vcf.gz"
#

runs/MAKE-GATK-REF/MAKE-GATK-REF.sh: $(REF_FNA) 
	echo $^ | \
	sed -r 's@'\
	'^(.+?).fasta'\
	'@'\
	'java -jar $(PIC) CreateSequenceDictionary '\
	'R=\1.fasta '\
	'O=\1.dict; '\
	'samtools faidx \1.fasta'\
	'@g' > $(RUNFILES)/make-gatk-ref.txt # end
	SLURM_Array -c $(RUNFILES)/make-gatk-ref.txt \
		--mail $(EMAIL) \
		-r runs/MAKE-GATK-REF \
		-l $(PICARD) $(SAMTOOLS) \
		--hold \
		-w $(ROOT_DIR)

$(REF_IDX): $(FASTA) runs/MAKE-GATK-REF/MAKE-GATK-REF.sh

# Pain points:
# 
# GATK is very picky as far as paths go. If it sees a relative
# path, it will use pwd. On this SLURM system, this results in
# a path that's not accessible.
#
# GATK assumes that you named your dict file with the basename
# of your file and not just appended dict on the end.
#
runs/MAKE-GVCF/MAKE-GVCF.sh: $(DUPMRK) | $(GVCF_DIR)
	echo $^ | \
	sed -r 's@'\
	'$(BAM_DIR)/([^ ]+?)_dupmrk.bam *'\
	'@'\
	'java -Djava.io.tmpdir=$(TMP) -jar $(gatk) '\
	'-T HaplotypeCaller '\
	'-R $(ROOT_DIR)/$(REF_FNA) '\
	'--emitRefConfidence GVCF '\
	'-ploidy 1 '\
	'-I $(BAM_DIR)/\1_dupmrk.bam '\
	'-o $(GVCF_DIR)/\1.g.vcf.gz\n'\
	'@g' > $(RUNFILES)/make-gvcf.txt
	SLURM_Array -c $(RUNFILES)/make-gvcf.txt \
		--mail $(EMAIL) \
		-r runs/MAKE-GVCF \
		-l $(GATK) \
		--hold \
		-m 25g \
		-w $(ROOT_DIR)

$(GVCF) : $(DUPMRK) runs/MAKE-GVCF/MAKE-GVCF.sh


# 
# Note for this step, memory matters more than the number of cores.
#
# For example, here I'm using 50g of memory by setting the -Xmx
# and the -m flag in the SLURM_Array command. Notice that for the
# Xmx flag, the number of corse must butt up against the flag. This
# is a Java thing.
#
# I also set the number of threads with -nt and -P flags, respectively
#
runs/MAKE-VCF/MAKE-VCF.sh: $(GVCF)
	printf "java -Xmx100g -Djava.io.tmpdir=$(TMP) "\
	"-jar $(gatk) "\
	"-nt 6 "\
	"-T GenotypeGVCFs "\
	"-R $(ROOT_DIR)/$(REF_FNA) "\
	"$(addprefix -V , $^) "\
	"-o $(GVCF_DIR)/res.\$$SLURM_ARRAY_TASK_ID.vcf.gz --intervals" | \
	./scripts/prepend-to-file.sh $(INTERVALS) $(RUNFILES)/make-vcf.txt
	SLURM_Array -c $(RUNFILES)/make-vcf.txt \
		--mail $(EMAIL) \
		-r runs/MAKE-VCF \
		-l $(GATK) \
		--hold \
		-m 100g \
		-t 24:00:00 \
		-P 6 \
		-w $(ROOT_DIR)

runs/CONCAT-VCF/CONCAT-VCF.sh: 
	echo 'vcf-concat $(shell ls $(GVCF_DIR)/res.*.vcf.gz | sort -t'.' -n -k2)'\
	' | gzip -c > $(GVCF_DIR)/res.vcf.gz' > \
	$(RUNFILES)/merge-vcf.txt
	SLURM_Array -c $(RUNFILES)/merge-vcf.txt \
		-r runs/CONCAT-VCF \
		-l vcftools/0.1 \
		-w $(ROOT_DIR)

$(VCF): $(GVCF) runs/MAKE-VCF/MAKE-VCF.sh

help :
	@echo
	@echo "COMMANDS"
	@echo "============"
	@echo
	@echo "all	almost all -- Make res.n.vcf.gz files"
	@echo "concat	concatenate res.n.vcf.gz files (to run after all)"
	@echo "index	generate the bowtie2 index"
	@echo "map	map reads and validate the SAM files"
	@echo "bam	convert sam to bam and filter"
	@echo "dup	deduplicate bam files and validate"
	@echo "vcf	create g.vcf and vcf files (this is the longest step)" 
	@echo
	@echo "help	show this message" 
	@echo "burn	REMOVE ALL GENERATED FILES"
	@echo "manifest	create a manifest of all generated files per read"
	@echo "runclean.JOB_NAME	clean runfiles"
	@echo
	@echo "PARAMETERS"
	@echo "============"
	@echo
	@echo "EMAIL     : " $(EMAIL)
	@echo "ROOT DIR  : " $(ROOT_DIR)
	@echo "TEMP DIR  : " $(TMP)
	@echo "INDEX DIR : " $(IDX_DIR)
	@echo "PREFIX    : " $(PREFIX)
	@echo "SAM DIR   : " $(SAM_DIR)
	@echo "BAM DIR   : " $(BAM_DIR)
	@echo "GENOME    : " $(FASTA)
	@echo "RUNFILES  : " $(RUNFILES)
	@echo "READS     : " $(READS)
	@echo

manifest :
	printf "READ1,"\
	"READ2,"\
	"SAM,"\
	"SAM VALIDATION,"\
	"SORTED BAM,"\
	"FIXED BAM,"\
	"FIXED BAM VALIDATION,"\
	"MARKED DUPLICATES BAM,"\
	"MARKED DUPLICATES BAM VALIDATION,"\
	"GVCF FILE\n" > manifest.csv
	printf "$(MANIFEST)" >> manifest.csv

runclean.%:
	$(RM) -r runs/$*

burn:
	$(RM) -r $(IDX_DIR) $(SAM_DIR) $(BAM_DIR) $(GVCF_DIR) runs $(REF_DIR) $(RUNFILES)



.PHONY: all index help map bam dup vcf clean burn manifest
