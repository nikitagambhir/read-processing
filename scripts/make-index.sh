#!/bin/bash

if [ $# -lt 2 ]; then
	echo 
	echo "Create index for bowtie2"
	echo
	echo
	echo "All arguments for this script are the arguments for bowtie2-index"
	echo
	echo "Usage:"
	echo
	echo "	bash make-index.sh <input_reference.fasta> <index_prefix>" 
	echo
	echo "	<input_reference.fasta> - fasta or gzipped fasta file"
	echo "	<index_prefix> - prefix for the output files"
	echo
	exit
fi

INPUT_REF=$1
INDEX_PREFIX=$2
CMD="bowtie2-build --seed 99"


if [[ $INPUT_REF == *gz ]]; then
	WRITE_TMP="zcat $INPUT_REF > \$TMPDIR/tmp.fa"
	CMD="$WRITE_TMP; $CMD \$TMPDIR/tmp.fa $INDEX_PREFIX"
else 
	CMD="$CMD $INPUT_REF $INDEX_PREFIX"
fi

printf "$CMD\n" > runfiles/make-index.txt

