#!/bin/bash

if [ $# -lt 3 ]; then
        echo
        echo "Align samples to indexed genome output."
        echo
        echo
        echo "This assumes that the data are paired-end reads."
        echo
        echo "Usage:"
        echo
        echo "  bash make-index.sh <bt2-idx> <dir> SAMPLES..."
        echo
        echo "  <bt2-idx>  - basename of reference (with path)"
        echo "  <dir>      - directory in which to place the SAM files"
	echo "	SAMPLES... - list of basenames for samples"
        echo
        exit
fi

IDX=$1
DIR=$2
# Slice the array from the 3rd position to the end
# http://stackoverflow.com/a/9057392
SAMPLES="${@:3}"
CMD="bowtie2 -x $IDX -S $DIR/"


printf "" > runfiles/make-alignment.txt

tmp=($SAMPLES)
tmp=${tmp[0]}"_1.fq.gz"
# Info Line: 
#          0	     1	          2	      3	      4	5	   6	       7	         8
# INSTRUMENT	RUN_NO	FLOWCELL_ID	LANE_NO	TILE_ID	X	Y	READ	FILTERED	CONTROL_NO
INFO=($(zcat $tmp | head -n 1 | sed 's/:/ /g'))

for s in $SAMPLES; do
	# Create named pipes because bowtie2 is stupid and can't handle gzipped files >:(
	MKFIFO="mkfifo \$TMPDIR/M1.fifo \$TMPDIR/M2.fifo"
	M1=$s"_1.fq.gz"
	M2=$s"_2.fq.gz"
	WFIFO="zcat $M1 > \$TMPDIR/M1.fifo & zcat $M2 > \$TMPDIR/M2.fifo &"

	# Setting arguments
	ARGS="-1 \$TMPDIR/M1.fifo -2 \$TMPDIR/M2.fifo --maxins 800 --fr"
	BASE=$(basename $s)
	
	# Setup read group information for GATK
	RGID="${INFO[2]}.${INFO[3]}"
	RG="--rg SM:$BASE PL:ILLUMINA LB:RUN.${INFO[1]}"

	printf "$MKFIFO; $WFIFO $CMD$BASE $ARGS $RGID $RG\n" >> runfiles/make-alignment.txt
done
