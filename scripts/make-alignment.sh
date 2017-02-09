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
# http://stackoverflow.com/a/1336245
SAMPLES="${@:3}"
CMD="bowtie2 -x $IDX -S $DIR/"


printf "" > runfiles/make-alignment.txt

for s in $SAMPLES; do
	MKFIFO="mkfifo \$TMPDIR/M1.fifo \$TMPDIR/M2.fifo"
	M1=$s"_1.fq.gz"
	M2=$s"_2.fq.gz"
	WFIFO="zcat $M1 > \$TMPDIR/M1.fifo & zcat $M2 > \$TMPDIR/M2.fifo &"
	ARGS="-1 \$TMPDIR/M1.fifo -2 \$TMPDIR/M2.fifo --maxins 800 --fr"
	BASE=$(basename $s)
	printf "$MKFIFO; $WFIFO $CMD$BASE $ARGS\n"
done
