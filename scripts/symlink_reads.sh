#!/bin/bash
#
# Script to symlink reads

if [ $# -eq 0 ]; then
	echo 
	echo "Symlink reads from BGI to the reads/ directory"
	echo
	echo "Usage:"
	echo "	$0 FILES..."
	echo
	echo "	- FILES... paths to pairs of read files from BGI in separate folders named SS.11.XX"
	echo "	  Where XX represents the sample number."
	echo
	exit
fi


FILES=$@

for i in $FILES; do
	LNK="$(echo $i | sed -r 's/^.+?(SS\.11\.[0-9]{2}).+?_([12].fq.gz)$/\1_\2/')"
	ln -s $i reads/$LNK
done
