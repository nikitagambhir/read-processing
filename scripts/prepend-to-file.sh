#!/usr/bin/env bash


if [ $# -lt 2 ]; then
	echo
	echo "Prepend stdin to each line of file"
	echo
	echo "Usage:"
	echo
	echo "	cat 'stdin to prepend' | prepend-to-file.sh INFILE OUTFILE"
	echo
	exit
fi 


input=$(cat)
file=$(cat $1)
out=$2

printf "" > $out

for f in $file; do
	printf "$input $f\n" >> $out
done



 
