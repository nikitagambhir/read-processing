#!/usr/bin/env python3.5


import sys
import io
import re
import argparse

def make_chromosomes(fasta):
	chrom_list = list()
	count = 0
	for line in fasta:
		line_stripped = line.strip()
		is_header = re.match(r'^>', line_stripped)
		if is_header:
			if count > 0:
				chrom_list.append(Chromosome(header, count))
				count = 0
			header = line_stripped[1:]
			next
		else:
			count = count + len(line_stripped)
	chrom_list.append(Chromosome(header, count))
	return(chrom_list)
	


def parse_input():
	parser = argparse.ArgumentParser(description = "Reads in a FASTA formatted file and returns a list of windows per chromosome for use with the GATK --intervals option")
	parser.add_argument("-f", "--fasta", required = False, dest = "fasta", default = "-", help = "a fasta formatted file (default: stdin)")
	parser.add_argument("-w", "--window", required = False, dest = "window", default = 100, type = int, help = "window size in bp (default: 100)")
	parser.add_argument("-o", "--outfile", required = False, dest = "outfile", default = None, help = "a file to print the windows (default: stdout)")
	args = parser.parse_args()
	
	fasta = None
	if args.fasta == "-":
		if sys.stdin.isatty():
			parser.print_help()
			quit()
		fasta = sys.stdin
	else:
		fasta = io.open(args.fasta, "r")
	args.fasta = fasta
	outfile = None
	if args.outfile is not None:
		outfile = io.open(args.outfile, "w")
	args.outfile = outfile
	return(args)


class Chromosome:
	def __init__(self, name, length):
		self.name = name
		self.length = length
		self.windows = list()
	def get_name(self):
		return(self.name)
	def get_length(self):
		return(self.length)
	def get_windows(self):
		return(self.windows)
	def print_windows(self, fh = None):
		if len(self.windows) == 0:
			return
		for w in self.windows:
			if fh is None:
				print(w)
			else:
				fh.write(w + "\n")
		if fh is not None:
			fh.close()
		return
	def make_windows(self, n):
		self.windows = list()
		# create a list of of coordinates in the GATK format
		# chromosome_1:1-25
		# chromosome_1:26-50
		windows = list(range(0, self.length, n))
		for w in range(1, len(windows)):
			start = str(windows[w - 1] + 1)
			end = str(windows[w])
			winstring = self.name + ":" + start + "-" + end
			self.windows.append(winstring)
		start = str(windows[len(windows) - 1] + 1)
		end = str(self.length)
		winstring = self.name + ":" + start + "-" + end
		self.windows.append(winstring)
		return(self.windows)

args = parse_input()
chroms = make_chromosomes(args.fasta)
for chrom in chroms:
	chrom.make_windows(args.window)
	chrom.print_windows(args.outfile)
