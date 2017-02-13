

import sys
import os
import argparse
from itertools import product


usage = '''This program will specify all possible primer sequences
		   (forward and reverse complement) for a given degenerate primer.
		   If a fasta file is specified as inputfile, the program will search
		   for primer matches and return sequence ids and position.'''

parser = argparse.ArgumentParser(description = usage)

parser.add_argument(
	'-v',
	action = 'version',
	version = '%(prog)s 1.0'
	)
parser.add_argument(
	'-p',
	dest = 'primer',
	metavar = 'PRIMER',
	required = True,
	help = 'sequence of degenerate primer'
	)
parser.add_argument(
	'-m',
	dest = 'mismatch',
	metavar = 'MISMATCHES',
	type = int,
	help = 'if mismatches allowed, specify number of mismatches'
	)
parser.add_argument(
	'-i',
	dest = 'infile',
	metavar = 'INFILE',
	help = 'fasta file to search for matching primer'
	)
parser.add_argument(
	'-o',
	dest = 'outfile',
	metavar = 'OUTFILE',
	type = argparse.FileType('w'),
	default = sys.stdout,
	help = 'output file, else write to STDOUT'
	)

args = parser.parse_args()


# File tests
if args.infile:
	if not os.path.isfile(args.infile):
		print('ERROR: Inputfile {} does not exist'.format(args.infile))
if args.outfile:
	if os.path.isfile(args.outfile.name):
		print('ERROR: Outputfile {} already exists'.format(args.outfile.name))
		overwrite = input('Do you want to overwrite the outputfile? (y/[n]): ')
		if not (overwrite == 'yes' or overwrite == 'y'):
			sys.exit()


def possible_primers(primer_seq):
	'''Input is a list of possible bases at each position in the primer.
	   Returns a list of all possible non-degenerate primers'''
	primers = ['']
	for i in range(len(primer_seq)):
		primers = [a+b for a,b in product(primers, primer_seq[i])]
	return(primers)


def reverse_complement(primers):
	'''Creates reverse complement of the primers in the input list'''
	rev_primers = []
	comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	for primer in primers:
		comp = ''.join([comp_dict[nuc] for nuc in primer])
		rev_primers.append(comp[::-1])
	return(rev_primers)


def search(primers, seq):
	'''Inputs a list of primers and a sequence to search for the primers and
	   prints sequence id, position and primer for each hit. If mismatches is 
	   specifies the function performs a fuzzy search.'''
	if not args.mismatch:
		positions = {}
		for primer in primers:
			if primer in seq:
				positions[seq.index(primer)+1] = primer
	else:
		positions = fuzzy_search(primers, seq)
	if positions:
		print_positions(positions)
		return(True)


def fuzzy_search(primers, seq):
	positions = {}
	for prime in primers:
			start = -1
			i = 0
			while i < len(seq)-1:
				mism = 0
				start += 1
				for i, j in zip(range(start,len(seq)), range(len(prime))):
					if seq[i] == prime[j]:
						if j == len(prime)-1 and (start+1) not in positions:
							positions[start+1] = prime
					else:
						mism += 1
						if mism < args.mismatch and j == len(prime)-1 and (start+1) not in positions:
							positions[start+1] = prime
						elif mism > args.mismatch:
							break
	return(positions)


def print_positions(positions):
	'''Print function for id line, position of hit and specific primer'''
	for key in sorted(positions):
		print('sequence id: {}, position: {}, primer: {}'.format(save_id, key, positions[key]), file = args.outfile)


# Dictionary to specify bases in degenerate primer
degen_dict = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': ['A', 'T'],
		 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'],
		 'Y': ['C', 'T'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
		 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']}


# Specify not-degenerate primer
primer_seq = [degen_dict[nuc] for nuc in args.primer.upper()]
# Create lists of forward and reverse primers
forw_primers = possible_primers(primer_seq)
rev_primers = reverse_complement(forw_primers)
primers = list(set(forw_primers + rev_primers))



# Prints possible primers
print('Forward primers: {}'.format(' '.join(forw_primers)))
print('Reverse primers: {}'.format(' '.join(rev_primers)))




# Run through sequences in fasta file and search for primers
if args.infile:
	with open(args.infile, 'r') as inf:
		seq = ''
		match = False
		for line in inf:
			if line[0] == '>':
				if seq:					     # when seq is complete
					match = search(primers, seq)     # search for primers
				save_id = line.rstrip()[1:]
				seq = ''				     # reset seq
			else:
				seq += line.rstrip().upper()   # extend sequence (spanning multiple lines in file)
		match = search(primers, seq)       # searches final sequence
		if not match:
			print('Primer not in sequence')











