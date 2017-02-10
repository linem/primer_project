

import sys
import itertools

#primer = sys.argv[1]
seq_file = sys.argv[1]
#seq = 'GGTTACCTTGTTACGACTT'
primer = 'CATTAGCGGCCAGGATGCT'
mismatch = 1




def possible_primers(primer_seq):
	'''Input is a list of possible bases at each position in the primer.
	   Returns a list of all possible non-degenerate primers'''
	sequences = ['']
	for i in range(len(primer_seq)):
		sequences = [a+b for a,b in itertools.product(sequences, primer_seq[i])]
	return(sequences)


def reverse_complement(sequences):
	'''Creates reverse complement of the primers in the input list'''
	rev_primers = []
	bases = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	for seq in sequences:
		comp = ''.join([bases[base] for base in seq])
		rev_primers.append(comp[::-1])
	return(rev_primers)


def search(primers, seq):
	'''Inputs a list of primers and a sequence to search for the primers.
	   Prints sequence id, position and primer for each hit'''
	positions = {}
	for prime in primers:
		if prime in seq:
			positions[seq.index(prime)+1] = prime
	if positions:
		print_positions(positions)
		return(True)


def fuzzy_search(primers, seq):
	'''Inputs a list of primers and a sequence to search for the primers, while allowing for
	   mismatches. Prints sequence id, position and primer for each hit'''
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
					if mism < mismatch and j == len(prime)-1 and (start+1) not in positions:
						positions[start+1] = prime
					elif mism > mismatch:
						break
	if positions:
		print_positions(positions)
		return(True)


def print_positions(positions):
	'''Print function for id line, position of hit and specific primer'''
	for key in sorted(positions):
		print('sequence id: {}, position: {}, primer: {}'.format(save_id, key, positions[key]))


# Dictionary to specify bases in degenerate primer
bases = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': ['A', 'T'],
		 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'],
		 'Y': ['C', 'T'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
		 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']}


# Specify not-degenerate primer
primer_seq = [bases[base] for base in primer]
# Create lists of forward and reverse primers
forw_primers = possible_primers(primer_seq)
rev_primers = reverse_complement(forw_primers)
primers = list(set(forw_primers + rev_primers))



# Prints possible primers
print('Forward primers: {}'.format(' '.join(forw_primers)))
print('Reverse primers: {}'.format(' '.join(rev_primers)))




# Run through sequences in fasta file and search for primers
with open(seq_file, 'r') as inf:
	seq = ''
	found = False
	for line in inf:
		if line[0] == '>':
			if seq:					     # when seq is complete
				found = search(primers, seq)     # search for primers
			save_id = line.rstrip()[1:]
			seq = ''				     # reset seq
		else:
			seq += line.rstrip()         # extend sequence (spanning multiple lines in file)
	found = search(primers, seq)       # searches final sequence
	if not found:
		print('Primer not in sequence')











