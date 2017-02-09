

import sys
import itertools

#primer = sys.argv[1]
seq_file = sys.argv[1]
#primer = 'GGTTACCTTGTTACGACTT'
primer = 'GWAS'

bases = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': ['A', 'T'],
		 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'],
		 'Y': ['C', 'T'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
		 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']}



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
	positions = set()
	for prime in primers:
		if prime in seq:
			positions.add(seq.index(prime)+1)
	for p in sorted(positions):
		print('sequence id: {}, position: {}, primer: {}'.format(save_id, p, seq[p-1:p+len(primer)-1]))
	if positions:
		return(True)


# Specify not-degenerate primer
primer_seq = [bases[base] for base in primer]
# Create list of forward and reverse primers
sequences = possible_primers(primer_seq)
rev_primers = reverse_complement(sequences)


# Prints possible primers
print('Forward primers: {}'.format(' '.join(sequences)))
print('Reverse primers: {}'.format(' '.join(rev_primers)))


# Run through sequences in fasta file and search for primers
with open(seq_file, 'r') as inf:
	sequences.extend(rev_primers)        # all possible primers 
	seq = ''
	found = False
	for line in inf:
		if line[0] == '>':
			if seq:					     # when seq is complete
				found = search(sequences, seq)     # search for primers
			save_id = line.rstrip()[1:]
			seq = ''				     # reset seq
		else:
			seq += line.rstrip()         # extend sequence (spanning multiple lines in file)
	found = search(sequences, seq)       # searches final sequence
	if not found:
		print('Primer not in sequence')











