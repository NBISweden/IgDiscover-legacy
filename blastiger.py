#!/usr/bin/env python3
"""
IgBLAST wrapper. Heavily inspired by igblastwrapper.

Desired output is a table with the following columns:
* read name
* V gene
* D gene
* J gene
* CDR3 sequence
"""
from sqt import HelpfulArgumentParser
import subprocess
import sys
import os
from collections import namedtuple

__author__ = "Marcel Martin"


IgblastRecord = namedtuple('IgblastRecord', 'query_name vdj')

def run_igblast():
	"igblastn -germline_db_V $IGDATA/database/rhesus_monkey_IG_H_V -germline_db_J $IGDATA/database/rhesus_monkey_IG_H_J -germline_db_D $IGDATA/database/rhesus_monkey_IG_H_D -auxiliary_data $IGDATA/optional_file/rhesus_monkey_gl.aux -organism rhesus_monkey -ig_seqtype Ig -num_threads 16 -domain_system imgt -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -out result3.txt -query head.fasta -outfmt '7 qseqid qstart qseq sstart sseq'"



def parse_igblast(path):
	"""
	This is how an IgBLAST "record" looks like:
	# IGBLASTN 2.2.29+
# Query: M00559:99:000000000-ACGRF:1:1101:5380:6946;size=515;
# Database: /proj/b2013006/sw/apps/igblastwrp-0.6/data/database/rhesus_monkey_IG_H_V /proj/b2013006/sw/apps/igblastwrp-0.6/data/database/rhesus_monkey_IG_H_D /proj/b2013006/sw/apps/igblastwrp
-0.6/data/database/rhesus_monkey_IG_H_J
# Domain classification requested: imgt

# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.
IGHV4-2*01  IGHD4-1*01  IGHJ4*01    VH  No  In-frame    Yes +

# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
TGCGA   AATA    TGACTACGGT  CCAAT   CTTTG

# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)
FR1-IMGT    101 175 75  72  3   0   96
CDR1-IMGT   176 199 24  18  6   0   75
FR2-IMGT    200 250 51  45  6   0   88.2
CDR2-IMGT   251 274 24  17  7   0   70.8
FR3-IMGT    275 388 114 106 8   0   93
CDR3-IMGT (germline)    389 392 4   4   0   0   100
Total   N/A N/A 292 262 30  0   89.7

# Hit table (the first field indicates the chain type of the hit)
# Fields: query id, q. start, query seq, s. start, subject seq
# 3 hits found
V   M00559:99:000000000-ACGRF:1:1101:5380:6946;size=515 101 CAGGTCCAGCTGCAGGAGTCGGGCCCAGGACTGCTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTCTAGTGGCTCTTTCAGCAGTTACTGGTGGAACTGGATCCGCCAGTCCCCAGGGAAGGGACTGGAGTGGATTGGGGAAATCCATGGTAATAATGAGATCACCAACTATAACCCCTCCCTCAAGAGTCGAGTCACCATTTCAAAAGACGCGTCCAAGAAGCAGGTCTCCCTGAAGCTGAGTTCTGTGACCGCCGCGGACACGGCCGTGTATTTCTGTGCGA    1   CAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTCTGGTGGCTCCATCAGCAGTAACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGGACGTATCTCTGGTAGTGGTGGGAGCACCGACTACAACCCCTCCCTCAAGAGTCGAGTCACCATTTCAACAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGA
D   M00559:99:000000000-ACGRF:1:1101:5380:6946;size=515 397 TGACTACGGT  1   TGACTACGGT
J   M00559:99:000000000-ACGRF:1:1101:5380:6946;size=515 412 CTTTGACTTGTGGGGCCAGGGAGTCCTGGTCACCGTCTCC    5   CTTTGACTACTGGGGCCAGGGAGTCCTGGTCACCGTCTCC

	last line in file:
	# BLAST processed 14 queries
	"""
	n = 0
	state = 'start'
	query_name = None
	with open(path) as f:
		for line in f:
			line = line.strip()
			print('state:', state, line)
			if state == 'start':
				if line.startswith('# IGBLASTN') or line.startswith('# BLAST processed'):
					if line.startswith('# IGBLASTN'):
						assert line == '# IGBLASTN 2.2.29+'
					if query_name is not None:
						yield IgblastRecord(query_name=query_name, vdj=vdj)
					else:
						assert not line.startswith('# BLAST processed')
					state = 'query'
			elif state == 'query':
				if line.startswith('# Query: '):
					query_name = line.split(': ')[1]
					print('found query name:', query_name)
					state = 'vdjheader'
			elif state == 'vdjheader':
				if line.startswith('# V-(D)-J rearrangement summary'):
					state = 'vdj'
			elif state == 'vdj':
				vdj = line
				state = 'start'



def get_argument_parser():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('igblastout', help='output file of igblastn')
	return parser


def main():
	parser = get_argument_parser()
	args = parser.parse_args()
	records = list(parse_igblast(args.igblastout))
	for r in records:
		print(r)
	print(len(records), 'records')

if __name__ == '__main__':
	main()
