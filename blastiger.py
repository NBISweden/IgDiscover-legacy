#!/usr/bin/env python3
"""
An IgBLAST wrapper. Heavily inspired by igblastwrapper.

Desired output is a table with the following columns:
* read name
* V gene
* D gene
* J gene
* CDR3 sequence
"""
import subprocess
import sys
import os
from collections import namedtuple
from sqt import HelpfulArgumentParser
from sqt import SequenceReader
from sqt.dna import GENETIC_CODE

__author__ = "Marcel Martin"


IgblastRecord = namedtuple('IgblastRecord',
	'query_name cdr3_start hits v_gene d_gene j_gene chain has_stop in_frame is_productive strand')
Hit = namedtuple('Hit', 'query_id query_start query_sequence subject_start subject_sequence')


def nt_to_aa(s):
	"""Translate nucleotide sequence to amino acid sequence"""
	return ''.join(GENETIC_CODE[s[i:i+3]] for i in range(0, len(s), 3))


def run_igblast():
	"igblastn -germline_db_V $IGDATA/database/rhesus_monkey_IG_H_V -germline_db_J $IGDATA/database/rhesus_monkey_IG_H_J -germline_db_D $IGDATA/database/rhesus_monkey_IG_H_D -auxiliary_data $IGDATA/optional_file/rhesus_monkey_gl.aux -organism rhesus_monkey -ig_seqtype Ig -num_threads 16 -domain_system imgt -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -out result3.txt -query head.fasta -outfmt '7 qseqid qstart qseq sstart sseq pident'"



def split_by_section(it, section_starts):
	"""
	Parse a stream of lines into chunks of sections. When one of the lines
	starts with a string given in section_starts, a new section is started, and
	a tuple (head, lines) is returned where head is the matching line and lines
	contains a list of the lines following the section header, up to (but
	excluding) the next section header.

	Works a bit like str.split(), but on lines.
	"""
	lines = None
	header = None
	for line in it:
		line = line.strip()
		for start in section_starts:
			if line.startswith(start):
				if header is not None:
					yield (header, lines)
				header = line
				lines = []
				break
		else:
			if header is None:
				raise ParseError("Expected a line starting with one of {}".format(', '.join(section_starts)))
			lines.append(line)
	if header is not None:
		yield (header, lines)


def parse_igblast_record(record_lines):
	BOOL = { 'Yes': True, 'No': False, 'N/A': None }
	FRAME = { 'In-frame': True, 'Out-of-frame': False, 'N/A': None }
	SECTIONS = set([
		'# Query:',
		'# V-(D)-J rearrangement summary',
		'# Alignment summary',
		'# Hit table',
	])
	hits = dict()
	# All of the sections are optional, so we need to set default values here.
	query_name = None
	cdr3_start = None
	v_gene, d_gene, j_gene, chain, has_stop, in_frame, is_productive, strand = [None] * 8
	for section, lines in split_by_section(record_lines, SECTIONS):
		if section.startswith('# Query: '):
			query_name = section.split(': ')[1]
		elif section.startswith('# V-(D)-J rearrangement summary'):
			fields = lines[0].split('\t')
			if len(fields) == 7:
				# No D assignment
				v_gene, j_gene, chain, has_stop, in_frame, is_productive, strand = fields
				d_gene = None
			else:
				v_gene, d_gene, j_gene, chain, has_stop, in_frame, is_productive, strand = fields
			v_gene = None if v_gene == 'N/A' else v_gene
			d_gene = None if d_gene == 'N/A' else d_gene
			j_gene = None if j_gene == 'N/A' else j_gene
			assert chain != 'N/A'
			has_stop = BOOL[has_stop]
			in_frame = FRAME[in_frame]
			is_productive = BOOL[is_productive]
			strand = strand if strand in '+-' else None
		elif section.startswith('# Alignment summary'):
			# This section is optional
			for line in lines:
				if line.startswith('CDR3-IMGT (germline)'):
					cdr3_start = int(line.split('\t')[1]) - 1
					break
		elif section.startswith('# Hit table'):
			for line in lines:
				if not line or line.startswith('#'):
					continue
				# TODO pident (last column)
				gene, query_id, query_start, query_sequence, subject_start, subject_sequence = line.split('\t')
				query_start = int(query_start) - 1
				subject_start = int(subject_start) - 1
				assert gene in ('V', 'D', 'J')
				assert gene not in hits, "Two hits for same gene found"
				hits[gene] = Hit(query_id, query_start, query_sequence, subject_start, subject_sequence)

	if __debug__:
		for gene, hit in hits.items():
			# IgBLAST removes the trailing semicolon (why, oh why??)
			qname = query_name[:-1] if query_name.endswith(';') else query_name
			qid = hit.query_id
			qid = qid[len('reversed|'):] if qid.startswith('reversed|') else qid
			assert hit.query_id.startswith('reversed|') == (strand == '-')
			assert qid == qname, (qid, qname)
			assert chain in (None, 'VL', 'VH', 'VK', 'NON'), chain
	return IgblastRecord(
		query_name=query_name,
		cdr3_start=cdr3_start,
		v_gene=v_gene,
		d_gene=d_gene,
		j_gene=j_gene,
		chain=chain,
		has_stop=has_stop,
		in_frame=in_frame,
		is_productive=is_productive,
		strand=strand,
		hits=hits)


def parse_igblast(path):
	"""
	Parse IgBLAST output, created with option -outfmt "7 qseqid qstart qseq sstart sseq"

	TODO add pident column
	"""
	with open(path) as f:
		for record_header, record_lines in split_by_section(f, ['# IGBLASTN']):
			assert record_header == '# IGBLASTN 2.2.29+'
			yield parse_igblast_record(record_lines)


def get_argument_parser():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('igblastout', help='output file of igblastn')
	parser.add_argument('fasta', help='File with original reads')
	return parser


def main():
	parser = get_argument_parser()
	args = parser.parse_args()
	n = 0
	for record in parse_igblast(args.igblastout):
		n += 1
		print(record)
		print()
	print(n, 'records')

if __name__ == '__main__':
	main()
