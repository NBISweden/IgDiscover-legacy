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
import re
from collections import namedtuple
from sqt import HelpfulArgumentParser
from sqt import SequenceReader
from sqt.dna import GENETIC_CODE, reverse_complement
from sqt.ansicolor import red, blue

__author__ = "Marcel Martin"


IgblastRecordNT = namedtuple('IgblastRecord',
	'full_sequence query_name cdr3_start hits v_gene d_gene j_gene chain has_stop in_frame is_productive strand')
Hit = namedtuple('Hit', 'query_id query_start query_sequence subject_start subject_sequence percent_identity')


class IgblastRecord(IgblastRecordNT):
	# TODO move computation of cdr3_span, cdr3_sequence, vdj_sequence into constructor
	# TODO maybe make all coordinates relative to full sequence

	cdr3regex = re.compile('(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC|CG)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCTA]')

	@property
	def vdj_sequence(self):
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		hit_v = self.hits['V']
		hit_j = self.hits['J']
		vdj_start = hit_v.query_start
		vdj_stop = hit_j.query_start + len(hit_j.query_sequence)
		return self.full_sequence[vdj_start:vdj_stop]

	def cdr3_span(self):
		"""
		Return (start, end) of CDR3 *relative to VDJ*.
		"""
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		hit_v = self.hits['V']
		hit_j = self.hits['J']

		match = self.cdr3regex.search(self.vdj_sequence)
		if not match:
			return None
		# The first three and the last two codons are not part of the CDR3.
		start = match.start() + 9
		end = match.end() - 6
		assert start < end
		# Make sure that the match starts within V and ends within J.
		if not (start <= len(hit_v.query_sequence) and end >= hit_j.query_start - hit_v.query_start):
			return None
		return (start, end)

	@property
	def cdr3_sequence(self):
		span = self.cdr3_span()
		if span is None:
			return None
		return self.vdj_sequence[span[0]:span[1]]

	def igblast_matches_regex(self):
		if self.cdr3_start is None:
			return True
		span = self.cdr3_span()
		if span is None:
			return True
		start, end = span
		return self.cdr3_start - self.hits['V'].query_start == start


def nt_to_aa(s):
	"""Translate nucleotide sequence to amino acid sequence"""
	return ''.join(GENETIC_CODE.get(s[i:i+3], '*') for i in range(0, len(s), 3))


def run_igblast(fasta, database, organism='rhesus_monkey'):
	"""
	fasta -- path to input FASTA file
	database -- directory that contains IgBLAST databases. Files in that
	directory must be databases created by the makeblastdb program and have
	names organism_gene, such as "rhesus_monkey_V".
	"""
	"""
	TODO
	Igblastwrapper has the databases in files named like this:
	$IGDATA/database/rhesus_monkey_IG_H_J
	Should we also include the _IG_H part?

	TODO
	When running igblastn, the IGDATA environment variable needs to point to
	the directory that contains the internal_data/ directory.
	Either require IGDATA to be set or set it here.
	"""

	arguments = ['igblastn']
	for gene in 'V', 'D', 'J':
		arguments += ['-germline_db_{gene}'.format(gene),
			'{database}/{organism}_{gene}'.format(database, organism, gene)]
	arguments += [
		#TODO '-auxiliary_data', '$IGDATA/optional_file/{organism}_gl.aux',
		'-organism', organism,
		'-ig_seqtype', 'Ig',
		'-num_threads' , '1',
		'-domain_system', 'imgt',
		'-num_alignments_V', '1',
		'-num_alignments_D', '1',
		'-num_alignments_J', '1',
		'-outfmt', "'7 qseqid qstart qseq sstart sseq pident'",
		'-out', 'result3.txt',
		'-query', fasta,
	]


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


def parse_igblast_record(record_lines, fasta_record):
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
				gene, query_id, query_start, query_sequence, subject_start, subject_sequence, percent_identity = line.split('\t')
				query_start = int(query_start) - 1
				subject_start = int(subject_start) - 1
				percent_identity = float(percent_identity)
				assert gene in ('V', 'D', 'J')
				assert gene not in hits, "Two hits for same gene found"
				hits[gene] = Hit(query_id, query_start, query_sequence.replace('-', ''), subject_start, subject_sequence, percent_identity)

	assert fasta_record.name == query_name
	full_sequence = fasta_record.sequence
	if strand == '-':
		full_sequence = reverse_complement(full_sequence)

	if __debug__:
		for gene in ('V', 'D', 'J'):
			if gene not in hits:
				continue
			hit = hits[gene]

			qsequence = hit.query_sequence
			#print(gene, hit.query_start, '-', hit.query_start + len(qsequence), end=' ')

			# IgBLAST removes the trailing semicolon (why, oh why??)
			qname = query_name[:-1] if query_name.endswith(';') else query_name
			qid = hit.query_id
			qid = qid[len('reversed|'):] if qid.startswith('reversed|') else qid
			assert hit.query_id.startswith('reversed|') == (strand == '-')
			assert qid == qname, (qid, qname)
			assert chain in (None, 'VL', 'VH', 'VK', 'NON'), chain
			assert qsequence == full_sequence[hit.query_start:hit.query_start+len(qsequence)]
		#print(len(full_sequence))

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
		hits=hits,
		full_sequence=full_sequence)


def parse_igblast(path, fasta_path):
	"""
	Parse IgBLAST output created with option -outfmt "7 qseqid qstart qseq sstart sseq pident"
	"""
	with SequenceReader(fasta_path) as fasta:
		with open(path) as f:
			for fasta_record, (record_header, record_lines) in zip(fasta, split_by_section(f, ['# IGBLASTN'])):
				assert record_header == '# IGBLASTN 2.2.29+'
				yield parse_igblast_record(record_lines, fasta_record)


def highlight(s, span):
	"""Highlight part of a string in red"""
	if span is None:
		return s
	start, stop = span
	return s[0:start] + red(s[start:stop]) + s[stop:]


def print_row(record):
	"""Print one tab-separated row
	# TODO use csv writer
	"""
	cdr3nt = record.cdr3_sequence
	cdr3aa = nt_to_aa(cdr3nt) if cdr3nt else None
	print(
		record.v_gene,
		record.d_gene,
		record.j_gene,
		record.query_name,
		cdr3nt,
		cdr3aa,
		sep='\t'
	)


def get_argument_parser():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('igblastout', help='output file of igblastn')
	parser.add_argument('fasta', help='File with original reads')
	return parser


def main():
	parser = get_argument_parser()
	args = parser.parse_args()
	n = 0
	for record in parse_igblast(args.igblastout, args.fasta):
		n += 1
		print_row(record)

		#print(record)
		#print('CDR3:', highlight(record.vdj_sequence, record.cdr3_span()))
		#cdr3s = record.cdr3_sequence
		#if cdr3s is not None:
			#aa = nt_to_aa(cdr3s)
			#print(cdr3s, aa)
			##assert ('*' in aa) == record.has_stop
		#else:
			#print('', '')
		#print()
	print(n, 'records')

if __name__ == '__main__':
	main()
