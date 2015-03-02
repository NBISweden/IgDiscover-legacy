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
import csv
import sys
import os
import re
import logging
import errno
from collections import namedtuple
from sqt import HelpfulArgumentParser
from sqt import SequenceReader
from sqt.dna import GENETIC_CODE, reverse_complement
from sqt.ansicolor import red, blue

__author__ = "Marcel Martin"

logger = logging.getLogger(__name__)

IgblastRecordNT = namedtuple('IgblastRecord',
	'full_sequence query_name cdr3_start hits v_gene d_gene j_gene chain has_stop in_frame is_productive strand size')
Hit = namedtuple('Hit', 'query_id query_start query_sequence subject_start subject_sequence percent_identity')

sizeregex = re.compile('(.*);size=(\d+);$')  # TODO move into class below

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

	m = sizeregex.match(query_name)
	if m:
		query_name = m.group(1)
		size = int(m.group(2))
	else:
		size = None
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
		full_sequence=full_sequence,
		size=size)


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


class TableWriter:
	def __init__(self, file):
		self._file = file
		self._writer = csv.writer(file, delimiter='\t')
		self._writer.writerow([
			"# count",
			"V_gene",
			"D_gene",
			"J_gene",
			"V_%identity",
			"stop",
			"CDR3_nt",
			"CDR3_aa",
			"name",
			"sequence",
		])

	def write(self, record):
		cdr3nt = record.cdr3_sequence
		cdr3aa = nt_to_aa(cdr3nt) if cdr3nt else None
		v_percent_identity = record.hits['V'].percent_identity if 'V' in record.hits else None
		self._writer.writerow([
			record.size,
			record.v_gene,
			record.d_gene,
			record.j_gene,
			v_percent_identity,
			"yes" if record.has_stop else "no",
			cdr3nt,
			cdr3aa,
			record.query_name,
			record.full_sequence,
		])


def main():
	logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('igblastout', help='output file of igblastn')
	parser.add_argument('fasta', help='File with original reads')
	args = parser.parse_args()
	n = 0
	writer = TableWriter(sys.stdout)
	for record in parse_igblast(args.igblastout, args.fasta):
		n += 1
		try:
			writer.write(record)
		except IOError as e:
			if e.errno == errno.EPIPE:
				sys.exit(1)
			raise
		#print('CDR3:', highlight(record.vdj_sequence, record.cdr3_span()))
	logger.info('%d records parsed and written', n)


if __name__ == '__main__':
	main()
