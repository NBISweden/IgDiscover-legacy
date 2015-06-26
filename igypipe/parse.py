"""
Parse IgBLAST output and write out a tab-separated table
"""
from collections import namedtuple
import re
import sys
import csv
import logging
import errno

from sqt import SequenceReader
from sqt.dna import GENETIC_CODE, reverse_complement

logger = logging.getLogger(__name__)


def add_subcommand(subparsers):
	subparser = subparsers.add_parser('parse', help=__doc__)
	subparser.set_defaults(func=parse_command)
	subparser.add_argument('--rename', default=None, metavar='PREFIX',
		help='Rename reads to PREFIXseqN (where N is a number starting at 1)')
	subparser.add_argument('igblast', help='IgBLAST output')
	subparser.add_argument('fasta', help='File with original reads')
	return subparser


def highlight(s, span):
	"""Highlight part of a string in red"""
	if span is None:
		return s
	start, stop = span
	return s[0:start] + red(s[start:stop]) + s[stop:]


IgblastRecordNT = namedtuple('IgblastRecordNT',
	'full_sequence query_name alignments ' \
	'hits v_gene d_gene j_gene chain has_stop in_frame is_productive ' \
	'strand size junction')
AlignmentSummary = namedtuple('AlignmentSummary', 'start stop length matches mismatches gaps percent_identity')
JunctionVDJ = namedtuple('JunctionVDJ', 'v_end vd_junction d_region dj_junction j_start')
JunctionVJ = namedtuple('JunctionVJ', 'v_end vj_junction j_start')


HitNT = namedtuple('HitNT', 'subject_id query_start query_sequence subject_start subject_sequence subject_length percent_identity evalue')
class Hit(HitNT):
	def covered(self):
		"""
		Return fraction of bases in the original subject sequence that are
		covered by this hit.
		"""
		return len(self.subject_sequence) / self.subject_length


sizeregex = re.compile('(.*);size=(\d+);$')  # TODO move into class below

class IgblastRecord(IgblastRecordNT):
	# TODO move computation of cdr3_span, cdr3_sequence, vdj_sequence into constructor
	# TODO maybe make all coordinates relative to full sequence

	# This is a slightly improved version of the regular expression by
	# D’Angelo et al.: The antibody mining toolbox.
	# http://dx.doi.org/10.4161/mabs.27105
	# The amino-acid version of the expression is:
	# [FY][FWVHY]C[ETNGASDRIKVM]X{5,32}W[GAV]
	cdr3regex = re.compile('(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC|CG)[AGCT]([ACGT]{3}){5,31}TGGG[GCT][GCTA]')

	def __init__(self, *args, **kwargs):
		if 'CDR3' in self.alignments:
			self.alignments['CDR3'] = self._fixed_cdr3_alignment()

	@property
	def vdj_sequence(self):
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		hit_v = self.hits['V']
		hit_j = self.hits['J']
		vdj_start = hit_v.query_start
		vdj_stop = hit_j.query_start + len(hit_j.query_sequence)
		return self.full_sequence[vdj_start:vdj_stop]

	def _fixed_cdr3_alignment(self):
		"""
		Return a repaired AlignmentSummary object for the CDR3 region, which
		does not use IgBLAST’s coordinates. IgBLAST does not determine the end
		of the CDR3 correctly.
		"""
		span = self._cdr3_span()
		if span is None:
			return None
		start, stop = span
		return AlignmentSummary(start=start, stop=stop, length=None, matches=None,
			mismatches=None, gaps=None, percent_identity=None)

	def _cdr3_span(self):
		"""
		Return (start, end) of CDR3 relative to query. The CDR3 is detected
		using a regular expression. Return None if no CDR3 detected.
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
		# Make coordinates relative to query
		return (start + hit_v.query_start, end + hit_v.query_start)

	def region_sequence(self, region):
		"""
		Return sequence of a named region. Allowed names are:
		CDR1, CDR2, CDR3, FR1, FR2, FR3. For all regions except CDR3, sequences
		are extracted from the full read using begin and end coordinates from
		IgBLAST’s "alignment summary" table.
		"""
		alignment = self.alignments.get(region, None)
		if alignment is None:
			return None
		if alignment.start is None or alignment.stop is None:
			return None
		return self.full_sequence[alignment.start:alignment.stop]

	"""
	def igblast_matches_regex(self):
		if self.cdr3_start is None:
			return True
		span = self.cdr3_span()
		if span is None:
			return True
		start, end = span
		return self.cdr3_start - self.hits['V'].query_start == start
	"""


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


def none_if_na(s):
	"""Return None if s == 'N/A'. Return s otherwise."""
	return None if s == 'N/A' else s


def parse_alignment_summary(fields):
	start, stop, length, matches, mismatches, gaps = (int(v) for v in fields[:6])
	percent_identity = float(fields[6])
	return AlignmentSummary(
		start=start - 1,
		stop=stop,
		length=length,
		matches=matches,
		mismatches=mismatches,
		gaps=gaps,
		percent_identity=percent_identity
	)


def parse_hit(line):
	"""
	Parse a line of the "Hit table" section and return a tuple (hit, gene) where
	hit is a Hit object.
	"""
	gene, subject_id, query_start, query_sequence, subject_start, subject_sequence, percent_identity, subject_length, evalue = line.split('\t')
	# subject_sequence and query_sequence describe the alignment:
	# They contain '-' characters for insertions and deletions.
	assert len(subject_sequence) == len(query_sequence)
	query_sequence = query_sequence.replace('-', '')
	subject_sequence = subject_sequence.replace('-', '')
	query_start = int(query_start) - 1
	subject_start = int(subject_start) - 1
	subject_length = int(subject_length)  # Length of original subject sequence
	#subject_end = subject_start + len(subject_sequence)
	percent_identity = float(percent_identity)
	evalue = float(evalue)
	hit = Hit(subject_id, query_start, query_sequence, subject_start, subject_sequence, subject_length, percent_identity, evalue)
	return hit, gene


def parse_igblast_record(record_lines, fasta_record):
	BOOL = { 'Yes': True, 'No': False, 'N/A': None }
	FRAME = { 'In-frame': True, 'Out-of-frame': False, 'N/A': None }
	SECTIONS = set([
		'# Query:',
		'# V-(D)-J rearrangement summary',
		'# V-(D)-J junction details',
		'# Alignment summary',
		'# Hit table',
	])
	hits = dict()
	# All of the sections are optional, so we need to set default values here.
	query_name = None
	junction = None
	v_gene, d_gene, j_gene, chain, has_stop, in_frame, is_productive, strand = [None] * 8
	alignments = dict()
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
			v_gene = none_if_na(v_gene)
			d_gene = none_if_na(d_gene)
			j_gene = none_if_na(j_gene)
			chain = none_if_na(chain)
			has_stop = BOOL[has_stop]
			in_frame = FRAME[in_frame]
			is_productive = BOOL[is_productive]
			strand = strand if strand in '+-' else None
		elif section.startswith('# V-(D)-J junction details'):
			fields = lines[0].split('\t')
			if len(fields) == 5:
				junction = JunctionVDJ(
					v_end=fields[0],
					vd_junction=fields[1],
					d_region=fields[2],
					dj_junction=fields[3],
					j_start=fields[4]
				)
			else:
				junction = JunctionVJ(
					v_end=fields[0],
					vj_junction=fields[1],
					j_start=fields[2])
		elif section.startswith('# Alignment summary'):
			for line in lines:
				fields = line.split('\t')
				if len(fields) == 8 and fields[0] != 'Total':
					summary = parse_alignment_summary(fields[1:])
					region_name, _, imgt = fields[0].partition('-')
					assert imgt in ('IMGT', 'IMGT (germline)')
					alignments[region_name] = summary
		elif section.startswith('# Hit table'):
			for line in lines:
				if not line or line.startswith('#'):
					continue
				hit, gene = parse_hit(line)
				assert gene in ('V', 'D', 'J')
				assert gene not in hits, "Two hits for same gene found"
				hits[gene] = hit

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
		alignments=alignments,
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
		size=size,
		junction=junction)


def parse_igblast(path, fasta_path):
	"""
	Parse IgBLAST output created with option -outfmt "7 sseqid qstart qseq sstart sseq pident slen"
	"""
	with SequenceReader(fasta_path) as fasta:
		with open(path) as f:
			for fasta_record, (record_header, record_lines) in zip(fasta, split_by_section(f, ['# IGBLASTN'])):
				assert record_header == '# IGBLASTN 2.2.29+'
				yield parse_igblast_record(record_lines, fasta_record)


def yesno(v):
	"""
	Return "yes", "no" or None for boolean value v, which may also be None.
	"""
	if v is None:
		return None
	return ["no", "yes"][v]


class TableWriter:
	def __init__(self, file, rename=None):
		"""
		If rename is not None, rename reads to {rename}seq{number}, where
		{number} is a sequential number starting from 1.
		"""
		self._file = file
		self._writer = csv.writer(file, delimiter='\t')
		self._writer.writerow([
			"count",
			"V_gene",
			"D_gene",
			"J_gene",
			"stop",
			"productive",
			"before_V",
			"V_covered",
			"D_covered",
			"J_covered",
			"V_evalue",
			"D_evalue",
			"J_evalue",
			"FR1_SHM",
			"CDR1_SHM",
			"FR2_SHM",
			"CDR2_SHM",
			"FR3_SHM",
			"V_SHM",
			"J_SHM",
			"CDR1_nt",
			"CDR1_aa",
			"CDR2_nt",
			"CDR2_aa",
			"CDR3_nt",
			"CDR3_aa",
			"V_nt",
			"V_aa",
			"V_end",
			"VD_junction",
			"D_region",
			"DJ_junction",
			"J_start",
			"name",
			"sequence",
		])
		self.read_number = 1
		self.rename = rename

	def write(self, record):
		cdr1nt = record.region_sequence('CDR1')
		cdr1aa = nt_to_aa(cdr1nt) if cdr1nt else None
		cdr2nt = record.region_sequence('CDR2')
		cdr2aa = nt_to_aa(cdr2nt) if cdr2nt else None
		cdr3nt = record.region_sequence('CDR3')
		cdr3aa = nt_to_aa(cdr3nt) if cdr3nt else None

		def shm(region):
			if region in record.alignments:
				rar = record.alignments[region]
				if rar is None or rar.percent_identity is None:
					return None
				return '{:.1f}'.format(100.0 - rar.percent_identity)
			else:
				return None

		if 'V' in record.hits:
			v_nt = record.hits['V'].query_sequence
			v_aa = nt_to_aa(v_nt)
			v_shm = '{:.1f}'.format(100.0 - record.hits['V'].percent_identity)
			before_v = record.full_sequence[:record.hits['V'].query_start]
			v_covered = '{:.1f}'.format(100*record.hits['V'].covered())
			v_evalue = '{:G}'.format(record.hits['V'].evalue)
		else:
			v_nt = None
			v_aa = None
			v_shm = None
			before_v = None
			v_covered = None
			v_evalue = None
		if 'D' in record.hits:
			d_covered = '{:.1f}'.format(100*record.hits['D'].covered())
			d_evalue = '{:G}'.format(record.hits['D'].evalue)
		else:
			d_covered = None
			d_evalue = None
		if 'J' in record.hits:
			j_shm = '{:1.1f}'.format(100.0 - record.hits['J'].percent_identity)
			j_covered = '{:.1f}'.format(100*record.hits['J'].covered())
			j_evalue = '{:G}'.format(record.hits['J'].evalue)
		else:
			j_shm = None
			j_covered = None
			j_evalue = None
		v_end = getattr(record.junction, 'v_end', None)
		vd_junction = getattr(record.junction, 'vd_junction', None)
		d_region = getattr(record.junction, 'd_region', None)
		dj_junction = getattr(record.junction, 'dj_junction', None)
		j_start = getattr(record.junction, 'j_start', None)
		if self.rename is not None:
			name = "{}seq{}".format(self.rename, self.read_number)
		else:
			name = record.query_name
		self.read_number += 1
		self._writer.writerow([
			record.size,
			record.v_gene,
			record.d_gene,
			record.j_gene,
			yesno(record.has_stop),
			yesno(record.is_productive),
			before_v,
			v_covered,
			d_covered,
			j_covered,
			v_evalue,
			d_evalue,
			j_evalue,
			shm('FR1'),
			shm('CDR1'),
			shm('FR2'),
			shm('CDR2'),
			shm('FR3'),
			v_shm,
			j_shm,
			cdr1nt,
			cdr1aa,
			cdr2nt,
			cdr2aa,
			cdr3nt,
			cdr3aa,
			v_nt,
			v_aa,
			v_end,
			vd_junction,
			d_region,
			dj_junction,
			j_start,
			name,
			record.full_sequence,
		])


def parse_command(args):
	"""
	Parse IgBLAST output
	"""
	n = 0
	writer = TableWriter(sys.stdout, args.rename)
	for record in parse_igblast(args.igblast, args.fasta):
		n += 1
		try:
			writer.write(record)
		except IOError as e:
			if e.errno == errno.EPIPE:
				sys.exit(1)
			raise
		#print('CDR3:', highlight(record.vdj_sequence, record.cdr3_span()))
	logger.info('%d records parsed and written', n)
