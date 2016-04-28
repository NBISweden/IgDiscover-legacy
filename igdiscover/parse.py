"""
Parse IgBLAST output and write out a tab-separated table.

IgBLAST must have been run with -outfmt "7 sseqid qstart qseq sstart sseq pident slen"

A few extra things are done in addition to parsing:
- The CDR3 is detected by using a regular expression
- The leader is detected within the sequence before the found V gene (by
  searching for the start codon).
- If the V sequence hit starts at base 2 in the reference, it is extended
  one to the left.
"""
import re
import sys
import csv
import errno
import logging
from collections import namedtuple

import pandas as pd
from sqt import SequenceReader, xopen
from sqt.dna import reverse_complement
from .utils import nt_to_aa
from .species import CDR3_REGEX, CDR3_SEARCH_START

logger = logging.getLogger(__name__)


def add_arguments(parser):
	parser.add_argument('--rename', default=None, metavar='PREFIX',
		help='Rename reads to PREFIXseqN (where N is a number starting at 1)')
	parser.add_argument('--vdatabase', '--vdb', metavar='FASTA',
		help="Path to FASTA file with V genes. Used to fix the 5' ends of V "
		"gene alignments. If not given, 'N' bases will be inserted instead.")
	parser.add_argument('--hdf5', metavar='FILE',
		help='Write table in HDF5 format to FILE')
	parser.add_argument('igblast', help='IgBLAST output')
	parser.add_argument('fasta', help='File with original reads')


def none_if_na(s):
	"""Return None if s == 'N/A'. Return s otherwise."""
	return None if s == 'N/A' else s


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


AlignmentSummary = namedtuple('AlignmentSummary', 'start stop length matches mismatches gaps percent_identity')
JunctionVDJ = namedtuple('JunctionVDJ', 'v_end vd_junction d_region dj_junction j_start')
JunctionVJ = namedtuple('JunctionVJ', 'v_end vj_junction j_start')


_Hit = namedtuple('_Hit', [
	'subject_id',  # name of database record, such as "VH4.11"
	'query_start',
	'query_sequence',  # aligned part of the query
	'subject_start',
	'subject_sequence',  # aligned part of the reference
	'subject_length',  # total length of reference, depends only on subject_id
	'errors',
	'percent_identity',
	'evalue',
])
class Hit(_Hit):
	# This avoids having a __dict__ attribute, which is necessary for namedtuple
	# subclasses that need _asdict() to work (http://bugs.python.org/issue24931)
	__slots__ = ()

	def covered(self):
		"""
		Return fraction of bases in the original subject sequence that are
		covered by this hit.
		"""
		return len(self.subject_sequence) / self.subject_length

	@property
	def query_end(self):
		return self.query_start + len(self.query_sequence)


def parse_header(header):
	"""Extract size= and barcode= fields from the FASTA/FASTQ header line"""
	fields = header.split(maxsplit=1)
	if len(fields) < 2:
		return header, None, None
	query_name, rest = fields
	fields = rest.split(';')
	size = barcode = None
	for field in rest.split(';'):
		if field == '':
			continue
		if '=' in field:
			key, value = field.split('=', maxsplit=1)
			if key == 'size':
				size = int(value)
			elif key == 'barcode':
				barcode = value
	return query_name, size, barcode


IgBlastRecord = namedtuple('IgBlastRecord',
	'full_sequence query_name alignments '
	'hits v_gene d_gene j_gene chain has_stop in_frame is_productive '
	'strand junction')


class ExtendedIgBlastRecord(IgBlastRecord):
	"""
	This extended record does a few extra things:
	- The CDR3 is detected by using a regular expression
	- The leader is detected within the sequence before the found V gene (by
	searching for the start codon).
	- The RACE-specific run of G in the beginning of the sequence is detected.
	- If the V sequence hit starts at base 2 in the reference, it is extended
	  one to the left.
	"""
	# TODO move computation of cdr3_sequence, vdj_sequence into constructor
	# TODO maybe make all coordinates relative to full sequence

	# Order of columns (use with asdict())
	columns = [
		'count',
		'V_gene',
		'D_gene',
		'J_gene',
		'chain',
		'stop',
		'productive',
		'V_covered',
		'D_covered',
		'J_covered',
		'V_evalue',
		'D_evalue',
		'J_evalue',
		'FR1_SHM',
		'CDR1_SHM',
		'FR2_SHM',
		'CDR2_SHM',
		'FR3_SHM',
		'V_SHM',
		'J_SHM',
		'V_errors',
		'J_errors',
		'UTR',
		'leader',
		'CDR1_nt',
		'CDR1_aa',
		'CDR2_nt',
		'CDR2_aa',
		'CDR3_nt',
		'CDR3_aa',
		'V_nt',
		'V_aa',
		'V_end',
		'VD_junction',
		'D_region',
		'DJ_junction',
		'J_nt',
		'name',
		'barcode',
		'race_G',
		'genomic_sequence',
	]
	def __new__(cls, record, v_database):
		d = record._asdict()
		query_name, size, barcode = parse_header(record.query_name)
		d['query_name'] = query_name
		obj = super().__new__(cls, **d)
		obj.size = size
		obj.barcode = barcode
		return obj

	def __init__(self, record, v_database):
		self.genomic_sequence = self.full_sequence
		self.race_g = None  # TODO since the parse script does not extract the race_G anymore, we don’t have this info available
		if 'V' in self.hits:
			self.hits['V'] = self._fixed_v_hit(v_database)
		self.utr, self.leader = self._utr_leader()
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

	def _utr_leader(self):
		"""
		Split the sequence before the V gene match into UTR and leader by
		searching for the start codon.
		"""
		if 'V' not in self.hits:
			return None, None
		before_v = self.full_sequence[:self.hits['V'].query_start]

		# Search for the start codon
		for offset in (0, 1, 2):
			for i in range(66, 42, -3):
				if before_v[-i + offset : -i + 3 + offset] == 'ATG':
					return before_v[:-i + offset], before_v[-i + offset:]
		return None, None

	def _fixed_cdr3_alignment(self):
		"""
		Return a repaired AlignmentSummary object for the CDR3 region which
		does not use IgBLAST’s coordinates. IgBLAST does not determine the end
		of the CDR3 correctly, at least when a custom database is used,

		Return (start, end) of CDR3 relative to query. The CDR3 is detected
		using a regular expression. Return None if no CDR3 detected.
		"""
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		if not self.chain in CDR3_REGEX:
			return None
		# Search in a window around the V(D)J junction for the CDR3
		if 'CDR3' in self.alignments:
			window_start = self.alignments['CDR3'].start - CDR3_SEARCH_START
		else:
			window_start = max(0, self.hits['V'].query_end - CDR3_SEARCH_START)
		window_end = self.hits['J'].query_end
		window = self.full_sequence[window_start:window_end]
		match = CDR3_REGEX[self.chain].search(window)
		if not match:
			return None
		start = match.start('cdr3') + window_start
		end = match.end('cdr3') + window_start
		assert start < end
		return AlignmentSummary(start=start, stop=end, length=None, matches=None,
			mismatches=None, gaps=None, percent_identity=None)

	def _fixed_v_hit(self, v_database):
		"""
		Extend the V hit by one base to the left if it starts at the second
		base in the V reference.
		"""
		hit = self.hits['V']
		if hit.subject_start != 1 or hit.query_start == 0:
			return hit
		d = hit._asdict()
		assert self.full_sequence[d['query_start']:d['query_start'] + len(d['query_sequence'])] == d['query_sequence']
		d['query_start'] -= 1
		d['subject_start'] -= 1
		d['query_sequence'] = self.full_sequence[d['query_start']:d['query_start'] + len(d['query_sequence']) + 1]
		if v_database:
			reference = v_database[hit.subject_id]
			preceding_base = reference[d['subject_start']]
		else:
			preceding_base = 'N'
		d['subject_sequence'] = preceding_base + d['subject_sequence']  # TODO maybe not needed
		# TODO one could adjust the no. of errors, percent_identity (and E-value)
		return Hit(**d)

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

	def asdict(self):
		"""
		Return a flattened representation of this record as a dictionary.
		The dictionary can then be used with e.g. a csv.DictWriter or
		pandas.DataFrame.from_items.
		"""
		cdr1nt = self.region_sequence('CDR1')
		cdr1aa = nt_to_aa(cdr1nt) if cdr1nt else None
		cdr2nt = self.region_sequence('CDR2')
		cdr2aa = nt_to_aa(cdr2nt) if cdr2nt else None
		cdr3nt = self.region_sequence('CDR3')
		cdr3aa = nt_to_aa(cdr3nt) if cdr3nt else None

		def shm(region):
			if region in self.alignments:
				rar = self.alignments[region]
				if rar is None or rar.percent_identity is None:
					return None
				return 100. - rar.percent_identity
			else:
				return None

		if 'V' in self.hits:
			v_nt = self.hits['V'].query_sequence
			v_aa = nt_to_aa(v_nt)
			v_shm = 100. - self.hits['V'].percent_identity
			v_errors = self.hits['V'].errors
			v_covered = 100. * self.hits['V'].covered()
			v_evalue = self.hits['V'].evalue
		else:
			v_nt = None
			v_aa = None
			v_shm = None
			v_errors = None
			v_covered = None
			v_evalue = None
		if 'D' in self.hits:
			d_covered = 100. * self.hits['D'].covered()
			d_evalue = self.hits['D'].evalue
		else:
			d_covered = None
			d_evalue = None
		if 'J' in self.hits:
			j_nt = self.hits['J'].query_sequence
			j_shm = 100. - self.hits['J'].percent_identity
			j_errors = self.hits['J'].errors
			j_covered = 100. * self.hits['J'].covered()
			j_evalue = self.hits['J'].evalue
		else:
			j_nt = None
			j_shm = None
			j_errors = None
			j_covered = None
			j_evalue = None
		v_end = getattr(self.junction, 'v_end', None)
		vd_junction = getattr(self.junction, 'vd_junction', None)
		d_region = getattr(self.junction, 'd_region', None)
		dj_junction = getattr(self.junction, 'dj_junction', None)
		return dict(
			count=self.size,
			V_gene=self.v_gene,
			D_gene=self.d_gene,
			J_gene=self.j_gene,
			chain=self.chain,
			stop=self.has_stop,
			productive=self.is_productive,
			V_covered=v_covered,
			D_covered=d_covered,
			J_covered=j_covered,
			V_evalue=v_evalue,
			D_evalue=d_evalue,
			J_evalue=j_evalue,
			FR1_SHM=shm('FR1'),
			CDR1_SHM=shm('CDR1'),
			FR2_SHM=shm('FR2'),
			CDR2_SHM=shm('CDR2'),
			FR3_SHM=shm('FR3'),
			V_SHM=v_shm,
			J_SHM=j_shm,
			V_errors=v_errors,
			J_errors=j_errors,
			UTR=self.utr,
			leader=self.leader,
			CDR1_nt=cdr1nt,
			CDR1_aa=cdr1aa,
			CDR2_nt=cdr2nt,
			CDR2_aa=cdr2aa,
			CDR3_nt=cdr3nt,
			CDR3_aa=cdr3aa,
			V_nt=v_nt,
			V_aa=v_aa,
			V_end=v_end,
			VD_junction=vd_junction,
			D_region=d_region,
			DJ_junction=dj_junction,
			J_nt=j_nt,
			name=self.query_name,
			barcode=self.barcode,
			race_G=self.race_g,
			genomic_sequence=self.genomic_sequence,
		)


class ParseError(Exception):
	pass


class IgBlastParser:
	"""
	Parser for IgBLAST results. Works only when IgBLAST was run with
	the option -outfmt "7 sseqid qstart qseq sstart sseq pident slen".
	"""
	BOOL = { 'Yes': True, 'No': False, 'N/A': None }
	FRAME = { 'In-frame': True, 'Out-of-frame': False, 'N/A': None }
	SECTIONS = frozenset([
		'# Query:',
		'# V-(D)-J rearrangement summary',
		'# V-(D)-J junction details',
		'# Alignment summary',
		'# Hit table',
	])

	def __init__(self, sequence_file, igblast_file):
		self._sequence_reader = SequenceReader(sequence_file)
		self._igblast_file = igblast_file

	def __iter__(self):
		"""
		Yield IgBlastRecord objects
		"""
		zipped = zip(self._sequence_reader, split_by_section(self._igblast_file, ['# IGBLASTN']))
		for fasta_record, (record_header, record_lines) in zipped:
			assert record_header == '# IGBLASTN 2.2.29+'
			yield self._parse_record(record_lines, fasta_record)

	def _parse_record(self, record_lines, fasta_record):
		"""
		Parse a single IgBLAST record
		"""
		hits = dict()
		# All of the sections are optional, so we need to set default values here.
		query_name = None
		junction = None
		v_gene, d_gene, j_gene, chain, has_stop, in_frame, is_productive, strand = [None] * 8
		alignments = dict()
		for section, lines in split_by_section(record_lines, self.SECTIONS):
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
				has_stop = self.BOOL[has_stop]
				in_frame = self.FRAME[in_frame]
				is_productive = self.BOOL[is_productive]
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
						summary = self._parse_alignment_summary(fields[1:])
						region_name, _, imgt = fields[0].partition('-')
						assert imgt in ('IMGT', 'IMGT (germline)')
						alignments[region_name] = summary
			elif section.startswith('# Hit table'):
				for line in lines:
					if not line or line.startswith('#'):
						continue
					hit, gene = self._parse_hit(line)
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

		return IgBlastRecord(
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
			junction=junction)

	def _parse_alignment_summary(self, fields):
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


	def _parse_hit(self, line):
		"""
		Parse a line of the "Hit table" section and return a tuple (hit, gene)
		where hit is a Hit object.
		"""
		gene, subject_id, query_start, query_sequence, subject_start, subject_sequence, percent_identity, subject_length, evalue = line.split('\t')
		# subject_sequence and query_sequence describe the alignment:
		# They contain '-' characters for insertions and deletions.
		assert len(subject_sequence) == len(query_sequence)
		alignment_length = len(subject_sequence)
		errors = sum(a != b for a,b in zip(subject_sequence, query_sequence))
		query_sequence = query_sequence.replace('-', '')
		subject_sequence = subject_sequence.replace('-', '')
		query_start = int(query_start) - 1
		subject_start = int(subject_start) - 1
		subject_length = int(subject_length)  # Length of original subject sequence
		#subject_end = subject_start + len(subject_sequence)
		# Percent identity is calculated by IgBLAST as
		# 100 - errors / alignment_length and then rounded to two decimal digits
		percent_identity = float(percent_identity)
		evalue = float(evalue)
		hit = Hit(subject_id, query_start, query_sequence, subject_start,
			subject_sequence, subject_length, errors, percent_identity, evalue)
		return hit, gene


class TableWriter:
	def __init__(self, file):
		"""
		If rename is not None, rename reads to {rename}seq{number}, where
		{number} is a sequential number starting from 1.
		"""
		self._file = file
		self._writer = csv.DictWriter(file, fieldnames=ExtendedIgBlastRecord.columns, delimiter='\t')
		self._writer.writeheader()

	@staticmethod
	def yesno(v):
		"""
		Return "yes", "no" or None for boolean value v, which may also be None.
		"""
		if v is None:
			return None
		return ["no", "yes"][v]

	def write(self, d):
		"""
		Write the IgBLAST record (must be given as dictionary) to the output
		file.
		"""
		d = d.copy()
		d['stop'] = self.yesno(d['stop'])
		d['productive'] = self.yesno(d['productive'])
		for name in ('V_covered', 'D_covered', 'J_covered',
				'FR1_SHM', 'CDR1_SHM', 'FR2_SHM', 'CDR2_SHM', 'FR3_SHM',
				'V_SHM', 'J_SHM'):
			if d[name] is not None:
				d[name] = '{:.1f}'.format(d[name])
		for name in ('V_evalue', 'D_evalue', 'J_evalue'):
			if d[name] is not None:
				d[name] = '{:G}'.format(d[name])
		self._writer.writerow(d)


def main(args):
	"""
	Parse IgBLAST output
	"""
	n = 0
	if args.hdf5:
		rows = []

	if args.vdatabase:
		with SequenceReader(args.vdatabase) as fr:
			v_database = { record.name: record.sequence.upper() for record in fr }
	else:
		v_database = None

	writer = TableWriter(sys.stdout)
	with xopen(args.fasta) as fasta, xopen(args.igblast) as igblast:
		parser = IgBlastParser(fasta, igblast)
		for record in parser:
			n += 1
			extended_record = ExtendedIgBlastRecord(record, v_database=v_database)
			d = extended_record.asdict()
			if args.rename is not None:
				d['name'] = "{}seq{}".format(args.rename, n)
			if args.hdf5:
				rows.append(d)
			try:
				writer.write(d)
			except IOError as e:
				if e.errno == errno.EPIPE:
					sys.exit(1)
				raise
	logger.info('%d records parsed and written', n)
	if args.hdf5:
		from igdiscover.table import STRING_COLUMNS, INTEGER_COLUMNS
		df = pd.DataFrame(rows, columns=ExtendedIgBlastRecord.columns)
		# TODO code below is copied from igdiscover.table
		# Convert all string columns to str to avoid a PerformanceWarning
		for col in STRING_COLUMNS:
			df[col].fillna('', inplace=True)
			df[col] = df[col].astype('str')
			# Empty strings have been set to NaN by read_csv. Replacing
			# by the empty string avoids problems with groupby, which
			# ignores NaN values.
			# Columns that have any NaN values in them cannot be converted to
			# int due to a numpy limitation.
			for col in INTEGER_COLUMNS:
				if all(df[col].notnull()):
					df[col] = df[col].astype(int)
		df.to_hdf(args.hdf5, 'table', mode='w', complevel=3, complib='zlib')
		logger.info('HDF5 file written')
