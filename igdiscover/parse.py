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
import csv
import logging
from collections import namedtuple
import functools

from sqt.dna import reverse_complement
from sqt.align import edit_distance, hamming_distance

from .utils import nt_to_aa
from .species import find_cdr3, CDR3_SEARCH_START

logger = logging.getLogger(__name__)


def none_if_na(s):
	"""Return None if s == 'N/A'. Return s otherwise."""
	return None if s == 'N/A' else s


def split_by_section(iterable, section_starts):
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
	for line in iterable:
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

# Each alignment summary describes a region in the V region (FR1, CDR1, etc. up to CDR3)
AlignmentSummary = namedtuple('AlignmentSummary', 'start stop length matches mismatches gaps percent_identity')
JunctionVDJ = namedtuple('JunctionVDJ', 'v_end vd_junction d_region dj_junction j_start')
JunctionVJ = namedtuple('JunctionVJ', 'v_end vj_junction j_start')


_Hit = namedtuple('_Hit', [
	'subject_id',  # name of database record, such as "VH4.11"
	'query_start',
	'query_alignment',  # aligned part of the query, with '-' for deletions
	'subject_start',
	'subject_alignment',  # aligned reference, with '-' for insertions
	'subject_length',  # total length of reference, depends only on subject_id
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

	@property
	def subject_end(self):
		return self.subject_start + len(self.subject_sequence)

	@property
	def query_sequence(self):
		return self.query_alignment.replace('-', '')

	@property
	def subject_sequence(self):
		return self.subject_alignment.replace('-', '')

	@property
	def errors(self):
		return sum(a != b for a, b in zip(self.subject_alignment, self.query_alignment))

	def query_position(self, reference_position):
		"""
		Given a position on the reference, return the same position but relative to
		the full query sequence.
		"""
		# Iterate over alignment columns
		ref_pos = self.subject_start
		query_pos = self.query_start
		if ref_pos == reference_position:
			return query_pos
		for ref_c, query_c in zip(self.subject_alignment, self.query_alignment):
			if ref_c != '-':
				ref_pos += 1
			if query_c != '-':
				query_pos += 1
			if ref_pos == reference_position:
				return query_pos
		return None


def parse_header(header):
	"""
	Extract size= and barcode= fields from the FASTA/FASTQ header line

	>>> parse_header("name;size=12;barcode=ACG;")
	('name', 12, 'ACG')
	>>> parse_header("another name;size=200;foo=bar;")
	('another name', 200, None)
	"""
	fields = header.split(';')
	query_name = fields[0]
	size = barcode = None
	for field in fields[1:]:
		if field == '':
			continue
		if '=' in field:
			key, value = field.split('=', maxsplit=1)
			if key == 'size':
				size = int(value)
			elif key == 'barcode':
				barcode = value
	return query_name, size, barcode


class IgBlastRecord:
	def __init__(
		self,
		full_sequence,
		query_name,
		alignments,
		hits,
		v_gene,
		d_gene,
		j_gene,
		chain,
		has_stop,
		in_frame,
		is_productive,
		strand,
		junction
	):
		self.full_sequence = full_sequence
		self.query_name = query_name
		self.alignments = alignments
		self.hits = hits
		self.v_gene = v_gene
		self.d_gene = d_gene
		self.j_gene = j_gene
		self.chain = chain
		self.has_stop = has_stop
		self.in_frame = in_frame
		self.is_productive = is_productive
		self.strand = strand
		self.junction = junction

	def region_sequence(self, region):
		"""
		Return the nucleotide sequence of a named region. Allowed names are:
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

	def __repr__(self):
		return 'IgBlastRecord(query_name={query_name!r}, ' \
			'v_gene={v_gene!r}, d_gene={d_gene!r}, j_gene={j_gene!r}, chain={chain!r}, ...)'.format(
				**vars(self))


class ExtendedIgBlastRecord(IgBlastRecord):
	"""
	This extended record does a few extra things:
	- The CDR3 is detected by using a regular expression
	- The leader is detected within the sequence before the found V gene (by
	searching for the start codon).
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
		'V_aa_mut',
		'J_aa_mut',
		'FR1_aa_mut',
		'CDR1_aa_mut',
		'FR2_aa_mut',
		'CDR2_aa_mut',
		'FR3_aa_mut',
		'V_errors',
		'D_errors',
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
		'V_CDR3_start',
		'VD_junction',
		'D_region',
		'DJ_junction',
		'J_nt',
		'VDJ_nt',
		'VDJ_aa',
		'name',
		'barcode',
		'genomic_sequence',
	]

	CHAINS = {
		'VH': 'heavy', 'VK': 'kappa', 'VL': 'lambda',
		'VA': 'alpha', 'VB': 'beta',
		'VG': 'gamma', 'VD': 'delta'
	}

	def __init__(self, database, **kwargs):
		super().__init__(**kwargs)
		self.query_name, self.size, self.barcode = parse_header(self.query_name)
		self.genomic_sequence = self.full_sequence
		self._database = database
		if 'V' in self.hits:
			self.hits['V'] = self._fixed_v_hit()
		self.utr, self.leader = self._utr_leader()
		self.alignments['CDR3'] = self._find_cdr3()

	@property
	def vdj_sequence(self):
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		hit_v = self.hits['V']
		hit_j = self.hits['J']
		vdj_start = hit_v.query_start
		vdj_stop = hit_j.query_start + len(hit_j.query_sequence)
		return self.full_sequence[vdj_start:vdj_stop]

	@property
	def v_cdr3_start(self):
		"""Start of CDR3 within V"""
		if 'V' not in self.hits or self.alignments['CDR3'] is None:
			return 0
		v_start = self.hits['V'].query_start
		cdr3_start = self.alignments['CDR3'].start
		return cdr3_start - v_start

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

	# TODO this is unused
	def _fixed_cdr3_alignment_by_regex(self):
		"""
		Return a repaired AlignmentSummary object for the CDR3 region which
		does not use IgBLAST’s coordinates. IgBLAST does not determine the end
		of the CDR3 correctly, at least when a custom database is used,

		Return (start, end) of CDR3 relative to query. The CDR3 is detected
		using a regular expression. Return None if no CDR3 detected.
		"""
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		# Search in a window around the V(D)J junction for the CDR3
		if 'CDR3' in self.alignments:
			window_start = self.alignments['CDR3'].start - CDR3_SEARCH_START
		else:
			window_start = max(0, self.hits['V'].query_end - CDR3_SEARCH_START)
		window_end = self.hits['J'].query_end
		window = self.full_sequence[window_start:window_end]
		match = find_cdr3(window, self.chain)
		if not match:
			return None
		start = match[0] + window_start
		end = match[1] + window_start
		assert start < end
		return AlignmentSummary(start=start, stop=end, length=None, matches=None,
			mismatches=None, gaps=None, percent_identity=None)

	def _find_cdr3(self):
		"""
		Return a repaired AlignmentSummary object that describes the CDR3 region.
		Return None if no CDR3 detected.
		"""
		if 'V' not in self.hits or 'J' not in self.hits:
			return None
		if self.chain not in self.CHAINS:
			return None

		# CDR3 start
		cdr3_ref_start = self._database.v_cdr3_start(self.hits['V'].subject_id, self.CHAINS[self.chain])
		if cdr3_ref_start is None:
			return None
		cdr3_query_start = self.hits['V'].query_position(reference_position=cdr3_ref_start)
		if cdr3_query_start is None:
			# Alignment is not long enough to cover CDR3 start position; try to rescue it
			# by assuming that the alignment would continue without indels.
			hit = self.hits['V']
			cdr3_query_start = hit.query_end + (cdr3_ref_start - hit.subject_end)

		# CDR3 end
		cdr3_ref_end = self._database.j_cdr3_end(self.hits['J'].subject_id, self.CHAINS[self.chain])
		if cdr3_ref_end is None:
			return None

		cdr3_query_end = self.hits['J'].query_position(reference_position=cdr3_ref_end)
		if cdr3_query_end is None:
			return None

		return AlignmentSummary(start=cdr3_query_start, stop=cdr3_query_end, length=None, matches=None,
			mismatches=None, gaps=None, percent_identity=None)

	def _fixed_v_hit(self):
		"""
		Extend the V hit to the left if it does not start at the first nucleotide of the V gene.
		"""
		hit = self.hits['V']
		d = hit._asdict()
		while d['subject_start'] > 0 and d['query_start'] > 0:
			d['query_start'] -= 1
			d['subject_start'] -= 1
			preceding_query_base = self.full_sequence[d['query_start']]
			d['query_alignment'] = preceding_query_base + d['query_alignment']
			if self._database.v:
				reference = self._database.v[hit.subject_id]
				preceding_base = reference[d['subject_start']]
			else:
				preceding_base = 'N'
			d['subject_alignment'] = preceding_base + d['subject_alignment']
		return Hit(**d)

	def asdict(self):
		"""
		Return a flattened representation of this record as a dictionary.
		The dictionary can then be used with e.g. a csv.DictWriter or
		pandas.DataFrame.from_items.
		"""
		nt_regions = dict()
		aa_regions = dict()
		for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3'):
			nt_seq = self.region_sequence(region)
			nt_regions[region] = nt_seq
			aa_regions[region] = nt_to_aa(nt_seq) if nt_seq else None

		vdj_nt = self.vdj_sequence
		vdj_aa = nt_to_aa(vdj_nt) if vdj_nt else None

		def nt_mutation_rate(region):
			"""Nucleotide-level mutation rate in percent"""
			if region in self.alignments:
				rar = self.alignments[region]
				if rar is None or rar.percent_identity is None:
					return None
				return 100. - rar.percent_identity
			else:
				return None

		def j_aa_mutation_rate():
			if 'J' not in self.hits:
				return None
			j_subject_id = self.hits['J'].subject_id
			if self.chain not in self.CHAINS:
				return None
			cdr3_ref_end = self._database.j_cdr3_end(j_subject_id, self.CHAINS[self.chain])
			if cdr3_ref_end is None:
				return None
			cdr3_query_end = self.hits['J'].query_position(reference_position=cdr3_ref_end)
			if cdr3_query_end is None:
				return None

			query = self.full_sequence[cdr3_query_end:self.hits['J'].query_end]
			try:
				query_aa = nt_to_aa(query)
			except ValueError:
				return None
			ref = self._database.j[j_subject_id][cdr3_ref_end:self.hits['J'].subject_end]
			try:
				ref_aa = nt_to_aa(ref)
			except ValueError:
				return None
			if not ref_aa:
				return None
			return 100. * edit_distance(ref_aa, query_aa) / len(ref_aa)

		def aa_mutation_rates():
			"""Amino-acid level mutation rates for all regions and V in percent"""
			rates = dict()
			rates['J'] = j_aa_mutation_rate()
			v_aa_mutations = 0
			v_aa_length = 0
			for region in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3'):
				if not aa_regions[region]:
					break
				try:
					reference_sequence = self._database.v_regions_aa[self.v_gene][region]
				except KeyError:
					break
				aa_sequence = aa_regions[region]

				# We previously computed edit distance, but some FR1 alignments are reported
				# with a frameshift by IgBLAST. By requiring that reference and query FR1
				# lengths are identical, we can filter out these cases (and use Hamming distance
				# to get a bit of a speedup)
				if len(nt_regions[region]) != len(self._database.v_regions_nt[self.v_gene][region]):
					break
				mutations = hamming_distance(reference_sequence, aa_sequence)
				length = len(reference_sequence)
				# If the mutation rate is still obviously too high, assume something went
				# wrong and skip mutation rate assignment
				if region == 'FR1' and mutations / length >= 0.8:
					break
				rates[region] = 100. * mutations / length
				v_aa_mutations += mutations
				v_aa_length += length
			else:
				rates['V'] = 100. * v_aa_mutations / v_aa_length
				return rates
			return dict(FR1=None, CDR1=None, FR2=None, CDR2=None, FR3=None, V=None, J=None)

		aa_rates = aa_mutation_rates()

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
			d_errors = self.hits['D'].errors
			d_covered = 100. * self.hits['D'].covered()
			d_evalue = self.hits['D'].evalue
		else:
			d_errors = None
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
			V_covered=v_covered,
			D_covered=d_covered,
			J_covered=j_covered,
			V_evalue=v_evalue,
			D_evalue=d_evalue,
			J_evalue=j_evalue,
			FR1_SHM=nt_mutation_rate('FR1'),
			CDR1_SHM=nt_mutation_rate('CDR1'),
			FR2_SHM=nt_mutation_rate('FR2'),
			CDR2_SHM=nt_mutation_rate('CDR2'),
			FR3_SHM=nt_mutation_rate('FR3'),
			V_SHM=v_shm,
			J_SHM=j_shm,
			V_aa_mut=aa_rates['V'],
			J_aa_mut=aa_rates['J'],
			FR1_aa_mut=aa_rates['FR1'],
			CDR1_aa_mut=aa_rates['CDR1'],
			FR2_aa_mut=aa_rates['FR2'],
			CDR2_aa_mut=aa_rates['CDR2'],
			FR3_aa_mut=aa_rates['FR3'],
			V_errors=v_errors,
			D_errors=d_errors,
			J_errors=j_errors,
			UTR=self.utr,
			leader=self.leader,
			CDR1_nt=nt_regions['CDR1'],
			CDR1_aa=aa_regions['CDR1'],
			CDR2_nt=nt_regions['CDR2'],
			CDR2_aa=aa_regions['CDR2'],
			CDR3_nt=nt_regions['CDR3'],
			CDR3_aa=aa_regions['CDR3'],
			V_nt=v_nt,
			V_aa=v_aa,
			V_end=v_end,
			V_CDR3_start=self.v_cdr3_start,
			VD_junction=vd_junction,
			D_region=d_region,
			DJ_junction=dj_junction,
			J_nt=j_nt,
			VDJ_nt=vdj_nt,
			VDJ_aa=vdj_aa,
			name=self.query_name,
			barcode=self.barcode,
			genomic_sequence=self.genomic_sequence,
		)


class ParseError(Exception):
	pass


class IgBlastParser:
	"""
	Parser for IgBLAST results. Works only when IgBLAST was run with
	the option -outfmt "7 sseqid qstart qseq sstart sseq pident slen".
	"""
	BOOL = {'Yes': True, 'No': False, 'N/A': None}
	FRAME = {'In-frame': True, 'Out-of-frame': False, 'N/A': None}
	SECTIONS = frozenset([
		'# Query:',
		'# V-(D)-J rearrangement summary',
		'# V-(D)-J junction details',
		'# Alignment summary',
		'# Hit table',
		'Total queries = ',
	])

	def __init__(self, sequences, igblast_lines, database=None):
		"""
		If a database is given, iterating over this object will
		yield ExtendedIgBlastRecord objects, otherwise 'normal' IgBlastRecord objects
		"""
		self._sequences = sequences
		self._igblast_lines = igblast_lines
		self._database = database
		if self._database is None:
			self._create_record = IgBlastRecord
		else:
			self._create_record = functools.partial(ExtendedIgBlastRecord, database=self._database)

	def __iter__(self):
		"""
		Yield (Extended-)IgBlastRecord objects
		"""
		zipped = zip(self._sequences, split_by_section(self._igblast_lines, ['# IGBLASTN']))
		for fasta_record, (record_header, record_lines) in zipped:
			# 'IGBLASTN 2.5.1+': IgBLAST 1.6.1
			assert record_header in {
				'# IGBLASTN 2.2.29+',  # IgBLAST 1.4.0
				'# IGBLASTN 2.3.1+',  # IgBLAST 1.5.0
				'# IGBLASTN 2.6.1+',  # IgBLAST 1.7.0
				'# IGBLASTN',  # IgBLAST 1.10
			}
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
			elif section.startswith('Total queries = '):
				continue

		assert fasta_record.name == query_name
		full_sequence = fasta_record.sequence.upper()
		if strand == '-':
			full_sequence = reverse_complement(full_sequence)

		if __debug__:
			for gene in ('V', 'D', 'J'):
				if gene not in hits:
					continue
				hit = hits[gene]
				qsequence = hit.query_sequence

				# IgBLAST removes the trailing semicolon (why, oh why??)
				qname = query_name[:-1] if query_name.endswith(';') else query_name
				assert chain in (None, 'VL', 'VH', 'VK', 'NON', 'VA', 'VB', 'VG', 'VD'), chain
				assert qsequence == full_sequence[hit.query_start:hit.query_start+len(qsequence)]

		return self._create_record(
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
		(gene, subject_id, query_start, query_alignment, subject_start, subject_alignment,
			percent_identity, subject_length, evalue) = line.split('\t')
		query_start = int(query_start) - 1
		subject_start = int(subject_start) - 1
		subject_length = int(subject_length)  # Length of original subject sequence
		# Percent identity is calculated by IgBLAST as
		# 100 - errors / alignment_length and then rounded to two decimal digits
		percent_identity = float(percent_identity)
		evalue = float(evalue)
		hit = Hit(subject_id, query_start, query_alignment, subject_start,
			subject_alignment, subject_length, percent_identity, evalue)
		return hit, gene


class TableWriter:
	def __init__(self, file):
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
		for name in ('V_covered', 'D_covered', 'J_covered',
				'FR1_SHM', 'CDR1_SHM', 'FR2_SHM', 'CDR2_SHM', 'FR3_SHM',
				'V_SHM', 'J_SHM', 'V_aa_mut', 'J_aa_mut',
				'FR1_aa_mut', 'CDR1_aa_mut', 'FR2_aa_mut', 'CDR2_aa_mut', 'FR3_aa_mut'):
			if d[name] is not None:
				d[name] = '{:.1f}'.format(d[name])
		for name in ('V_evalue', 'D_evalue', 'J_evalue'):
			if d[name] is not None:
				d[name] = '{:G}'.format(d[name])
		self._writer.writerow(d)
