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

from tinyalign import edit_distance, hamming_distance

from .utils import nt_to_aa
from .dna import reverse_complement


logger = logging.getLogger(__name__)


def none_if_na(s):
    """Return None if s == 'N/A'. Return s otherwise."""
    return None if s == 'N/A' else s


def gene_without_prefix(s):
    if s == "N/A":
        return None
    else:
        assert s.startswith("%")
        return s[1:]


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
                raise ParseError("Expected a line starting with one of {}".format(
                    ', '.join(section_starts)))
            lines.append(line)
    if header is not None:
        yield (header, lines)


JunctionVDJ = namedtuple('JunctionVDJ', 'v_end vd_junction d_region dj_junction j_start')
JunctionVJ = namedtuple('JunctionVJ', 'v_end vj_junction j_start')


class AlignmentSummary:
    """An alignment summary describes a framework region or complementarity-determining region
    (FR1/2/3, CDR1/2/3)"""

    def __init__(self, start, stop, length, matches, mismatches, gaps, percent_identity):
        self.start = start
        self.stop = stop
        self.length = length
        self.matches = matches
        self.mismatches = mismatches
        self.gaps = gaps
        if matches is not None:
            assert matches + mismatches + gaps == length
        if percent_identity is not None:
            assert matches is not None and length is not None
            assert abs(100. * matches / length - percent_identity) < 0.1
        self.percent_identity = percent_identity


class Hit:
    def __init__(
        self,
        subject_id: str,  # name of database record, such as "VH4.11"
        query_start: int,
        query_alignment: str,  # aligned part of the query, with '-' for deletions
        subject_start: int,
        subject_alignment: str,  # aligned part of reference, with '-' for insertions
        subject_length: int,  # total length of reference, depends only on subject_id
        percent_identity: float,
        evalue: float,
    ):
        assert len(subject_alignment) == len(query_alignment)
        self.subject_id = subject_id
        self.query_start = query_start
        self.query_alignment = query_alignment
        self.subject_start = subject_start
        self.subject_alignment = subject_alignment
        self.subject_length = subject_length
        self.percent_identity = percent_identity
        self.evalue = evalue

        # Derived attributes
        self.errors = self._errors(self.subject_alignment, self.query_alignment)
        self.query_sequence = self.query_alignment.replace('-', '')
        self.subject_sequence = self.subject_alignment.replace('-', '')

        assert abs(self.percent_identity - self._percent_identity()) < 0.01

    def extend_left_ungapped(self, query_sequence, subject_sequence):
        """
        Extend this hit to the left until it reaches the first nucleotide of the subject sequence.
        Used for extending V hits to the 5' end.
        """
        query_bases = []
        subject_bases = []

        while self.subject_start > 0 and self.query_start > 0:
            self.query_start -= 1
            self.subject_start -= 1
            query_base = query_sequence[self.query_start]
            query_bases.append(query_base)
            if subject_sequence is not None:
                subject_base = subject_sequence[self.subject_start]
            else:
                subject_base = 'N'
            subject_bases.append(subject_base)

        query_bases = ''.join(query_bases[::-1])
        subject_bases = ''.join(subject_bases[::-1])
        self.query_alignment = query_bases + self.query_alignment
        self.subject_alignment = subject_bases + self.subject_alignment
        self.query_sequence = query_bases + self.query_sequence
        self.subject_sequence = subject_bases + self.subject_sequence
        self.errors += self._errors(query_bases, subject_bases)
        self.percent_identity = self._percent_identity()

    def _percent_identity(self):
        """This is how IgBLAST computes percent identity"""
        matches = len(self.subject_alignment) - self.errors
        return 100. * matches / len(self.subject_alignment)

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

    @staticmethod
    def _errors(alignment1, alignment2):
        return sum(a != b for a, b in zip(alignment1, alignment2))

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
        CDR1, CDR2, CDR3, FR1, FR2, FR3. Sequences are extracted from the full read
        using begin and end coordinates from IgBLAST’s "alignment summary" table.
        """
        if region not in ("CDR1", "CDR2", "CDR3", "FR1", "FR2", "FR3"):
            raise KeyError(f"Region {region!r} not allowed")
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


class Region:
    """A CDR or FR region in a V(D)J rearranged sequence (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4)"""

    def __init__(self, nt_sequence, nt_reference, aa_reference=None, percent_identity=None):
        self.nt_sequence = nt_sequence
        self.aa_sequence = nt_to_aa(nt_sequence) if nt_sequence else None
        self.nt_reference = nt_reference
        if aa_reference is None and nt_reference is not None:
            aa_reference = nt_to_aa(nt_reference)
        self.aa_reference = aa_reference
        self.aa_mutations = self._aa_mutations()
        self.percent_identity = percent_identity
        if percent_identity is None:
            self.percent_identity = self._percent_identity()

    def _percent_identity(self):
        # FIXME This is not quite how IgBLAST computes percent identity
        if not self.nt_reference or not self.nt_sequence:
            return None
        dist = edit_distance(self.nt_reference, self.nt_sequence)
        return 100. - 100. * dist / len(self.nt_reference)

    def _aa_mutations(self):
        # Earlier versions of this code used edit distance to compute the number of mutations,
        # but some FR1 alignments are reported with a frameshift by IgBLAST. By requiring that
        # reference and query lengths are identical, we can filter out these cases (and use
        # Hamming distance to get some speedup)
        if (
            not self.aa_reference
            or not self.aa_sequence
            or len(self.nt_sequence) != len(self.nt_reference)
        ):
            return None
        dist = hamming_distance(self.aa_reference, self.aa_sequence)

        # If the mutation rate is still obviously too high, assume something went
        # wrong and ignore the computed value
        if dist / len(self.aa_reference) >= 0.8:
            return None
        return dist

    def aa_mutation_rate(self):
        if self.aa_mutations is None or not self.aa_reference:
            return None
        return 100. * self.aa_mutations / len(self.aa_reference)

    def nt_mutation_rate(self):
        """Return nucleotide-level mutation rate in percent"""
        if self.percent_identity is None:
            return None
        return 100. - self.percent_identity


class ExtendedIgBlastRecord(IgBlastRecord):
    """
    This extended record does a few extra things:
    - The CDR3 is detected
    - The leader is detected within the sequence before the found V gene (by
    searching for the start codon).
    - If the V sequence hit starts not at base 1 in the reference, it is extended
    to the left.
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
        'FR4_SHM',
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
            subject_sequence = self._database.v[self.hits['V'].subject_id]
            self.hits['V'].extend_left_ungapped(self.full_sequence, subject_sequence)
        self.utr, self.leader = self._utr_leader()
        self.alignments['CDR3'] = self._find_cdr3()

        self.regions = {
            name: self._make_region(name) for name in
            ('FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'CDR3')}
        self.regions['FR4'] = self._make_fr4_region()
        self.vdj_sequence = self._make_vdj_sequence()

    def _make_region(self, name: str):
        nt_sequence = self.region_sequence(name)
        if self.v_gene in self._database.v_regions_nt:
            nt_reference = self._database.v_regions_nt[self.v_gene].get(name)
        else:
            nt_reference = None
        if self.v_gene in self._database.v_regions_aa:
            aa_reference = self._database.v_regions_aa[self.v_gene].get(name)
        else:
            aa_reference = None
        if self.alignments.get(name, None) is not None:
            percent_identity = self.alignments[name].percent_identity
        else:
            percent_identity = None
        return Region(nt_sequence, nt_reference, aa_reference, percent_identity)

    def _make_fr4_region(self):
        if 'J' not in self.hits:
            return Region(None, None)
        j_subject_id = self.hits['J'].subject_id
        if self.chain not in self.CHAINS:
            return Region(None, None)
        cdr3_ref_end = self._database.j_cdr3_end(j_subject_id, self.CHAINS[self.chain])
        if cdr3_ref_end is None:
            return Region(None, None)
        cdr3_query_end = self.hits['J'].query_position(reference_position=cdr3_ref_end)
        if cdr3_query_end is None:
            return Region(None, None)

        query = self.full_sequence[cdr3_query_end:self.hits['J'].query_end]
        ref = self._database.j[j_subject_id][cdr3_ref_end:self.hits['J'].subject_end]

        return Region(query, ref)

    def _make_vdj_sequence(self):
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
                if before_v[-i + offset:-i + 3 + offset] == 'ATG':
                    return before_v[:-i + offset], before_v[-i + offset:]
        return None, None

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

        return AlignmentSummary(start=cdr3_query_start, stop=cdr3_query_end, length=None,
            matches=None, mismatches=None, gaps=None, percent_identity=None)

    def fr4_aa_mutation_rate(self):
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

    def v_aa_mutation_rate(self):
        """
        TODO This returns actually the total mutation rate of the FR1, CDR1, FR2, CDR2, FR3 regions
        (The CDR3 part of V is excluded.)
        """
        mutations = 0
        length = 0
        for name in ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3'):
            region = self.regions.get(name)
            if region is None:
                return None
            if region.aa_reference is None or region.aa_mutations is None:
                return None
            mutations += region.aa_mutations
            length += len(region.aa_reference)

        return 100. * mutations / length

    def asdict(self):
        """
        Return a flattened representation of this record as a dictionary.
        The dictionary can then be used with e.g. a csv.DictWriter or
        pandas.DataFrame.from_items.
        """
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
            FR1_SHM=self.regions['FR1'].nt_mutation_rate(),
            CDR1_SHM=self.regions['CDR1'].nt_mutation_rate(),
            FR2_SHM=self.regions['FR2'].nt_mutation_rate(),
            CDR2_SHM=self.regions['CDR2'].nt_mutation_rate(),
            FR3_SHM=self.regions['FR3'].nt_mutation_rate(),
            FR4_SHM=self.regions['FR4'].nt_mutation_rate(),
            V_SHM=v_shm,
            J_SHM=j_shm,
            V_aa_mut=self.v_aa_mutation_rate(),
            J_aa_mut=self.regions['FR4'].aa_mutation_rate(),  # TODO J vs FR4
            FR1_aa_mut=self.regions['FR1'].aa_mutation_rate(),
            CDR1_aa_mut=self.regions['CDR1'].aa_mutation_rate(),
            FR2_aa_mut=self.regions['FR2'].aa_mutation_rate(),
            CDR2_aa_mut=self.regions['CDR2'].aa_mutation_rate(),
            FR3_aa_mut=self.regions['FR3'].aa_mutation_rate(),
            # FR4_aa_mut=aa_rates['FR4'],  # TODO
            V_errors=v_errors,
            D_errors=d_errors,
            J_errors=j_errors,
            UTR=self.utr,
            leader=self.leader,
            CDR1_nt=self.regions['CDR1'].nt_sequence,
            CDR1_aa=self.regions['CDR1'].aa_sequence,
            CDR2_nt=self.regions['CDR2'].nt_sequence,
            CDR2_aa=self.regions['CDR2'].aa_sequence,
            CDR3_nt=self.regions['CDR3'].nt_sequence,
            CDR3_aa=self.regions['CDR3'].aa_sequence,
            V_nt=v_nt,
            V_aa=v_aa,
            V_end=v_end,
            V_CDR3_start=self.v_cdr3_start,
            VD_junction=vd_junction,
            D_region=d_region,
            DJ_junction=dj_junction,
            J_nt=j_nt,
            VDJ_nt=self.vdj_sequence,
            VDJ_aa=nt_to_aa(self.vdj_sequence) if self.vdj_sequence else None,
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
                    d_gene = gene_without_prefix(d_gene)
                v_gene = gene_without_prefix(v_gene)
                j_gene = gene_without_prefix(j_gene)
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
        assert abs(percent_identity - 100. * matches / length) < 0.1
        # Note length is not necessarily equal to stop - start. Not sure why.
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
        # Names have been mangled by adding a '%' to the beginning. Otherwise, IgBLAST
        # may recognize sequence ids such as AB123456 as being an accession and mangle
        # them to conform to BLAST’s naming scheme.
        # Undo this here.
        assert subject_id.startswith('%')
        subject_id = subject_id[1:]
        query_start = int(query_start) - 1
        subject_start = int(subject_start) - 1
        subject_length = int(subject_length)  # Length of original subject sequence
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
                'FR1_SHM', 'CDR1_SHM', 'FR2_SHM', 'CDR2_SHM', 'FR3_SHM', 'FR4_SHM',
                'V_SHM', 'J_SHM', 'V_aa_mut', 'J_aa_mut',
                'FR1_aa_mut', 'CDR1_aa_mut', 'FR2_aa_mut', 'CDR2_aa_mut', 'FR3_aa_mut'):
            if d[name] is not None:
                d[name] = '{:.1f}'.format(d[name])
        for name in ('V_evalue', 'D_evalue', 'J_evalue'):
            if d[name] is not None:
                d[name] = '{:G}'.format(d[name])
        self._writer.writerow(d)
