"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from sqt.dna import amino_acid_regex
from .utils import nt_to_aa


# Regular expressions for CDR3 detection
#
# The idea comes from D’Angelo et al.: The antibody mining toolbox.
# http://dx.doi.org/10.4161/mabs.27105
# The heavy-chain regex was taken directly from there, but the difference
# is that we express everything in terms of amino acids, not nucleotides.
# This simplifies the expressions and makes them more readable.
#
_CDR3_REGEX = {
	# Heavy chain
	'VH': re.compile("""
		[FY] [FHVWY] C
		(?P<cdr3>
			[ADEGIKMNRSTV] .{5,31}
		)
		W[GAV]
		""", re.VERBOSE),

	# Light chain, kappa
	'VK': re.compile("""
		[FSVY] [CFHNVY] [CDFGLSW]
		(?P<cdr3>
			.{5,15}
		)
		[FLV][GRV]
		""", re.VERBOSE),

	# Light chain, lambda
	'VL': re.compile("""
		# the negative lookahead assertion ensures that the rightmost start is found
		[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]
		(?P<cdr3>
			.{5,15}
		)
		[FS]G
		""", re.VERBOSE)
}

_CDR3_VH_ALTERNATIVE_REGEX = re.compile("""
		C
		(?P<cdr3> . [RK] .{3,30})
		[WF]G.G
""", re.VERBOSE)


def find_cdr3(sequence, chain):
	"""
	Find the CDR3 in the given sequence, assuming it comes from the given chain ('VH', 'VK', 'VL').
	If the chain is not one of 'VH', 'VK', 'VL', return None.

	Return a tuple (start, stop) if found, None otherwise.
	"""
	try:
		regex = _CDR3_REGEX[chain]
	except KeyError:
		return None
	matches = []
	for offset in 0, 1, 2:
		aa = nt_to_aa(sequence[offset:])
		match = regex.search(aa)
		if not match and chain == 'VH':
			match = _CDR3_VH_ALTERNATIVE_REGEX.search(aa)
		if match:
			start, stop = match.span('cdr3')
			matches.append((start * 3 + offset, stop * 3 + offset))
	return min(matches, default=None)


# When searching for the CDR3, start this many bases to the left of the end of
# the V match.
CDR3_SEARCH_START = 30


def _build_V_gene_regex():
	r = '('
	r += '|'.join(amino_acid_regex(aa) for aa in 'DV EA EM EV LV QV QL QE VE'.split())
	r += ')' + amino_acid_regex('Q')
	r += '([ACGT]{3})*'

	#r += '(' + amino_acid_regex('F') + '|' + amino_acid_regex('Y') + ')'
	#r += '([ACGT]{3})*'  # any codon
	# beginning of the CDR3 expression
	#r += '(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC|CG)[AGCT]'
	#r += '[ACGT]{4,6}$'
	return re.compile(r)


# TODO this is unused
_V_GENE_REGEX = _build_V_gene_regex()


# TODO this is work in progress
def looks_like_V_gene(s, chain):
	"""
	Check whether the given sequence matches our expectations of how a V gene
	should look like.
	"""
	s = s.upper()
	for start in START_MOTIFS[chain]:
		if s.startswith(start):
			return True #break
	else:
		return False

	# TODO unusued
	for end in 'TATTACTGT TTTTACTGT TATTACTGC TATTACTGC TATTGTGCA TATTACTGC TATTATTGT'.split():
		if s[-len(end)-13:].find(end) != -1:
			return True
	return False

	#return bool(_V_GENE_REGEX.match(s))


# A newly discovered V gene must start with one of these motifs
START_MOTIFS = {
	'VH': [
		'CAAAT',
		'CAGAT',
		'CAGCT',
		'CAGGA',
		'CAGGT',
		'CGGCT',
		'CGGGT',
		'CTGGT',
		'GAAAT',
		'GAAGT',
		'GACGT',
		'GAGAT',
		'GAGGC',
		'GAGGT',
		'GTGGA',
	],
	'VK': [
		'AACATC',
		'AATATT',
		'AATTTC',
		'CAAGTT',
		'CAGACT',
		'GAAACA',
		'GAAACG',
		'GAAATA',
		'GAAATG',
		'GAAATT',
		'GACATC',
		'GACATG',
		'GACATT',
		'GAGATT',
		'GATACT',
		'GATATC',
		'GATATT',
		'GATGCT',
		'GATGTT',
		'GATTTC',
		'GCCATC',
		'GCTGTT',
		'GTTATT',
	],
	'VL': [
		'AAGCCT',
		'AATTTT',
		'ACATAT',
		'AGGTCT',
		'CAAGCT',
		'CAATCT',
		'CACCCT',
		'CAGACT',
		'CAGCCT',
		'CAGGCT',
		'CAGTCT',
		'CGGTCT',
		'CTGACT',
		'CTGTCT',
		'GAGGCT',
		'GAGGTT',
		'GCATCT',
		'TCCACA',
		'TCCTAT',
	]
}

for _motif_list in START_MOTIFS.values():
	assert len(set(_motif_list)) == len(_motif_list)
