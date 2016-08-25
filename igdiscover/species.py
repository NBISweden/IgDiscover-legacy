"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from sqt.dna import amino_acid_regex
from .utils import nt_to_aa


# Regular expressions for CDR3 detection
_CDR3_REGEX = {
	# Heavy chain.
	#
	# This is a slightly improved version of the regular expression by
	# D’Angelo et al.: The antibody mining toolbox.
	# http://dx.doi.org/10.4161/mabs.27105
	# The amino-acid version of the expression is:
	# [FY][FWVHY]C[ETNGASDRIKVM]X{5,32}W[GAV]
	'VH': re.compile("""
		(TT[CT] | TA[CT])                            # F or Y
		(TT[CT] | TA[CT] | CA[CT] | GT[ACGT] | TGG)  # any of F, Y, H, V, W
		(TG[CT])                                     # C
		(?P<cdr3>                                    # actual CDR3 starts here
			(([GA][AGCT]) | TC | CG) [ACGT]          # any of ETNGASDRIKVM
			([ACGT]{3}){5,31}                        # between five and 31 codons
		)                                            # end of CDR3
		TGG                                          # W
		G[CGT][ACGT]                                 # G, A or V
		""", re.VERBOSE),


	# Light chain, kappa.
	# The amino-acid version is: [YFVS][YNHFCV][CWFGDLS]X{5,15}[FLV][GRV]
	'VK': re.compile("""
		(TA[CT] | TT[CT] | GT[ACGT] | TC[ACGT]|AG[CT] )           # Y, F, V, S
		(TA[CT] | AA[CT] | CA[CT] | TT[CT] | TG[CT] | GT[ACGT] )  # Y, N, H, F, C, V
		(TG[CT] | TGG | TT[CT] | GG[ACGT] | GA[CT] |              # C, W, F, G, D ...
			TC[ACGT] | AG[CT] |                                   # ... or S
			CT[ACGT] | TT[AG] )                                   # ... or L,
		(?P<cdr3>
			([ACGT]{3}){5,15}   # between five and fifteen codons
		)
		(
			TT[CT]   | GT[ACGT] |                                # F, V ...
			CT[ACGT] | TT[AG]                                    # ... or L
		)
		(
			GG[ACGT] | GT[ACGT] |                                # G, V ...
			CG[ACGT] | AG[AG]T                                   # ... or R
		)
		""", re.VERBOSE),


	# Light chain, lambda.
	# The amino-acid version is: [CDY][CFHSY][CFGW]X{5,15}[FS]G
	'VL': re.compile("""
		(TG[CT] | GA[CT] | TA[CT] )                             # C, D, Y
		(TG[CT] | TT[CT] | CA[CT] | TC[ACGT]|AG[CT] | TA[CT] )  # C, F, H, S, Y
		(TG[CT] | TT[CT] | GG[ACGT] | TGG )                     # C, F, G, W
		(?P<cdr3>
			([ACGT]{3}){5,15}   # between five and fifteen codons
		)
		( TT[CT] | TC[ACGT]|AG[CT] )                   # F, S
		GG[ACGT]                                       # G
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
	if chain not in _CDR3_REGEX:
		return None
	match = _CDR3_REGEX[chain].search(sequence)
	if match:
		return match.span('cdr3')
	if chain == 'VH':
		for offset in 0, 1, 2:
			aa = nt_to_aa(sequence[offset:])
			match = _CDR3_VH_ALTERNATIVE_REGEX.search(aa)
			if match:
				start, stop = match.span('cdr3')
				start = start * 3 + offset
				stop = stop * 3 + offset
				return start, stop
	return None


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
