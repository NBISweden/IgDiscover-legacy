"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from sqt.dna import amino_acid_regex


# Regular expressions for CDR3 detection
CDR3_REGEX = {
	# Heavy chain.
	#
	# This is an extended version of the regular expression by
	# D’Angelo et al.: The antibody mining toolbox.
	# http://dx.doi.org/10.4161/mabs.27105
	#
	# The amino-acid version of this regular expression is:
	# [FYHED][FYHVWDS][CRGW][ETNGASDRIKVM]X{5,31}[WG][GAVDC]
	'VH': re.compile("""
		(TT[CT] | TA[CT] | CA[CT] | GA[AG] | GA[CT] )          # F, Y, H, E, D
		(TT[CT] | TA[CT] | CA[CT] | GT[ACGT] | TGG | GA[CT] |  # F, Y, H, V, W, D,
			TC[ACGT]|AG[CT])                         # ... S
		(TG[CT] | CG[ACGT]|AG[AG] | GG[ACGT] | TGG ) # C, R, G, W
		(?P<cdr3>                                    # actual CDR3 starts here
			(([GA][AGCT]) | TC | CG) [ACGT]          # any of ETNGASDRIKVM
			([ACGT]{3}){5,31}                        # between five and 31 codons
		)                                            # end of CDR3
		(TGG | GG[ACGT])                             # W, G
		(G[CGT][ACGT] | GA[CT] | TG[CT])             # G, A, V, D, C
		""", re.VERBOSE),

	# Light chain, kappa.
	# The amino-acid version is: [YFV][YNHFC][CWFGLS]X{5,15}[FL][GR]
	'VK': re.compile("""
		(TA[CT] | TT[CT] | GT[ACGT )                   # Y, F, V
		(TA[CT] | AA[CT] | CA[CT] | TT[CT] | TG[CT] )  # Y, N, H, F, C
		(TG[CT] | TGG | TT[CT] | GG[ACGT] |            # C, W, F, G ...
			TC[ACGT] | AG[CT] |                        # ... or S
			CT[ACGT] | TT[AG])                         # ... or L,
		(?P<cdr3>
			([ACGT]{3}){5,15}   # between five and fifteen codons
		)
		(
			TT[CT]   |                                 # F ...
			CT[ACGT] | TT[AG]                          # ... or L
		)
		(
			GG[ACGT] |                                 # G ...
			CG[ACGT] | AG[AG]T                         # ... or R
		)
		""", re.VERBOSE),

	# Light chain, lambda.
	# The amino-acid version is: YYCX{5,15}FG
	'VL': re.compile("""
		TA[CT]    # Y
		TA[CT]    # Y
		TG[CT]    # C
		(?P<cdr3>
			([ACGT]{3}){5,15}   # between five and fifteen codons
		)
		TT[CT]    # F
		GG[ACGT]  # G
		""", re.VERBOSE)
}

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
		'CAGCT',
		'CAGGA',
		'CAGGT',
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
