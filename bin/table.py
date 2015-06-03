"""
Function for reading the table created by the 'parse' subcommand.
"""
import logging
import pandas as pd

logger = logging.getLogger(__name__)

# Thresholds for filtering
V_GENE_COVERAGE = 90  # at least
J_GENE_COVERAGE = 60  # at least
V_GENE_EVALUE = 1E-3  # at most


def read_table_and_filter(path, log=False):
	d = pd.read_csv(path, sep='\t')
	if log: logger.info('%s rows in input table', len(d))

	# Both V and J must be assigned
	filtered = d.dropna(subset=('V_gene', 'J_gene'))[:]
	if log: logger.info('%s rows have both V and J assignment', len(filtered))
	filtered['V_gene'] = pd.Categorical(filtered['V_gene'])

	# Filter out sequences that have a stop codon
	filtered = filtered[filtered.stop == 'no']
	if log: logger.info('%s of those do not have a stop codon', len(filtered))

	# Filter out sequences with a too low V gene hit E-value
	filtered = filtered[filtered.V_evalue <= V_GENE_EVALUE]
	if log: logger.info('%s of those have an E-value of at most %s', len(filtered), V_GENE_EVALUE)

	# Filter out sequences with too low V gene coverage
	filtered = filtered[filtered.V_covered >= V_GENE_COVERAGE]
	if log: logger.info('%s of those cover the V gene by at least %s%%', len(filtered), V_GENE_COVERAGE)

	# Filter out sequences with too low J gene coverage
	filtered = filtered[filtered.J_covered >= J_GENE_COVERAGE]
	if log: logger.info('%s of those cover the J gene by at least %s%%', len(filtered), J_GENE_COVERAGE)

	return filtered
