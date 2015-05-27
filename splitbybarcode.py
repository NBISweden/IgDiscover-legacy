#!/usr/bin/env python3
"""
Split sequences in a FASTA file by barcode.
"""
__author__ = "Marcel Martin"

import sys
import random
from sqt import HelpfulArgumentParser, FastaReader


def main():
	parser = HelpfulArgumentParser(description=__doc__)
	parser.add_argument('--barcode-length', '-l', metavar='N', type=int, default=12,
		help="Barcode length (default: %(default)s)")
	parser.add_argument("fasta", metavar='FASTA', help="Input FASTA file")
	parser.add_argument("prefix", metavar='PREFIX', nargs='?', default='split',
		help='Prefix for output file names. The output files will be named <PREFIX><BARCODE>.fasta')
	args = parser.parse_args()
	n = 0
	barcodes = set()
	for record in FastaReader(args.fasta):
		barcode = record.sequence[:args.barcode_length]
		barcodes.add(barcode)
		# Open/close files every time to avoid too many open files
		with open(args.prefix + barcode + '.fasta', mode='a') as f:
			f.write('>{}\n{}\n'.format(record.name, record.sequence))
		n += 1
	print('Wrote {} sequences in total to {} files.'.format(n, len(barcodes)))


if __name__ == '__main__':
        main()
