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
	parser.add_argument('--sequential', '-s', action='store_true', default=False,
		help="Number output files sequentially, do not use BARCODE in the file name.")
	parser.add_argument("fasta", metavar='FASTA', help="Input FASTA file")
	parser.add_argument("prefix", metavar='PREFIX', nargs='?', default='split',
		help='Prefix for output file names. The output files will be named '
		'<PREFIX><BARCODE>.fasta or <PREFIX><NUMBER>.fasta.')
	args = parser.parse_args()
	n = 0
	barcodes = dict()
	n_barcodes = 1
	for record in FastaReader(args.fasta):
		barcode = record.sequence[:args.barcode_length]
		if barcode in barcodes:
			index = barcodes[barcode]
		else:
			index = n_barcodes
			barcodes[barcode] = index
			n_barcodes += 1
		if args.sequential:
			path = args.prefix + str(index) + '.fasta'
		else:
			path = args.prefix + barcode + '.fasta'
		# Open/close files every time to avoid too many open files
		with open(path, mode='a') as f:
			f.write('>{}\n{}\n'.format(record.name, record.sequence))
		n += 1
	print('Wrote {} sequences to {} files.'.format(n, len(barcodes)))


if __name__ == '__main__':
        main()
