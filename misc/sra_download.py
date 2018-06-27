#!/usr/bin/env python3
import sys
import subprocess


def download_sra(sra_id):
	if not sra_id.startswith('SRR') and not sra_id.startswith('ERR'):
		raise ValueError('only ids starting with SRR or ERR supported')
	number = sra_id[3:]

	if len(number) == 6:
		url = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{sra_id[:6]}/{sra_id}/'
	elif len(number) == 7:
		url = f'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{sra_id[:6]}/00{number[6]}/{sra_id}/'
	else:
		raise ValueError('number not supported')

	single_end_url = f'{url}{sra_id}.fastq.gz'
	paired_end_urls = [f'{url}{sra_id}_{i}.fastq.gz' for i in (1, 2)]

	result = subprocess.run(['wget', single_end_url])
	if result.returncode != 0:
		for url in paired_end_urls:
			subprocess.run(['wget', url], check=True)


if __name__ == '__main__':
	download_sra(sys.argv[1])
