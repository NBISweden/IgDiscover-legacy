#!/usr/bin/env python3
"""
Create a new IgBLAST database.
"""

"""A Python re-write of PrepareDatabaseNew.groovy (part of igblastwrapper).

Original copyright of the groovy version:

 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
"""

import shutil
import os
import subprocess
from argparse import ArgumentParser


def main():
	p = ArgumentParser(description=__doc__)
	p.add_argument("-o", "--output", default=".", help="Output directory (default: current directory)")
	p.add_argument("segments", help="segment descriptions table")
	args = p.parse_args()

	inputFileName = args.segments
	database_path = os.path.join(args.output, 'database')
	jref_path = os.path.join(args.output, 'jref.txt')

	if os.path.exists(database_path):
		shutil.rmtree(database_path)
	os.mkdir(database_path)
	try:
		os.remove("jref.txt")
	except:
		pass

	speciesGeneHash = set()

	speciesAliasMap = {
		"HomoSapiens"         : "human",
		"MusMusculus"         : "mouse",
		"RattusNorvegicus"    : "rat",
		"OryctolagusCuniculus": "rabbit",
		"MacacaMulatta"       : "rhesus_monkey"
	}

	for line in open(inputFileName):
		it = line.split("\t")
		(species, geneFull, segment, segmentFull, refPoint, seq) = it
		seq = seq.strip()

		if 'N' in seq:
			continue
		# Change names to IgBlast semantics
		species = speciesAliasMap.get(species, None)
		segment = segment[0]
		segmentFull = segmentFull.replace("/", "_")  # otherwise makeblastdb will crash
		_prefix = "{database_path}/{species}_{gene1}_{gene2}".format(
			database_path=database_path, species=species, gene1=geneFull[:2],
			gene2=geneFull[2])

		if species:  # TODO???
			majorAllele = segmentFull.endswith("*01")

			if majorAllele:
				refs = [_prefix]
			else:
				refs = []
			refs.append(_prefix + "_all")
			for prefix in refs:
				speciesGeneHash.add(prefix)

				with open("{prefix}_{segment}.fa".format(prefix=prefix, segment=segment), "at") as writer:
					print(">{segmentFull}\n{seq}".format(segmentFull=segmentFull, seq=seq), file=writer)

				if segment == "J":
					with open(jref_path, "a") as writer:
						print(species, geneFull, segmentFull, refPoint, seq, sep="\t", file=writer)
		else:
			pass
			# Ignore those species, we won't be able to use them without framework markup anyway

	# Add dummy references.
	# This is absolutely necessary as IgBlast requires D segment references
	# and is not aware whether a certain chain has D segment or not

	for it in speciesGeneHash:
		fileName = it + "_D.fa"
		if not os.path.exists(fileName):
			with open(fileName, 'wt') as pw:
				print(">.\nGGGGGGGGGGGGGGGG", file=pw)

	for path in os.listdir(database_path):
		subprocess.call(["makeblastdb", "-parse_seqids", "-dbtype", "nucl", "-in", database_path + "/" + path, "-out", database_path + "/" + path[:-3]])


if __name__ == '__main__':
	main()
