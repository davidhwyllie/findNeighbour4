""" scans a COMPASS maskfile and extracts the 'always N' bases """

import glob
import sys
import csv
from fastaMasker import fastaMasker


if __name__=='__main__':
	
	print("Loading positions from fastamasker")
	# now add the positions from the fastaMasker
	fr=fastaMasker(gene2positionFile="../reference/refgenome.txt")    # load positions of the genes
	fr.proscribeGenes(['rrl', 'rrs', 'rpoB', 'Rv2082'])
	Ns = fr._proscribedPositions
	
	print("Loading positions from rep array")
	ix = 0
	filename = "../COMPASS_reference/R39/R00000039_repmask.array.orig"
	with open(filename) as f:
		reader = csv.reader(f)
		for row in reader:
			if row == ['1']:
				Ns.add(ix)
			ix += 1
				
	print("There are {0} masked bases.  Writing to file ...".format(len(Ns)))


	outputfile = "../reference/TB-exclude.txt"
	with open(outputfile, 'wt') as f:
		for item in sorted(Ns):
			f.write("{0}\n".format(item))
			
	

	