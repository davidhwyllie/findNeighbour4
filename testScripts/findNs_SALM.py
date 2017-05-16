""" scans a directory of fasta files and identifies locations at which
there are always Ns. Used for generating 'exclude.txt' files which
require this information. """

import glob
import sys
from Bio import SeqIO
import gzip

if __name__=='__main__':
	
	# command line usage.  Pass the pass to glob as a single argument.
	print("Startup ..")
	globpath = "/mnt/md0/dwyllie/data/relatednesstest/SE_FASTA/*.fa.gz"
	
	fn = glob.glob(globpath)
	print("Found {0} files.".format(len(fn)))
	
	nScanned = 0
	for filename in fn:
		print(filename)
		if filename.endswith('.gz'):
			with gzip.open(filename) as f:
				for r in SeqIO.parse(f,'fasta'):
					s = str(r.seq)
					nScanned += 1
					if nScanned == 1:
						# first one.
						expected_len = len(s)
						notNs = set(list(range(0, expected_len)))
						
					if not expected_len == len(s):
						raise ValueError("Sequence is the wrong length ({0} vs {1})".format(len(s), expected_len))
					
					bases = zip(range(0,expected_len),list(s))
					for (pos, base) in bases:
						if  not base == 'N':
							if pos in notNs:
								notNs.remove(pos)
					
					print(nScanned, len(notNs))
	
		if nScanned > 100:
			break
	outputfile = "/mnt/md0/dwyllie/dev/ELEPHANTWALK2/reference/SALM-exclude.txt"
	with open(outputfile, 'wt') as f:
		for item in sorted(notNs):
			f.write("{0}\n".format(item))