""" scans a directory of fasta files and identifies locations at which
there are always Ns. Used for generating 'exclude.txt' files which
require this information. """

import glob
import sys
from Bio import SeqIO
import gzip

if __name__=='__main__':
	
	print("Startup ..")
	globpath = "/mnt/md0/dwyllie/data/relatednesstest/TB_FASTA/*.fasta"
	
	fn = glob.glob(globpath)
	print("Found {0} files.".format(len(fn)))
	if len(fn) == 0:
		exit()
		
	nScanned = 0
	for filename in fn:
		print(filename)
		if filename.endswith('.gz'):
			f = gzip.open(filename)
		else:
			f= open(filename)
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
					if  not base in ('N','-'):
						if pos in notNs:
							notNs.remove(pos)
			f.close()	
			print(nScanned, len(notNs))

		if nScanned > 1850:
			break
	outputfile = "/mnt/md0/dwyllie/dev/ELEPHANTWALK2/reference/tb-exclude-new.txt"
	with open(outputfile, 'wt') as f:
		for item in sorted(notNs):
			f.write("{0}\n".format(item))