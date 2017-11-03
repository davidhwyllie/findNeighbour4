#!/usr/bin/env python3

# calibrate.py

""" assists with calibration of MAXN_STORAGE parameter,
which is the most important parameter regarding server performance.
"""

import argparse
import os
import glob
import gzip
from Bio import SeqIO
import numpy as np
import scipy.special

from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Calibrate MAXN_STORAGE by analysing the distribution of Ns in fasta files.  Example:  python3 calibrate.py /home/dwyllie/dev/findNeighbour2/reference/TB-exclude.txt "/home/dwyllie/data/relatednesstest/TB_FASTA_MINIMASKTEST/fastas/*fasta.gz"')
	parser.add_argument('EXCLUDEFILE', help='A file containing a list of bases to be excluded.')
	parser.add_argument('GLOBPATH', help='A directory containing FASTA files to examine')
	
	args = parser.parse_args()
	
	excluded = set()
	if os.path.exists(args.EXCLUDEFILE):
			with open(args.EXCLUDEFILE,'rt') as f:
				rows=f.readlines()
				for row in rows:
					excluded.add(int(row))
			print("Excluded {0} positions.".format(len(excluded)))


	# check that the files exist
	fastaFiles = glob.glob(args.GLOBPATH)
	if len(fastaFiles)==0:
		print("Scan finished.  No files found on {0}".format(CONFIG['GLOBPATH']))
	else:
		print("Initiating scan. Found {0} files ".format(len(fastaFiles)))
	
	nTested = 0
	residual_ns = set()
	
	for fastaFileGz in fastaFiles:
	
		#Convert fasta.gz to fasta to a tmp file
		fastaFile=fastaFileGz
		fastaFile=fastaFile.replace('.gz','')
		fo=open(fastaFile,'wb')
		# it appears that BioPython can't cope with reading the gzip file on the fly
		with gzip.open(fastaFileGz,'rb') as fi:
			fileContent=fi.read()
			fo.write(fileContent)     # so we decompress it
		fo.close()

		nInFile = 0
		with open(fastaFile, 'rt') as f:
			for seq_record in SeqIO.parse(f, 'fasta'):
				guid = str(os.path.basename(fastaFile)[0:36])

				nTested += 1
				nInFile += 1
				if nInFile > 1:
					raise ValueError("Multi-fasta file detected.  These are not supported.  Filename = {0}".format(fastaFile))
				
				seq = str(seq_record.seq)
				ns = seq.count('N')
				seql = list(seq)
				
				# subtract the ns from the excluded position
				for i in excluded:
					if seql[i]=='N':
						ns = ns -1		# doesn't count
				residual_ns.add(ns)		# add this entry
				
		print(nTested, '/', len(fastaFiles), guid, ns)

		#Delete the tmp file
		os.remove(fastaFile)
		residual_ns.add(ns)
		hist, edges = np.histogram(list(residual_ns), density=True, bins='fd')


	print("Scan finished.  Added {0}".format(nTested))
	p1 = figure(
		title="Histogram of residuals Ns once exclusions applied from {0} nFiles".format(len(fastaFiles)),
		tools="save")
	p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="#036564", line_color="#033649")
	p1.xaxis.axis_label = 'Number of residual Ns (as required by MAXN_STORAGE)'
	p1.yaxis.axis_label = 'Pr(x)'
	show(p1)	
	
	print("Scan finished.  Added {0}.  Histogram 'calibrate.html' should be in current working directory.".format(nTested))