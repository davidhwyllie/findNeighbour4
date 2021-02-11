#!/usr/bin/env python3
""" removes fasta deflines from a set of fasta files, replacing them with a code number.

Some fasta deflines in the test data sets we provide include laboratory identifiers.
Under GDPR, these are considered patient identifiable; therefore,
this utility provides a new version of the fasta files, gzipped, with the defline replaced by the
first 36 characters of the filename, which are assumed to be a guid.

"""

import os
import glob
import gzip
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
     
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Produce gzipped fasta files with sanitised deflines
inputdir    inputdirectory in which .fasta files are sought
outputdir   output directory into which .fasta files are written.
""")
    parser.add_argument('inputdir', type=str, nargs=1,
                        help='Input directory')
    parser.add_argument('outputdir', type=str, nargs=1, 
                        help='Output directory')
    
    args = parser.parse_args()
    inputdir = os.path.abspath(args.inputdir[0])
    outputdir = os.path.abspath(args.outputdir[0])
  
    if inputdir == outputdir:
        print("inputdir and outputdir must be different; both are {0}".format(inputdir))
        
    inputfiles = glob.glob(os.path.join(inputdir, "*.fasta"))
    if len(inputfiles)==0:
        print("no .fasta files found in {0}".format(inputdir))
        exit()
        
    # if the output directory does not exist, create it
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        print("Created output directory {0}".format(outputdir))
        
    for inputfile in inputfiles:
        guid = os.path.basename(inputfile)[0:36]
        outputfile = os.path.join(outputdir, "{0}.fasta".format(guid))
        if not os.path.exists(outputfile):
            # don't overwrite
            with open(inputfile, 'rt') as f:
                    for record in SeqIO.parse(f,'fasta'):               
                        record.id = guid
                        record.name = guid
                        print("Reheadering",guid)
                        record.description = "Test fasta file derived from consensus base calling of TB id = {0}".format(guid)
                        with open(outputfile, 'wt') as fo:
                            SeqIO.write(record, fo, 'fasta')
        