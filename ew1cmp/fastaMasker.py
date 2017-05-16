#!/usr/bin/env python

import os
import unittest
import csv
import sys
import gzip
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

class fastaMasker():
    """ turns particular regions within multi-fasta files to 'N'
    
    example usage #1: you specify the proscribed genes in code
    ================
    fr=fastaMasker()                                # instantiate
                                                    # by default, it looks for a file 'refgenome.txt' in the current directory.
                                                    # you can override this by supplying a gene2positionFile argument to the constructor
                                                    # fr=fastaMasker(gene2positionFile='/some/other/file.txt')
                                                    # the gene2position file is a tab delimited file with a particular structure.  The constructor checks the file is valid.
                                                    
                                                   
                                                    
    fr.proscribeGenes(['rrl'])                      # mask the rrl gene     - not ncessary if a proscribedGenesFile is provided containing the relevant gene(s)
    fr.mask(inputFastaFile='NC_000962.fasta')       # in this sequence
    fr.write(outputFile='test.fasta')               # writing masked output to this file.


    example usage #2: you specify the proscribed genes in a file which is loaded on startup
    =================
    fr=fastaMasker(proscribedGenesFile='proscribedGenes.txt')
                                                    # instantiate
                                                    # by default, it looks for a file 'refgenome.txt' in the current directory.
                                                    # you can override this by supplying a gene2positionFile argument to the constructor
                                                    # fr=fastaMasker(gene2positionFile='/some/other/file.txt')
    
                                                    # the gene2position file is a tab delimited file with a particular structure.  The constructor checks the file is valid.   
                                                    # you can also supply a file containing, one row per entry, the genes to be masked.
                                                    # by default, does not expect one.
    
    """
    def __init__(self, gene2positionFile="../reference/refgenome.txt", proscribedGenesFile=None):
        # loads the gene to position table into memory.
        self._proscribedPositions=set()
        self._proscribedGenes=set()
        self._geneLocations={}
        self._gene2positionFile=gene2positionFile
        self._refresh()
        
        # will throw a file not found error if the file isn't present.
        with open(gene2positionFile,'rt') as f:
            reader=csv.DictReader(f, dialect='excel-tab')
            
            # the file read must lookup start position, end position, and gene.
            # it must have headers including 'name','start','end'.
            # we throw an error if it doesn't exist
            headers=set(reader.fieldnames)
            nBases=0
            if len(set(['name','start','end']).intersection(headers))==3:
                for row in reader:
                    nBases=int(row['end'])
                    if not row['name'] in self._geneLocations.keys():
                        self._geneLocations[row['name']]=[int(row['start']),int(row['end'])]
            else:       # the input file is invalid
                raise TypeError('The gene to position lookup table provided to fasta masker is %s.  The headers are: %s.   The file must be tab separated, and have columns including name, start and end reflecting the gene name, start and end' % (self._gene2positionFile, ','.join(str(headers))))
        self._refLength=nBases
        
        ## optionally, the constructor can read a list of proscribed genes from a file, and proscribe them.
        if proscribedGenesFile is not None:
            # will try to read it
            # will throw an error if it doesn't exist.
            with open(proscribedGenesFile,'rt') as f:
                reader=csv.reader(f, dialect='excel-tab')
                for row in reader:
                    if not len(row)==1:
                        raise TypeError('The proscribeGenesFile must have one entry per row.  It does not.  The offending row is :  %s.' % (row))
                    self.proscribeGenes(row)
    def _refresh(self):
        """ clears out any masked records currently stored """
        self.maskedRecords={}       
    def proscribeGenes(self, genes):
        """ specify, by gene name, which genes are to be masked """
        for gene in genes:
            if gene in self._geneLocations.keys():
                self._proscribedGenes.add(gene)
                start=self._geneLocations[gene][0]
                end=self._geneLocations[gene][1]
                for proscribedPosition in range(start,end+1):
                    self._proscribedPositions.add(proscribedPosition)
            else:
                # it's an error: the gene must be defined.
                raise KeyError("The gene %s is not included in the name field in the file %s.  Proscribed genes must be defined." % (gene, self._gene2positionFile))                 
    def export_proscribedPositions(self, outputfile):
        """ write the proscribed positions to an output file """
        with open(outputfile, 'wt') as f:
            for i in sorted(self._proscribedPositions):
                f.write('{0}\n'.format(i))
                
    def mask(self, inputFastaFile):
        """ reads and masks a fasta file.  Stores the results in an internal dictionary."""
        nRead=0
        self._refresh()
        
        # test whether it is a gzip file
        if inputFastaFile.endswith('.gz'):
            f=gzip.open(inputFastaFile)
        else:
            f=open(inputFastaFile,'rt')
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
            nRead+=1
            # mask the record                
            self.maskedRecords[nRead]=self._maskRecord(record)
        f.close()
    def _maskRecord(self, record):
        """ mask a particular SeqIO.seqRecord """
        if not(type(record)==SeqRecord):
            raise TypeError('maskRecord requires a SeqIO.SeqRecord object.')
    
        seqstr=self.maskString(str(record.seq))
        newRecord=record
        newRecord.seq=Seq(seqstr)

        return(newRecord)

    def maskString(self, seq):
        """ mask a string containing the sequence.
        The string, seq,  must be of the same length as the reference.
        """
        if not(type(seq)==str):
            raise TypeError('maskString requires a string as a sequence.  class is: %s To analyse a SeqIO.SeqRecord objects, use _maskRecord(); to analyse a file, use mask()' % type(seq))
          
        seq=list(seq)       # convert to string
        if not len(seq)>=self._refLength:
            raise TypeError( "Sequence has length %i but reference has length %i" % (len(seq),self._refLength))
        
        for proscribedPosition in self._proscribedPositions:
            seq[proscribedPosition-1]='N'
        return(str(''.join(seq)))

    def write(self, outputFile):
        with open(outputFile,'wt') as f:
            for key in self.maskedRecords.keys():
                record=self.maskedRecords[key]
                SeqIO.write(record,f,'fasta')
              
# unit tests
class Test_fastaMasker1(unittest.TestCase):
    """ integration test, checking all functions in fasta masker """
    def runTest(self):
        
        # remove any target file output
        targetfile='../unitTest_tmp/output.txt'
        try:
            os.remove(targetfile)
        except:
            pass    # do nothing

        fr=fastaMasker()
        fr.proscribeGenes(['rrl'])
        fr.mask(inputFastaFile='../reference/NC_000962.fasta')
        fr.write(outputFile=targetfile)

        # generate output
        if not os.path.isfile(targetfile):
            self.fail("target file not created.")

        # check there are no Ns in the inputfile
        with open('../reference/NC_000962.fasta', 'rt') as f:
            nRead=0
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    nRead+=1
        self.assertEqual(1,nRead)       # one item output
        Ns=record.seq.count('N')
        self.assertEqual(0,Ns)
        
        # check there are  Ns in the output file.
        with open(targetfile, 'rt') as f:
            nRead=0
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    nRead+=1
        self.assertEqual(1,nRead)       # one item output
        Ns=record.seq.count('N')
        self.assertEqual(3138,Ns)
class Test_fastaMasker2(unittest.TestCase):
    """ checks construction fails if no gene2positionFile is available"""
    def runTest(self):        
       with  self.assertRaises(IOError):
            fr=fastaMasker(gene2positionFile='None available')
class Test_fastaMasker3(unittest.TestCase):
    """ checks loading of the default (TB) gene list succeeded """
    def runTest(self):      
        fr=fastaMasker()
        self.assertEqual(fr._refLength,4410929)
class Test_fastaMasker4(unittest.TestCase):
    """ checks proscribing genes  """
    def runTest(self):      
        fr=fastaMasker()
        fr.proscribeGenes([])
        self.assertEqual(len(fr._proscribedPositions),0)
class Test_fastaMasker5(unittest.TestCase):
    """ checks proscribing genes  """
    def runTest(self):      
        fr=fastaMasker()
        fr.proscribeGenes(['rrs','rrl','inhA'])
        self.assertEqual(len(fr._proscribedPositions),5485)
class Test_fastaMasker6(unittest.TestCase):
    """ checks proscribing genes  """
    def runTest(self):      
        fr=fastaMasker()
        with self.assertRaises(KeyError):
            fr.proscribeGenes(['not_here'])
class Test_fastaMasker7(unittest.TestCase):
    """ checks proscribing genes and masking the right base  """
    def runTest(self):      
        fr=fastaMasker()
        fr.proscribeGenes(['rrs'])
        # from NCBI entry:
        #  gene            1471846..1473382
        #             /gene="rrs"
        #             /locus_tag="Rvnr01"
        #             /db_xref="GeneID:2700429"
        fr.mask(inputFastaFile='../reference/NC_000962.fasta')
        self.assertEqual(fr.maskedRecords[1].seq[1471845],'N')      # due to zero indexing, this position should be N
        self.assertEqual(fr.maskedRecords[1].seq[1473381],'N')      # due to zero indexing, this position should be N
        self.assertEqual(fr.maskedRecords[1].seq[1473382],'A')      # due to zero indexing, this position should be A
class Test_fastaMasker8(unittest.TestCase):
    """ tests the function mask string  """
    def runTest(self):    

        fr=fastaMasker()
        fr.proscribeGenes(['rrl'])
        with self.assertRaises(TypeError):
            fr.maskString("AAAAA")  # should fail, is not the right length
class Test_fastaMasker9(unittest.TestCase):
    """ tests the function mask string  """
    def runTest(self):    

        fr=fastaMasker()
        fr.proscribeGenes(['rrl'])
        testString=''.join('A'*fr._refLength)
        result=fr.maskString(testString)  # should fail, is not the right length
        Ns=result.count('N')
        self.assertEqual(3138,Ns) 
class Test_fastaMasker10(unittest.TestCase):
    """ tests the use of a proscribedgenes file with the initiator  """
    def runTest(self):
        
        # remove any target file output
        targetfile='../unitTest_tmp/output.txt'
        try:
            os.remove(targetfile)
        except:
            pass    # do nothing

        fr=fastaMasker(proscribedGenesFile='../reference/proscribedGenes.txt')
        fr.mask(inputFastaFile='../reference/NC_000962.fasta')
        fr.write(outputFile=targetfile)

        # generate output
        if not os.path.isfile(targetfile):
            self.fail("target file not created.")

        # check there are no Ns in the inputfile
        with open('../reference/NC_000962.fasta', 'rt') as f:
            nRead=0
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    nRead+=1
        self.assertEqual(1,nRead)       # one item output
        Ns=record.seq.count('N')
        self.assertEqual(0,Ns)
        
        # check there are  Ns in the output file.
        with open(targetfile, 'rt') as f:
            nRead=0
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    nRead+=1
        self.assertEqual(1,nRead)       # one item output
        Ns=record.seq.count('N')
        self.assertEqual(10792,Ns)                                     
# testing
# testing
# run unittests
if __name__ == '__main__':
    logging.getLogger().setLevel(logging.WARNING)

    # make a list of excluded bases
    excluded_genes = ('rrl','rrs','rpoB', 'Rv2082')
    fr=fastaMasker(proscribedGenesFile='../reference/proscribedGenes.txt')
    fr.proscribeGenes(excluded_genes)
    fr.export_proscribedPositions(outputfile= "proscribed.txt")
    
    unittest.main()            ## test everything
    

    
