""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from collections import Counter
    from fn4client import fn4Client

    # instantiate client
    fn4c = fn4Client("http://findneighbours04.unix.phe.gov.uk:5023")      # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    print("There are {0} existing guids".format(len(existing_guids)))

    my_special_samples = set(['ALDP-12ABA31','ALDP-12AA63E'])           # pretend these are the samples of interest

    missing_special_samples = my_special_samples - existing_guids
    if not len(missing_special_samples) == 0:
        print("Not all of my special samples are present; missing are: {0}".format(missing_special_samples))

    # get the neighbours of my_special_samples
    for_msa = my_special_samples.copy()

    for my_special_sample in my_special_samples:
        res = fn4c.guid2neighbours(my_special_sample, threshold = 2)        # find neighbours within 2 SNV
        for related_sample, distance in res:
                for_msa.add(related_sample)

    # OTPIONAL : if you want an outgroup / ancestor to root your tree with, there is a special sample called --Wuhan-Reference--
    for_msa.add('--Wuhan-Reference--')

    # build an MSA
    # various other kinds of output are possible including
    #    json 
    #    json-records
    #    html
    #    json-fasta
    #    fasta
        
    # note that the for_msa call can return
    print("Building MSA with {0} sequences.".format(len(for_msa)))
    msa_df = fn4c.msa(for_msa, output_format='json-records', what='N_or_M')
  
    # export to excel
    excel_outputfile = 'msa.xlsx' 
    msa_df.to_excel(excel_outputfile)

    # export to fasta
    msa_creation_timestamp = datetime.datetime.now().isoformat()
    seqs = []
    for ix in msa_df.index:
        guid = msa_df.at[ix,'guid']
        seq = msa_df.at[ix, 'aligned_mseq']         # should be aligned_seq, but looks like aligned_seq and aligned_mseq outputs are the wrong way round
        sr = SeqRecord(     Seq(seq), 
                            id= guid,
                            description="| Variant sites only shown within MSA created {0}".format(msa_creation_timestamp))
        seqs.append(sr)

    fasta_outputfile =  'msa.fasta' 
    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
            

    print("Complete.  Output is in {0} and {1}".format(fasta_outputfile, excel_outputfile))

 