""" illustrates use of findNeighbour4 with covid samples,
in particular the generation of multisequence alignments
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':

    import os
    import glob
    import datetime
    import json
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from collections import Counter
    from fn4client import fn4Client

    # instantiate client
    fn4c = fn4Client("http://findneighbours04.unix.phe.gov.uk:5023")      # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    runtime = datetime.datetime.now().isoformat()
    print("Analysis run at {1}; there are {0} existing guids".format(len(existing_guids), runtime))

    # LOAD THE LIST OF SAMPLES YOU WANT AN MSA FROM HERE
    my_special_samples = set(['ALDP-12ABA31','ALDP-12AA63E'])           # pretend these are the samples of interest

    print("Analysing {0} samples, obtaining neighbourhood, and constructing msa.".format(len(my_special_samples)))
    missing_special_samples = my_special_samples - existing_guids
    if not len(missing_special_samples) == 0:
        warning("Not all of my special samples are present in the server; missing are: {0}".format(missing_special_samples))

    # get the neighbours of my_special_samples
    for_msa = my_special_samples.copy()

    for my_special_sample in my_special_samples:
        res = fn4c.guid2neighbours(my_special_sample, threshold = 1)        # find neighbours within 1 SNV
        for related_sample, distance in res['neighbours']:
                for_msa.add(related_sample)

    # build an MSA
    # various other kinds of output are possible including
    #    json           # includes json-records but also positions of variation
    #    json-records
    #    html
    #    json-fasta
    #    fasta
        
    # note that the for_msa call can return
    print("Building MSA with {0} sequences.".format(len(for_msa)))
    msa_json= fn4c.msa(for_msa, output_format='json', what='N_or_M') 
    print(msa_json)
    msa_dict = json.loads(msa_json)

    print(msa_dict.keys())

    exit()

    # to just get the MSA
    msa_df = fn4c.msa(for_msa, output_format='json-records', what='N_or_M') 

    # get all data
    #res = json.loads(fn4c.msa(for_msa, output_format='json', what='N_or_M'))
    
    # export to excel
    excel_outputfile = 'msa.xlsx' 
    msa_df.to_excel(excel_outputfile)

    # export to fasta
    msa_creation_timestamp = datetime.datetime.now().isoformat()
    seqs = []
    for ix in msa_df.index:
        guid = msa_df.at[ix,'guid']
        seq = msa_df.at[ix, 'aligned_seq'].replace('M','N')         # replace non  ACGT with N
        sr = SeqRecord(     Seq(seq), 
                            id= guid,
                            description="" )        # | Variant sites only shown within MSA created {0}".format(msa_creation_timestamp)
        seqs.append(sr)

    fasta_outputfile =  'msa.fasta' 
    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
            
    # test whether an environment variable, IQTREE, is present.  if it is, build the MSA
    if os.environ.get("IQTREE") is not None:
        iqtree_exe = os.environ.get("IQTREE")
        print("Path to iqtree executable found in environment variable, IQTREE")
        cmd = '{0} -m GTR+G -mem 20G -o "{1}" -s {2} --redo'.format(iqtree_exe, outgroup_name, fasta_outputfile)
        print(cmd)
    else:
        print("No environment variable IQTREE found.  If using virtual environment, set it in the .env file.")


    print("Complete.  Output is in {0} and {1}".format(fasta_outputfile, excel_outputfile))
    print("iqtree build command is \n{0}".format(cmd))

 