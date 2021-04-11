""" illustrates use of findNeighbour4 with covid samples,
in particular the generation of multisequence alignments
assumes a findNeighbour4 server is running
"""

if __name__ == "__main__":

    import os
    import datetime
    import pandas as pd
    import warnings
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from fn4client import fn4Client

    ############################################# INPUT #######################################################################
    # LOAD THE LIST OF SAMPLES YOU WANT AN MSA FROM HERE
    my_special_samples = set(["ALDP-12ABA31", "ALDP-12AA63E"])  # pretend these are the samples of interest
    snp_cutoff = 1  # msa of samples within 1 snp of either of the above
    output_file_stem = "ALDP_31_3E"
    ######################################### End of input ####################################################################

    # instantiate client
    fn4c = fn4Client("http://findneighbours04.unix.phe.gov.uk:5023")  # this is the experimental covid-19 server

    existing_guids = set(fn4c.guids())
    runtime = datetime.datetime.now().isoformat()
    print("Analysis run at {1}; there are {0} existing guids".format(len(existing_guids), runtime))

    print("Analysing {0} samples, obtaining neighbourhood, and constructing msa.".format(len(my_special_samples)))
    missing_special_samples = my_special_samples - existing_guids
    if not len(missing_special_samples) == 0:
        warnings.warning(
            "Not all of my special samples are present in the server; missing are: {0}".format(missing_special_samples)
        )

    # get the neighbours of my_special_samples
    for_msa = my_special_samples.copy()

    for my_special_sample in my_special_samples:
        res = fn4c.guid2neighbours(my_special_sample, threshold=1)  # find neighbours within 1 SNV
        for related_sample, distance in res["neighbours"]:
            for_msa.add(related_sample)

    # optional - add outgroup/root
    outgroup_name = "--Wuhan-Reference--"  # the reference sequence, which is in the server
    for_msa.add(outgroup_name)

    # build an MSA
    # various kinds of output are possible including
    #    json           # includes all info in json-records but also positions of variation, constant site numbers (for iqtree) etc
    #    json-records   # a table
    #    html
    #    json-fasta
    #    fasta

    # get a dictionary containing lots of details about the msa
    print("Building MSA with {0} sequences.".format(len(for_msa)))
    msa_dict = fn4c.msa(for_msa, output_format="json", what="N_or_M")

    # returned object has keys as follows
    # variant_positions - which positions the msa items correspond to
    # valid/invalid_guids  - those that pass/fail %N cutoff
    # expected_p1 - ignore, to do with %N in the sequences overall
    # df_dict - convertable to a pandas dataframe
    # what_tested in mixpore computation, N, M, or N_or_M
    # outgroup - if specified
    # fconst - constant sites outside the msa

    print(msa_dict.keys())  # as above

    print(msa_dict["fconst"])  # if needed can construct iqtree command with these
    msa_df = pd.DataFrame.from_dict(msa_dict["df_dict"], orient="index")
    print(msa_df)

    bonferroni_corrected_p = 0.05 / len(msa_df.index)
    print(
        "The following items failed a compositional (MIXPORE) test with more Ns/Ms than expected in the alignment; their positions in the ML tree are likely unstable"
    )
    print(msa_df.loc[msa_df["p_value3"] <= bonferroni_corrected_p])

    # export to excel
    excel_outputfile = "{0}.xlsx".format(output_file_stem)
    msa_df.to_excel(
        excel_outputfile
    )  
    #  study this carefully.  The p-values, if significant, are indicative of various types of mixture.
    # see https://github.com/davidhwyllie/findNeighbour4/blob/master/doc/mixtureTesting.md

    # note: the column aligned_seq may contain M characters if there is an IUPAC code
    # the column aligned_mseq will contain actual IUPAC code
    # in general for iqTree etc you just want Ns - see below

    # export to fasta
    msa_creation_timestamp = datetime.datetime.now().isoformat()
    seqs = []
    for guid in msa_df.index:
        seq = msa_df.at[guid, "aligned_seq"].replace("M", "N")  # replace non  ACGT mixed characters 'M' with N
        sr = SeqRecord(
            Seq(seq), id=guid, description=""
        )  # | Variant sites only shown within MSA created {0}".format(msa_creation_timestamp)
        seqs.append(sr)

    fasta_outputfile = "{0}.fasta".format(output_file_stem)
    with open(fasta_outputfile, "w") as f:
        SeqIO.write(seqs, f, "fasta")

    # test whether an environment variable, IQTREE, is present.  if it is, build the MSA
    if os.environ.get("IQTREE") is not None:
        iqtree_exe = os.environ.get("IQTREE")
        print("Path to iqtree executable found in environment variable, IQTREE")
        run_ok = True
    else:
        iqtree_exe = "<<INSERT PATH TO IQTREE EXECUTABLE HERE>>"
        print("No environment variable IQTREE found.  If using virtual environment, set it in the .env file.")
        run_ok = False

    cmd = '{0} -m GTR+G -mem 20G -o "{1}" -s {2} -fconst {3},{4},{5},{6} --redo'.format(
        iqtree_exe,
        outgroup_name,
        fasta_outputfile,
        msa_dict["fconst"]["A"],
        msa_dict["fconst"]["C"],
        msa_dict["fconst"]["G"],
        msa_dict["fconst"]["T"],
    )

    print("Complete.  Output is in {0} and {1}".format(fasta_outputfile, excel_outputfile))
    print("iqtree build command is \n{0}".format(cmd))
    if run_ok:
        os.system(cmd)
        pass
