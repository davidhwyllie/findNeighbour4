""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import Bio
    import random
    from fn4client import fn4Client


    ### Modify this line to reflect where the fasta files are 
    # define directory where the fastas are
    fastadir = os.path.join('/srv','data','covid')
    fastafile = os.path.join(fastadir, 'elan.consensus.fasta')
    utputfile = os.path.join(fastadir, 'milk_nano.fas')

    nSkipped = 0
    i = 0
    seqs = list()
    guids = list()
    for record in Bio.SeqIO.parse(fastafile, 'fasta'):
        if len(record.seq) == 29903:
            guids.append(record.id)
    print("There are {0}".format(len(guids)))
    exit()
    for record in Bio.SeqIO.parse(fastafile, 'fasta'):

        guid = record.id
        seq = str(record.seq)

        if guid in g2000:
            i = i + 1
            seqs.append(record)

    else:

        nSkipped +=1

    with open(outputfile, "wt") as fw:
        Bio.SeqIO.write(seqs, fw, 'fasta')

    print("Complete; exported ",i)


# ./iqtree2 -s /srv/data/mixfiles/covid/milk_micro.fas -T 6 -m GTR+I+G --redo-tree  ## GTR is faster
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4620419/