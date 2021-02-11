""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import pandas as pd
    import Bio
    from collections import Counter
    from fn4client import fn4Client

    ### Modify this line to reflect where the fasta files are 
    # define directory where the fastas are
    fastadir = os.path.join('/srv','data','mixfiles','covid')
    fastafile = os.path.join(fastadir, 'milk_all.fas')

    # instantiate client
    fn4c = fn4Client("http://localhost:5020")      # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    clustering_created = False
    print("There are {0} existing guids".format(len(existing_guids)))
    nSkipped = 0
    nBad = 0
    i = 0
    for record in Bio.SeqIO.parse(fastafile, 'fasta'):
        i = i + 1
        t1 = datetime.datetime.now()
        guid = record.id
        if not guid in existing_guids:
            seq = str(record.seq)
            counter = Counter(list(seq))

            print(guid, counter, len(seq), seq[:10], seq[-10:])
         
            if len(seq)== 29903:
                res = fn4c.insert(guid=guid,seq=seq)
            else:
                nBad +=1
                
            t2 = datetime.datetime.now()
            i1 = t2-t1
            s = i1.total_seconds()
            
            #print(i, guid, len(seq), res, datetime.datetime.now(),  "in", s, "seconds")

        else:

            nSkipped +=1

    print("Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(nSkipped, nBad))


