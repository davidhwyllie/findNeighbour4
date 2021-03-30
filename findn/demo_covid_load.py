""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import pandas as pd
    import shutil
    import Bio
    import logging
    from collections import Counter
    from fn4client import fn4Client
    import sentry_sdk


  
    ### Modify this line to reflect where the fasta files are 
    # define directory where the fastas are
    fastadir = os.path.join('/srv','data','covid')
    globpath  = "COVID_MSA_*.fasta"
   
    # instantiate client
    fn4c = fn4Client("http://localhost:5023")      # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    clustering_created = False
    print("There are {0} existing guids".format(len(existing_guids)))

    # add the reference sequence as the root if not already present
    ref_guid = '--Wuhan-Reference--'
    ref_guid_present = fn4c.guid_exists(ref_guid)
    if not ref_guid_present:
        print("Adding reference")
        for record in Bio.SeqIO.parse("../reference/nc_045512.fasta", 'fasta'):

            seq = str(record.seq).upper()
            res = fn4c.insert(guid=ref_guid,seq=seq)
    else:
        print("Reference already present")

    for fastafile in glob.glob(os.path.join(fastadir, globpath)):
        print("Scanning", fastafile)

        nSkipped = 0
        nBad = 0
        nGood = 0
        failed = []
        i = 0
        for record in Bio.SeqIO.parse(fastafile, 'fasta'):
            i = i + 1
            t1 = datetime.datetime.now()
            guid = record.id
            guid = guid.replace('/','-')
            guid = guid.replace(':','-')

            seq = str(record.seq).upper()
            seq = seq.replace('?','N')
            seq = seq.replace(' ','N')
            counter = Counter(list(seq))
            
            res = {'record_id':record.id,'guid':guid, 'seqlen':len(seq), **counter}
            if not guid in existing_guids:
                res = fn4c.insert(guid=guid,seq=seq)
                if not res.status_code == 200:
                    # failed to add
                    failed.append((i,guid))
                    msg = "** FAILED **"

                else:
                    msg = "Succeeeded"
                    nGood +=1

                t2 = datetime.datetime.now()
                i1 = t2-t1
                s = i1.total_seconds()
                print(i, nGood, guid, msg, datetime.datetime.now(),  "in", s, "seconds", len(guid),  len(seq), seq[:10], seq[-10:],counter)
                
            else:
                nSkipped +=1

            if i % 1 == 0:
                print("Examined ",i,"skipped",nSkipped)   
            #if i > 5:
            #    break       
        print("Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(nSkipped, nBad))

        print("Failed samples are:")
        print(failed)
