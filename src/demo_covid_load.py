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
    import sentry_sdk


    CONFIG={}
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")

    else:
        print("No sentry key provided")
   
    sentry_sdk.init(
    "https://9a1bf5a4e715450aa731e34866f611b9@o528497.ingest.sentry.io/5645887"
    
    )
    
    ### Modify this line to reflect where the fasta files are 
    # define directory where the fastas are
    fastadir = os.path.join('/srv','data','covid')
    fastafile = os.path.join(fastadir, 'naive_msa.fasta')

    # instantiate client
    fn4c = fn4Client("http://localhost:5023")      # expects operation on local host; pass baseurl if somewhere else.


    # add the reference sequence as the root
    for record in Bio.SeqIO.parse("../reference/nc_045512.fasta", 'fasta'):
        guid = '--Wuhan-Reference--'
        seq = str(record.seq).upper()
        res = fn4c.insert(guid=guid,seq=seq)

    existing_guids = set(fn4c.guids())
    clustering_created = False
    print("There are {0} existing guids".format(len(existing_guids)))

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

            t2 = datetime.datetime.now()
            i1 = t2-t1
            s = i1.total_seconds()
            print(i, nGood, guid, msg, datetime.datetime.now(),  "in", s, "seconds", len(guid),  len(seq), seq[:10], seq[-10:],counter)
            
        else:
            nBad +=1
            

       
    print("Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(nSkipped, nBad))

    print("Failed samples are:")
    print(failed)
