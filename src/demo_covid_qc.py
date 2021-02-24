""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import pandas as pd
    from Bio import SeqIO
    from collections import Counter
    import sentry_sdk

    CONFIG={}
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")

    else:
        print("No sentry key provided")
       
    ### Modify this line to reflect where the fasta files are 
    # define directory where the fastas are
    fastadir = os.path.join('/srv','data','covid')
    fastafile = os.path.join(fastadir, 'naive_msa.fasta')

    results = []
    i = 0
    for record in SeqIO.parse(fastafile, 'fasta'):
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
        #print(res)
        results.append(res)  
        if i % 5000 == 0:
            print("Analysed", i)

    targetfile1 = os.path.join(fastadir, 'naive_msa.consensus.fasta.csv')
    df = pd.DataFrame.from_records(results)
    df.to_csv(targetfile1)
   
    cnts = df.groupby(['seqlen']).size().reset_index(name='counts')
    targetfile2 = os.path.join(fastadir, 'naive_msa.consensus.seqlen_count.csv')
    cnts.to_csv(targetfile2)
    print("Complete.  Details of {0} sequences are in {1} and {2}".format(i, targetfile1, targetfile2))
    
