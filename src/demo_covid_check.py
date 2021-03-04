""" checks that catwalk and findneighbour stored neighbours are the same
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import pandas as pd
    import Bio
    from collections import Counter
    import requests
    from fn4client import fn4Client

    import sentry_sdk
  
    # instantiate client
    fn4c = fn4Client("http://localhost:5025")      # expects operation on local host; pass baseurl if somewhere else. ## 5023 for covid ## 5027 for TB

    cw_url = "http://localhost:5026"        # 5024 for covid, 5007 for TB
     
    existing_guids = set(fn4c.guids())
    clustering_created = False
    print("There are {0} existing guids".format(len(existing_guids)))

    nSkipped = 0
    nBad = 0
    nGood = 0
    failed = []
    i = 0
    for guid in existing_guids:
        t1 = datetime.datetime.now()
        url = "{0}/neighbours/{1}/{2}".format(cw_url, guid,5)
        r = requests.get(url)
        if r.status_code == 200:
            cw_reply = r.json()
        
        t2 = datetime.datetime.now()
        i1 = t2-t1
        s1 = i1.total_seconds()

        t1 = datetime.datetime.now()            
        fn4_reply = fn4c.guid2neighbours(guid,5,0.5,None)
        #print(guid, "  ---> ", fn4_reply)
        t2 = datetime.datetime.now()
        i1 = t2-t1
        s2 = i1.total_seconds()
        diff =set()
        category='SUCCESS'
        if not len(cw_reply) == len(fn4_reply):
            cw_set = set()
            fn4_set = set()
            for item in cw_reply:
                cw_set.add(item[0])
            for item in fn4_reply:
                fn4_set.add(item[0])
            category='FAILURE'
                
            diff = cw_set ^ fn4_set

        print("Replies cw/fn: timings", s1,s2 ,' for ',guid, 'neighbours', len(cw_reply), len(fn4_reply), len(diff), category)
           

      
    
        