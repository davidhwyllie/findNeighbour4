""" compares output from two different findNeighbour instances.
"""

import os
import glob
import datetime
import pandas as pd
from fn4client import fn4Client

if __name__ == '__main__':
    # instantiate clients; should contain same data
    fn4c_1 = fn4Client("http://localhost:5032")      
    fn4c_2 = fn4Client("http://localhost:5042")      

    print("Server 1 is : ",fn4c_1.server_time())
    print("Server 2 is : ",fn4c_2.server_time())

    guids1 = set(fn4c_1.guids())
    guids2 = set(fn4c_2.guids())

    # are they different?
    different = guids1 ^ guids2
    print("there are {0} differences between samples in server 1 vs 2".format(len(different)))
    same = guids1.intersection(guids2)
    print("there are {0} samples in both server 1 vs 2".format(len(same)))

    # get neighbours of all
    snpcmp_1 = {}
    snpcmp_2 = {}
    distrib_1 = list(range(20+1))
    distrib_2 = list(range(20+1))
    for i,guid in enumerate(same):
        for guid2,dist_1 in fn4c_1.guid2neighbours(guid, threshold=20):
            if guid2 in same:
                snpcmp_1["{0}|{1}".format(guid,guid2)]=dist_1 
                distrib_1[dist_1]+=1

    for i,guid in enumerate(same):
        for guid2,dist_2 in fn4c_2.guid2neighbours(guid, threshold=20):
            if guid2 in same:
               snpcmp_2["{0}|{1}".format(guid,guid2)]=dist_2 
               distrib_2[dist_2]+=1

    print(1,distrib_1)
    print(2,distrib_2)
    if distrib_1==distrib_2:
        print("Distributions are the same")
    else:
        print("Fail: distributions differ")

    if not set(snpcmp_1.keys())==set(snpcmp_2.keys()):
        print("FAIL: pairs identified differ")
    print("Examining {0} pairs, comparing both methods; will report any discrepancies".format(len(snpcmp_1)))
    failures =0
    for key in snpcmp_1.keys():
        try:
            if not snpcmp_1[key] == snpcmp_2[key]:
                print("FAIL: Distances differ for ",key, snpcmp_1[key] , snpcmp_2[key])
                failures +=1
        except KeyError:
                print("FAIL: Pair ",key, "is not present in snpcmp_2")
                failures +=1
    print('Finished, failures = {0}'.format(failures))
    #exit()


    # get clusters
    cl1 = fn4c_1.clusters("SNV12_include")
    cl2 = fn4c_2.clusters("SNV12_include")

    # get guids mentioned in clusters
    g1 = set(cl1['guid'])
    g2 = set(cl2['guid'])
    if not g1==g2:
        print("Guids differ in the comparators {0} vs {1}")
        print(len(g1),len(g2))
    if not len(cl1.index)==len(cl2.index):
        print("clusters differ: see cluster sizes")
        cc1 = cl1['cluster_id'].value_counts()
        cc2 = cl2['cluster_id'].value_counts()
        print("Max cluster size with #1 is {0} : top 5 {1}".format(max(cc1), cc1[:5]))
        print("Max cluster size with #2 is {0} : top 5 {1}".format(max(cc2), cc2[:5]))

    exit()
    print(cl1['cluster_id'].value_counts())
