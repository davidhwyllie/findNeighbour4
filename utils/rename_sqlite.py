# renames sqlite files removing any leaving 0- from the filename
import os
import glob

searchdir = "/data/data/pca/cp_subsets_output_400_plus/"
globpath = os.path.join(searchdir, "0-*.sqlite")

for filename in glob.glob(globpath):
    bn = os.path.basename(filename)
    tn = bn.replace("0-202","202")
    print(bn, tn)
    targetfile = os.path.join(searchdir, tn)
    os.rename(filename, targetfile)