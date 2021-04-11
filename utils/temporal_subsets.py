""" generates lists of SARS-CoV-2 samples which occurred before a particular date
Also generates a dictionary of reference compressed sequences
And a subset of these

Together, these can be passed to a ram_persistence object which 
can be used instead of an fn3persistence object to test the performance of PCA, or for other 
unit testing purposes.

Also useful for investigating how PCA detected the ingress of new strains over time.
Uses public cog metadata downloaded from COG-UK 7/4/2021, saved in
testdata/pca/cog_metadata.csv.gz, and requires access to an fn3persistence object containing the same data.
"""

import os
import pandas as pd
import datetime
import gzip
import pickle
import progressbar
import random
from findn.mongoStore import fn3persistence
from findn.common_utils import ConfigManager

# open connection to existing covid datastore
config_file = os.path.join("..", "demos", "covid", "covid_config_v3.json")
cfm = ConfigManager(config_file)
CONFIG = cfm.read_config()
PERSIST = fn3persistence(dbname=CONFIG["SERVERNAME"], connString=CONFIG["FNPERSISTENCE_CONNSTRING"], debug=CONFIG["DEBUGMODE"])
inputfile = "/data/software/fn4dev/testdata/pca/cog_metadata.csv.gz"
outputdir = "/data/data/pca/subsets"  # or wherever

# read metadata file into pandas
with gzip.open(inputfile, "rt") as f:
    df = pd.read_csv(f)

# we are using the middle part of the cog_id as the sample name as the sample_id; extract this.
sample_ids = df["sequence_name"].to_list()
df["sample_id"] = [x.split("/")[1] for x in sample_ids]

# load a small subset of the reference compressed sequences, for testing purposes
# load the reference compressed sequences
storage_dict = {}
sampled = random.sample(df["sample_id"].to_list(), 5000)
bar = progressbar.ProgressBar(max_value=len(sampled))

for i, sample_id in enumerate(sampled):
    res = PERSIST.refcompressedsequence_read(sample_id)
    bar.update(i)
    storage_dict[sample_id] = res
bar.finish()
# write out the dictionary
outputfile = "/data/software/fn4dev/testdata/pca/seqs_5000test.pickle"
with open(outputfile, "wb") as f:
    pickle.dump(storage_dict, f)
outputfile = "/data/software/fn4dev/testdata/pca/seqs_5000test_ids.pickle"
with open(outputfile, "wb") as f:
    pickle.dump(sampled, f)
exit()
# load the reference compressed sequences
storage_dict = {}
bar = progressbar.ProgressBar(max_value=len(df.index))
for i, sample_id in enumerate(df["sample_id"]):
    res = PERSIST.refcompressedsequence_read(sample_id)
    bar.update(i)
    storage_dict[sample_id] = res
bar.finish()

# write out the dictionary
outputfile = os.path.join(outputdir, "seqs_20210401.pickle.gz")
with gzip.open(outputfile, "wb") as f:
    pickle.dump(storage_dict, f)

# construct counts between 1 June 2020 and end March 2021
cnts = df.groupby(["sample_date"]).size()
cnts = cnts[cnts.index >= "2020-06-01"]
cnts = cnts[cnts.index < "2021-04-01"]
cnts = pd.DataFrame(cnts)
cnts.columns = ["count"]
cnts["dow"] = [datetime.date.fromisoformat(item).weekday() for item in cnts.index]
cnts["isodate"] = [datetime.date.fromisoformat(item) for item in cnts.index]

# write samples to consider in the PCA into a series of json files in the output directory
for cutoff_date in cnts.index:
    dow = cnts.loc[cutoff_date, "dow"]
    df_subset = df[df["sample_date"] < cutoff_date]
    sample_ids = df_subset["sample_id"].to_list()

    outputfile = os.path.join(outputdir, "{0}-{1}.picklen.gz".format(dow, cutoff_date))
    with gzip.open(outputfile, "wb") as f:
        pickle.dump(sample_ids, f)
        print(outputfile)
