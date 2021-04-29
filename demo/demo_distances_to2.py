""" illustrates use of findNeighbour4

Illustrates the SNP distribution of all samples in the server.
communicates directly with mongodb; does not use the REST API
"""

if __name__ == "__main__":
    import pandas as pd
    import gzip
    import json
    import progressbar

    from findn.common_utils import ConfigManager
    from findn.mongoStore import fn3persistence

    # instantiate storage class
    config_file = "demos/covid/covid_config_v3.json"
    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()
    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"], connString=CONFIG["FNPERSISTENCE_CONNSTRING"], debug=0
        )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception:
        raise

    existing_guids = PERSIST.guids()
    print("There are {0} existing guids".format(len(existing_guids)))

    inputfile = "/data/software/fn4dev/testdata/pca/cog_metadata.csv.gz"  # COG-UK data
    # read metadata file into pandas
    with gzip.open(inputfile, "rt") as f:
        df = pd.read_csv(f)

    # we are using the middle part of the cog_id as the sample name as the sample_id; extract this.
    sample_ids = df["sequence_name"].to_list()
    df["sample_id"] = [x.split("/")[1] for x in sample_ids]
    # make sample_id the index
    df.set_index("sample_id", inplace=True)
    regions = df["adm1"].unique()
    starting_counts = {}
    for region in regions:
        starting_counts[region] = {0: 0, 1: 0, 2: 0, 3: 0}
    starting_counts["ALL"] = {0: 0, 1: 0, 2: 0, 3: 0}

    print("There are {0} samples in the COG-UK list".format(len(df.index)))

    both = set(df.index).intersection(existing_guids)
    print("There are {0} samples in the server which are also annotated".format(len(both)))
    bar = progressbar.ProgressBar(max_value=len(both))
    all_dists = dict()
    for i, guid in enumerate(sorted(both)):
        bar.update(i)
        neighbours = PERSIST.guid2neighbours(guid, returned_format=1, cutoff=3)
        all_dists[guid] = starting_counts.copy()

        for (neighbouring_guid, distance) in neighbours["neighbours"]:
            if neighbouring_guid in both:
                adm1 = df.at[neighbouring_guid, "adm1"]
                all_dists[guid][adm1][distance] += 1
                all_dists[guid]["ALL"][distance] += 1

        # if i> 200:
        #    break
    outputfile = "/data/data/pca/distribs.json"
    with open(outputfile, "w") as f:
        json.dump(all_dists, f)
    print("Complete. wrote json to ", outputfile)
