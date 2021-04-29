""" illustrates use of findNeighbour4
assumes a findNeighbour4 server is running on port 5025.
Illustrates the SNP distribution of all samples in the server.
"""

if __name__ == "__main__":
    import pandas as pd
    import gzip
    import progressbar
    from fn4client import fn4Client

    # instantiate client
    fn4c = fn4Client("http://localhost:5025")  # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
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

    print("There are {0} samples in the COG-UK list".format(len(df.index)))

    both = set(df.index).intersection(existing_guids)
    print("There are {0} samples in the server which are also annotated".format(len(both)))
    bar = progressbar.ProgressBar(max_value=len(both))
    all_dists = dict()
    for i, guid in enumerate(sorted(both)):
        bar.update(i)
        neighbours = fn4c.guid2neighbours(guid, threshold=3, timeout=10)

        if False:
            distances = {"ALL": {0: 0, 1: 0, 2: 0, 3: 0}}
            for (neighbouring_guid, distance) in neighbours:
                if neighbouring_guid in both:
                    adm1 = df.at[neighbouring_guid, "adm1"]

                    if adm1 not in distances.keys():
                        distances[adm1] = {0: 0, 1: 0, 2: 0, 3: 0}

                    distances[adm1][distance] = distances[adm1][distance] + 1
                    distances["ALL"][distance] = distances["ALL"][distance] + 1

            # print(distances)
            all_dists[guid] = distances

    print("Complete.")
