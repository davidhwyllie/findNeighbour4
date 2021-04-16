""" 
Analyses subsets of SARS-CoV-2 samples which occurred before a particular date
Uses data sets collated by make_temporal_subset
"""

import os
import glob

from pca.pca import PersistenceTest, VariantMatrix, PCARunner
from findn.common_utils import ConfigManager

# read configuration
config_file_name = "demos/covid/covid_config_v3.json"
cfm = ConfigManager(config_file_name)
CONFIG = cfm.read_config()

# read stored data
inputdir = "/data/data/pca/subsets"  # or wherever
seqfile = os.path.join(inputdir, "seqs_20210401.pickle")
globpath = os.path.join(inputdir, "0-*.pickle")
sample_id_files = glob.glob(globpath)
outputdir = "/data/data/pca/subsets_output"  # or wherever
# for each Sunday (sequence file begins with 0_), run pca.
for i, sample_id_file in enumerate(sample_id_files):

    analysis_name = os.path.basename(sample_id_file).replace('.pickle', '')
    print(i, sample_id_file, analysis_name)

    if not os.path.exists(os.path.join(outputdir, analysis_name)):
        TPERSIST = PersistenceTest()
        TPERSIST.load_data(
            sample_id_file,
            seqfile
        )

        v = VariantMatrix(CONFIG, TPERSIST, show_bar=True)
        v.build()       # build matrix
        pcr = PCARunner(v, show_bar=True)  # run pca
        pcr.run(n_components=200, pca_parameters={})    # extract pcs
        v = pcr.cluster()   # cluster pcs
        v.to_sqlite(outputdir = outputdir, analysis_name = analysis_name)

    else:
        print("Skipped, output exists")

print("Finished")
