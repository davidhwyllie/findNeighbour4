""" 
A component of a findNeighbour4 server which provides relatedness information for bacterial genomes.
It does so using PCA, and supports PCA based cluster generation.

The associated classes compute a variation model for samples in a findNeighbour4 server.
Computation uses data in MongoDb, and is not memory intensive, using configuration information in a 
config file. 


Functionality is provided in following classes:
* VariationModel - stores results of producing variant matrix and running PCA
* VariantMatrix - computes sample x variant matrix (requires: PERSIST object for mongodb access; server configuration file) 
* PCARunner - runs PCA on VariantMatrix

Unit testing is facilitated by a 
* PersistenceTest class.  This exposes a small subset of the fn3persist object's methods, sufficient to test PCA.  It can be used to store subsets of data for testing purposes 
without then need to access a real fn3persistence data store.

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

"""

# import libraries
import os
import logging
import warnings
import datetime
import random
from typing import Tuple, Set
from collections import defaultdict
import pickle
import hashlib

import pandas as pd
import numpy as np
from scipy.stats import poisson
import progressbar
import sqlalchemy

from scipy.stats import binom_test, median_abs_deviation
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


class PersistenceTest:
    """a class which mimics some methods available in an fn3persistence object, sufficient to unit test PCA generation.

    Only these methods are implemented:
    __init__
    refcompressedsequence_guids
    refcompressedsequence_read

    Additionally, a load_data method is provided which loads data into the object.
    Data in the correct format can be generated by utils/temporal_subsets.py
    """

    def __init__(self, **kwargs):
        """ constructs the object.  Any parameters are accepted, and none have any effect """
        self.seqs = {}
        self.sample_ids = set([])

    def load_data(self, sample_ids_file, sequences_file):

        with open(sequences_file, "rb") as f:
            self.seqs = pickle.load(f)
        with open(sample_ids_file, "rb") as f:
            self.sample_ids = set(pickle.load(f))
        # sanity check
        # check there are no samples in sample_ids which are not present in seqs
        # sample_ids are allowed to be a subset of seqs, but
        # no samples should exists in sample_ids which aren't in seqs
        missing = self.sample_ids - set(self.seqs.keys())
        if len(missing) > 0:
            raise KeyError(
                "Provided with sample_ids which are not in seqs.  There are {0} such sequences.  Examples are: {1}".format(
                    len(missing), missing
                )
            )

    def refcompressedsequence_guids(self):
        return self.sample_ids

    def refcompressedsequence_read(self, guid):
        """ read a single sequence """
        if guid not in self.sample_ids:
            return None
        return self.seqs[guid]


class MNStats:
    """ computes the number of M and N bases in a reference compressed object """

    def __init__(self, select_positions, analysed_reference_length):
        """input:
        select_positions: the positions contributing to the pca model, as generated by ModelBuilder.
        analysed_reference_length: the number of reference bases analysed."""
        self.select_positions = select_positions
        self.analysed_reference_length = analysed_reference_length

    def examine(self, obj):
        """examines the reference compressed object obj,
        reporting
        * number of Ns and Ms in the sequence, subdivided by whether they are in
          select_positions
        * "Test 2" (binomial test, as per findNeighbour3 paper) testing whether the frequency of Ns/Ms in the selected_positions exceed those elsewhere,
                   indicative of a mixture."""

        missing = {
            "M_in_model": 0,
            "N_in_model": 0,
            "model_positions": len(self.select_positions),
            "reference_positions": self.analysed_reference_length,
        }

        for base in ["M", "N"]:  # compute missingness

            # record total numbers of N and M per guid
            try:
                missing["{0}_total".format(base)] = len(obj[base])
            except KeyError:
                missing["{0}_total".format(base)] = 0

            # examine all missing (N/M) sites, adding to a missingness model
            try:
                for pos in obj[base]:
                    if pos in self.select_positions:
                        try:
                            missing["{0}_in_model".format(base)] += 1
                        except KeyError:
                            pass

            except KeyError:
                pass  # if there are no M,N then we can ignore these

            ## do binomial test
            not_model = self.analysed_reference_length - len(self.select_positions)
            p_expected = (missing["{0}_total".format(base)] - missing["{0}_in_model".format(base)]) / not_model
            missing["{0}_expected_proportion".format(base)] = p_expected
            p_observed = missing["{0}_in_model".format(base)] / len(self.select_positions)
            missing["{0}_observed_proportion".format(base)] = p_observed
            p_val = binom_test(
                missing["{0}_in_model".format(base)], len(self.select_positions), p_expected, alternative="greater"
            )

            missing["{0}_p_value".format(base)] = p_val

        return missing


class VariationModel:
    """Stores a VariantMatrix, the output of a PCA of the matrix, and (optionally) a clustering of the principal components.
    You should not normally have to call this class directly to create a VariationModel - the VariantMatrix class would do this for you.

    - You might wish to instantiate this class directly if you are restoring a previously serialised VariationModel - see constructor"""

    def __init__(self):
        """
        creates a new Variation model.

        """
        self.model = {"built": False}

        return

    def __getitem__(self, key):
        if key not in self.model:
            raise KeyError(f"Key {key} not found")
        return self.model[key]

    def __setitem__(self, key, value):
        """ adds a key-value pair to the model """
        if key in self.model.keys():
            raise KeyError(f"Cannot replace key {key}")
        else:
            self.model[key] = value

    def _coefficients_hash(self):
        """computes a hash on the coefficients in the variant model.
        This is useful for version tracking & storing patterns of masking."""
        h = hashlib.md5()
        h.update(self.model["eigenvectors"].to_csv().encode("utf-8"))
        md5_l = h.hexdigest()
        return "{0}".format(md5_l)

    def finish(self):
        """ completes construction of the VariationModel """
        self.model["build_time"] = datetime.datetime.now().isoformat()
        self.model["coefficients_hash"] = self._coefficients_hash()
        self.model["built"] = True

    def to_sqlite(self, outputdir="", analysis_name="pca_output", rebuild_databases_if_present=True):
        """write output to sqlite database

        Inputs
        =======
        persistdir                         the directory the SQLite database goes into
        analysis_name                      name of the analysis.  Will become the first part of the file name
        rebuild_databases_if_present       delete any existing SQLite database

        Returns
        =======
        path to sqlite database
        """

        # configure sqlite file for output.
        sqlite_file = os.path.join(outputdir, "{0}.sqlite".format(analysis_name))
        engine = sqlalchemy.create_engine("sqlite:///{0}".format(sqlite_file), echo=True)

        # run checks on sqlite file
        if rebuild_databases_if_present:
            try:
                os.unlink(sqlite_file)
            except FileNotFoundError:
                pass

        # open connection
        conn = engine.connect()

        metadata = []
        for key in self.model:

            if not (
                key == "variant_matrix" or isinstance(self.model[key], PCA)
            ):  # we don't serialise these; one is massive and the other can't be serialised
                logging.info("Writing {0}".format(key))
                native_type = type(self.model[key])
                if native_type in [bool, int, float, str]:
                    metadata.append({"variable": key, "value": str(self.model[key]), "native_type": str(native_type)})
                elif type(self.model[key]) in [set, list]:
                    if type(self.model[key]) == set:
                        list_data = sorted(list(self.model[key]))
                    else:
                        list_data = self.model[key]
                    tmp = pd.DataFrame(list_data, columns=[key])

                    tmp.to_sql(key, conn, if_exists="fail")  # we don't serialise these at present
                elif type(self.model[key]) in [dict, defaultdict]:
                    output_records = []
                    for this_key in self.model[key].keys():
                        item = self.model[key][this_key]
                        if type(item) in [float, bool, int, str]:
                            item = [item]
                        if not type(item) == list:
                            raise TypeError(
                                "Can only export dictionaries which are of key:list or key:scalar format; the list element is of type {0} : {1}".format(
                                    type(item), item
                                )
                            )
                        for list_element in item:
                            output_records.append({key: this_key, "value": list_element})
                    tmp = pd.DataFrame.from_records(output_records)

                elif type(self.model[key]) == np.int64:
                    metadata.append({"variable": key, "value": str(int(self.model[key]))})
                elif type(self.model[key]) == pd.core.frame.DataFrame:
                    self.model[key].to_sql(key, conn, if_exists="fail")
                else:
                    warnings.warn("Not handled {0} with class {1}".format(key, type(self.model[key])))

        metadata_df = pd.DataFrame.from_records(metadata)
        metadata_df.to_sql("Metadata", conn, if_exists="fail")
        conn.close()
        return sqlite_file


class VariantMatrix:
    """In charge of producing a sample x SNP matrix"""

    def __init__(self, CONFIG, PERSIST, show_bar=True):
        """Construct a variant matrix

        Parameters:
        CONFIG: a configuration dictionary, as produced by findn.common_utils.ConfigManager.read_config()
        PERSIST: a persistence object providing access to stored sequence data.
                Either a findn.mongoStore.fn3persistence object, or a PersistenceTest object, the latter being useful for unit testing.
        show_bar: whether or not to show a progress bar

        """
        # store the persistence object as part of the object
        self.PERSIST = PERSIST

        # set easy to read properties from the config
        self.analysed_reference_length = len(CONFIG["reference"]) - len(set(CONFIG["excludePositions"]))

        # store whether we're using bars for display
        self.show_bar = show_bar

        # we start without any variation model
        self._reset()

    def _reset(self):
        """ clears existing variation model and pca result """
        self.vm = VariationModel()
        self._invalid = set()  # invalid guids for which we can't compute pcs
        self.model = {"built": False, "built_with_guids": []}
        self.validation_data = None

    def guids(self):
        """ returns list of guids currently in the findNeighbour4 database"""
        return sorted(
            self.PERSIST.refcompressedsequence_guids()
        )  # sorting is not essential, but makes more deterministic for testing

    def _column_name(self, pos, base):
        """ given a base at a position, returns a position:base string suitable for use as a pandas column name """
        return f"{pos}:{base}"

    def get_position_counts(self, guids=None) -> Tuple[set, dict, dict]:
        """returns positions of variation across the genome

        Parameters:
        guids : a set of sequence identifiers to analyse. If None, all samples in self.PERSIST are analysed

        Returns:
        a tuple consisting of:
        the sample ids (guids) analysed (set)
        a dictionary consisting of the positions of variation, and the numbers of samples at each position; example: {28281: 2692, 5387: 2704, 23603: 2722, 23270: 2708, 6953: 2685, 16175: 2689, 24913: 2707, ...}
        a dictionary consisting of the positions where missingness/gaps (either N, - or IUPAC mixture codes) are present, and the numbers of samples at each position. Format as above"""
        vmodel = defaultdict(int)  # variants
        mmodel = defaultdict(int)  # missingness

        # set default values
        if guids is None:
            guids = self.guids()

        guids_analysed = set()
        if self.show_bar:
            bar = progressbar.ProgressBar(max_value=len(guids))
        for num_loaded, guid in enumerate(guids):
            if self.show_bar:
                bar.update(num_loaded)
            refcompressed_sample = self.PERSIST.refcompressedsequence_read(guid)  # ref compressed sequence

            if refcompressed_sample["invalid"] != 1:
                guids_analysed.add(guid)
                # for definite calls, compute variation at each position
                for base in ["A", "C", "G", "T"]:
                    var_positions = refcompressed_sample.get(base, [])
                    for var_pos in var_positions:
                        vmodel[var_pos] += 1
                # compute missingness/gaps if it's mixed (M) or N
                for base in ["M", "N"]:  # compute missingness/gaps if it's mixed or N
                    missingness_positions = refcompressed_sample.get(base, [])
                    for missingness_pos in missingness_positions:
                        mmodel[missingness_pos] += 1

        if self.show_bar:
            bar.finish()
        return guids_analysed, vmodel, mmodel

    def get_missingness_cutoff(self, positions: Set[int], mmodel: dict) -> int:
        """computes a missingness cutoff, applicable at a per-sequence level.

        Samples which have high levels of missingness (i.e. N, -, or IUPAC mixture codes)
        may be unsuitable for incorporation into PCA models.  Some level of missingness is expected,
        but samples with very high missingness may compromise modelling.

        Parameters:
        positions: a set of integers, representing the positions across which missingness is to be estimated
        mmodel:    a missingness model, which is a dictionary of the form {28281: 2692, 5387: 2704, 23603: 2722}
                   where 28281  is a position, and 2692 is the number of sequences with missingness at that position.
                   this kind of dictionary is generated by .get_position_counts

        Returns:
                   a estimate of how many missing positions are unexpected.  This is based on the idea that the number of missing positions
                   for the majority of samples is estimated by Poisson process; the cutoff returned approximates the 99.9% upper confidence
                   interval on the expected number of missing positions if this is true, and mu = 2 * the median missingness.


        Note:      This criterion is somewhat arbitrary.
                   The impact of this approximation has not been systematically evaluated, and could be the subject of further work.
        """

        missingness = list(map(lambda pos: mmodel.get(pos, 0), positions))
        median_missingness = np.median(
            missingness
        )  # study median, as this will be relatively insensitive to samples with high missingness
        # missingness_distribution = Counter(missingness)

        # use Poisson cdf with mu = 2x median_missingness; find 99% CI
        upper_cutoff = poisson.ppf(0.999, 2 * median_missingness)  # crude approximation to upper CI
        return upper_cutoff

    def build(self, min_variant_freq=None, num_train_on=None, deterministic=True):
        """
        input:
            min_variant_freq: the minimum proportion of samples with variation at that site for the site to be included.  If none, is set to 3/train_on, i.e. each variant has to appear 3 times to be considered
            num_train_on: only compute PCA on a subset of train_on samples.  Set to None for all samples.
            deterministic:  if num_train on is not None, setting deterministic = True (default) ensures the same samples are analysed each time.  If num_train_on is None, has no effect.
        """
        # determine guids there in the database
        guids = self.guids()

        ########################################################################################################
        # if we have been told to use a subset of these to build the model, construct that subset.
        if num_train_on is None:  # if we are not told how many to use then we use
            num_train_on = len(guids)  # all samples
        else:
            # randomise order for model training purposes if required
            if not deterministic:
                random.shuffle(guids)  # makes pipeline non-deterministic if not all samples are analysed
            else:
                guids = sorted(guids)  # keep order constant between runs
        if num_train_on < len(guids):
            guids = guids[:num_train_on]
        self.vm["num_train_on"] = num_train_on  # persist parameters used
        #########################################################################################################

        #########################################################################################################
        # if minimum variation is not set, only analyse variants seen at least 3 times.
        if min_variant_freq is None:
            min_variant_freq = 3 / num_train_on
        #########################################################################################################

        #########################################################################################################
        # determine the variation model.  In the first stage, we analyse by position
        # positions with unexpectedly high levels of missingness are excluded, as these may be hard to call.
        logging.info("Assessing per-base variation from {0} samples".format(num_train_on))
        guids_analysed_stage1, vmodel, mmodel = self.get_position_counts(guids)

        # store variant model
        self.vm["variant_frequencies"] = vmodel
        self.vm["min_variant_freq"] = min_variant_freq
        self.vm["analysed_reference_length"] = self.analysed_reference_length

        # from the totality of the variation, select positions with > cutoff % variation
        cutoff_variant_number = num_train_on * min_variant_freq
        self.vm["cutoff_variant_number"] = cutoff_variant_number

        select_positions = set()
        for pos, variant_count in vmodel.items():
            if variant_count >= cutoff_variant_number:
                select_positions.add(pos)
        logging.info(
            "Found {0} positions which vary at frequencies more than {1}.".format(len(select_positions), min_variant_freq)
        )

        if len(select_positions) == 0:
            raise ValueError(
                "No variation found above cutoff. normally this is because you ran the PCA operation against an empty database; this is what happens if you omit a config file parameter, when a test database is examined by default.   Cannot continue"
            )

        self.model["variant_positions_gt_cutoff_variant_number"] = len(select_positions)  # store result
        upper_cutoff = self.get_missingness_cutoff(select_positions, mmodel)  # define cutoff
        self.vm["max_ok_missingness"] = float(upper_cutoff)
        self.vm["max_ok_missingness_pc"] = int(100 * upper_cutoff / num_train_on)

        # remove any positions with high levels of missingness from the variation model to be used in the PCA
        num_removed = 0
        for pos in mmodel.keys():  # positions
            if mmodel[pos] > upper_cutoff and pos in select_positions:
                num_removed += 1
                select_positions.remove(pos)
        self.select_positions = select_positions
        pc_removed = int(100 * (num_removed / (len(select_positions) + num_removed)))
        logging.info(
            f"Removed {num_removed} positions with missingness > cutoff {upper_cutoff} sequences.  Represents removal of {pc_removed} %"
        )
        logging.info(
            f"There remain {len(select_positions)} positions which vary at frequencies more than {min_variant_freq} and pass the missingness cutoff."
        )
        self.vm["variant_positions_ok_missingness"] = len(select_positions)
        #########################################################################################################

        #########################################################################################################
        # determine the variation model.  In the second stage, we analyse by sample.
        # find any samples which have a high number of missing or uncertain bases in the positions of variation
        # vs in other positions.  Such samples may be mixed.
        self.mns = MNStats(select_positions, self.vm.model["analysed_reference_length"])
        logging.info(
            "Scanning for samples with unexpectedly high missingness (N), or likely to be mixed (M) in variation model"
        )
        if self.show_bar:
            bar = progressbar.ProgressBar(max_value=len(guids_analysed_stage1))

        guid2missing = {}
        for nLoaded, guid in enumerate(guids_analysed_stage1):
            if self.show_bar:
                bar.update(nLoaded)
            obj = self.PERSIST.refcompressedsequence_read(guid)  # ref compressed sequence
            for base in ["M", "N"]:  # compute how many bases in this position are either M or N
                # examine all missing (N/M) sites, adding to a missingness model
                try:
                    for pos in obj[base]:
                        if pos in select_positions:
                            try:
                                mmodel[pos] = mmodel[pos] + 1
                            except KeyError:
                                if pos not in vmodel.keys():
                                    mmodel[pos] = 1  # first occurrence at this position
                except KeyError:
                    pass  # if there are no M,N then we can ignore these

            ## do binomial test for unexpectedly high missingness in the variant sites, as well as an n/m count
            guid2missing[guid] = self.mns.examine(obj)
        if self.show_bar:
            bar.finish()

        #########################################################################################################
        # reject samples with higher missingness in the variation model vs. other sites.
        logging.info("Scan complete.  Collating information.")

        # collate mixture quality information, and identify low quality (mixed, \
        # as judged by high Ns or Ms in the variant sites)
        mix_quality_info = pd.DataFrame.from_dict(guid2missing, orient="index")
        self.vm["mix_quality_info"] = mix_quality_info

        # identify any mixed samples.  we don't build the model from these.
        # mixed are defined as having significantly more N or M in the variant
        # positions than in other bases.
        # Note: this impact of this step has been evaluated in TB, but not as extensively in SARS-COV-2
        mix_quality_cutoff = 0.01 / len(
            mix_quality_info.index
        )  # 0.01 divided by the number of samples analysed;  Bonferroni adj.
        suspect_quality = mix_quality_info.query("M_p_value < {0} or N_p_value < {0}".format(mix_quality_cutoff))
        self.vm["suspect_quality_seqs"] = suspect_quality
        n_suspect = len(set(suspect_quality.index.to_list()))
        pc_suspect = int(100 * n_suspect / len(mix_quality_info.index))
        logging.info(
            "Identified {0} ({1}%) sequences with higher N/M in variant vs. non-variant bases (composition p cutoff {2}); excluded from model as may be mixed.".format(
                n_suspect, pc_suspect, mix_quality_cutoff
            )
        )
        guids_analysed_stage2 = guids_analysed_stage1 - set(suspect_quality.index.to_list())

        #########################################################################################################
        # build a variation matrix for variant sites which pass, and samples which pass
        vmodel = {}
        nLoaded = 0
        logging.info(
            "Gathering variation for matrix construction from {0} unmixed samples into dictionary ".format(
                len(guids_analysed_stage2)
            )
        )

        if self.show_bar:
            bar = progressbar.ProgressBar(max_value=len(guids_analysed_stage2))

        self.model["built_with_guids"] = []
        for guid in guids_analysed_stage2:
            nLoaded += 1

            if self.show_bar:
                bar.update(nLoaded)

            obj = self.PERSIST.refcompressedsequence_read(guid)  # ref compressed sequence
            # for definite calls, compute variation at each position

            # for invalid samples, compute nothing
            if obj["invalid"] == 1:
                self._invalid.add(guid)
            else:
                # compute a variation model - a list of bases and variants where variation occurs
                variants = {}  # variation for this guid

                # for definite calls, compute variation at each position
                # positions of variation where a call was made
                for base in set(["A", "C", "G", "T"]).intersection(obj.keys()):
                    target_positions = select_positions.intersection(obj[base])
                    called_positions = dict((self._column_name(pos, base), 1) for pos in target_positions)

                    variants = {**variants, **called_positions}

                vmodel[guid] = variants
        if self.show_bar:
            bar.finish()

        #########################################################################################################
        # build a variation matrix for variant sites using pandas - may take several minutes for giant matrices
        logging.info("Building variant matrix as pandas DataFrame.  May take several minutes for huge matrices.")
        vmodel = pd.DataFrame.from_dict(vmodel, orient="index")
        vmodel.fillna(value=0, inplace=True)  # if not completed, then it's reference
        # unless it's null, which we are ignoring at present- we have preselected sites as having low null frequencies
        logging.info("Matrix construction complete.  There are {0} sequences in the variation model".format(len(vmodel.index)))
        self.vm["variant_matrix"] = vmodel
        return None


class PCARunner:
    """ Performs PCA on a VariantMatrix """

    def __init__(self, snp_matrix: VariantMatrix):
        self.vm = snp_matrix.vm
        self.eigenvalues = None

    def run(self, n_components, pca_parameters={}, deterministic=True) -> VariationModel:
        """conducts pca on a snp_matrix, storing the results in the snp_matrix's VariantModel object.

        input:
            n_components: the maximum number of components to extract.
            pca_parameters: a dictionary of parameters passed to the scikit-learn PCA command
                            The contents of the dictionary are passed as-is to the PCA command, without any checking.
        """
        logging.info("Performing pca, extracting {0} components".format(n_components))
        self.vm["n_pca_components"] = n_components

        # if necessary, can perform incremental PCA see https://stackoverflow.com/questions/31428581/incremental-pca-on-big-data
        pca = PCA(n_components=n_components, **pca_parameters)
        variant_matrix = self.vm["variant_matrix"]
        pca.fit(variant_matrix)
        contribs = []
        # summarise the positions and variants responsible for each pc
        pc2contributing_pos = {}
        contributing_basepos = set()
        contributing_pos = set()
        for i, row in enumerate(pca.components_, 0):
            # mark values far from the median, which is close to zero
            # this information is not used, but it is retained for depiction purposes if needed
            row_median = np.median(row)
            row_mad = median_abs_deviation(row)
            row_upper_ci = row_median + 3 * row_mad
            row_lower_ci = row_median - 3 * row_mad

            pc2contributing_pos[i] = set()
            for j, cell in enumerate(row, 0):
                if cell > row_upper_ci or cell < row_lower_ci:
                    outside_3mad = True
                else:
                    outside_3mad = False
                pos = int(variant_matrix.columns[j].split(":")[0])
                allele = variant_matrix.columns[j].split(":")[1]

                # indicate whether positions are strongly weighted
                if outside_3mad:
                    pc2contributing_pos[i].add(pos)
                    contributing_basepos.add(variant_matrix.columns[j])
                    contributing_pos.add(pos)
                contribs.append(
                    {
                        "pc": i,
                        "pos": pos,
                        "allele": allele,
                        "col": variant_matrix.columns[j],
                        "weight": cell,
                        "outside_3mad": outside_3mad,
                    }
                )
            pc2contributing_pos[i] = sorted(list(pc2contributing_pos[i]))  # can be json serialised, unlike set

        # report eigenvectors which are different from median +- 3 median absolute deviations
        self.eigenvectors = pd.DataFrame.from_records(contribs)

        if len(self.eigenvectors.index) == 0:
            raise KeyError(
                "PCA problem.  No eigenvectors found.  Contributions found are as follows: {0}.  This usually means there is insufficient data to build PCs.  Try increasing the sample number".format(
                    contribs
                )
            )

        # compute eigenvalues for the samples on which the fit was performed.
        logging.info("Computing eigenvalues")
        eigenvalues_dict = {}
        for guid, evs in zip(variant_matrix.index, pca.transform(variant_matrix)):
            eigenvalues_dict[guid] = evs
        self.eigenvalues = pd.DataFrame.from_dict(eigenvalues_dict, orient="index")
        self.eigenvalues.columns = range(n_components)

        self.vm["pca"] = pca
        self.vm["eigenvalues"] = self.eigenvalues
        self.vm["eigenvectors"] = self.eigenvectors
        self.vm["explained_variance_ratio"] = list(pca.explained_variance_ratio_)
        self.vm["n_contributing_positions"] = len(contributing_pos)
        self.vm["pc2_contributing_positions"] = pc2contributing_pos
        self.vm["n_contributing_variants"] = len(contributing_basepos)
        self.vm["contributing_basepos"] = contributing_basepos
        self.vm["contributing_pos"] = contributing_pos
        self.vm["built_with_guids"] = variant_matrix.index.tolist()
        self.vm["pos_per_pc"] = [len(x) for x in self.vm.model["pc2_contributing_positions"].values()]

        print("PCA completed, identified {0} strongly contributing base/positions".format(len(contributing_basepos)))
        self.vm.finish()

        return self.vm

    def cluster(self, initial_cats_per_unit=8):
        """clusters the eigenvalues obtained by run()

        Categorises eigenvalues.  Uses kmeans clustering, and uses a crude approximation to estimate the number of clusters.

        The technique used operates per pc; we bin eigenvalues into bins 1/initial_cats_per_unit wide, and count the non-zero bins.  This is used as an estimate of
        the numbers of clusters, and the pca is provided with the bin centres as a set of starting values.

        Empirically, the technique was found to provide better discrimination of emerging SARS-CoV-2 genomes than approaches based on model fitting,
        such as Gaussian mixture modelling, although it undoubtedly splits some PCs arbitrarily.

        Parameters:
            initial_cats_per_unit (default 8).  Used to estimate the number of k-means clusters to generate.

        Outputs:
            None

        Sets:
            self.eigenvalue_categories: a data frame containing cluster names for each cluster
            self.eigenvalue_category_meta: a data frame containging centroids of each cluster and other metadata

        """

        # check there is a model
        if self.eigenvalues is None:
            raise NotImplementedError("No eigenvalues.  You must call .run() before calling .cluster()")

        # prepare data for clustering

        ev = self.eigenvalues.copy()  # eigenvalues.  option to drop PCs of technical origin could be dropped.

        for col in ev.columns:
            logging.info("Clustering eigenvalue {0}".format(col))
            this_ev = ev[col].to_frame()
            this_ev.columns = ["eigenvalue"]
            this_ev["pc"] = col
            this_ev["initial_cat"] = [int(x * initial_cats_per_unit) for x in this_ev["eigenvalue"]]

            # how many non-zero categories
            cats = this_ev.groupby(["initial_cat"])["eigenvalue"].describe()

            # convert to arrays to fit
            to_fit = this_ev["eigenvalue"].to_numpy().reshape(-1, 1)
            centres = cats["mean"].to_numpy().reshape(-1, 1)
            km = KMeans(n_clusters=len(cats.index), n_init=1, init=centres).fit(to_fit)
            this_ev["cat"] = km.labels_
            this_ev["sample_id"] = this_ev.index
            cats = this_ev.groupby(["cat"])["eigenvalue"].describe()
            this_ev.drop(["initial_cat"], axis=1)
            if col == 0:
                evs = this_ev
                categories = cats
            else:
                evs = evs.append(this_ev, ignore_index=True)
                categories = categories.append(cats, ignore_index=True)

        self.vm["eigenvalue_categories"] = evs
        self.vm["eigenvalue_category_meta"] = categories
        self.vm.finish()
        return self.vm
