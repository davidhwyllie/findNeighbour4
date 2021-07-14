""" contains a class to store the a dictionary of properties,
    as produced by the VariantMatrix class .model property,
    as json.

    This is required because of 
    a) the need to store these dictionaries for unittesting
    b) the insecurity of pickle
    c) the fact that pickled data doesn't unpickle properly if the pickle & pandas versions aren't sync'd
     which cannot be assured in CI environments.
"""

import unittest
import json
import pickle
import pandas as pd
import collections
import datetime


class DictStorage:
    """stores and recovers models as produced by pca.VariantMatrix

    The model is a dictionary with the following keys:

    contains a class to store the a dictionary of properties,
    as produced by the VariantMatrix class .model property,
    as json.

    This is required because of
    a) the need to store these dictionaries for unittesting
    b) the insecurity of pickle
    c) the fact that pickled data doesn't unpickle properly if the pickle & pandas versions aren't sync'd
     which cannot be assured in CI environments.

    An example of the data stored is a dictionary with following key/values:

    built <class 'bool'>
    num_train_on <class 'int'>
    variant_frequencies <class 'collections.defaultdict'>
    min_variant_freq <class 'float'>
    analysed_reference_length <class 'int'>
    cutoff_variant_number <class 'float'>
    max_ok_missingness <class 'float'>
    max_ok_missingness_pc <class 'int'>
    variant_positions_ok_missingness <class 'int'>
    mix_quality_info <class 'pandas.core.frame.DataFrame'>
    suspect_quality_seqs <class 'pandas.core.frame.DataFrame'>
    variant_matrix <class 'pandas.core.frame.DataFrame'>
    n_pca_components <class 'int'>
    pca <class 'sklearn.decomposition._pca.PCA'>
    transformed_coordinates <class 'pandas.core.frame.DataFrame'>
    eigenvectors <class 'pandas.core.frame.DataFrame'>
    explained_variance_ratio <class 'list'>
    n_contributing_positions <class 'int'>
    pc2_contributing_positions <class 'dict'>
    n_contributing_variants <class 'int'>
    contributing_basepos <class 'set'>
    contributing_pos <class 'set'>
    sample_id <class 'list'>
    pos_per_pc <class 'list'>
    build_time <class 'str'>
    coefficients_hash <class 'str'>
    transformed_coordinate_categories <class 'pandas.core.frame.DataFrame'>

    """

    def __init__(self):
        pass

    def to_json(self, obj):
        """converts the dictionary to json"""

        if not isinstance(obj, dict):
            raise TypeError("Must pass a dictionary, not a {0}".format(type(obj)))

        # output variables
        outputdict = {}
        nativetype = {}

        # remove the pca entry, as it can't be serialised & isn't needed for testing.
        try:
            del obj["pca"]  # this is an object which can't be serialised
        except KeyError:
            pass  # no pca key

        for key in obj.keys():

            if isinstance(obj[key], collections.defaultdict):
                outputdict[key] = dict(obj[key])
                nativetype[key] = "dict"
            elif isinstance(obj[key], pd.core.frame.DataFrame):
                outputdict[key] = obj[key].to_json(date_format="iso")
                nativetype[key] = "pd_df"
            elif isinstance(obj[key], set):
                outputdict[key] = list(obj[key])
                nativetype[key] = "set"
            elif isinstance(obj[key], datetime.date):
                outputdict[key] = obj[key].isoformat()
                nativetype[key] = "datetime.date"
            else:
                nativetype[key] = str(type(obj[key]))
                pass
                outputdict[key] = obj[key]

        outputdict["nativetype"] = nativetype
        return json.dumps(outputdict)

    def from_json(self, jsonobj):
        """converts the output of to_json back to native types"""
        outputdict = json.loads(jsonobj)
        nativetype = outputdict["nativetype"]
        del outputdict["nativetype"]
        for key in nativetype.keys():
            if nativetype[key] == "pd_df":
                str_df = outputdict[key]
                outputdict[key] = pd.read_json(str_df)
            elif nativetype[key] == "<class 'set'>":
                outputdict[key] = set(outputdict[key])
            elif nativetype[key] == "datetime.date":
                outputdict[key] = datetime.date.fromisoformat(outputdict[key])
        return outputdict


class Test_jsonstore(unittest.TestCase):
    """tests storage of json data for testing"""

    def runTest(self):

        # load a variation model for testiong
        inputfile = "testdata/pca/vm.pickle"
        with open(inputfile, "rb") as f:
            vm = pickle.load(f)
        m = DictStorage()
        j1 = m.to_json(vm.model)
        self.assertIsInstance(j1, str)
        m.from_json(j1)
