""" contains a class to store reference compressed sequences as json
"""

import json


class SeqDictConverter:
    """stores and recovers reference compressed sequences as json.

    The reference compressed sequence is a dictionary similar to the below:
                {
                "G": set([]),
                "A": set([]),
                "C": set([]),
                "T": set([]),
                "N": set([0]),
                "U": set([0, 1, 3]),
                "M": {1: "Y", 3: "Q"},
                "invalid": 0,
            }

    contains a class to store the a dictionary of properties,
    as produced by the py_seqComparer.compress() method,
    as json.

    This is required because of
    a) the need to serialise store these dictionaries for storage
    b) the insecurity of pickle

    The key barrier here is that 
    - internally differences from the reference are stored as sets, but json only support lists
    - on back from json, in some case we expect dictionary keys as integers

    """

    def __init__(self):
        self.bases = ["A", "C", "G", "T", "N"]

    def to_json(self, obj):
        """converts the dictionary to a json string"""

        if not isinstance(obj, dict):
            raise TypeError("Must pass a dictionary, not a {0}".format(type(obj)))

        # output variables
        outputdict = {}

        for key in obj.keys():
            if isinstance(obj[key], set):
                outputdict[key] = list(obj[key])
            else:
                outputdict[key] = obj[key]

        return json.dumps(outputdict)

    def from_json(self, jsonobj):
        """converts the output of to_json ( a string ) back to native types"""
        outputdict = json.loads(jsonobj)
        for key in outputdict.keys():
            if key in self.bases:
                outputdict[key] = set(outputdict[key])

        # the keys of 'M' are converted to strings
        if "M" in outputdict.keys():
            M_dict = {}
            for base in outputdict["M"].keys():
                M_dict[int(base)] = outputdict["M"][base]
            outputdict["M"] = M_dict

        return outputdict
