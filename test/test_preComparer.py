""" tests preComparer.py

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

import unittest
import os
import io
import glob
import gzip
from Bio import SeqIO

from findn.preComparer import preComparer


class test_preComparer_1(unittest.TestCase):
    """tests __init__ method"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            over_selection_cutoff_ignore_factor=5,
            uncertain_base="M",
        )
        self.assertIsInstance(sc, preComparer)


class test_preComparer_1b(unittest.TestCase):
    """tests check_operating_parameters method"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            over_selection_cutoff_ignore_factor=5,
            uncertain_base="M",
        )

        sc.set_operating_parameters(
            selection_cutoff=20,
            over_selection_cutoff_ignore_factor=5,
            uncertain_base="M",
        )

        self.assertTrue(
            sc.check_operating_parameters(
                selection_cutoff=20,
                over_selection_cutoff_ignore_factor=5,
                uncertain_base="M",
            )
        )

        sc.set_operating_parameters(
            selection_cutoff=10,
            over_selection_cutoff_ignore_factor=5,
            uncertain_base="M",
        )

        self.assertFalse(
            sc.check_operating_parameters(
                selection_cutoff=20,
                over_selection_cutoff_ignore_factor=5,
                uncertain_base="M",
            )
        )


class test_preComparer_2(unittest.TestCase):
    """tests storage"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )
        obj = {"A": set([1, 2, 3, 4])}
        sc.persist(obj, "guid1")
        self.assertEqual(sc.composition["guid1"]["A"], 4)
        self.assertEqual(sc.composition["guid1"]["C"], 0)
        self.assertEqual(sc.composition["guid1"]["T"], 0)
        self.assertEqual(sc.composition["guid1"]["G"], 0)
        self.assertEqual(sc.composition["guid1"]["N"], 0)
        self.assertEqual(sc.composition["guid1"]["M"], 0)
        self.assertEqual(sc.composition["guid1"]["invalid"], 0)


class test_preComparer_2a(unittest.TestCase):
    """tests storage"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
            catWalk_parameters={
                "cw_binary_filepath": None,
                "reference_name": "h37rv",
                "reference_filepath": "reference/TB-ref.fasta",
                "mask_filepath": "reference/TB-exclude-adaptive.txt",
                "bind_host": "127.0.0.1",
                "bind_port": 5999,
            },
        )
        obj = {"A": set([1, 2, 3, 4])}
        sc.persist(obj, "guid1")
        self.assertEqual(sc.composition["guid1"]["A"], 4)
        self.assertEqual(sc.composition["guid1"]["C"], 0)
        self.assertEqual(sc.composition["guid1"]["T"], 0)
        self.assertEqual(sc.composition["guid1"]["G"], 0)
        self.assertEqual(sc.composition["guid1"]["N"], 0)
        self.assertEqual(sc.composition["guid1"]["M"], 0)
        self.assertEqual(sc.composition["guid1"]["invalid"], 0)


class test_preComparer_3(unittest.TestCase):
    """tests storage of invalid samples"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )
        obj = {"A": set([1, 2, 3, 4]), "invalid": 1}
        sc.persist(obj, "guid1")
        self.assertEqual(sc.composition["guid1"]["A"], 4)
        self.assertEqual(sc.composition["guid1"]["C"], 0)
        self.assertEqual(sc.composition["guid1"]["T"], 0)
        self.assertEqual(sc.composition["guid1"]["G"], 0)
        self.assertEqual(sc.composition["guid1"]["N"], 0)
        self.assertEqual(sc.composition["guid1"]["M"], 0)
        self.assertEqual(sc.composition["guid1"]["invalid"], 1)

        self.assertEqual(
            set(sc.seqProfile["guid1"].keys()),
            set(["invalid", "A", "T", "C", "U", "G"]),
        )
        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")
        self.assertEqual(sc.composition["guid2"]["A"], 4)
        self.assertEqual(sc.composition["guid2"]["C"], 0)
        self.assertEqual(sc.composition["guid2"]["T"], 0)
        self.assertEqual(sc.composition["guid2"]["G"], 0)
        self.assertEqual(sc.composition["guid2"]["N"], 0)
        self.assertEqual(sc.composition["guid2"]["M"], 0)
        self.assertEqual(sc.composition["guid2"]["invalid"], 0)

        obj = {"A": set([1, 2, 3, 4]), "N": set([10, 11]), "invalid": 0}
        sc.persist(obj, "guid3")
        self.assertEqual(sc.composition["guid3"]["A"], 4)
        self.assertEqual(sc.composition["guid3"]["C"], 0)
        self.assertEqual(sc.composition["guid3"]["T"], 0)
        self.assertEqual(sc.composition["guid3"]["G"], 0)
        self.assertEqual(sc.composition["guid3"]["N"], 2)
        self.assertEqual(sc.composition["guid3"]["M"], 0)
        self.assertEqual(sc.composition["guid3"]["invalid"], 0)

        obj = {
            "A": set([1, 2, 3, 4]),
            "N": set([10, 11, 12]),
            "M": {13: "Y"},
            "invalid": 0,
        }
        sc.persist(obj, "guid4")
        self.assertEqual(sc.composition["guid4"]["A"], 4)
        self.assertEqual(sc.composition["guid4"]["C"], 0)
        self.assertEqual(sc.composition["guid4"]["T"], 0)
        self.assertEqual(sc.composition["guid4"]["G"], 0)
        self.assertEqual(sc.composition["guid4"]["N"], 3)
        self.assertEqual(sc.composition["guid4"]["M"], 1)
        self.assertEqual(sc.composition["guid4"]["invalid"], 0)

        self.assertEqual(
            sc.guidscachedinram(), set(["guid1", "guid2", "guid3", "guid4"])
        )

        self.assertEqual(
            set(sc.seqProfile["guid2"].keys()),
            set(["A", "C", "T", "G", "U", "invalid"]),
        )
        self.assertEqual(
            set(sc.seqProfile["guid1"].keys()),
            set(["A", "C", "T", "G", "U", "invalid"]),
        )


class test_preComparer_3a(unittest.TestCase):
    """tests storage of invalid samples"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
            catWalk_parameters={
                "cw_binary_filepath": None,
                "reference_name": "h37rv",
                "reference_filepath": "reference/TB-ref.fasta",
                "mask_filepath": "reference/TB-exclude-adaptive.txt",
                "bind_host": "127.0.0.1",
                "bind_port": 5999,
            },
        )
        obj = {"A": set([1, 2, 3, 4]), "invalid": 1}
        sc.persist(obj, "guid1")
        self.assertEqual(sc.composition["guid1"]["A"], 4)
        self.assertEqual(sc.composition["guid1"]["C"], 0)
        self.assertEqual(sc.composition["guid1"]["T"], 0)
        self.assertEqual(sc.composition["guid1"]["G"], 0)
        self.assertEqual(sc.composition["guid1"]["N"], 0)
        self.assertEqual(sc.composition["guid1"]["M"], 0)
        self.assertEqual(sc.composition["guid1"]["invalid"], 1)

        self.assertEqual(set(sc.seqProfile["guid1"].keys()), set(["invalid"]))

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")
        self.assertEqual(sc.composition["guid2"]["A"], 4)
        self.assertEqual(sc.composition["guid2"]["C"], 0)
        self.assertEqual(sc.composition["guid2"]["T"], 0)
        self.assertEqual(sc.composition["guid2"]["G"], 0)
        self.assertEqual(sc.composition["guid2"]["N"], 0)
        self.assertEqual(sc.composition["guid2"]["M"], 0)
        self.assertEqual(sc.composition["guid2"]["invalid"], 0)

        obj = {"A": set([1, 2, 3, 4]), "N": set([10, 11]), "invalid": 0}
        sc.persist(obj, "guid3")
        self.assertEqual(sc.composition["guid3"]["A"], 4)
        self.assertEqual(sc.composition["guid3"]["C"], 0)
        self.assertEqual(sc.composition["guid3"]["T"], 0)
        self.assertEqual(sc.composition["guid3"]["G"], 0)
        self.assertEqual(sc.composition["guid3"]["N"], 2)
        self.assertEqual(sc.composition["guid3"]["M"], 0)
        self.assertEqual(sc.composition["guid3"]["invalid"], 0)

        obj = {
            "A": set([1, 2, 3, 4]),
            "N": set([10, 11, 12]),
            "M": {13: "Y"},
            "invalid": 0,
        }
        sc.persist(obj, "guid4")
        self.assertEqual(sc.composition["guid4"]["A"], 4)
        self.assertEqual(sc.composition["guid4"]["C"], 0)
        self.assertEqual(sc.composition["guid4"]["T"], 0)
        self.assertEqual(sc.composition["guid4"]["G"], 0)
        self.assertEqual(sc.composition["guid4"]["N"], 3)
        self.assertEqual(sc.composition["guid4"]["M"], 1)
        self.assertEqual(sc.composition["guid4"]["invalid"], 0)

        self.assertEqual(
            sc.guidscachedinram(), set(["guid1", "guid2", "guid3", "guid4"])
        )
        self.assertEqual(set(sc.seqProfile["guid2"].keys()), set(["invalid"]))
        self.assertEqual(set(sc.seqProfile["guid1"].keys()), set(["invalid"]))


class test_preComparer_4(unittest.TestCase):
    """tests reporting of server status"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )
        res1 = sc.summarise_stored_items()
        self.assertEqual({"server|pcstat|nSeqs": 0}, res1)

        obj = {"A": set([1, 2, 3, 4])}
        sc.persist(obj, "guid1")
        res2 = sc.summarise_stored_items()
        self.assertEqual({"server|pcstat|nSeqs": 1}, res2)


class test_preComparer_4a(unittest.TestCase):
    """tests reporting of server status"""

    def runTest(self):
        # initialise comparer
        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
            catWalk_parameters={
                "cw_binary_filepath": None,
                "reference_name": "h37rv",
                "reference_filepath": "reference/TB-ref.fasta",
                "mask_filepath": "reference/TB-exclude-adaptive.txt",
                "bind_host": "127.0.0.1",
                "bind_port": 5999,
                "unittesting": True,
            },
        )
        res1 = sc.summarise_stored_items()
        self.assertEqual(res1["server|pcstat|nSeqs"], 0)
        self.assertEqual(
            res1["server|catwalk|mask_name"], "reference/TB-exclude-adaptive.txt"
        )


class test_preComparer_5(unittest.TestCase):
    """tests comparison"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")

        obj = {"A": set([1, 2, 3, 4]), "N": set([10, 11]), "invalid": 0}
        sc.persist(obj, "guid3")

        res = sc.compare("guid2", "guid3")
        self.assertTrue(res is not None)


class test_preComparer_6(unittest.TestCase):
    """tests comparison"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")

        obj = {"A": set([3, 4, 5, 6]), "N": set([10, 11]), "invalid": 0}
        sc.persist(obj, "guid3")

        res = sc.compare("guid2", "guid3")
        self.assertTrue(res is not None)


class test_preComparer_7(unittest.TestCase):
    """tests comparison"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")

        obj = {"A": set([3, 4, 5, 6]), "N": set([10, 11]), "invalid": 1}
        sc.persist(obj, "guid3")

        # check invalid set
        self.assertEqual(sc.seqProfile["guid3"]["invalid"], 1)
        self.assertEqual(sc.seqProfile["guid2"]["invalid"], 0)

        res = sc.compare("guid2", "guid3")
        self.assertIsNone(res["dist"])


class test_preComparer_8(unittest.TestCase):
    """tests comparison"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")

        # check all keys are present
        self.assertEqual(
            set(sc.seqProfile["guid2"].keys()),
            set(["U", "A", "C", "G", "T", "invalid"]),
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 1}
        sc.persist(obj, "guid3")

        # check all keys are present
        self.assertEqual(
            set(sc.seqProfile["guid3"].keys()),
            set(["U", "A", "C", "G", "T", "invalid"]),
        )

        # check invalid set
        self.assertEqual(sc.seqProfile["guid3"]["invalid"], 1)


class test_preComparer_9(unittest.TestCase):
    """tests mcompare"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        sc.persist(obj, "guid2")

        # check all keys are present
        self.assertEqual(
            set(sc.seqProfile["guid2"].keys()),
            set(["U", "A", "C", "G", "T", "invalid"]),
        )

        obj = {"A": set([1, 2, 3, 4, 5]), "invalid": 0}
        sc.persist(obj, "guid3")

        res = sc.mcompare("guid2")
        self.assertEqual(len(res), 1)


class test_preComparer_9a(unittest.TestCase):
    """tests mcompare"""

    def runTest(self):
        # print("Starting precomparer")
        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
            catWalk_parameters={
                "cw_binary_filepath": None,
                "reference_name": "h37rv",
                "reference_filepath": "reference/TB-ref.fasta",
                "mask_filepath": "reference/TB-exclude-adaptive.txt",
                "bind_host": "127.0.0.1",
                "bind_port": 5999,
            },
        )

        obj = {"A": set([1, 2, 3, 4]), "invalid": 0}
        # print("persisting")
        sc.persist(obj, "guid2")

        # check only invalid key is present
        self.assertEqual(set(sc.seqProfile["guid2"].keys()), set(["invalid"]))

        obj = {"A": set([1, 2, 3, 4, 5]), "invalid": 0}
        sc.persist(obj, "guid3")

        res = sc.mcompare("guid2")
        self.assertEqual(len(res), 1)


class test_preComparer_10(unittest.TestCase):
    """tests comparison"""

    def runTest(self):

        sc = preComparer(
            selection_cutoff=20,
            uncertain_base="M",
            over_selection_cutoff_ignore_factor=5,
        )

        obj = {"invalid": 1}
        sc.persist(obj, "guid2")

        obj = {"A": set([1, 2, 3, 4]), "N": set([10, 11]), "invalid": 0}
        sc.persist(obj, "guid3")

        obj = {"A": set([1, 2, 3, 4]), "N": set([10, 11]), "invalid": 0}
        res = sc.persist(
            obj, "guid2"
        )  # should return results for the previously stored guid2

        self.assertEqual(res, 1)


class test_preComparer_11(unittest.TestCase):
    """compares catwalk vs python comparisons with real data"""

    def compress(self, sequence, reference):
        """reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases.
        does not use a mask, in this toy example.

        """
        if not len(sequence) == len(reference):
            raise TypeError(
                "sequence must of the same length as reference; seq is {0} and ref is {1}".format(
                    len(sequence), len(reference)
                )
            )
        if len(reference) == 0:
            raise TypeError("reference cannot be of zero length")

        # we consider - characters to be the same as N
        sequence = sequence.replace("-", "N")

        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        # a dictionary, M, records the mixed base calls.
        diffDict = {
            "A": set([]),
            "C": set([]),
            "T": set([]),
            "G": set([]),
            "N": set([]),
            "M": {},
        }
        for i in range(len(sequence)):  # no mask: consider all sequences
            if not sequence[i] == reference[i]:  # if it's not reference
                if sequence[i] in ["A", "C", "T", "G", "N"]:
                    diffDict[sequence[i]].add(i)  # if it's a definitively called base
                else:
                    # we regard it as a code representing a mixed base.  we store the results in a dictionary
                    diffDict["M"][i] = sequence[i]

        diffDict["invalid"] = 0
        return diffDict

    def read_fasta_file(self, fastafile):
        """reads the content of a fasta file into memory.
        returns a dictionary {seqid:(first part of defline), seq:(nucleic acid), content:(entire file content)}.
        Supports both .gz and uncompressed files transparently.
        Does not support multi-fasta files.  Will raise an error if such are detected.
        """
        # first determine whether it is a .gz file or not; read into RAM.
        if fastafile.endswith(".gz"):
            # we decompress it on the fly.
            with gzip.open(fastafile, "r") as f:
                content = f.read().decode("utf-8")
        else:
            with open(fastafile, "rt") as f:
                content = f.read()

        # use BioPython3 SeqIO library to read the file.
        nFiles = 0
        with io.StringIO(content) as f:
            for record in SeqIO.parse(f, "fasta"):
                nFiles += 1
                if nFiles > 1:  # that's a multifasta, and we don't support that
                    raise ValueError(
                        "Multifasta file is present in {0}.  Multifasta files are not supported".format(
                            fastafile
                        )
                    )
                else:
                    res = {
                        "seq": str(record.seq),
                        "seqid": str(record.id),
                        "content": content,
                    }
                    return res
        raise IOError(
            "no content parsed from result of length {0}".format(len(content))
        )

    def runTest(self):
        """do comparison between cw and standard snv computation methods"""

        uncertain_base = "N_or_M"  # consider Ns or Ms as unknown. (cw calls these 'N')
        selection_cutoff = 5000000
        print(
            "#1 is running conventional python based comparisons with {0} snv cutoff".format(
                selection_cutoff
            )
        )
        print("#1 is using catwalk with {0} snv cutoff".format(selection_cutoff))
        sc1 = preComparer(
            selection_cutoff=selection_cutoff,
            over_selection_cutoff_ignore_factor=1,
            uncertain_base=uncertain_base,
        )

        sc2 = preComparer(
            selection_cutoff=selection_cutoff,
            over_selection_cutoff_ignore_factor=1,
            uncertain_base=uncertain_base,
            catWalk_parameters={
                "cw_binary_filepath": None,
                "reference_name": "h37rv",
                "reference_filepath": "reference/TB-ref.fasta",
                "mask_filepath": "reference/nil.txt",
                "bind_host": "127.0.0.1",
                "bind_port": 5999,
                "unittesting": True,
            },
        )

        # define directory where the fastas are
        fastadir = os.path.join("..", "demos", "AC587", "fasta")

        reference = self.read_fasta_file("reference/TB-ref.fasta")["seq"]

        # we load randomly selected guids
        guids = list()
        for i, fastafile in enumerate(
            sorted(glob.glob(os.path.join(fastadir, "test", "*.mfasta.gz")))
        ):
            guid = os.path.basename(fastafile).replace(".mfasta.gz", "")
            guids.append(guid)
            seq = self.read_fasta_file(fastafile)["seq"]
            rc = self.compress(seq, reference)

            # add to both
            sc1.persist(rc, guid)
            sc2.persist(rc, guid)

            if i > 5:
                break

        # get neighbours of all
        snpcmp_1 = {}
        snpcmp_2 = {}
        distrib_1 = list(range(20 + 1))
        distrib_2 = list(range(20 + 1))
        for i, guid in enumerate(guids):
            for res in sc1.mcompare(guid):
                if res["dist"] <= selection_cutoff:
                    snpcmp_1["{0} vs {1}".format(guid, res["guid2"])] = res["dist"]
                if res["dist"] <= 20:
                    distrib_1[res["dist"]] += 1

        for i, guid in enumerate(guids):
            for res in sc2.mcompare(guid):
                if res["dist"] <= selection_cutoff:
                    snpcmp_2["{0} vs {1}".format(guid, res["guid2"])] = res["dist"]
                if res["dist"] <= 20:
                    distrib_2[res["dist"]] += 1

        print(1, distrib_1)
        print(2, distrib_2)
        if distrib_1 == distrib_2:
            print("Distributions are the same")
        else:
            print("Fail: distributions differ")

        if not set(snpcmp_1.keys()) == set(snpcmp_2.keys()):
            print("FAIL: pairs identified differ")
        print(
            "Examining {0} pairs, comparing both methods; will report any discrepancies".format(
                len(snpcmp_1)
            )
        )
        failures = 0
        for key in sorted(snpcmp_1.keys()):  # compare results for both methods
            try:
                if not snpcmp_1[key] == snpcmp_2[key]:
                    print(
                        "FAIL: Distances differ for ", key, snpcmp_1[key], snpcmp_2[key]
                    )
                    failures += 1
            except KeyError:
                print(
                    "FAIL: Pair ",
                    key,
                    "is not present (likely >20) in snpcmp_2.  Python distance is ",
                    snpcmp_1[key],
                )
                failures += 1
        print("Finished, failures = {0}".format(failures))
        self.assertEqual(failures, 0)
