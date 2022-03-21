""" test localstore system for local json storage """
import unittest
import datetime
import os
from localstore.localstoreutils import LocalStore


class Test_LS_1(unittest.TestCase):
    """tests storage of json data for testing"""

    def runTest(self):

        # test whether a reference compressed data structure can be recovered from json
        input = {
            "A": set([1, 2, 3]),
            "C": set([6]),
            "T": set([4]),
            "G": set([13, 14]),
            "invalid": 1,
        }

        js = LocalStore("unitTest_tmp/test.tar")

        res = js._compress(input)
        self.assertIsInstance(res, bytes)
        recycled = js._decompress(res)

        self.assertEqual(input, recycled)


class Test_LS_2(unittest.TestCase):
    """tests creation of a tarfile,
    including addition of items and reading them"""

    def runTest(self):

        for compression_method in ["lzma", "gzip", "pickle"]:
            # delete any existing test
            tarfile_name = "unitTest_tmp/test.tar"
            try:
                os.unlink(tarfile_name)
            except FileNotFoundError:
                pass

            self.assertEqual(False, os.path.exists(tarfile_name))
            # create the tar file
            js = LocalStore(
                "unitTest_tmp/test.tar",
                write_batch_size=3,
                compression_method=compression_method,
            )

            self.assertTrue(True, os.path.exists(tarfile_name))

            self.assertEqual([], js.sequence_ids())

            input = {
                "A": set([]),
                "C": set([6]),
                "T": set([4]),
                "G": set([13, 14]),
                "invalid": 1,
            }

            expected_sequence_ids = []
            for i in range(11):
                input["A"] = set([100 + i])
                sequence_id = "SEQ_{0}".format(i)
                js.store(sequence_id, input)
                expected_sequence_ids.append(sequence_id)
            js.flush()

            self.assertEqual(expected_sequence_ids, js.sequence_ids())

            # check we can read back the right thing by name
            for i, sequence_id in enumerate(expected_sequence_ids):
                read_sequence_id, res = js.read(sequence_id)
                self.assertEqual(read_sequence_id, sequence_id)
                self.assertIsInstance(res, dict)
                self.assertEqual(res["A"], set([i + 100]))

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.read_many():
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.read_many(
                select_sequence_ids=["SEQ_0", "SEQ_1"]
            ):
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(["SEQ_0", "SEQ_1"], observed_sequence_ids)

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.read_all():
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)


class Test_LS_2a(unittest.TestCase):
    """tests creation of a tarfile,
    including addition of items and reading them
    uses functions which alias sequence_ids(), read_many() and read_all()"""

    def runTest(self):

        for compression_method in ["lzma", "gzip", "pickle"]:
            # delete any existing test
            tarfile_name = "unitTest_tmp/test.tar"
            try:
                os.unlink(tarfile_name)
            except FileNotFoundError:
                pass

            self.assertEqual(False, os.path.exists(tarfile_name))
            # create the tar file
            js = LocalStore(
                "unitTest_tmp/test.tar",
                write_batch_size=3,
                compression_method=compression_method,
            )

            self.assertTrue(True, os.path.exists(tarfile_name))

            self.assertEqual(set([]), js.refcompressedsequence_guids())

            input = {
                "A": set([]),
                "C": set([6]),
                "T": set([4]),
                "G": set([13, 14]),
                "invalid": 1,
            }

            expected_sequence_ids = []
            for i in range(11):
                input["A"] = set([100 + i])
                sequence_id = "SEQ_{0}".format(i)
                js.store(sequence_id, input)
                expected_sequence_ids.append(sequence_id)
            js.flush()

            self.assertEqual(set(expected_sequence_ids), js.refcompressedsequence_guids())

            # check we can read back the right thing by name
            for i, sequence_id in enumerate(expected_sequence_ids):
                read_sequence_id, res = js.read(sequence_id)
                self.assertEqual(read_sequence_id, sequence_id)
                self.assertIsInstance(res, dict)
                self.assertEqual(res["A"], set([i + 100]))

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.refcompressedsequence_read_many():
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.refcompressedsequence_read_many(
                select_sequence_ids=["SEQ_0", "SEQ_1"]
            ):
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(["SEQ_0", "SEQ_1"], observed_sequence_ids)

            i = 0
            observed_sequence_ids = []
            for (sequence_id, res) in js.refcompressedsequence_read_all():
                if sequence_id is None:
                    break
                else:
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    self.assertEqual(res["A"], set([i + 100]))
                    observed_sequence_ids.append(sequence_id)
                    i = i + 1

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)


class Test_LS_3(unittest.TestCase):
    """tests read_benchmark, which measures read speeds"""

    def runTest(self):

        for compression_method in ["lzma", "gzip", "pickle"]:
            # delete any existing test
            tarfile_name = "unitTest_tmp/test.tar"
            try:
                os.unlink(tarfile_name)
            except FileNotFoundError:
                pass

            self.assertEqual(False, os.path.exists(tarfile_name))
            # create the tar file
            js = LocalStore(
                "unitTest_tmp/test.tar",
                write_batch_size=3,
                compression_method=compression_method,
            )

            self.assertTrue(True, os.path.exists(tarfile_name))

            self.assertEqual([], js.sequence_ids())

            input = {
                "A": set([]),
                "C": set([6]),
                "T": set([4]),
                "G": set([13, 14]),
                "invalid": 1,
            }

            expected_sequence_ids = []
            for i in range(50):
                sequence_id = "SEQ_{0}".format(i)
                js.store(sequence_id, input)
                expected_sequence_ids.append(sequence_id)
            js.flush()

            self.assertEqual(expected_sequence_ids, js.sequence_ids())

            # check we can read back the right thing by name
            res = js.read_benchmark()
            self.assertIsInstance(res, dict)
            res = js.read_benchmark(restrict_to_first_n=10)
            self.assertIsInstance(res, dict)


class Test_LS_4(unittest.TestCase):
    """tests functions operating on an empty tar file"""

    def runTest(self):

        for compression_method in ["lzma", "gzip", "pickle"]:
            # delete any existing test
            tarfile_name = "unitTest_tmp/test.tar"
            try:
                os.unlink(tarfile_name)
            except FileNotFoundError:
                pass

            self.assertEqual(False, os.path.exists(tarfile_name))
            # create the tar file
            js = LocalStore(
                "unitTest_tmp/test.tar",
                write_batch_size=3,
                compression_method=compression_method,
            )

            self.assertTrue(True, os.path.exists(tarfile_name))

            self.assertEqual([], js.sequence_ids())

            expected_sequence_ids = []
            self.assertEqual(expected_sequence_ids, js.sequence_ids())

            observed_sequence_ids = []
            for (sequence_id, res) in js.read_many():
                if sequence_id is None:
                    break
                else:
                    observed_sequence_ids.append(sequence_id)

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)

            observed_sequence_ids = []
            for (sequence_id, res) in js.read_many(
                select_sequence_ids=["SEQ_0", "SEQ_1"]
            ):
                if sequence_id is None:
                    break
                else:
                    observed_sequence_ids.append(sequence_id)

            self.assertEqual([], observed_sequence_ids)

            observed_sequence_ids = []
            for (sequence_id, res) in js.read_all():
                if sequence_id is None:
                    break
                else:
                    observed_sequence_ids.append(sequence_id)

            self.assertEqual(expected_sequence_ids, observed_sequence_ids)


@unittest.skip(reason="benchmark")
class Test_LS_benchmark(unittest.TestCase):
    """benchmark creation of a tarfile,
    including addition of items and reading them
    based on a synthetic reference sequence with 10002 variants from reference"""

    def runTest(self):

        # delete any existing test

        for compression_method in ["lzma", "gzip", "pickle"]:
            # delete any existing test

            tarfile_name = "unitTest_tmp/test.tar"
            try:
                os.unlink(tarfile_name)
            except FileNotFoundError:
                pass

            self.assertEqual(False, os.path.exists(tarfile_name))
            # create the tar file
            js = LocalStore(tarfile_name, compression_method=compression_method)

            self.assertTrue(True, os.path.exists(tarfile_name))

            self.assertEqual([], js.sequence_ids())

            input = {
                "A": set([]),
                "C": set([6]),
                "T": set([4]),
                "G": set(range(10000)),
                "invalid": 1,
            }

            expected_sequence_ids = []
            t0 = datetime.datetime.now()
            for i in range(5000):
                # if i % 1000 == 0:
                #    print("storing", i)
                input["A"] = set([100 + i])
                sequence_id = "SEQ_{0}".format(i)
                js.store(sequence_id, input)
                expected_sequence_ids.append(sequence_id)
            js.flush()

            fsize = os.path.getsize(tarfile_name)

            t1 = datetime.datetime.now()
            i = 0
            for (sequence_id, res) in js.read_all():
                if sequence_id is None:
                    break
                else:
                    # if i % 1000== 0:
                    #    print("reading", i)
                    self.assertEqual(sequence_id, "SEQ_{0}".format(i))
                    self.assertIsInstance(res, dict)
                    # self.assertEqual(res['A'], set([i+100]))
                    i = i + 1

            t2 = datetime.datetime.now()
            d2 = (t2 - t1).total_seconds()
            d1 = (t1 - t0).total_seconds()
            print(compression_method, d1, d2, fsize)
