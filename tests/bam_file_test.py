import unittest
from lib import BamFile

class BamFileTest(unittest.TestCase):
    def test_workflow(self):
        with BamFile('samples/bam/test.bam') as bamfile:
            for alignment in bamfile.get_alignments_matching_cigarstring("100M"):
                self.assertEqual("100M", alignment.cigarstring)
                break
            self.assertEqual("BamFile(path = \"samples/bam/test.bam\")", str(bamfile))