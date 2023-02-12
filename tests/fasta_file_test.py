import unittest
from lib import FastaFile

class FastaFileTest(unittest.TestCase):
    def test_workflow(self):
        with FastaFile('samples/fasta/test.fasta') as fastafile:
            for reference_name in fastafile.get_reference_names():
                reference = fastafile.get_reference(reference_name)
                reference.get_slice(0, 10)
                reference.get_slice(0)
                reference.get_slice()
                self.assertEqual("Reference(name = \"ENm006\", sequence = \"ACGT\")", str(reference))
            self.assertEqual("FastaFile(path = \"samples/fasta/test.fasta\")", str(fastafile))