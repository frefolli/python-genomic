import unittest
from lib import FastaFile

class FastaFileTest(unittest.TestCase):
    def test_workflow(self):
        with FastaFile('samples/fasta/test.fasta') as fastafile:
            for reference_name in fastafile.get_reference_names():
                reference = fastafile.get_reference(reference_name)
                self.assertEqual("CGTCA", reference.get_slice(5, 10))
                self.assertEqual("ATTCGTCAGAAATGAGCTAAACAAATTTAAATCATTAAATGCGAGCGGCGAATCCGGA", reference.get_slice(2))
                self.assertEqual("GAATTCGTCAGAAATGAGCTAAACAAATTTAAATCATTAAATGCGAGCGGCGAATCCGGA", reference.get_slice())
                self.assertEqual("GAATTCGTCA", reference.get_slice(end_index=10))
                self.assertEqual("AATT", reference.get_slice(1, 5))
                self.assertEqual("AATT", reference.get_slice(1, 5, True))
                self.assertEqual("TCCGGATTCGCCGCTCGCATTTAATGATTTAAATTTGTTTAGCTCATTTCTGACGAATTC", reference.get_slice(reverse=True))
                self.assertEqual("Reference(name = \"ENm006\", sequence = \"GAATTCGTCAGAAATGAGCTAAACAAATTTAAATCATTAAATGCGAGCGGCGAATCCGGA\")", str(reference))
            self.assertEqual("FastaFile(path = \"samples/fasta/test.fasta\")", str(fastafile))