"""
    @does test FastaFile
"""

import unittest

from genomic import FastaFile


class FastaFileTest(unittest.TestCase):
    """
        @does test FastaFile
    """

    def test_workflow(self):
        """
            @does test
        """
        with FastaFile('samples/fasta/test.fasta') as fastafile:
            for reference_name in fastafile.get_reference_names():
                reference = fastafile.get_reference(reference_name)
                self.assertEqual("CGTCA", reference.get_slice(5, 10))

                self.assertEqual("ATTCGTCAGAAATGAGCTAAACAAATTTA" +
                                 "AATCATTAAATGCGAGCGGCGAATCCGGA",
                                 reference.get_slice(2))

                self.assertEqual("GAATTCGTCAGAAATGAGCTAAACAAATT" +
                                 "TAAATCATTAAATGCGAGCGGCGAATCCGGA",
                                 reference.get_slice())

                self.assertEqual("GAATTCGTCA",
                                 reference.get_slice(end_index=10))

                self.assertEqual("AATT", reference.get_slice(1, 5))
                self.assertEqual("AATT", reference.get_slice(1, 5, True))

                self.assertEqual("TCCGGATTCGCCGCTCGCATTTAATGATTT" +
                                 "AAATTTGTTTAGCTCATTTCTGACGAATTC",
                                 reference.get_slice(reverse=True))

                self.assertEqual("Reference(name = \"ENm006\", " +
                                 "sequence = \"" +
                                 "GAATTCGTCAGAAATGAGCTAAACAAATTT" +
                                 "AAATCATTAAATGCGAGCGGCGAATCCGGA\")",
                                 str(reference))
            self.assertEqual("FastaFile(path = \"samples/fasta/test.fasta\")",
                             str(fastafile))
