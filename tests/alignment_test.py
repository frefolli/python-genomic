import unittest
from lib import Alignment
from lib import BamFile
from lib import FastaFile

class AlignmentTest(unittest.TestCase):
    def test_workflow(self):
        with FastaFile('../BDGP6.X.fasta') as fastafile:
            with BamFile('samples/bam/test2.bam') as bamfile:
                matching = [_ for _ in enumerate(bamfile.get_alignments())]
                align = None
                match_rule = "([0-9]+[M])*[0-9]+N([0-9]+[M])*"
                matching = [_ for _ in enumerate(bamfile.get_alignments_matching_cigarstring(match_rule))]
                reference = fastafile.get_reference('X')
                for (id, alignment) in matching:
                    align = Alignment(id, alignment)
                    ris = align.process_alignment(reference)
                    self.assertEqual(str(ris[3][:len(ris[4])]), ris[4])
                    self.assertEqual(ris[3][len(ris[4]):len(ris[4])+ris[7]], ris[5])
                    self.assertEqual(ris[3][-len(ris[6]):], ris[6])
