import unittest
from lib import BamFile, AlignedSegment

class BamFileTest(unittest.TestCase):
    def test_workflow(self):
        with BamFile('samples/bam/test.bam') as bamfile:
            for aligned_segment in bamfile.get_alignments_matching_cigarstring("100M"):
                self.assertEqual("100M", aligned_segment.cigarstring)

                result = AlignedSegment.from_pysam_alignment_segment(aligned_segment)
                self.assertEqual(aligned_segment.cigartuples, result.get_cigar_tuples())
                self.assertEqual(aligned_segment.cigarstring, result.get_cigar_string())
                self.assertEqual(aligned_segment.query_alignment_start, result.get_query_alignment_start())
                self.assertEqual(aligned_segment.query_alignment_end, result.get_query_alignment_end())
                self.assertEqual(len(aligned_segment.query_sequence), result.get_query_sequence_length())
                self.assertEqual(aligned_segment.query_sequence, result.get_query_sequence())
                self.assertEqual(aligned_segment.query_alignment_sequence, result.get_query_alignment_sequence())
                self.assertEqual(aligned_segment.reference_start, result.get_reference_start())
                self.assertEqual(aligned_segment.reference_end, result.get_reference_end())
                self.assertEqual(aligned_segment.is_reverse, result.get_is_reverse())
                self.assertEqual(aligned_segment.query_name, result.get_read_id())
                break
            self.assertEqual("BamFile(path = \"samples/bam/test.bam\")", str(bamfile))