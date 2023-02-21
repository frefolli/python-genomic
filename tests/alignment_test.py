import unittest
from lib import Alignment
from lib import AlignedSegment
from lib import CigarOperation
from lib import Reference
from Bio import SeqRecord

class AlignmentTest(unittest.TestCase):
    def test_get_cigar_tuples(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
            (i, i) for i in range(9)
        ])
        alignment = Alignment(0, aligned_segment)
        for (op, length) in alignment.get_cigar_tuples():
            self.assertEqual(CigarOperation(length), op)
    
    def test_get_intronic_site(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(1, alignment.get_intronic_site())
    
    def test_get_first_aligned_exon_length(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(5, alignment.get_first_aligned_exon_length())
    
    def test_get_first_aligned_exon_length_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(7, alignment.get_first_aligned_exon_length())
    
    def test_get_first_aligned_exon_range(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual((0, 5), alignment.get_first_aligned_exon_range())

    def test_get_first_aligned_exon_range_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2)
        alignment = Alignment(0, aligned_segment)
        self.assertEqual((2, 7), alignment.get_first_aligned_exon_range())
    
    def test_get_second_aligned_exon_length(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(6, alignment.get_second_aligned_exon_length())
    
    def test_get_second_aligned_exon_length_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(8, alignment.get_second_aligned_exon_length())
    
    def test_get_second_aligned_exon_range(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual((10, 16), alignment.get_second_aligned_exon_range())
    
    def test_get_second_aligned_exon_range_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 1)
            ],
            query_alignment_end=16)
        alignment = Alignment(0, aligned_segment)
        self.assertEqual((10, 16), alignment.get_second_aligned_exon_range())
    
    def test_get_first_aligned_exon(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCAAAAACCGGCC")
        alignment = Alignment(0, aligned_segment)
        self.assertEqual("ACGTC", alignment.get_first_aligned_exon())
    
    def test_get_first_aligned_exon_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2,
            query_sequence="CCACGTCAAAAACCGGCC")
        alignment = Alignment(0, aligned_segment)
        self.assertEqual("ACGTC", alignment.get_first_aligned_exon())
    
    def test_get_second_aligned_exon(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCAAAAACCGGCC")
        alignment = Alignment(0, aligned_segment)
        self.assertEqual("CCGGCC", alignment.get_second_aligned_exon())
    
    def test_get_second_aligned_exon_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ],
            query_alignment_end=16,
            query_sequence="ACGTCAAAAACCGGCCAA")
        alignment = Alignment(0, aligned_segment)
        self.assertEqual("CCGGCC", alignment.get_second_aligned_exon())
    
    def test_get_first_exon_length(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(5, alignment.get_first_exon_length())
    
    def test_get_first_exon_length_with_soft_clipping(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
        ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(5, alignment.get_first_exon_length())
    
    def test_get_intron_length(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 7),
                (CigarOperation.M, 6)
            ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual(7, alignment.get_intron_length())
    
    def test_get_intron_range(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = Alignment(0, aligned_segment)
        self.assertEqual((5, 10), alignment.get_intron_range())
    
    def test_get_intron(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        self.assertEqual("TCCCC", alignment.get_intron(reference))
    
    def test_get_intron_with_head_sc(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        self.assertEqual("TCCCC", alignment.get_intron(reference))
    
    def test_get_intron_with_tail_sc(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            query_alignment_end=16,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        self.assertEqual("TCCCC", alignment.get_intron(reference))
    
    def test_get_intron_with_both_sc(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            query_alignment_end=18,
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        self.assertEqual("TCCCC", alignment.get_intron(reference))

    def test_get_intron_with_reverse(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            is_reverse = True,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        self.assertEqual("CCCCG", alignment.get_intron(reference))
    
    def test_get_intron_with_head_sc_and_reverse(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            is_reverse = True,
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        self.assertEqual("CCCCG", alignment.get_intron(reference))
    
    def test_get_intron_with_tail_sc_and_reverse(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            is_reverse = True,
            query_alignment_end=16,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        self.assertEqual("CCCCG", alignment.get_intron(reference))
    
    def test_get_intron_with_both_sc_and_reverse(self):
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples = [
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            is_reverse = True,
            query_alignment_end=18,
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        alignment = Alignment(0, aligned_segment)
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        self.assertEqual("CCCCG", alignment.get_intron(reference))