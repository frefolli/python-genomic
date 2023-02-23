"""
    @define test class for Alignment and AlignedSegment
"""

import unittest

from lib import AlignmentWorker
from lib import AlignedSegment
from lib import CigarOperation
from lib import Reference


def craft_alignment_worker(aligned_segment: AlignedSegment,
                           reference: Reference = None):
    """
        @does craft alignment worker for AlignedSegment and Reference
    """
    alignment_worker = AlignmentWorker(0)
    alignment_worker.set_aligned_segment(aligned_segment)
    alignment_worker.set_reference(reference)
    return alignment_worker


class AlignmentGenericTest(unittest.TestCase):
    """
        @define test class for Alignment and AlignedSegment
    """

    def test_get_cigar_tuples(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (i, i) for i in range(9)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        for (operator, length) in alignment.get_cigar_tuples():
            self.assertEqual(CigarOperation(length), operator)

    def test_get_intronic_site(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(1, alignment.get_intronic_site())


class AlignmentLeftGetExonTest(unittest.TestCase):
    """
        @define test class for getting first exon
    """

    def test_get_left_aligned_exon_length(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(5, alignment.get_left_aligned_exon_length())

    def test_get_left_aligned_exon_length_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(7, alignment.get_left_aligned_exon_length())

    def test_get_left_aligned_exon_range(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual((0, 5), alignment.get_left_aligned_exon_range())

    def test_get_left_aligned_exon_range_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2)
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual((2, 7), alignment.get_left_aligned_exon_range())

    def test_get_left_aligned_exon(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual("ACGTC", alignment.get_left_aligned_exon())

    def test_get_left_aligned_exon_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2,
            query_sequence="CCACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual("ACGTC", alignment.get_left_aligned_exon())

    def test_get_left_exon_length(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(5, alignment.get_left_exon_length())

    def test_get_left_exon_length_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(5, alignment.get_left_exon_length())


class AlignmentGetRightExonTest(unittest.TestCase):
    """
        @define test class for getting second exon
    """

    def test_get_right_aligned_exon_length(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(6, alignment.get_right_aligned_exon_length())

    def test_get_right_aligned_exon_length_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(8, alignment.get_right_aligned_exon_length())

    def test_get_right_aligned_exon_range(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual((5, 11), alignment.get_right_aligned_exon_range())

    def test_get_right_aligned_exon_range_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(11, (
            alignment.get_aligned_segment().get_query_alignment_end()))
        self.assertEqual((5, 11), alignment.get_right_aligned_exon_range())

    def test_get_right_aligned_exon(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCCCGGCC")
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual("CCGGCC", alignment.get_right_aligned_exon())

    def test_get_right_aligned_exon_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ],
            query_sequence="ACGTCCCGGCCAA")
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(0, (
            alignment.get_aligned_segment().get_query_alignment_start()))
        self.assertEqual(13, (
            alignment.get_aligned_segment().get_query_sequence_length()))
        self.assertEqual(11, (
            alignment.get_aligned_segment().get_query_alignment_end()))
        self.assertEqual("CCGGCC", alignment.get_right_aligned_exon())


class AlignmentGetIntronTest(unittest.TestCase):
    """
        @define test class for getting intron
    """

    def test_get_intron_length(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 7),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual(7, alignment.get_intron_length())

    def test_get_intron_range(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ])
        alignment = craft_alignment_worker(aligned_segment)
        self.assertEqual((5, 10), alignment.get_intron_range())

    def test_get_intron(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("TCCCC", alignment.get_intron())

    def test_get_intron_with_head_sc(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("TCCCC", alignment.get_intron())

    def test_get_intron_with_tail_sc(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            query_alignment_end=16,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("TCCCC", alignment.get_intron())

    def test_get_intron_with_both_sc(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
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
        reference = Reference.from_name_and_sequence("X", "AAAAATCCCCGGGGGG")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("TCCCC", alignment.get_intron())

    def test_get_intron_with_reverse(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            is_reverse=True,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("CCCCG", alignment.get_intron())

    def test_get_intron_with_head_sc_and_reverse(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            is_reverse=True,
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("CCCCG", alignment.get_intron())

    def test_get_intron_with_tail_sc_and_reverse(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            is_reverse=True,
            query_alignment_end=16,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("CCCCG", alignment.get_intron())

    def test_get_intron_with_both_sc_and_reverse(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 3)
            ],
            is_reverse=True,
            query_alignment_end=18,
            query_alignment_start=2,
            reference_start=0,
            reference_end=16,
            query_sequence="AAAAAGGGGGG")
        reference = Reference.from_name_and_sequence("X", "CCCCCCGGGGATTTTT")
        alignment = craft_alignment_worker(aligned_segment, reference)
        self.assertEqual("CCCCG", alignment.get_intron())


class AlignmentGetPartsOnForwardStrandTest(unittest.TestCase):
    """
        @define test class for getting parts on forward strand
    """

    def test_get_parts_from_alignment(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCCCGGCC")
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("ACGTC", first_exon)
        self.assertEqual("AAAAA", str(intron))
        self.assertEqual("CCGGCC", second_exon)

    def test_get_parts_from_alignment_with_soft_clipping(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="CCACGTCCCGGCC")
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("ACGTC", first_exon)
        self.assertEqual("AAAAA", str(intron))
        self.assertEqual("CCGGCC", second_exon)

    def test_get_parts_from_alignment_with_soft_clipping_2(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ],
            query_sequence="ACGTCCCGGCCCC")
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("ACGTC", first_exon)
        self.assertEqual("AAAAA", str(intron))
        self.assertEqual("CCGGCC", second_exon)


class AlignmentGetPartsOnReverseStrandTest(unittest.TestCase):
    """
        @define test class for getting parts on reverse strand
    """

    def test_get_parts_from_alignment_on_reverse_strand(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="ACGTCCCGGCC",
            is_reverse=True)
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("GGCCGG", first_exon)
        self.assertEqual("TTTTT", str(intron))
        self.assertEqual("GACGT", second_exon)

    def test_get_parts_from_alignment_with_soft_clipping_right(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.S, 2),
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6)
            ],
            query_sequence="CCACGTCCCGGCC",
            is_reverse=True)
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("GGCCGG", first_exon)
        self.assertEqual("TTTTT", str(intron))
        self.assertEqual("GACGT", second_exon)

    def test_get_parts_from_alignment_with_soft_clipping_left(self):
        """
            @does test
        """
        aligned_segment = AlignedSegment.from_properties(
            cigar_tuples=[
                (CigarOperation.M, 5),
                (CigarOperation.N, 5),
                (CigarOperation.M, 6),
                (CigarOperation.S, 2)
            ],
            query_sequence="ACGTCCCGGCCCC",
            is_reverse=True)
        alignment = craft_alignment_worker(aligned_segment)
        reference = Reference.from_name_and_sequence("X", "ACGTCAAAAACCGGCC")
        alignment = craft_alignment_worker(aligned_segment, reference)
        (first_exon, intron, second_exon) = alignment.split_read()
        self.assertEqual("GGCCGG", first_exon)
        self.assertEqual("TTTTT", str(intron))
        self.assertEqual("GACGT", second_exon)
