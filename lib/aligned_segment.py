"""
    @represents a useful abstraction for pysam.AlignedSegment
"""

import pysam
from lib import cigar_tuples_to_cigar_string


class AlignedSegment:
    """
        @represents a useful abstraction for pysam.AlignedSegment
    """

    def __init__(self) -> None:
        """
            @does define an empty AlignedSegment
        """
        self.cigar_tuples = None
        self.cigar_string = None
        self.query_alignment_start = None
        self.query_alignment_end = None
        self.query_sequence_length = None
        self.query_sequence = None
        self.query_alignment_sequence = None
        self.reference_start = None
        self.reference_end = None
        self.is_reverse = None
        self.read_id = None

    @staticmethod
    def from_properties(
            query_sequence=None,
            cigar_tuples=None,
            cigar_string=None,
            query_sequence_length=None,
            query_alignment_start=None,
            query_alignment_end=None,
            reference_start=None,
            reference_end=None,
            is_reverse=None,
            read_id=None):
        """
            @does build an AlignedSegment using properties
        """

        result = AlignedSegment()

        if query_sequence is None:
            query_sequence = ""
        if cigar_tuples is None:
            cigar_tuples = [(0, len(query_sequence))]
        if cigar_string is None:
            cigar_string = cigar_tuples_to_cigar_string(cigar_tuples)
        if query_sequence_length is None:
            query_sequence_length = sum([op[1] for op in cigar_tuples])
        if query_alignment_start is None:
            query_alignment_start = 0
        if query_alignment_end is None:
            query_alignment_end = query_sequence_length
        if reference_start is None:
            reference_start = 0
        if reference_end is None:
            reference_end = 0
        if is_reverse is None:
            is_reverse = False
        if read_id is None:
            read_id = 0

        result.cigar_tuples = cigar_tuples
        result.cigar_string = cigar_string
        result.query_alignment_start = query_alignment_start
        result.query_alignment_end = query_alignment_end
        result.query_sequence_length = query_sequence_length
        result.query_sequence = query_sequence
        result.query_alignment_sequence = (
            query_sequence[query_alignment_start: query_alignment_end])
        result.reference_start = reference_start
        result.reference_end = reference_end
        result.is_reverse = is_reverse
        result.read_id = read_id

        return result

    @staticmethod
    def from_pysam_alignment_segment(aligned_segment: pysam.AlignedSegment):
        """
            @does build an AlignedSegment using pysam.AlignedSegment
        """

        result = AlignedSegment()

        result.cigar_tuples = aligned_segment.cigartuples
        result.cigar_string = aligned_segment.cigarstring
        result.query_alignment_start = aligned_segment.query_alignment_start
        result.query_alignment_end = aligned_segment.query_alignment_end
        result.query_sequence_length = len(aligned_segment.query_sequence)
        result.query_sequence = aligned_segment.query_sequence
        result.query_alignment_sequence = (
            aligned_segment.query_alignment_sequence)
        result.reference_start = aligned_segment.reference_start
        result.reference_end = aligned_segment.reference_end
        result.is_reverse = aligned_segment.is_reverse
        result.read_id = aligned_segment.query_name

        return result

    def get_cigar_string(self) -> str:
        """
            @returns cigar_string
        """
        return self.cigar_string

    def get_cigar_tuples(self) -> list[tuple[int, int]]:
        """
            @returns cigar_tuples
        """
        return self.cigar_tuples

    def get_query_alignment_start(self) -> int:
        """
            @returns query_alignment_start
        """
        return self.query_alignment_start

    def get_query_alignment_end(self) -> int:
        """
            @returns query_alignment_end
        """
        return self.query_alignment_end

    def get_query_sequence_length(self) -> int:
        """
            @returns query_sequence_length
        """
        return self.query_sequence_length

    def get_query_sequence(self) -> str:
        """
            @returns query_sequence
        """
        return self.query_sequence

    def get_query_alignment_sequence(self) -> str:
        """
            @returns query_alignment_sequence
        """
        return self.query_alignment_sequence

    def get_reference_start(self) -> int:
        """
            @returns reference_start
        """
        return self.reference_start

    def get_reference_end(self) -> int:
        """
            @returns reference_end
        """
        return self.reference_end

    def get_is_reverse(self) -> int:
        """
            @returns is_reverse
        """
        return self.is_reverse

    def get_read_id(self) -> str:
        """
            @returns read_id
        """
        return str(self.read_id)
