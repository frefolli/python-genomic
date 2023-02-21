import pysam

"""
    @represents a useful abstraction for pysam.AlignedSegment
"""
class AlignedSegment:
    @staticmethod
    def from_properties(
        query_sequence = None,
        cigartuples = None,
        query_sequence_length = None,
        query_alignment_start = None,
        query_alignment_end = None,
        reference_start = None,
        reference_end = None,
        is_reverse = None):
        result = AlignedSegment()

        if (not query_sequence): query_sequence = ""
        if (not cigartuples): cigartuples = [(0, len(query_sequence))]
        if (not query_sequence_length): query_sequence_length = sum([op[1] for op in cigartuples])
        if (not query_alignment_start): query_alignment_start = 0
        if (not query_alignment_end): query_alignment_end = query_sequence_length
        if (not reference_start): reference_start = 0
        if (not reference_end): reference_end = 0
        if (not is_reverse): is_reverse = False

        result.cigartuples = cigartuples
        result.query_alignment_start = query_alignment_start
        result.query_alignment_end = query_alignment_end
        result.query_sequence_length = query_sequence_length
        result.query_sequence = query_sequence
        result.query_alignment_sequence = query_sequence[query_alignment_start: query_alignment_end]
        result.reference_start = reference_start
        result.reference_end = reference_end
        result.is_reverse = is_reverse

        return result

    @staticmethod
    def from_pysam_alignment_segment(aligned_segment : pysam.AlignedSegment):
        result = AlignedSegment()

        result.cigartuples = aligned_segment.cigartuples
        result.query_alignment_start = aligned_segment.query_alignment_start
        result.query_alignment_end = aligned_segment.query_alignment_end
        result.query_sequence_length = len(aligned_segment.query_sequence)
        result.query_sequence = aligned_segment.query_sequence
        result.query_alignment_sequence = aligned_segment.query_alignment_sequence
        result.reference_start = aligned_segment.reference_start
        result.reference_end = aligned_segment.reference_end
        result.is_reverse = aligned_segment.is_reverse

        return result

    def get_cigar_tuples(self) -> list[tuple[int, int]]:
        return self.cigartuples

    def get_query_alignment_start(self) -> int:
        return self.query_alignment_start

    def get_query_alignment_end(self) -> int:
        return self.query_alignment_end
    
    def get_query_sequence_length(self) -> int:
        return self.query_sequence_length
    
    def get_query_sequence(self) -> str:
        return self.query_sequence
    
    def get_query_alignment_sequence(self) -> str:
        return self.query_alignment_sequence
    
    def get_reference_start(self) -> int:
        return self.reference_start
    
    def get_reference_end(self) -> int:
        return self.reference_end
    
    def get_is_reverse(self) -> int:
        return self.is_reverse