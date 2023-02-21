from lib import AlignedSegment
from lib import Intron
from lib import Reference
from lib import CigarOperation, CIGAR_OPERATIONS_THAT_CONSUME_QUERY, CIGAR_OPERATIONS_THAT_CONSUME_REFERENCE

"""
    @does parse introns and exons from AlignedSegments
"""
class Alignment:
    def __init__(self, id : int, aligned_segment : AlignedSegment):
        self.id = id
        self.aligned_segment = aligned_segment
    
    """
        @uses pysam.cigar_tuples
        @returns Alignment.cigar_tuples
    """
    def get_cigar_tuples(self) -> list[tuple[CigarOperation, int]]:
        return [(CigarOperation(op[0]), op[1]) for op in self.aligned_segment.get_cigar_tuples()]

    """
        @uses Alignment.cigar_tuples
        @returns intronic_site
    """
    def get_intronic_site(self) -> int:
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = 0
        while (cigar_tuples[intronic_site][0] != CigarOperation.N):
            intronic_site += 1
        return intronic_site
    
    """
        @uses Alignment.cigar_tuples
        @uses intronic_site
        @returns first_aligned_exon_length
    """
    def get_first_aligned_exon_length(self) -> int:
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[:intronic_site]])

    """
        @uses first_aligned_exon_length
        @uses query_alignment_start
        @returns first_aligned_exon_range
    """
    def get_first_aligned_exon_range(self) -> tuple[int, int]:
        first_aligned_exon_length = self.get_first_aligned_exon_length()
        head_soft_clipping_length = self.aligned_segment.get_query_alignment_start()
        return (head_soft_clipping_length, first_aligned_exon_length)

    """
        @uses Alignment.cigar_tuples
        @uses intronic_site
        @returns second_aligned_exon_length
    """
    def get_second_aligned_exon_length(self) -> int:
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[intronic_site+1:]])

    """
        @uses second_aligned_exon_length
        @uses query_sequence_length
        @uses query_alignment_end
        @returns second_aligned_exon_range
    """
    def get_second_aligned_exon_range(self) -> tuple[int, int]:
        second_aligned_exon_length = self.get_second_aligned_exon_length()
        query_sequence_length = self.aligned_segment.get_query_sequence_length()
        query_alignment_end = self.aligned_segment.get_query_alignment_end()
        return (query_sequence_length - second_aligned_exon_length, query_alignment_end)
    
    """
        @uses first_aligned_exon_range
        @returns first_aligned_exon
    """
    def get_first_aligned_exon(self):
        first_aligned_exon_range = self.get_first_aligned_exon_range()
        return self.aligned_segment.get_query_sequence()[first_aligned_exon_range[0]:first_aligned_exon_range[1]]
    
    """
        @uses second_aligned_exon_range
        @returns second_aligned_exon
    """
    def get_second_aligned_exon(self):
        second_aligned_exon_range = self.get_second_aligned_exon_range()
        return self.aligned_segment.get_query_sequence()[second_aligned_exon_range[0]:second_aligned_exon_range[1]]

    """
        @uses Alignment.cigar_tuples
        @uses intronic_site
        @returns first_exon_length
    """
    def get_first_exon_length(self) -> int:
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_REFERENCE else 0 for op in cigar_tuples[:intronic_site]])

    """
        @uses
        @returns intron_range
    """
    def get_intron_length(self) -> int:
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return cigar_tuples[intronic_site][1]

    """
        @uses reference_start
        @uses first_exon_length
        @uses intron_length
        @returns intron_range
    """
    def get_intron_range(self) -> tuple[int, int]:
        reference_start = self.aligned_segment.get_reference_start()
        first_exon_length = self.get_first_exon_length()
        intron_length = self.get_intron_length()
        return (reference_start + first_exon_length, reference_start + first_exon_length + intron_length)

    """
        @uses intron_range
        @uses is_reverse
        @returns intron
    """
    def get_intron(self, reference : Reference) -> str:
        intron_range = self.get_intron_range()
        is_reverse = self.aligned_segment.get_is_reverse()
        return reference.get_slice(intron_range[0], intron_range[1], is_reverse)
    
    """
        @uses reference
        @uses first_aligned_exon
        @uses second_aligned_exon
        @uses intron
    """
    def split_read(self, reference : Reference) -> tuple[str, str, str]:
        first_aligned_exon = self.get_first_aligned_exon()
        second_aligned_exon = self.get_second_aligned_exon()
        intron = self.get_intron(reference)
        return (first_aligned_exon, intron, second_aligned_exon)

    def process_alignment(self, reference : Reference) -> list:
        (first_aligned_exon, intron, second_aligned_exon) = self.split_read(reference)
        return [
            self.id,
            self.aligned_segment.get_cigar_string(),
            self.aligned_segment.get_query_sequence(),
            reference.get_slice(self.aligned_segment.get_reference_start(), self.aligned_segment.get_reference_end()),
            first_aligned_exon,
            intron,
            second_aligned_exon,
            len(intron),
            Intron(intron).is_canonic()
        ]
    
    def __str__(self) -> str:
        return f"Alignment(id = {self.id}, alignment = {self.aligned_segment})"