"""
    boh
"""

from lib import AlignedSegment
from lib import Intron
from lib import Reference
from lib import CigarOperation
from lib import CIGAR_OPERATIONS_THAT_CONSUME_QUERY
from lib import CIGAR_OPERATIONS_THAT_CONSUME_REFERENCE


class AlignmentWorker:
    """
        @does parse introns and exons from AlignedSegments
    """

    def __init__(self, worker_id: int):
        """
            @does initialize AlignmentWorker
        """
        self.worker_id = worker_id
        self.aligned_segment = None
        self.reference = None

    def set_aligned_segment(self, aligned_segment: AlignedSegment):
        """
            @does assigns aligned_segment
        """
        self.aligned_segment = aligned_segment

    def set_reference(self, reference: Reference):
        """
            @does assigns reference
        """
        self.reference = reference

    def get_cigar_tuples(self) -> list[tuple[CigarOperation, int]]:
        """
        @uses pysam.cigar_tuples
        @returns Alignment.cigar_tuples
        """
        return [(CigarOperation(op[0]), op[1])
                for op in self.aligned_segment.get_cigar_tuples()]

    def get_intronic_site(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @returns intronic_site
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = 0
        while cigar_tuples[intronic_site][0] != CigarOperation.N:
            intronic_site += 1
        return intronic_site

    def get_first_aligned_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns first_aligned_exon_length
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY
                    else 0 for op in cigar_tuples[:intronic_site]])

    def get_first_aligned_exon_range(self) -> tuple[int, int]:
        """
            @uses first_aligned_exon_length
            @uses query_alignment_start
            @returns first_aligned_exon_range
        """
        first_aligned_exon_length = self.get_first_aligned_exon_length()
        head_soft_clipping_length = (
            self.aligned_segment.get_query_alignment_start())
        return (head_soft_clipping_length, first_aligned_exon_length)

    def get_second_aligned_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns second_aligned_exon_length
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY
                    else 0 for op in cigar_tuples[intronic_site+1:]])

    def get_second_aligned_exon_range(self) -> tuple[int, int]:
        """
            @uses second_aligned_exon_length
            @uses query_sequence_length
            @uses query_alignment_end
            @returns second_aligned_exon_range
        """
        second_aligned_exon_length = self.get_second_aligned_exon_length()
        query_sequence_length = (
            self.aligned_segment.get_query_sequence_length())
        query_alignment_end = self.aligned_segment.get_query_alignment_end()
        return (query_sequence_length - second_aligned_exon_length,
                query_alignment_end)

    def get_first_aligned_exon(self):
        """
            @uses first_aligned_exon_range
            @returns first_aligned_exon
        """
        first_aligned_exon_range = self.get_first_aligned_exon_range()
        return (self.aligned_segment.get_query_sequence()
                [first_aligned_exon_range[0]:first_aligned_exon_range[1]])

    def get_second_aligned_exon(self):
        """
            @uses second_aligned_exon_range
            @returns second_aligned_exon
        """
        second_aligned_exon_range = self.get_second_aligned_exon_range()
        return (self.aligned_segment.get_query_sequence()
                [second_aligned_exon_range[0]:second_aligned_exon_range[1]])

    def get_first_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns first_exon_length
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_REFERENCE
                    else 0 for op in cigar_tuples[:intronic_site]])

    def get_intron_length(self) -> int:
        """
            @uses
            @returns intron_range
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return cigar_tuples[intronic_site][1]

    def get_intron_range(self) -> tuple[int, int]:
        """
            @uses reference_start
            @uses first_exon_length
            @uses intron_length
            @returns intron_range
        """
        reference_start = self.aligned_segment.get_reference_start()
        first_exon_length = self.get_first_exon_length()
        intron_length = self.get_intron_length()
        return (reference_start + first_exon_length,
                reference_start + first_exon_length + intron_length)

    def get_intron(self) -> str:
        """
            @uses intron_range
            @uses is_reverse
            @returns intron
        """
        intron_range = self.get_intron_range()
        is_reverse = self.aligned_segment.get_is_reverse()
        return self.reference.get_slice(intron_range[0],
                                        intron_range[1],
                                        is_reverse)

    def split_read(self) -> tuple[str, str, str]:
        """
            @uses reference
            @uses first_aligned_exon
            @uses second_aligned_exon
            @uses intron
        """
        first_aligned_exon = self.get_first_aligned_exon()
        second_aligned_exon = self.get_second_aligned_exon()
        intron = self.get_intron()
        return (first_aligned_exon, intron, second_aligned_exon)

    def process_alignment(self) -> list:
        """"
            @does process alignment and produce csv record
        """
        (first_aligned_exon, intron, second_aligned_exon) = self.split_read()
        return [
            self.aligned_segment.get_read_id(),
            self.aligned_segment.get_cigar_string(),
            self.aligned_segment.get_query_sequence(),
            self.reference.get_slice(
                self.aligned_segment.get_reference_start(),
                self.aligned_segment.get_reference_end()),
            first_aligned_exon,
            intron,
            second_aligned_exon,
            len(intron),
            Intron(intron).is_canonical()
        ]

    def __str__(self) -> str:
        """
            @does return string representation
        """
        return (f"Alignment(id = {self.worker_id}, " +
                f"alignment = {self.aligned_segment})")
