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
        self.__worker_id = worker_id
        self.__aligned_segment = None
        self.__reference = None

    @staticmethod
    def reverse_complement(seq : str) -> str:
        '''
            @does compute R&C of the seq
        '''
        reverse = seq[::-1]
        return (reverse.translate(str.maketrans('ATCG', "TAGC")))
    
    def get_worker_id(self):
        """
            @does return the worker id
        """
        return self.__worker_id
    
    def get_aligned_segment(self):
        """
            @does return the aligned_segment obj of the worker
        """
        return self.__aligned_segment
    
    def get_reference(self):
        """
            @does return the reference obj of the worker
        """
        return self.__reference
    
    def set_aligned_segment(self, aligned_segment: AlignedSegment):
        """
            @does assigns aligned_segment
        """
        self.__aligned_segment = aligned_segment

    def set_reference(self, reference: Reference):
        """
            @does assigns reference
        """
        self.__reference = reference

    def get_cigar_tuples(self) -> list[tuple[CigarOperation, int]]:
        """
        @uses pysam.cigar_tuples
        @returns Alignment.cigar_tuples
        """
        return [(CigarOperation(op[0]), op[1])
                for op in self.__aligned_segment.get_cigar_tuples()]

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

    def get_left_aligned_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns left_aligned_exon_length
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY
                    else 0 for op in cigar_tuples[:intronic_site]])

    def get_left_aligned_exon_range(self) -> tuple[int, int]:
        """
            @uses left_aligned_exon_length
            @uses query_alignment_start
            @returns left_aligned_exon_range
        """
        left_aligned_exon_length = self.get_left_aligned_exon_length()
        head_soft_clipping_length = (
            self.__aligned_segment.get_query_alignment_start())
        return (head_soft_clipping_length, left_aligned_exon_length)

    def get_right_aligned_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns right_aligned_exon_length
        """
        cigar_tuples = self.get_cigar_tuples()
        intronic_site = self.get_intronic_site()
        return sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY
                    else 0 for op in cigar_tuples[intronic_site+1:]])

    def get_right_aligned_exon_range(self) -> tuple[int, int]:
        """
            @uses right_aligned_exon_length
            @uses query_sequence_length
            @uses query_alignment_end
            @returns right_aligned_exon_range
        """
        right_aligned_exon_length = self.get_right_aligned_exon_length()
        query_alignment_end = self.__aligned_segment.get_query_alignment_end()
        query_alignment_length = (
            self.__aligned_segment.get_query_sequence_length())
        return (query_alignment_length - right_aligned_exon_length,
                query_alignment_end)

    def get_left_aligned_exon(self):
        """
            @uses left_aligned_exon_range
            @returns left_aligned_exon
        """
        left_aligned_exon_range = self.get_left_aligned_exon_range()
        return (self.__aligned_segment.get_query_sequence()
                [left_aligned_exon_range[0]:left_aligned_exon_range[1]])

    def get_right_aligned_exon(self):
        """
            @uses right_aligned_exon_range
            @returns right_aligned_exon
        """
        right_aligned_exon_range = self.get_right_aligned_exon_range()
        return (self.__aligned_segment.get_query_sequence()
                [right_aligned_exon_range[0]:right_aligned_exon_range[1]])

    def get_left_exon_length(self) -> int:
        """
            @uses Alignment.cigar_tuples
            @uses intronic_site
            @returns left_exon_length
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
            @uses left_exon_length
            @uses intron_length
            @returns intron_range
        """
        reference_start = self.__aligned_segment.get_reference_start()
        left_exon_length = self.get_left_exon_length()
        intron_length = self.get_intron_length()
        return (reference_start + left_exon_length,
                reference_start + left_exon_length + intron_length)

    def get_intron(self) -> str:
        """
            @uses intron_range
            @uses is_reverse
            @returns intron
        """
        intron_range = self.get_intron_range()
        is_reverse = self.__aligned_segment.get_is_reverse()
        return self.__reference.get_slice(intron_range[0],
                                        intron_range[1],
                                        is_reverse)

    def split_read(self) -> tuple[str, str, str]:
        """
            @uses reference
            @uses left_aligned_exon
            @uses right_aligned_exon
            @uses Seq
            @uses intron
        """
        left_aligned_exon = self.get_left_aligned_exon()
        right_aligned_exon = self.get_right_aligned_exon() 
        intron = self.get_intron()

        is_reverse = self.__aligned_segment.get_is_reverse()
        
        first_aligned_exon = (self.reverse_complement(right_aligned_exon) 
                              if is_reverse else left_aligned_exon)
        second_aligned_exon = (self.reverse_complement(left_aligned_exon) 
                              if is_reverse else right_aligned_exon)
        return (first_aligned_exon, intron, second_aligned_exon)
            

    def process_alignment(self) -> list:
        """"
            @does process alignment and produce csv record
        """
        (first_aligned_exon, intron, second_aligned_exon) = self.split_read()
        return [
            self.__aligned_segment.get_read_id(),
            self.__aligned_segment.get_cigar_string(),
            self.__aligned_segment.get_query_sequence(),
            self.__reference.get_slice(
                self.__aligned_segment.get_reference_start(),
                self.__aligned_segment.get_reference_end()),
            self.__aligned_segment.get_is_reverse(),
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
        return (f"Alignment(id = {self.__worker_id}, " +
                f"alignment = {self.__aligned_segment})")
