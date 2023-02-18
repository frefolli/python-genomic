from enum import Enum
from pysam import AlignedSegment
from lib import Intron
from lib import Reference



class Alignment:
    class CigarOperation(Enum):
        M = 0
        I = 1
        D = 2
        N = 3
        S = 4
        H = 5
        P = 6
        E = 7
        X = 8

    CIGAR_OPERATIONS_THAT_CONSUME_QUERY = [
        CigarOperation.M,
        CigarOperation.I,
        CigarOperation.S,
        CigarOperation.E,
        CigarOperation.X
    ]
    
    def __init__(self, id : int, alignment : AlignedSegment):
        self.id = id
        self.alignment = alignment

    def split_read(self, reference_sequence : str) -> tuple[str, str, Intron, str, str]:
        cigar_tuples = [(Alignment.CigarOperation(op[0]), op[1]) for op in self.alignment.cigartuples]
        intron_op_index = 0
        while (cigar_tuples[intron_op_index][0] != Alignment.CigarOperation.N):
            intron_op_index += 1
        first_exon_length = sum([op[1] for op in cigar_tuples[:intron_op_index]])
        first_aligned_exon_length = sum([op[1] if op[0] in Alignment.CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[:intron_op_index]])
        second_exon_length = sum([op[1] for op in cigar_tuples[intron_op_index+1:]])
        second_aligned_exon_length = sum([op[1] if op[0] in Alignment.CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[intron_op_index+1:]])
        intron_length = cigar_tuples[intron_op_index][1]
        first_exon = reference_sequence[:first_exon_length]
        first_aligned_exon = self.alignment.query_sequence[:first_aligned_exon_length]
        intron = reference_sequence[first_exon_length:first_exon_length+intron_length]
        second_exon = reference_sequence[-second_exon_length:]
        second_aligned_exon = self.alignment.query_sequence[-second_aligned_exon_length:]
        return (first_exon, first_aligned_exon, Intron(intron), second_exon, second_aligned_exon)

    def process_alignment(self, reference : Reference) -> list[int, str, str, str, str, str, str, int, bool]:
        reference_slice = reference.get_slice(self.alignment.reference_start, self.alignment.reference_end, self.alignment.is_reverse)
        ID = self.id
        CIGAR_STRING = self.alignment.cigarstring
        QUERY_SEQUENCE = self.alignment.query_alignment_sequence # remove SC
        REFERENCE_SEQUENCE = reference_slice
        (first_exon, first_aligned_exon, intron, second_exon, second_aligned_exon) = self.split_read(REFERENCE_SEQUENCE)
        FIRST_ALIGNED_EXON = first_aligned_exon
        INTRON = intron.sequence
        SECOND_ALIGNED_EXON = second_aligned_exon
        INTRON_LENGTH = len(INTRON)
        INTRON_IS_CANONIC = intron.is_canonic()
        return [ID, CIGAR_STRING, QUERY_SEQUENCE, REFERENCE_SEQUENCE, FIRST_ALIGNED_EXON, INTRON, SECOND_ALIGNED_EXON, INTRON_LENGTH, INTRON_IS_CANONIC]
    
    def __str__(self) -> str:
        return f"Alignment(id = {self.id}, alignment = {self.alignment})"