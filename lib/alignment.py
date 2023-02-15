from pysam import AlignedSegment
from lib import Intron
from lib import Reference

class Alignment:
    def __init__(self, id : int, alignment : AlignedSegment, reference : Seq):
        self.id = id
        self.alignment = alignment

def split_read(self, query_sequence : str, reference_sequence : str, cigar_tuples : list):
    cigar_tuples = [(Alignment.CigarOperation(op[0]), op[1]) for op in cigar_tuples]
    intron_op_index = 0
    while (cigar_tuples[intron_op_index][0] != Alignment.CigarOperation.N):
        intron_op_index += 1
    first_exon_length = sum([op[1] for op in cigar_tuples[:intron_op_index]])
    first_aligned_exon_length = sum([op[1] if op[0] in Alignment.CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[:intron_op_index]])
    second_exon_length = sum([op[1] for op in cigar_tuples[intron_op_index+1:]])
    second_aligned_exon_length = sum([op[1] if op[0] in Alignment.CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[intron_op_index+1:]])
    intron_length = cigar_tuples[intron_op_index][1]
    first_exon = reference_sequence[:first_exon_length]
    first_aligned_exon = query_sequence[:first_aligned_exon_length]
    intron = reference_sequence[first_exon_length:first_exon_length+intron_length]
    second_exon = reference_sequence[-second_exon_length:]
    second_aligned_exon = query_sequence[-second_aligned_exon_length:]
    return (first_exon, first_aligned_exon, Intron(intron), second_exon, second_aligned_exon)

def process_alignment(self, id : int, reference : Reference, alignment : AlignedSegment):
    reference_slice = reference.get_slice(alignment.reference_start, alignment.reference_end, alignment.is_reverse)
    ID = id
    CIGAR_STRING = alignment.cigarstring
    QUERY_SEQUENCE = alignment.query_sequence
    REFERENCE_SEQUENCE = reference_slice
    (first_exon, first_aligned_exon, intron, second_exon, second_aligned_exon) = split_read(QUERY_SEQUENCE, REFERENCE_SEQUENCE, alignment.cigartuples)
    FIRST_ALIGNED_EXON = first_aligned_exon
    INTRON = intron.sequence
    SECOND_ALIGNED_EXON = second_aligned_exon
    INTRON_LENGTH = len(INTRON)
    INTRON_IS_CANONIC = intron.is_canonic()
    return [ID, CIGAR_STRING, QUERY_SEQUENCE, REFERENCE_SEQUENCE, FIRST_ALIGNED_EXON, INTRON, SECOND_ALIGNED_EXON, INTRON_LENGTH, INTRON_IS_CANONIC]