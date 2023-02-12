#!/usr/bin/env python3
from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from tqdm import tqdm
from enum import Enum

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

def split_read(query_sequence : str, reference_sequence : str, cigar_tuples : list):
    cigar_tuples = [(CigarOperation(op[0]), op[1]) for op in cigar_tuples]
    intron_op_index = 0
    while (cigar_tuples[intron_op_index][0] != CigarOperation.N):
        intron_op_index += 1
    first_exon_length = sum([op[1] for op in cigar_tuples[:intron_op_index]])
    first_aligned_exon_length = sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[:intron_op_index]])
    second_exon_length = sum([op[1] for op in cigar_tuples[intron_op_index+1:]])
    second_aligned_exon_length = sum([op[1] if op[0] in CIGAR_OPERATIONS_THAT_CONSUME_QUERY else 0 for op in cigar_tuples[intron_op_index+1:]])
    intron_length = cigar_tuples[intron_op_index][1]
    first_exon = reference_sequence[first_exon_length:]
    first_aligned_exon = query_sequence[first_aligned_exon_length:]
    intron = reference_sequence[first_exon_length:first_exon_length+intron_length]
    second_exon = reference_sequence[:-second_exon_length]
    second_aligned_exon = query_sequence[:-second_aligned_exon_length]
    return (first_exon, first_aligned_exon, intron, second_exon, second_aligned_exon)

def process_alignment(id, reference, alignment):
    reference_slice = reference.get_slice(alignment.reference_start, alignment.reference_end)
    ID = id
    CIGAR_STRING = alignment.cigarstring
    QUERY_SEQUENCE = alignment.query_sequence
    REFERENCE_SEQUENCE = reference_slice.seq
    (first_exon, first_aligned_exon, intron, second_exon, second_aligned_exon) = split_read(QUERY_SEQUENCE, REFERENCE_SEQUENCE, alignment.cigartuples)
    FIRST_ALIGNED_EXON = first_aligned_exon
    INTRON = intron
    SECOND_ALIGNED_EXON = second_aligned_exon
    INTRON_LENGTH = len(INTRON)
    return [ID, CIGAR_STRING, QUERY_SEQUENCE, REFERENCE_SEQUENCE, FIRST_ALIGNED_EXON, INTRON, SECOND_ALIGNED_EXON, INTRON_LENGTH]

# sample of workflow
def workflow(bamfile_path : str, fastafile_path : str, reference_name : str):
    bamfile : BamFile
    with BamFile(bamfile_path) as bamfile:
        fastafile : FastaFile
        with FastaFile(fastafile_path) as fastafile:
            reference = fastafile.get_reference(reference_name)
            csvfile : CsvFile
            with CsvFile("reads.csv") as csvfile:
                csvfile.rebase(["id", "cigar_string", "query_sequence", "reference_sequence", "first_aligned_exon", "intron", "second_aligned_exon", "intron_length"])
                match_rule = "([0-9]+[MIDSHP=X])*[0-9]+N([0-9]+[MIDSHP=X])*"
                matching = [_ for _ in enumerate(bamfile.get_alignments_matching_cigarstring(match_rule))]
                with tqdm(total=len(matching), desc=f"processing alignments") as progress:
                    for (id, alignment) in matching:
                        progress.update(1)
                        csvfile.add_line(process_alignment(id, reference, alignment))
                csvfile.save()

if __name__ == "__main__":
    workflow('sample.bam', 'BDGP6.X.fasta', 'X')
