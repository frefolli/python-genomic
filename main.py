#!/usr/bin/env python3
# import pysam
import pysam, re
from Bio import SeqIO

# define reference import function
def open_reference(path):
    reference = SeqIO.index(path, "fasta")
    return reference['X']

# define samfile import function
def open_samfile(path):
    pysam.index(path)
    samfile = pysam.AlignmentFile(path, "rb")
    return samfile

# extract n-containing alignments
def extract_candidates(samfile):
    for read in samfile.fetch():
        if 'N' in read.cigarstring:
            yield read
 
# workflow
if __name__ == "__main__":
    samfile = open_samfile('./sample.bam')
    candidates = [_ for _ in extract_candidates(samfile)]
    reference = open_reference('GRCh37.X.fasta')

    c = candidates[0]
    print(c.cigarstring)
    print(c.query_sequence.upper())
    print(reference[c.reference_start: c.reference_end].seq.upper())


