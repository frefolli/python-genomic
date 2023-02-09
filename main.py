#!/usr/bin/env python3
# import pysam
import pysam, re

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
    print(f"Table Length: {len(candidates)}")
    for i in range(20):
        print(candidates[i])
    print("... ... ... ...")
