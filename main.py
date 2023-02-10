#!/usr/bin/env python3
from lib import FastaFile
from lib import BamFile

# sample of workflow
def workflow():
    bamfile : BamFile
    with BamFile('./sample.bam') as bamfile:
        reference : FastaFile
        with FastaFile('GRCh37.X.fasta') as reference:
            chrX = reference.get_reference('X')
            print(chrX[10304:170701])

if __name__ == "__main__":
    workflow()