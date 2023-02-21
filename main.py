#!/usr/bin/env python3
from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from tqdm import tqdm
from lib import AlignmentWorker
from lib import AlignedSegment
import argparse, sys

# sample of workflow
def workflow(bamfile_path : str, fastafile_path : str, reference_name : str):
    bamfile : BamFile
    with BamFile(bamfile_path) as bamfile:
        fastafile : FastaFile
        with FastaFile(fastafile_path) as fastafile:
            reference = fastafile.get_reference(reference_name)
            csvfile : CsvFile
            with CsvFile("reads.csv") as csvfile:
                csvfile.rebase([
                    "id",
                    "cigar_string",
                    "query_sequence",
                    "reference_sequence",
                    "first_aligned_exon",
                    "intron",
                    "second_aligned_exon",
                    "intron_length",
                    "intron_is_canonic"
                ])
                match_rule = "^([0-9]+[MIDSHP=X])*[0-9]+N([0-9]+[MIDSHP=X])*$"
                matching = [_ for _ in bamfile.get_alignments_matching_cigarstring(match_rule)]
                worker = AlignmentWorker(0)
                with tqdm(total=len(matching), desc=f"processing alignments") as progress:
                    for alignment in matching:
                        progress.update(1)
                        worker.set_aligned_segment(AlignedSegment.from_pysam_alignment_segment(alignment))
                        worker.set_reference(reference)
                        csvfile.add_line(worker.process_alignment())
                csvfile.save()

def command_line_interface():
    argument_parser = argparse.ArgumentParser(
        prog = "main",
        description = "this python module identifies introns and exons in reads when given FASTA and BAM")
    argument_parser.add_argument("-b", "--bam", type = str, default = "sample.bam", help = "supply BAM file path, default = 'sample.bam'")
    argument_parser.add_argument("-f", "--fasta", type = str, default = "BDGP6.X.fasta", help = "supply FASTA file path, default = 'BDGP6.X.fasta'")
    argument_parser.add_argument("-c", "--chromosome", type = str, default = "X", help = "supply chromosome name, default = 'X'")
    config = argument_parser.parse_args()
    
    bam = config.bam
    fasta = config.fasta
    chromosome = config.chromosome
    workflow(bam, fasta, chromosome)

if __name__ == "__main__":
    command_line_interface()