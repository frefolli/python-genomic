#!/usr/bin/env python3
from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from tqdm import tqdm
from lib import Alignment
from lib import AlignedSegment
import pysam

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
                matching = [_ for _ in enumerate(bamfile.get_alignments_matching_cigarstring(match_rule))]
                with tqdm(total=len(matching), desc=f"processing alignments") as progress:
                    align = None
                    for (id, alignment) in matching:
                        align = Alignment(id, AlignedSegment.from_pysam_alignment_segment(alignment))
                        progress.update(1)
                        csvfile.add_line(align.process_alignment(reference))
                csvfile.save()

if __name__ == "__main__":
    workflow('sample.bam', 'BDGP6.X.fasta', 'X')