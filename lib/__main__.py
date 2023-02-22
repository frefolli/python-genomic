#!/usr/bin/env python3
import argparse
from multiprocessing.pool import ThreadPool
import logging
import os
import datetime

from tqdm import tqdm

from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from lib import AlignmentWorker
from lib import AlignedSegment
from lib import Reference

def work_on_alignment(reference : Reference, alignment : AlignedSegment):
    worker = AlignmentWorker(0)
    worker.set_aligned_segment(AlignedSegment.from_pysam_alignment_segment(alignment))
    worker.set_reference(reference)
    return worker.process_alignment()

# sample of workflow
def workflow(
        bamfile_path : str, fastafile_path : str,
        reference_name : str, csvfile_path : str,
        jobs : int):
    logging.info(f"STARTED AT {datetime.datetime.now()}")
    bamfile : BamFile
    with BamFile(bamfile_path) as bamfile:
        logging.info(f"OPENED bamfile {bamfile_path}")
        fastafile : FastaFile
        with FastaFile(fastafile_path) as fastafile:
            logging.info(f"OPENED fastafile {fastafile_path}")
            reference = fastafile.get_reference(reference_name)
            csvfile : CsvFile
            with CsvFile(csvfile_path) as csvfile:
                logging.info(f"OPENED csvfile {csvfile_path}")
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
                logging.info(f"GOT {len(matching)} matching_alignments")
                worker = (lambda alignment : work_on_alignment(reference, alignment))
                with tqdm(total=len(matching), desc=f"processing alignments") as progress:
                    logging.info(f"OPENED tqdm")
                    with ThreadPool(jobs) as pool:
                        logging.info(f"OPENED thread_pool WITH {jobs} jobs")
                        for result in pool.imap(worker, matching):
                            csvfile.add_line(result)
                            progress.update(1)
                csvfile.save()
                logging.info(f"SAVED csvfile {csvfile_path}")
            logging.info(f"CLOSED csvfile {csvfile_path}")
        logging.info(f"CLOSED fastafile {fastafile_path}")
    logging.info(f"CLOSED bamfile {bamfile_path}")
    logging.info(f"FINISHED AT {datetime.datetime.now()}")

def command_line_interface():
    argument_parser = argparse.ArgumentParser(
        prog = "main",
        description = "this python module identifies introns and exons in reads when given FASTA and BAM",
        epilog = "log level is applied using the priority queue : --verbose => --silent => --loglevel")
    
    argument_parser.add_argument(
        "-j", "--jobs",
        type = int, default = None,
        help = "supply number of running jobs, default = # of CPU Threads")
    
    argument_parser.add_argument(
        "-b", "--bam",
        type = str, default = "sample.bam",
        help = "supply BAM file path, default = 'sample.bam'")
    
    argument_parser.add_argument(
        "-f", "--fasta",
        type = str, default = "BDGP6.X.fasta",
        help = "supply FASTA file path, default = 'BDGP6.X.fasta'")
    
    argument_parser.add_argument(
        "-c", "--chromosome",
        type = str, default = "X",
        help = "supply chromosome name, default = 'X'")
    
    argument_parser.add_argument(
        "-o", "--output",
        type = str, default = "reads.csv",
        help = "supply output csv file path, default = 'reads.csv'")
    
    argument_parser.add_argument(
        "-r", "--logfile",
        type = str, default = None,
        help = "supply log redirection file path, default = None; (= to stdout)")
    
    argument_parser.add_argument(
        "-l", "--loglevel",
        type = str, default = "WARNING",
        choices = ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help = "supply log level, default = 'WARNING'")
    
    argument_parser.add_argument(
        "-s", "--silent",
        action='store_true',
        help = "set log level to ERROR")
    
    argument_parser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help = "set log level to DEBUG")

    config = argument_parser.parse_args()
    apply_config(config)

def apply_config(config):
    bam = config.bam
    fasta = config.fasta
    chromosome = config.chromosome
    csv = config.output

    loglevel = config.loglevel
    if config.silent:
        loglevel = "ERROR"
    if config.verbose:
        loglevel = "INFO"

    logfile = config.logfile
    logging.basicConfig(filename=logfile, level=loglevel)
    
    jobs = config.jobs
    if (jobs is None):
        jobs = os.cpu_count()
    elif (jobs < 1):
        jobs = os.cpu_count()
        logging.warning("DEFAULTED jobs")

    workflow(bam, fasta, chromosome, csv, jobs)

if __name__ == "__main__":
    command_line_interface()