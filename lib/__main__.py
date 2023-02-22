#!/usr/bin/env python3
"""
    @uses FastaFile
    @uses BamFile
    @uses CsvFile
    @uses AlignmentWorker
    @uses AlignedSegment
    @uses Reference
    @does application
"""

import argparse
from multiprocessing.pool import ThreadPool
import logging
import os
import datetime
from functools import partial

from tqdm import tqdm
import pysam

from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from lib import AlignmentWorker
from lib import AlignedSegment
from lib import Reference


def work_on_alignment(reference: Reference, alignment: pysam.AlignedSegment):
    """
        @does setup worker and apply work force
    """
    worker = AlignmentWorker(0)
    aligned_segment = AlignedSegment.from_pysam_alignment_segment(alignment)
    worker.set_aligned_segment(aligned_segment)
    worker.set_reference(reference)
    return worker.process_alignment()


def unpack(iterator) -> list:
    """
        @does unpack iterator
    """
    return [_ for _ in iterator]


def workflow(
        bamfile_path: str, fastafile_path: str,
        reference_name: str, csvfile_path: str,
        jobs: int):
    """
        @does default workflow
    """
    logging.info("STARTED AT %s", datetime.datetime.now())
    bamfile: BamFile
    with BamFile(bamfile_path) as bamfile:
        logging.info("OPENED bamfile %s", bamfile_path)

        fastafile: FastaFile
        with FastaFile(fastafile_path) as fastafile:
            logging.info("OPENED fastafile %s", fastafile_path)
            reference = fastafile.get_reference(reference_name)

            csvfile: CsvFile
            with CsvFile(csvfile_path) as csvfile:
                logging.info("OPENED csvfile %s", csvfile_path)
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
                matching_iterator = (
                    bamfile.get_alignments_matching_cigarstring(match_rule))
                matching = unpack(matching_iterator)
                logging.info("GOT %s matching_alignments", len(matching))
                worker = partial(work_on_alignment, reference)

                with tqdm(total=len(matching),
                          desc="processing alignments") as progress:
                    logging.info("OPENED tqdm")

                    with ThreadPool(jobs) as pool:
                        logging.info("OPENED thread_pool WITH %s jobs", jobs)
                        for result in pool.imap(worker, matching):
                            csvfile.add_line(result)
                            progress.update(1)
                csvfile.save()
                logging.info("SAVED csvfile %s", csvfile_path)
            logging.info("CLOSED csvfile %s", csvfile_path)
        logging.info("CLOSED fastafile %s", fastafile_path)
    logging.info("CLOSED bamfile %s", bamfile_path)
    logging.info("FINISHED AT %s", datetime.datetime.now())


def command_line_interface():
    """
        @does parse command line arguments
    """
    argument_parser = argparse.ArgumentParser(
        prog="main",
        description=("this python module identifies" +
                     "introns and exons in reads when given FASTA and BAM"),
        epilog=("log level is applied using" +
                "the priority queue : --verbose => --silent => --loglevel"))

    argument_parser.add_argument(
        "-j", "--jobs",
        type=int, default=None,
        help="supply number of running jobs, default = # of CPU Threads")

    argument_parser.add_argument(
        "-b", "--bam",
        type=str, default="sample.bam",
        help="supply BAM file path, default = 'sample.bam'")

    argument_parser.add_argument(
        "-f", "--fasta",
        type=str, default="BDGP6.X.fasta",
        help="supply FASTA file path, default = 'BDGP6.X.fasta'")

    argument_parser.add_argument(
        "-c", "--chromosome",
        type=str, default="X",
        help="supply chromosome name, default = 'X'")

    argument_parser.add_argument(
        "-o", "--output",
        type=str, default="reads.csv",
        help="supply output csv file path, default = 'reads.csv'")

    argument_parser.add_argument(
        "-r", "--logfile",
        type=str, default=None,
        help="supply log redirection file path, default = None; (= to stdout)")

    argument_parser.add_argument(
        "-l", "--loglevel",
        type=str, default="WARNING",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        help="supply log level, default='WARNING'")

    argument_parser.add_argument(
        "-s", "--silent",
        action='store_true',
        help="set log level to ERROR")

    argument_parser.add_argument(
        "-v", "--verbose",
        action='store_true',
        help="set log level to DEBUG")

    config = argument_parser.parse_args()
    apply_config(config)


def apply_config(config):
    """
        @uses cli config
        @does start workflow
    """
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
    if jobs is None:
        jobs = os.cpu_count()
    elif jobs < 1:
        jobs = os.cpu_count()
        logging.warning("DEFAULTED jobs")

    workflow(bam, fasta, chromosome, csv, jobs)


if __name__ == "__main__":
    command_line_interface()
