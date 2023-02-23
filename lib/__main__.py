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

from multiprocessing.pool import ThreadPool
import logging
import datetime
from functools import partial
import argparse

from tqdm import tqdm
import pysam

from lib import FastaFile
from lib import BamFile
from lib import CsvFile
from lib import AlignmentWorker
from lib import AlignedSegment
from lib import Reference
from lib import CLI


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
    return list(iterator)


CIGAR_STRING_REGEX = "^([0-9]+[MIDSHP=X])*[0-9]+N([0-9]+[MIDSHP=X])*$"


def workflow(
        bamfile_path: str, fastafile_path: str,
        reference_name: str, csvfile_path: str,
        intermediate_path: str, jobs: int):
    """
        @does default workflow
    """
    logging.info("STARTED AT %s", datetime.datetime.now())
    bamfile: BamFile
    with BamFile(bamfile_path) as bamfile:

        fastafile: FastaFile
        with FastaFile(fastafile_path) as fastafile:
            reference = fastafile.get_reference(reference_name)

            csvfile: CsvFile
            with CsvFile(intermediate_path) as csvfile:
                csvfile.rebase([
                    "id",
                    "cigar_string",
                    "query_sequence",
                    "reference_sequence",
                    "is_reverse",
                    "first_aligned_exon",
                    "intron",
                    "second_aligned_exon",
                    "intron_length",
                    "intron_is_canonic"
                ])
                matching = unpack(
                    bamfile.get_alignments_matching_cigarstring(
                        CIGAR_STRING_REGEX))
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
    with CsvFile(intermediate_path) as csvfile:
        body = csvfile.body
        logging.info("GOT %s alignments", len(body))

        with CsvFile(csvfile_path) as outputfile:
            outputfile.rebase([
                    "id",
                    "first_exon",
                    "first_bases_of_intron",
                    "last_bases_of_intron",
                    "second_exon",
                    "intron_is_canonic"
                ])

            with tqdm(total=len(body),
                      desc="saving results") as progress:
                logging.info("OPENED tqdm")
                for alignment in body:
                    outputfile.add_line([
                        alignment[0],
                        alignment[5],
                        alignment[6][:20],
                        alignment[6][-20:],
                        alignment[7],
                        alignment[9],
                    ])
                    progress.update(1)
            outputfile.save()
            logging.info("SAVED %s formatted alignments", len(body))
    logging.info("FINISHED AT %s", datetime.datetime.now())


def command_line_interface(args: list[str] = None):
    """
        @does parse command line arguments
    """
    cli = CLI()
    config = cli.parse_args(args)
    apply_config(config)


def apply_config(config: argparse.Namespace):
    """
        @uses cli config
        @does start workflow
    """
    logging.basicConfig(
        filename=config.logfile,
        level=config.loglevel)
    workflow(
        config.bam,
        config.fasta,
        config.chromosome,
        config.output,
        config.intermediate_path,
        config.jobs)


def main():
    """
        @does define entry point
    """
    command_line_interface()


if __name__ == "__main__":
    main()
