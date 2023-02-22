"""
    @does define BamFile abstraction
"""

from typing import Generator
import re
import logging

import pysam


class BamFile:
    """
    @does define BamFile abstraction
    """

    def __init__(self, path: str):
        """
            @does initialize BamFile
        """
        self.path = path
        pysam.index(path)
        self.alignment_file = pysam.AlignmentFile(path, "rb")
        self.alignment_header = self.alignment_file.header

    def get_alignments(self) -> Generator[pysam.AlignedSegment, None, None]:
        """
            @does returns a generator for getting all alignments
        """
        for read in self.alignment_file.fetch():
            yield read

    def get_alignments_matching_cigarstring(
            self,
            match_rule: str) -> Generator[pysam.AlignedSegment, None, None]:
        """
            @does filter alignments with a regex
        """
        match_rule_engine = re.compile(match_rule)
        for alignment in self.get_alignments():
            if match_rule_engine.match(alignment.cigarstring):
                yield alignment

    def __enter__(self):
        """
            @does enter for with-as
        """
        logging.info("OPENED bamfile %s", self.path)
        return self

    def __exit__(self, exception_type,
                 exception_value, exception_traceback):
        """
            @does exit for with-as
        """
        self.alignment_file.close()
        logging.info("CLOSED bamfile %s", self.path)

    def __str__(self) -> str:
        """
            @returns string representation of BamFile
        """
        return f"BamFile(path = \"{self.path}\")"
