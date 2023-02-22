"""
    @does define Command Line Interface abstraction
"""

import argparse
import os


class CLI:
    """
        @does define Command Line Interface abstraction
    """

    def __init__(self) -> None:
        """
            @does initialize CLI
        """
        self.argument_parser = argparse.ArgumentParser(
            prog="main",
            description=("this python module identifies " +
                         "introns and exons in reads " +
                         "when given FASTA and BAM"),
            epilog=("log level is applied using" +
                    "the priority queue : " +
                    "--verbose => --silent => --loglevel"))

        self.argument_parser.add_argument(
            "-j", "--jobs",
            type=int, default=None,
            help="supply number of running jobs, default = # of CPU Threads")

        self.argument_parser.add_argument(
            "-b", "--bam",
            type=str, default="sample.bam",
            help="supply BAM file path, default = 'sample.bam'")

        self.argument_parser.add_argument(
            "-f", "--fasta",
            type=str, default="BDGP6.X.fasta",
            help="supply FASTA file path, default = 'BDGP6.X.fasta'")

        self.argument_parser.add_argument(
            "-c", "--chromosome",
            type=str, default="X",
            help="supply chromosome name, default = 'X'")

        self.argument_parser.add_argument(
            "-o", "--output",
            type=str, default="reads.csv",
            help="supply output csv file path, default = 'reads.csv'")

        self.argument_parser.add_argument(
            "-r", "--logfile",
            type=str, default=None,
            help="supply log redirection file path, " +
                 "default = None; (= to stdout)")

        self.argument_parser.add_argument(
            "-l", "--loglevel",
            type=str, default="WARNING",
            choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
            help="supply log level, default='WARNING'")

        self.argument_parser.add_argument(
            "-s", "--silent",
            action='store_true',
            help="set log level to ERROR")

        self.argument_parser.add_argument(
            "-v", "--verbose",
            action='store_true',
            help="set log level to DEBUG")

    def parse_args(self, args: list[str] = None) -> argparse.Namespace:
        """
            @does parse arguments with given arg list
            if args = None, uses sys.argv[1:]
        """
        config = self.argument_parser.parse_args(args)

        if config.verbose:
            config.loglevel = "INFO"
        elif config.silent:
            config.loglevel = "ERROR"
        if config.jobs is None or config.jobs < 1:
            config.jobs = os.cpu_count()
        config.csv = config.output

        return config
