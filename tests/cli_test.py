"""
    @does test CLI
"""

import unittest
import os

from genomic import CLI


class CLITest(unittest.TestCase):
    """
        @does test CLI
    """

    def test_workflow(self):
        """
            @does test
        """
        cli = CLI()
        config = cli.parse_args(["--verbose"])
        self.assertEqual("INFO", config.loglevel)

        config = cli.parse_args(["--silent"])
        self.assertEqual("ERROR", config.loglevel)

        config = cli.parse_args([])
        self.assertEqual("WARNING", config.loglevel)

        config = cli.parse_args(["-j", "1"])
        self.assertEqual(1, config.jobs)

        config = cli.parse_args([])
        self.assertEqual(os.cpu_count(), config.jobs)
