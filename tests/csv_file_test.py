"""
    @does test CsvFile
"""

import unittest
import tempfile
import os

from lib import CsvFile


class CsvFileTest(unittest.TestCase):
    """
        @does test CsvFile
    """

    def test_workflow(self):
        """
            @does test
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "some.csv")

            with CsvFile(path) as csvfile:
                csvfile.add_column('uno')
                csvfile.remove_column('uno')
                csvfile.add_column('due')
                csvfile.add_line(['1'])
                csvfile.add_line(['2'])
                csvfile.add_line(['3'])
                csvfile.add_column('tre')
                csvfile.remove_line(0)
                csvfile.save()
                self.assertEqual(f"CsvFile(path = \"{path}\")", str(csvfile))

            with CsvFile(path) as csvfile:
                self.assertEqual(['3', ''], csvfile.get_line(len(csvfile)-1))
                csvfile.rebase(['altra'])
                csvfile.change_path(path)
                csvfile.save_as(path)
