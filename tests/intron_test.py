import unittest
from lib import Intron

class CsvFileTest(unittest.TestCase):
    def test_workflow(self):
        self.assertFalse(Intron("ggag").is_canonical())
        self.assertFalse(Intron("gag").is_canonical())
        self.assertTrue(Intron("gtag").is_canonical())
        self.assertEqual(f"Intron(sequence = \"acgt\")", str(Intron("acgt")))