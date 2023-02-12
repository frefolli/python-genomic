import unittest
from lib import Intron

class CsvFileTest(unittest.TestCase):
    def test_workflow(self):
        self.assertFalse(Intron("ggag").is_canonic())
        self.assertFalse(Intron("gag").is_canonic())
        self.assertTrue(Intron("gtag").is_canonic())
        self.assertEqual(f"Intron(sequence = \"acgt\")", str(Intron("acgt")))