"""
    @does test Intron
"""

import unittest
from lib import Intron


class IntronTest(unittest.TestCase):
    """
        @does test Intron
    """

    def test_workflow(self):
        """
            @does test
        """
        self.assertFalse(Intron("ggag").is_canonical())
        self.assertFalse(Intron("gag").is_canonical())
        self.assertTrue(Intron("gtag").is_canonical())
        self.assertEqual("Intron(sequence = \"acgt\")",
                         str(Intron("acgt")))
