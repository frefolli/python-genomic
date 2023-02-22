"""
    @does define Reference abstraction
"""

from Bio import SeqRecord
from Bio import Seq


class Reference:
    """
        @does define Reference abstraction
    """

    @staticmethod
    def from_name_and_sequence(name: str, sequence: str):
        """
            @does build Refernce from string sequence and name
        """
        return Reference(
            name=name,
            sequence=SeqRecord.SeqRecord(Seq.Seq(sequence)))

    def __init__(self, name: str, sequence: SeqRecord.SeqRecord):
        """
            @does initialize Reference
        """
        self.name = name
        self.sequence = sequence

    def get_slice(self,
                  start_index: int = None,
                  end_index: int = None,
                  reverse: bool = False) -> SeqRecord.SeqRecord:
        """
            @returns slice of DNA and Reverse Complement is if reverse = True
        """
        seq = self.sequence.seq
        if start_index is None:
            start_index = 0
        if end_index is None:
            seq = seq[start_index:]
        else:
            seq = seq[start_index: end_index]
        return seq.reverse_complement() if reverse else seq

    def __str__(self) -> str:
        """
            @returns string representation of Reference
        """
        return (f"Reference(name = \"{self.name}\", " +
                f"sequence = \"{self.sequence.seq}\")")
