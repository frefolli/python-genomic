"""
    @does define FastaFile abstraction
"""

from Bio import SeqIO
from .reference import Reference


class FastaFile:
    """
        @does define FastaFile abstraction
    """

    def __init__(self, path: str):
        """
            @does initialize FastaFile
        """
        self.path = path
        self.index = SeqIO.index(path, "fasta")

    def get_reference_names(self) -> list[str]:
        """
            @returns list of reference names
        """
        return [_ for _ in self.index]

    def get_reference(self, name: str) -> Reference:
        """
            @returns reference from name
        """
        return Reference(name, self.index[name])

    def __enter__(self):
        """
            @does enter for with-as
        """
        return self

    def __exit__(self, exception_type,
                 exception_value, exception_traceback):
        """
            @does exit for with-as
        """
        self.index.close()

    def __str__(self) -> str:
        """
            @does
        """
        return f"FastaFile(path = \"{self.path}\")"
