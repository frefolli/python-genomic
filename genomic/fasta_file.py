"""
    @does define FastaFile abstraction
"""

import logging

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
        self.__path = path
        self.__index = SeqIO.index(path, "fasta")

    def get_reference_names(self) -> list[str]:
        """
            @returns list of reference names
        """
        return list(self.__index)

    def get_reference(self, name: str) -> Reference:
        """
            @returns reference from name
        """
        return Reference(name, self.__index[name])

    def __enter__(self):
        """
            @does enter for with-as
        """
        logging.info("OPENED fastafile %s", self.__path)
        return self

    def __exit__(self, exception_type,
                 exception_value, exception_traceback):
        """
            @does exit for with-as
        """
        self.__index.close()
        logging.info("CLOSED fastafile %s", self.__path)

    def __str__(self) -> str:
        """
            @does
        """
        return f"FastaFile(path = \"{self.__path}\")"
