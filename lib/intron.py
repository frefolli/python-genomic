"""
    @does fine Intron abstraction
"""


class Intron:
    """
        @does fine Intron abstraction
    """

    def __init__(self, sequence: str) -> None:
        """
            @does initialize Intron
        """
        self.__sequence = sequence

    def is_canonical(self) -> bool:
        """
            @returns True/False if Intron is Canonical or not
        """

        if len(self.__sequence) >= 4:
            return (self.__sequence[:2].lower() == 'gt' and
                    self.__sequence[-2:].lower() == 'ag')
        return False

    def __str__(self) -> str:
        """
            @returns string representation of Intron
        """
        return f"Intron(sequence = \"{self.__sequence}\")"
