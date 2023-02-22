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
        self.sequence = sequence

    def is_canonical(self) -> bool:
        """
            @returns True/False if Intron is Canonical or not
        """

        if len(self.sequence) >= 4:
            return (self.sequence[:2].lower() == 'gt' and
                    self.sequence[-2:].lower() == 'ag')
        return False

    def __str__(self) -> str:
        """
            @returns string representation of Intron
        """
        return f"Intron(sequence = \"{self.sequence}\")"
