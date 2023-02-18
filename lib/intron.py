class Intron:
    def __init__(self, sequence : str) -> None:
        self.sequence = sequence

    def is_canonic(self) -> bool:
        if len(self.sequence) >= 4:
            return self.sequence[:2].lower() == 'gt' and self.sequence[-2:].lower() == 'ag'
        return False
    
    def __str__(self) -> str:
        return f"Intron(sequence = \"{self.sequence}\")"