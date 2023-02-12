from Bio import SeqRecord

class Reference:
    def __init__(self, name : str, sequence : SeqRecord.SeqRecord):
        self.name = name
        self.sequence = sequence

    def get_slice(self, start_index : int = None, end_index : int = None) -> SeqRecord.SeqRecord:
        if (start_index == None):
            return self.get_slice(self, 0, end_index)
        if (end_index == None):
            return self.sequence[start_index:]
        return self.sequence[start_index: end_index]
    
    def __str__(self) -> str:
        return f"Reference(name = \"{self.name}\", sequence = {self.sequence.seq})"