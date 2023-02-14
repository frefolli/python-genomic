from Bio import SeqRecord

class Reference:
    def __init__(self, name : str, sequence : SeqRecord.SeqRecord):
        self.name = name
        self.sequence = sequence

    def get_slice(self, start_index : int = None, end_index : int = None, reverse : bool = False) -> SeqRecord.SeqRecord:
        seq = self.sequence.seq
        if (start_index == None): start_index=0
        if (end_index == None):
            seq = seq[start_index:]
        else:
            seq = seq[start_index: end_index]
        return seq.reverse_complement() if reverse else seq
    
    def __str__(self) -> str:
        return f"Reference(name = \"{self.name}\", sequence = \"{self.sequence.seq}\")"