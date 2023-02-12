from Bio import SeqIO
from .reference import Reference

class FastaFile:
    def __init__(self, path : str):
        self.path = path
        self.index = SeqIO.index(path, "fasta")
    
    def get_reference_names(self) -> list[str]:
        return [_ for _ in self.index]
    
    def get_reference(self, name : str) -> Reference:
        return Reference(name, self.index[name])
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, tb):
        self.index.close()
    
    def __str__(self) -> str:
        return f"FastaFile(path = \"{self.path}\")"