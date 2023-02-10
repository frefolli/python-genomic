from typing import Generator
import pysam

class BamFile:
    def __init__(self, path : str):
        self.path = path
        pysam.index(path)
        self.alignment_file = pysam.AlignmentFile(path, "rb")
        self.alignment_header = self.alignment_file.header
    
    def get_alignments(self) -> Generator[pysam.AlignedSegment, None, None]:
        for read in self.alignment_file.fetch():
            yield read
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, tb):
        self.alignment_file.close()
    
    def __str__(self) -> str:
        return f"BamFile(path = '{self.path}')"