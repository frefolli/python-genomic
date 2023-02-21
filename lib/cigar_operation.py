from enum import Enum

class CigarOperation(Enum):
    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    E = 7
    X = 8
    
CIGAR_OPERATIONS_THAT_CONSUME_QUERY = [
    CigarOperation.M,
    CigarOperation.I,
    CigarOperation.S,
    CigarOperation.E,
    CigarOperation.X
]
    
CIGAR_OPERATIONS_THAT_CONSUME_REFERENCE = [
    CigarOperation.M,
    CigarOperation.D,
    CigarOperation.N,
    CigarOperation.E,
    CigarOperation.X
]

def cigar_tuple_to_cigar_string(cigar_tuple : tuple[int, int]):
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.M): return f"M{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.I): return f"I{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.D): return f"D{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.N): return f"N{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.S): return f"S{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.H): return f"H{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.P): return f"P{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.E): return f"E{cigar_tuple[1]}"
    if (CigarOperation(cigar_tuple[0]) == CigarOperation.X): return f"X{cigar_tuple[1]}"

def cigar_tuples_to_cigar_string(cigar_tuples : list[tuple[int, int]]):
    return "".join([cigar_tuple_to_cigar_string(_) for _ in cigar_tuples])