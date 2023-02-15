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