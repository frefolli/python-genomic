"""
    @does define CigarOperation abstraction
"""

from enum import Enum


class CigarOperation(Enum):
    """
        @does define CigarOperation abstraction
    """
    M = 0
    I = 1 # noqa
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


def cigar_tuple_to_cigar_string(cigar_tuple: tuple[int, int]):
    """
        @does convert cigar tuple to cigar string
    """
    cigar_string = str(cigar_tuple[1])
    match CigarOperation(cigar_tuple[0]):
        case CigarOperation.M:
            cigar_string = "M" + cigar_string
        case CigarOperation.I:
            cigar_string = "I" + cigar_string
        case CigarOperation.D:
            cigar_string = "D" + cigar_string
        case CigarOperation.N:
            cigar_string = "N" + cigar_string
        case CigarOperation.S:
            cigar_string = "S" + cigar_string
        case CigarOperation.H:
            cigar_string = "H" + cigar_string
        case CigarOperation.P:
            cigar_string = "P" + cigar_string
        case CigarOperation.E:
            cigar_string = "E" + cigar_string
        case CigarOperation.X:
            cigar_string = "X" + cigar_string
    return cigar_string


def cigar_tuples_to_cigar_string(cigar_tuples: list[tuple[int, int]]):
    """
        @does convert cigar tuples to cigar string
    """
    return "".join([cigar_tuple_to_cigar_string(_) for _ in cigar_tuples])
