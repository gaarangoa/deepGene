"""
d_iupac_nt is a dictionary where the keys are upac notation and values possible standard nucleotides (A, C, G, T)
"""

__version__ = '1.1.0'
__author__ = 'Dhoha Abid'

from json import dump


d_iupac_nt = {
    'A':{'A'}
    , 'C':{'C'}
    , 'G':{'G'}
    , 'T':{'T'}
    , 'R':{'A', 'G'}
    , 'Y':{'C', 'T'}
    , 'M':{'A', 'C'}
    , 'K':{'G', 'T'}
    , 'S':{'G', 'C'}
    , 'W':{'A', 'T'}
    , 'H':{'A', 'C', 'T'}
    , 'B':{'C', 'G', 'T'}
    , 'V':{'A', 'C', 'G'}
    , 'D':{'A', 'G', 'T'}
    , 'N':{'A', 'C', 'G', 'T'}   
}


def dump_dict(f_name):
    with open(f_name, 'w') as f:
        dump(d_iupac_nt, f)


def get_dict():
    return d_iupac_nt