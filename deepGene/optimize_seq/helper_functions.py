"""
This module include helper functions for blast_iupac and optimize_seq
"""

__version__ = '1.0.0'
__author__ = 'Dhoha Abid'


from Bio.Seq import Seq


def index_aa_seq(seq_aa):
    """
    index an amino acids(aa) sequence. indices start from zero.
    :param seq_aa: aa sequence
    :return: a dictionary where the keys are the aa and values are indices of sequence.
    """
    d_aa_idx = {}
    for i, aa in enumerate(seq_aa):
        if aa != '*':
            if aa not in d_aa_idx.keys():
                d_aa_idx[aa] = {i}
            else:
                d_aa_idx[aa].add(i)
    return d_aa_idx


def transform_reg_exp_to_list(string):
    """
    given a string that follow a regular expression annotation, this method transforms it to a list of character.
    :param string: string of regular expression 'AT[GT]A'
    :return: list of characters ['A', 'T', ['G', 'T'], 'A']
    """
    l = []
    flag = False
    for s in string:
        if s == '[':
            flag = True
            l_sub = []
        elif s == ']':
            flag = False
            l.append(l_sub)
        elif flag:
            l_sub.append(s)
        elif s != '[' and s != ']' and not flag:
            l.append(s)
    return l


def translate_nt_sequence(seq_nt):
        seq_nt = Seq(seq_nt)
        return str(seq_nt.translate())


