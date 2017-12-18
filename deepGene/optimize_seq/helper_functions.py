"""
This file include the statistics methods
for the
"""

__version__ = '1.0.0'
__author__ = 'Dhoha Abid'


import reverse_translate as rt
import encode_iupac_to_nt as iupac

D_IUPAC_NT = iupac.d_dg_nt
D_AA_CODON = rt.d_aa_nt


def index_aa_seq(seq_aa):
    """
    index the amino acids(aa) of the input sequence. index starts from zero
    :param seq_aa: aa sequence
    :return: a dictionary where the keys are the aa and values are index in sequence.
    """
    d_aa_idx = {}
    for i, aa in enumerate(seq_aa):
        if aa != '*':
            if aa not in d_aa_idx.keys():
                d_aa_idx[aa] = {i}
            else:
                d_aa_idx[aa].add(i)
    return d_aa_idx


def attach_motif_to_seq(seq_motif, seq_aa):
    """
    attach the motif to the sequence. Attachment is based on the codons (aa) of the sequence and a trimer
    of the motif. Attach an index of the motif to the index of the sequence and save the results in a dictionary.
    :param seq_motif: motif sequence
    :param d_aa_idx: dictionary where keys are aa and values are the index of aa in the sequence
    :return: a dictionary where the keys are motif index and values are the indices of the sequence
    """
    d_aa_idx = index_aa_seq(seq_aa)
    d = {}
    for motif_idx, motif_trimer in enumerate([seq_motif[i:i + 3] for i in range(len(seq_motif) - 2)]):
        for aa, seq_idx in d_aa_idx.items():
            for codon in D_AA_CODON[aa]:
                if (codon[0] & D_IUPAC_NT[motif_trimer[0]]) and (codon[1] & D_IUPAC_NT[motif_trimer[1]]) \
                        and (codon[2] & D_IUPAC_NT[motif_trimer[2]]):
                    if motif_idx not in d.keys():
                        d[motif_idx] = seq_idx
                    else:
                        d[motif_idx] |= seq_idx
    return d
