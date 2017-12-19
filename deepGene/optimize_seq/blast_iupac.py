"""
finds possible motifs
algo similat to blast
finds the spots that can be converted to motifs

"""

__version__ = '1.1.0'
__author__ = 'Dhoha Abid'

import reverse_translate as rt
import encode_iupac_to_nt as iupac
import helper_functions as helper
import re
from Bio import SeqIO
import sys
import traceback
from functools import reduce
from Bio.Seq import Seq


D_IUPAC_NT = iupac.d_dg_nt
D_AA_CODON = rt.d_aa_nt


class ExtensionOutOfRangeException(Exception):
    def __init__(self, msg):
        self.msg = msg


class BlastIUPAC:
    """
    an implementation of blast to support the iupac alignments
    """
    def __init__(self, seq_nt):
        """
        initialize the input variables and define output variables
        :param seq_nt: nucleotide sequence (of type string)
        :param seq_aa: amino acid sequence (of type string)
        """
        self.seq_nt = seq_nt
        self.seq_aa = self.translate_seq_nt(seq_nt)
        self.len_aa = len(self.seq_aa)
        self.len_nt = len(seq_nt)
        self.seq_motif = None
        self.res_motifs_exist = {}  # will contain the motifs that already exists in the nt sequence
        self.res_motif_optimize = {}  # ill contain the motifs that can optimize (match) the nt sequence

    @staticmethod
    def translate_seq_nt(seq_nt):
        seq_nt = Seq(seq_nt)
        return str(seq_nt.translate())

    def run_batch(self, l_motif):
        for m in l_motif:
            self.run(seq_motif=m)

    def run(self, seq_motif):
        """
        after aligning all the trimer of the motif sequence to the aa sequence using the helper function
        attach_motif_to_aa_seq, we loop over all these seeds to extend them using extend_right and extend_left.
        The results are the set of motifs that can be aligned to the sequence
        :param seq_motif: motif sequence, iupac is supported
        :return:
        """
        self.seq_motif = seq_motif
        d_motif_idx__seq_aa_idx = helper.attach_motif_to_aa_seq(seq_motif, self.seq_aa)
        for motif_idx, l_aa_idx in d_motif_idx__seq_aa_idx.items():
            for aa_idx in l_aa_idx:
                try:
                    right_extension, right_nt_align, aa_upper_bound = self.extend_right(motif_idx, aa_idx)
                    left_extension, left_nt_align, aa_lower_bound = self.extend_left(motif_idx, aa_idx)
                    reg_exp_seq_optimized = left_extension + right_extension
                    wild_align = left_nt_align + right_nt_align
                    self.update_d_results(reg_exp_seq_optimized, wild_align, aa_lower_bound, aa_upper_bound)
                except ValueError as e:
                    if str(e) == 'too many values to unpack': # if extend_righ and extend_left return 'none', we won't be able to unpack the values and thus go to the next value in the loop
                        continue
                    else:
                        print e
                except:
                    traceback.print_exc()
                    sys.exit(1)

    def extend_right(self, motif_idx, aa_idx):
        """
        extend the alignment to the right. get all the trimers of the motif, and corresponding aa in the sequence
        using get_trimers_right. Then align them using the helper align_aa_trimer. Concatenate the resultant alignment,
        then return it. Regular expressions are used to return the right extension
        :param motif_idx: start index of the motif sequence
        :param aa_idx: the aa index in the sequence
        :return: regular expression of the right extension
        """
        try:
            extension = ''  # this one is the intersection
            wild = ''
            l, aa_upper_bound = self.get_trimers_right(motif_idx, aa_idx)
            for motif_trimer, wild_trimer, aa in l:
                reg_exp_trimer = helper.align_aa_trimer(aa, motif_trimer)
                if reg_exp_trimer != 'none':
                    extension += reg_exp_trimer
                    wild += wild_trimer
                else:
                    return 'none'
            return extension, wild, aa_upper_bound
        except ExtensionOutOfRangeException, e:
            return 'none'

    def get_trimers_right(self, motif_idx, aa_idx):
        """
        generate a list of trimers for right extension. The list would have tuple of motif trimers, nt trimer, and aa.
        The items in the tuple can be aligned in the method extend_right to see if we can match the motifs to the aa
        codons. In case the extension out range the nt sequence, we throw an exception ExtensionOutOfRangeException
        :param motif_idx: start index in the motif sequence (the nail)
        :param aa_idx: start index in the aa sequence
        :return: list of tuples (motif_trimer, nt_trimer, aa)
        """
        l_m = [self.seq_motif[i:i + 3].ljust(3, 'N')
               for i in range(motif_idx, len(self.seq_motif), 3)] # trimers of motif for right extension
        l_n = []  # trimers of the nt sequence for the right extension
        l_a = []  # aa acids for right extension
        aa_upper_bound = None  # inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx, aa_idx + len(l_m))]:
            if i_nt + 3 > self.len_nt:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_n.append(self.seq_nt[i_nt:i_nt + 3])
                l_a.append(self.seq_aa[i_aa])
                aa_upper_bound = i_aa
        if not aa_upper_bound:
            raise ExtensionOutOfRangeException(l_m)
        return [(m, n, a) for m, n, a in zip(l_m, l_n, l_a)], aa_upper_bound

    def extend_left(self, motif_idx, aa_idx):
        """
        ame description as extend_right, just this method is doing it for the left extension
        :param motif_idx: motif_idx: start index of the motif sequence
        :param aa_idx: aa_idx: the aa index in the sequence
        :return: regular expression of the left extension
        """
        try:
            extension = ''
            wild = ''
            l, aa_lower_bound = self.get_trimers_left(motif_idx, aa_idx)
            for motif_timer, wild_trimer, aa in l:
                reg_exp_trimer = helper.align_aa_trimer(aa, motif_timer)
                if reg_exp_trimer != 'none':
                    extension = reg_exp_trimer + extension
                    wild = wild_trimer + wild
                else:
                    return 'none'
            return extension, wild, aa_lower_bound
        except ExtensionOutOfRangeException, e:
            return 'none'

    def get_trimers_left(self, motif_idx, aa_idx):
        """
        same description as for get_trimers_right, but this method is doing it for the left extension
        :param motif_idx: start index in the motif sequence (the nail)
        :param aa_idx: start index in the aa sequence
        :return: list of tuples (motif_trimer, nt_trimer, aa)
        """
        l_m = [self.seq_motif[y - 3:y].rjust(3, 'N') if y - 3 >= 0 else self.seq_motif[0:y].rjust(3, 'N')
               for y in range(motif_idx, 0, - 3)]
        l_n = []
        l_a = []
        aa_lower_bound = None  # inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx - 1, aa_idx - len(l_m) - 1, -1)]:
            if i_nt < 0:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_n.append(self.seq_nt[i_nt: i_nt + 3])
                l_a.append(self.seq_aa[i_aa])
                aa_lower_bound = i_aa
        if not aa_lower_bound:
            raise ExtensionOutOfRangeException(l_m)
        return [(m, n, a) for m, n, a in zip(l_m, l_n, l_a)], aa_lower_bound

    def update_d_results(self, reg_exp_seq_optimized, wild_align, aa_lower_bound, aa_upper_bound):
        """
        aa
        :param reg_exp_seq_optimized:
        :param wild_align:
        :param aa_lower_bound: inclusive
        :param aa_upper_bound: inclusive
        :return:
        """
        res = (reg_exp_seq_optimized, wild_align, (aa_lower_bound, aa_upper_bound)
               , (aa_lower_bound * 3, aa_upper_bound * 3 + 3))
        if re.match(reg_exp_seq_optimized, wild_align):
            if self.seq_motif in self.res_motifs_exist.keys():
                if not res in self.res_motifs_exist[self.seq_motif]:
                    self.res_motifs_exist[self.seq_motif].append(res)
            else:
                self.res_motifs_exist[self.seq_motif] = [res]
        else:
            if self.seq_motif in self.res_motif_optimize.keys():
                if not res in self.res_motif_optimize[self.seq_motif]:
                    self.res_motif_optimize[self.seq_motif].append(res)
            else:
                self.res_motif_optimize[self.seq_motif] = [res]


def run_with_one_motif():
    seq_motif = 'WNNVYTAATTARYYNNN'
    with open('data/f8.orig.fasta', 'r') as f_seq_nt:
        seq_nt = str(SeqIO.read(f_seq_nt, 'fasta').seq)
    b = BlastIUPAC(seq_nt)
    b.run(seq_motif)
    print b.res_motif_optimize


def run_with_multiple_motif():
    with open('data/f8.orig.fasta', 'r') as f_seq_nt\
            , open('data/liver_motifs.fasta') as f_motifs:
        seq_nt = str(SeqIO.read(f_seq_nt, 'fasta').seq)
        l_motif = [str(record.seq) for record in SeqIO.parse(f_motifs, 'fasta')]
    b = BlastIUPAC(seq_nt)
    b.run_batch(l_motif)
    print b.res_motif_optimize

    pass


if __name__ == '__main__':
    # run_with_one_motif()
    run_with_multiple_motif()