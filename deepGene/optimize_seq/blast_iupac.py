"""
finds possible motifs
"""

__version__ = '1.1.0'
__author__ = 'Dhoha Abid'

import reverse_translate as rt
import encode_iupac_to_nt as iupac
import helper_functions as helper
import re
import math
from Bio import SeqIO
import sys
# import traceback
from functools import reduce
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

D_IUPAC_NT = iupac.d_dg_nt
D_AA_CODON = rt.d_aa_nt


class ExtensionOutOfRangeException(Exception):
    def __init__(self, msg='broken'):
        self.msg = msg


class BlastIUPAC:
    """

    """
    def __init__(self, seq_nt, seq_aa):
        """
        initialize the input variables and define output variables
        :param seq_nt: nucleotide sequence
        :param seq_aa: amino acid sequence
        """
        # TODO should be string or Seq (to decide about the input of the constructor)l
        self.seq_nt = seq_nt
        self.seq_aa = seq_aa
        self.len_aa = len(seq_aa) # TODO I need to remove the star from the sequence
        self.len_nt = len(seq_nt)
        self.seq_motif = None
        self.res_motifs_exist = {}  # will contain the motifs that already exists in the nt sequence
        self.res_motif_optimize = {}  # ill contain the motifs that can optimize (match) the nt sequence

    def run_batch(self, l_motif):
        for m in l_motif:
            self.run(seq_motif=m)

    def run(self, seq_motif):
        self.seq_motif = seq_motif
        d_motif_idx__seq_idx = helper.attach_motif_to_seq(seq_motif, self.seq_aa)
        for motif_idx, l_loc_aa in d_motif_idx__seq_idx.items():
            for loc_aa in l_loc_aa:  # extend of enh anchor, we loop over these gn location of the previous dict
                try:
                    right_ext = self.extend_right(seq_motif, motif_idx, loc_aa)
                    if right_ext != 'none':
                        left_ext = self.extend_left(seq_motif, motif_idx, loc_aa - 1)
                        if left_ext != 'none':
                            seq_opt = left_ext + right_ext
                            self.update_res(seq_opt, seq_motif, motif_idx, loc_aa)
                except KeyError, e:
                    # print 'I got a KeyError - reason "%s"' % str(e)
                    continue
                except IndexError, e:
                    # print 'I got an IndexError - reason "%s"' % str(e)
                    continue
                # except:
                #     pass
                # traceback.print_exc()

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
            extension = '' # this one is the intersection
            l, aa_upper_bound = self.get_trimers_right(motif_idx, aa_idx)
            for motif_trimer, nt_trimer, aa in l:
                reg_exp_trimer = helper.align_aa_trimer(aa, motif_trimer)
                if reg_exp_trimer != 'none':
                    extension += reg_exp_trimer
                else:
                    return 'none'
            nt_trimers = reduce(lambda x, y: x + y, [t[1] for t in l])
            return extension, nt_trimers, aa_upper_bound
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
        l_n = [] # trimers of the nt sequence for the right extension
        l_a = []  # aa acids for right extension
        aa_upper_bound = None #inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx, aa_idx + len(l_m))]:
            if i_nt + 3 > self.len_nt:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_n.append(self.seq_nt[i_nt:i_nt + 3])
                l_a.append(self.seq_aa[i_aa])
                aa_upper_bound = i_aa

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
            l, aa_lower_bound = self.get_trimers_left(motif_idx, aa_idx)
            for motif_timer, nt_trimer, aa in l:
                reg_exp_trimer = helper.align_aa_trimer(aa, motif_timer)
                if reg_exp_trimer != 'none':
                    extension = reg_exp_trimer + extension
                else:
                    return 'none'
            nt_trimers = reduce(lambda x, y: y + x, [t[1] for t in l])
            return extension, nt_trimers, aa_lower_bound
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
               for y in range (motif_idx, 0, -3)]
        l_n = []
        l_a = []
        aa_lower_bound = None #inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx - 1, aa_idx - len(l_m) -1, -1)]:
            if i_nt < 0:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_n.append(self.seq_nt[i_nt: i_nt + 3])
                l_a.append(self.seq_aa[i_aa])
                aa_lower_bound = i_aa
        return [(m, n, a) for m, n, a in zip(l_m, l_n, l_a)], aa_lower_bound

    def update_res(self, seq_opt, seq_motif, idx_motif, idx_aa):
        # get the indices of the aligned motif in the aa and nt seq
        left_idx_nt = idx_aa * 3 - len(seq_motif[:idx_motif]) - (0 if len(seq_motif[:idx_motif]) % 3 == 0 else (3 - (len(seq_motif[:idx_motif]) % 3)))
        right_idx_nt = idx_aa * 3 + len(seq_motif[idx_motif:]) + (0 if len(seq_motif[idx_motif:]) % 3 == 0 else (3 - (len(seq_motif[idx_motif:]) % 3)))
        left_idx_aa = idx_aa - math.ceil(len(seq_motif[:idx_motif]) / 3)
        right_idx_aa = idx_aa + math.ceil(len(seq_motif[idx_motif:]) / 3)

        loc_aa = (left_idx_aa, right_idx_aa)
        loc_nt = (left_idx_nt, right_idx_nt)
        seq_wild = str(self.seq_nt[left_idx_nt:right_idx_nt])
        res = (seq_opt, seq_wild, loc_aa, loc_nt)

        if re.match(seq_opt, seq_wild):
            if seq_motif in self.res_motifs_exist.keys():
                if not res in self.res_motifs_exist[seq_motif]:
                    self.res_motifs_exist[seq_motif].append(res)
            else:
                self.res_motifs_exist[seq_motif] = [res]
        else:
            if seq_motif in self.res_motif_optimize.keys():
                if not res in self.res_motif_optimize[seq_motif]:
                    self.res_motif_optimize[seq_motif].append(res)
            else:
                self.res_motif_optimize[seq_motif] = [res]

    @staticmethod
    def reg_exp_helper(string):
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
            if s != '[' and s != ']' and not flag:
                l.append(s)
        return l

    # def print_res(self, ):
    #     print self.res
    def process_res(self):
        seq_empty = ['N' for i in range(len(self.seq_nt))]
        for k, v in self.res_motif_optimize.items():
            for pm in v:
                idx_start = pm[3][0]
                idx_end = pm[3][1]
                if idx_start > 7044:
                    continue
                if len(pm[0]) != len(pm[1]):
                    continue
                if len(set(seq_empty[idx_start:idx_end])) == 1:
                    for i, s in enumerate(self.reg_exp_helper(pm[0]), idx_start):
                        try:
                            if type(s).__name__ == 'str':
                                seq_empty[i] = s
                            else:
                                seq_empty[i] = s[0]
                        except IndexError:
                            print k
                            print pm
                            sys.exit(0)
        j = 0
        l = [None] * len(self.seq_nt)
        for s1, s2 in zip(seq_empty, self.seq_nt):
            if s1 == 'N':
                l[j] = s2
            else:
                l[j] = s1
            j += 1

        return reduce(lambda x, y: x + y, l)

    def check_res(self, string):
        seq_res = Seq(string)
        translated = seq_res.translate()
        if translated != self.seq_aa:
            print 'we got serious pb'
        with open('res/res.liver.motif.v2.fasta', 'w') as f:
            SeqIO.write(SeqRecord(seq_res, id='liver.motif.v1'), f, 'fasta')

    def check_how_much_different(self, string):
        if string != self.seq_nt:
            print 'it ''s diffrent'
        j = 0
        for s1, s2 in zip(string, self.seq_nt):
            if s1 != s2:
                print j
                print s1
                print s2
            j += 1


if __name__ == '__main__':
    # motif = 'YCYMMAYW'
    with open('data/f8.orig.fasta', 'r') as f_seq\
            , open('data/liver_motifs.fasta', 'r') as f_liver_motifs\
            , open('data/enhancer_all.fasta', 'r') as f_motifs:
        nt = SeqIO.read(f_seq, 'fasta').seq
        aa = str(nt.translate())
        ll_motif = [str(record_motif.seq) for record_motif in SeqIO.parse(f_motifs, 'fasta')]
        l_liver_motif = [str(record_motif.seq) for record_motif in SeqIO.parse(f_liver_motifs, 'fasta')]
        # seq_motif = 'NRCGTGNNN'
        sseq_motif = 'NGDBCA'

    b = BlastIUPAC(seq_nt=nt, seq_aa=aa)
    b.run_batch(l_liver_motif)
    # print b.possible_motif
    # b.run(sseq_motif)
    res = b.process_res()
    print res
    b.check_res(res)
    b.check_how_much_different(res)


    # b.run(seq_motif)
    # print b.possible_motif

