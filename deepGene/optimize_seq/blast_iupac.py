"""
finds possible motifs
"""

__version__ = '1.0.0'
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


class BlastIUPAC:
    """

    """
    def __init__(self, seq_nt, seq_aa):
        """
        initialize the input variables and define output variables
        :param seq_nt: nucleotide sequence
        :param seq_aa: amino acid sequence
        """
        self.seq_nt = seq_nt
        self.seq_aa = seq_aa
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
        intersection of the motif and aa. right extension
        :param motif_idx: motif start index
        :param aa_idx: aa index in the sequence
        :return:
        """
        reg_exp_extension = ''
        for aa, motif_trimer in self.get_trimers_right(motif_idx, aa_idx):
            reg_exp_trimer = self.align_aa_trimer(aa, motif_trimer)
            if reg_exp_trimer != 'none':
                reg_exp_extension += reg_exp_trimer
            else:
                return 'none'
        return reg_exp_extension

    def get_trimers_right(self, motif_idx, aa_idx):
        """
        generate the list of the trimers for right extension. The extension can results in a longer right extension
        compared to the motif as we extend by three to preserve the aa sequence
        :param motif_idx: motif start index
        :param aa_idx: aa start index
        :return: list of tuples (aa, motif_trimer)
        """
        try:
            l_trimer = [(self.seq_aa[x], self.seq_motif[y:y + 3].ljust(3, 'N'))
               for x, y in enumerate([i for i in range(motif_idx, len(self.seq_motif), 3)], aa_idx)]
        except KeyError, e:
            # TODO I want to catch the key and if = to *, pass, else exit and show the stack
            pass
            # TODO if right extension depasse the length of the sequence o handle and control
        except IndexError, e:
            pass
        return l_trimer


    def extend_left(self, seq_motif, idx_motif, loc_aa):
        left_ext = ''
        for aa, e in self.get_trimers_left(seq_motif, idx_motif, loc_aa):
            reg_trimer= self.align_aa_trimer(aa, e)
            if reg_trimer != 'none':
                left_ext = reg_trimer + left_ext
            else:
                return 'none'
        return left_ext

    @staticmethod
    def align_aa_trimer(aa, trimer):
        for nt in D_AA_CODON[aa]:
            nt1 = nt[0] & D_IUPAC_NT[trimer[0]]
            nt2 = nt[1] & D_IUPAC_NT[trimer[1]]
            nt3 = nt[2] & D_IUPAC_NT[trimer[2]]
            if nt1 and nt2 and nt3:
                reg_trimer = reduce(lambda x, y: x + y, [reduce(lambda x, y: x + y, i) if len(i) == 1
                                                         else '[' + reduce(lambda x, y: x + y, i) + ']'
                                                         for i in [nt1, nt2, nt3]])
                return reg_trimer
        return 'none'



    def get_trimers_left(self, seq_motif, idx_motif, loc_aa):
        l_trimer = [seq_motif[y - 3:y].rjust(3, 'N') if y - 3 >= 0 else seq_motif[0:y].rjust(3, 'N')
                for xx, y in enumerate([i for i in range(idx_motif, 0, -3)], loc_aa)]
        l_idx_aa = [self.seq_aa[i] for i in range(loc_aa, loc_aa - int(math.ceil(idx_motif/float(3))), -1)]
        return zip(l_idx_aa, l_trimer)

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

