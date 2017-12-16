import reverse_translate as rt
import encode_iupac_to_nt as iupac
import re
import math
from Bio import SeqIO


class BlastIUPAC:
    """

    """
    d_dg_nt = iupac.d_dg_nt
    d_aa_nt = rt.d_aa_nt

    def __init__(self, seq_nt, seq_aa):
        self.seq_nt = seq_nt
        self.seq_aa = seq_aa
        self.d_idxenh_idxgn = {}
        self.d_aa_idx = {}
        for i, aa in enumerate(seq_aa):
            if aa != '*':
                if aa not in self.d_aa_idx.keys():
                    self.d_aa_idx[aa] = {i}
                else:
                    self.d_aa_idx[aa].add(i)
        self.res = {}

    def run_batch(self, l_motif):
        for seq_motif in l_motif:
            self.run(seq_motif)

    def run(self, seq_motif):
        self.anchor_motif(seq_motif)
        for idx_motif, l_loc_aa in self.d_idxenh_idxgn.items():
            for loc_aa in l_loc_aa:  # extend of enh anchor, we loop over these gn location of the previous dict
                try:
                    right_ext = self.extend_right(seq_motif, idx_motif, loc_aa)
                    if right_ext != 'none':
                        left_ext = self.extend_left(seq_motif, idx_motif, loc_aa - 1)
                        if left_ext != 'none':
                            seq_opt = left_ext + right_ext
                            self.update_res(seq_opt)
                except:
                    # TODO output the error
                    continue

    def anchor_motif(self, seq_motif):
        for idxenh, enh3nt in enumerate([seq_motif[i:i + 3] for i in range(len(seq_motif) - 2)]):
            for aa, idx in self.d_aa_idx.items():
                for nt in self.d_aa_nt[aa]:
                    if (nt[0] & self.d_dg_nt[enh3nt[0]]) and (nt[1] & self.d_dg_nt[enh3nt[1]]) \
                            and (nt[2] & self.d_dg_nt[enh3nt[2]]):
                        if idxenh not in self.d_idxenh_idxgn.keys():
                            self.d_idxenh_idxgn[idxenh] = idx
                        else:
                            self.d_idxenh_idxgn[idxenh] |= idx

    def extend_right(self, seq_motif, idx_motif, loc_aa):
        right_ext = ''
        for aa, e in self.get_trimers_right(idx_motif=idx_motif, loc_aa=loc_aa, seq_motif=seq_motif):
            reg_trimer = self.match_aa_trimer(aa, e)
            if reg_trimer != 'none':
                right_ext += reg_trimer
            else:
                return 'none'
        return right_ext

    def extend_left(self, seq_motif, idx_motif, loc_aa):
        left_ext = ''
        for aa, e in self.get_trimers_left(seq_motif, idx_motif, loc_aa):
            reg_trimer= self.match_aa_trimer(aa, e)
            if reg_trimer != 'none':
                left_ext = reg_trimer + left_ext
            else:
                return 'none'
        return left_ext

    @staticmethod
    def match_aa_trimer(aa, trimer):
        for nt in BlastIUPAC.d_aa_nt[aa]:
            nt1 = nt[0] & BlastIUPAC.d_dg_nt[trimer[0]]
            nt2 = nt[1] & BlastIUPAC.d_dg_nt[trimer[1]]
            nt3 = nt[2] & BlastIUPAC.d_dg_nt[trimer[2]]
            if nt1 and nt2 and nt3:
                reg_trimer = reduce(lambda x, y: x + y, [reduce(lambda x, y: x + y, i) if len(i) == 1
                                                         else '[' + reduce(lambda x, y: x + y, i) + ']'
                                                         for i in [nt1, nt2, nt3]])
                return reg_trimer
        return 'none'

    def get_trimers_right(self, seq_motif, idx_motif, loc_aa):
        return[(self.seq_aa[x], seq_motif[y:y + 3].ljust(3, 'N'))
               for x, y in enumerate([i for i in range(idx_motif, len(seq_motif), 3)], loc_aa)]

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
        seq_wild = self.seq_nt[left_idx_nt:right_idx_nt]
        res = (seq_opt, seq_motif, loc_aa, loc_nt)

        if not re.match(seq_opt, seq_wild):
            if seq_motif in self.res.keys():
                if not res in self.res[seq_motif]:
                    self.res[seq_motif].append(res)
            else:
                self.res[seq_motif] = [res]

    def print_res(self):
        print self.res


if __name__ == '__main__':
    motif = 'YCYMMAYW'
    with open('test/f8.orig.fasta', 'r') as seq_f:
        nt = SeqIO.read(seq_f, 'fasta').seq
        aa = str(nt.translate())

    v = BlastIUPAC(seq_nt=nt, seq_aa=aa)
    v.run(motif)
