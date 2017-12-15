import reverse_translate as rt
import encode_iupac_to_nt as iupac
import re
import math
from Bio import SeqIO

class BlastIUPAC:
    """

    """
    def __init__(self, seq_nt, seq_aa):
        self.seq_nt = seq_nt
        self.seq_aa = seq_aa
        self.d_aa_nt = rt.d_aa_nt
        self.d_dg_nt = iupac.d_dg_nt
        self.d_idxenh_idxgn = {}
        self.d_aa_idx = {}
        for i, aa in enumerate(seq_aa):
            if aa != '*':
                if aa not in self.d_aa_idx.keys():
                    self.d_aa_idx[aa] = {i}
                else:
                    self.d_aa_idx[aa].add(i)

    def run(self, seq_motif):
        self.anchor_motif(seq_motif)
        self.extend_left_right(seq_motif)

    def run_batch(self, l_motif):
        for seq_motif in l_motif:
            self.run(seq_motif)

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

    def extend_left_right(self, seq_motif):
        for idxenh, s_idxgn in self.d_idxenh_idxgn.items():
            for idxgn in s_idxgn:  # extend of enh anchor, we loop over these gn location of the previous dict
                # Right extension
                rext = ''
                enhrpad = ''
                try:
                    for aa, e in [(self.seq_aa[x], seq_motif[y:y + 3].ljust(3, 'N')) for x, y in
                                  enumerate([i for i in range(idxenh, len(seq_motif), 3)], idxgn)]:
                        rreg = 'none'
                        for nt in self.d_aa_nt[aa]:
                            nt1 = nt[0] & self.d_dg_nt[e[0]]
                            nt2 = nt[1] & self.d_dg_nt[e[1]]
                            nt3 = nt[2] & self.d_dg_nt[e[2]]
                            if nt1 and nt2 and nt3:
                                rreg = reduce(lambda x, y: x + y, [
                                    reduce(lambda x, y: x + y, i) if len(i) == 1 else '[' + reduce(lambda x, y: x + y,
                                                                                                   i) + ']' for i in
                                    [nt1, nt2, nt3]])
                        if rreg != 'none':
                            rext = rext + rreg
                            enhrpad = enhrpad + e
                        else:
                            break
                    # Left extension: if the right extension got the end of the enhancer
                    if rreg != 'none':
                        lext = ''
                        enhlpad = ''
                        d = {}
                        x = idxgn - 1
                        for e in [seq_motif[y - 3:y].rjust(3, 'N') if y - 3 >= 0
                                  else seq_motif[0:y].rjust(3, 'N') for xx, y in
                                  enumerate([i for i in range(idxenh, 0, -3)], idxgn - 1)]:
                            lreg = 'none'
                            aa = self.seq_aa[x]
                            for nt in self.d_aa_nt[aa]:
                                nt1 = nt[0] & self.d_dg_nt[e[0]]
                                nt2 = nt[1] & self.d_dg_nt[e[1]]
                                nt3 = nt[2] & self.d_dg_nt[e[2]]
                                if nt1 and nt2 and nt3:
                                    lreg = reduce(lambda x, y: x + y, [
                                        reduce(lambda x, y: x + y, i) if len(i) == 1 else '[' + reduce(lambda x, y: x + y,
                                                                                                       i) + ']' for i in
                                        [nt1, nt2, nt3]])
                            if lreg != 'none':
                                lext = lreg + lext
                                enhlpad = e + enhlpad
                                x = x - 1
                            else:
                                break
                        if lreg != 'none':
                            # populate the dictionary to save it in the file
                            #                             d['enh.orig'] = enh
                            d['enh.orig'] = reduce(lambda x, y: x + y, [
                                reduce(lambda x, y: x + y, self.d_dg_nt[dg]) if len(self.d_dg_nt[dg]) == 1 else '[' + reduce(
                                    lambda x, y: x + y, self.d_dg_nt[dg]) + ']' for dg in seq_motif])
                            d['res'] = lext + rext
                            enhpad = enhlpad + enhrpad
                            d['enh.pad'] = enhpad
                            lidxgn = idxgn * 3 - len(seq_motif[:idxenh]) - (
                            0 if len(seq_motif[:idxenh]) % 3 == 0 else (3 - (len(seq_motif[:idxenh]) % 3)))
                            ridxgn = idxgn * 3 + len(seq_motif[idxenh:]) + (
                            0 if len(seq_motif[idxenh:]) % 3 == 0 else (3 - (len(seq_motif[idxenh:]) % 3)))
                            d['gn.nt'] = str(self.seq_nt[lidxgn:ridxgn])
                            d['gn.nt.loc'] = (lidxgn, ridxgn)
                            lidxaa = idxgn - math.ceil(len(seq_motif[:idxenh]) / 3)
                            ridxaa = idxgn + math.ceil(len(seq_motif[idxenh:]) / 3)
                            d['gn.aa'] = self.seq_aa[lidxaa:ridxaa]
                            d['gn.aa.loc'] = (lidxaa, ridxaa)
                            if re.match(d['res'], d['gn.nt']):
                                d['enh.exists'] = 1
                                counter = counter + 1
                            else:
                                d['enh.exists'] = 0
                                # here I should store the results
                except:
                    continue

        print(d)

if __name__ == '__main__':
    motif = 'YCYMMAYW'
    with open('test/f8.orig.fasta', 'r') as seq_f:
        nt = SeqIO.read(seq_f, 'fasta').seq
        aa = str(nt.translate())

    v = BlastIUPAC(seq_nt=nt, seq_aa=aa)
    v.run(motif)
