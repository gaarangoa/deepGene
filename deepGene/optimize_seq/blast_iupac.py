"""
Implementation of blast that supports the iupac notations of nucleotides. The ultimate results of this class
when fed a list of motifs is to output two dictionaries existent_motifs that contains the motifs that already
exist in the nucleotide sequence, and motifs_to_integrate that has the motifs that have the motifs that can
optimize the nucleotide sequence
"""

__version__ = '1.1.0'
__author__ = 'Dhoha Abid'


import helper_functions as helper
import re
from Bio import SeqIO
import sys
import traceback
from dict_aa_codon import D_AA_CODON
from dict_iupac_nt import D_IUPAC_NT


class ExtensionOutOfRangeException(Exception):
    def __init__(self, msg):
        self.msg = msg


class BlastIUPAC:
    def __init__(self, seq_wild_nt):
        """
        :param seq_wild_nt: nucleotide sequence (of type string)
        """
        self.seq_wild_nt = seq_wild_nt
        self.seq_wild_aa = helper.translate_nt_sequence(seq_wild_nt)
        self.len_seq_wild_aa = len(self.seq_wild_aa)
        self.len_seq_wild_nt = len(seq_wild_nt)
        self.seq_motif = None
        self.existent_motifs = {}  # will contain the motifs that already exists in the nt sequence
        self.motifs_to_integrate = {}  # will contain the motifs that can optimize the nt sequence

    def run_batch(self, l_motif):
        for m in l_motif:
            self.run(seq_motif=m)

    def run(self, seq_motif):
        """
        This method joins all the pieces of this class together to run this version of blast that supports IUPAC.
        It generates the seeds - alignment of a trimer of the motif to the aa index. Then, it extends these seeds
        to the right and left side. Then, if the right and left extension result in an exact match, we update
        the results. The results include two dictionaries the existing motif in the sequence and the motifs that
        can be integrated in the wild sequences with the corresponding start and end locations.
        :param seq_motif: motif sequence
        :return: update self.existent_motifs and motifs_to_integrate
        """
        self.seq_motif = seq_motif
        d_motif_idx__aa_idx = self.attach_motif_to_aa_seq()
        for motif_idx, l_aa_idx in d_motif_idx__aa_idx.items():
            for aa_idx in l_aa_idx:
                try:
                    right_align_extend, right_align_wild, aa_idx_upper_bound = self.extend_right(motif_idx, aa_idx)
                    left_align_extend, left_align_wild, aa_idx_lower_bound = self.extend_left(motif_idx, aa_idx)
                    align_optimized = left_align_extend + right_align_extend
                    align_wild = left_align_wild + right_align_wild
                    self.update_results(align_optimized, align_wild, aa_idx_lower_bound, aa_idx_upper_bound)
                except ValueError as e:
                    if str(e) == 'too many values to unpack':
                        continue
                    else:
                        print e
                except Exception:
                    traceback.print_exc()
                    sys.exit(1)

    def attach_motif_to_aa_seq(self):
        """
        This is the seed implementation of blast: attach the motif to the sequence. Attachment is based on the codons
        (aa) of the sequence and a trimer of the motif. Attach an index of the motif to the index of the sequence
        and save the results in a dictionary.
        :return: a dictionary where the keys are motif index and values are the indices of aa of the sequence
        """
        d_aa_idx = helper.index_aa_seq(self.seq_wild_aa)
        d = {}
        for motif_idx, motif_trimer in enumerate([self.seq_motif[i:i + 3] for i in range(len(self.seq_motif) - 2)]):
            for aa, seq_idx in d_aa_idx.items():
                for codon in D_AA_CODON[aa]:
                    if (codon[0] & D_IUPAC_NT[motif_trimer[0]]) and (codon[1] & D_IUPAC_NT[motif_trimer[1]]) \
                            and (codon[2] & D_IUPAC_NT[motif_trimer[2]]):
                        if motif_idx not in d.keys():
                            d[motif_idx] = seq_idx
                        else:
                            d[motif_idx] |= seq_idx
        return d

    def extend_right(self, motif_idx, aa_idx):
        """
        extend the alignment to the right. get all the trimers of the motif starting from motif_idx, and corresponding
        aa in the nt sequence via get_trimers_right. Then, align them via align_aa_trimer. Concatenate the resulting
        alignment and return it.
        :param motif_idx: start index in the motif sequence
        :param aa_idx: start index in the aa sequence
        :return: a tuple: (regular expression of the right extension (intersection of the alignment of motif and aa)
                            , corresponding alignment of the wild sequence, upper bound index of aa -inclusive)
        """
        try:
            align_extend = ''
            align_wild = ''
            l, aa_idx_upper_bound = self.get_trimers_right(motif_idx, aa_idx)
            for motif_trimer, wild_trimer, aa in l:
                align_trimer = self.align_aa_trimer(aa, motif_trimer)
                if align_trimer != 'none':
                    align_extend += align_trimer
                    align_wild += wild_trimer
                else:
                    return 'none'
            return align_extend, align_wild, aa_idx_upper_bound
        except ExtensionOutOfRangeException:
            return 'none'

    def get_trimers_right(self, motif_idx, aa_idx):
        """
        generate the list of trimers of the right extension. This is a  list of tuples: (motif trimer, wild trimer, aa)
        The trimer of the tuple will be aligned by the method extend_right. We can notice that these trimers are aligned
        to codons of aa. In case the extension out ranges the nt sequence, an exception ExtensionOutOfRangeException
        is thrown
        :param motif_idx: start index in the motif sequence (the nail)
        :param aa_idx: start index in the aa sequence
        :return: list of tuples (motif_trimer, nt_trimer, aa)
        """
        l_m = [self.seq_motif[i:i + 3].ljust(3, 'N')
               for i in range(motif_idx, len(self.seq_motif), 3)]  # trimers of motif for right extension
        l_w = []  # trimers of the wild nt sequence for the right extension
        l_a = []  # aa acids for right extension
        aa_idx_upper_bound = None  # inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx, aa_idx + len(l_m))]:
            if i_nt + 3 > self.len_seq_wild_nt:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_w.append(self.seq_wild_nt[i_nt:i_nt + 3])
                l_a.append(self.seq_wild_aa[i_aa])
                aa_idx_upper_bound = i_aa
        if aa_idx_upper_bound is None:
            raise ExtensionOutOfRangeException(l_m)
        return [(m, n, a) for m, n, a in zip(l_m, l_w, l_a)], aa_idx_upper_bound  # here is the alignment

    def extend_left(self, motif_idx, aa_idx):
        """
        same description as extend_right, just this method is doing it for the left side
        :param motif_idx: start index in the motif sequence
        :param aa_idx: aa_idx: start index in the aa sequence
        :return: a tuple: (regular expression of the right extension (intersection of the alignment of motif and aa)
                            , corresponding alignment of the wild sequence, upper bound index of aa -inclusive)
        """
        try:
            align_extend = ''
            align_wild = ''
            l, aa_idx_lower_bound = self.get_trimers_left(motif_idx, aa_idx)
            for motif_timer, wild_trimer, aa in l:
                align_trimer = self.align_aa_trimer(aa, motif_timer)
                if align_trimer != 'none':
                    align_extend = align_trimer + align_extend
                    align_wild = wild_trimer + align_wild
                else:
                    return 'none'
            return align_extend, align_wild, aa_idx_lower_bound
        except ExtensionOutOfRangeException:
            return 'none'

    def get_trimers_left(self, motif_idx, aa_idx):
        """
        same description as for get_trimers_right, just this method is doing it for the left side
        :param motif_idx: start index in the motif sequence (the nail)
        :param aa_idx: start index in the aa sequence
        :return: list of tuples (motif_trimer, nt_trimer, aa)
        """
        l_m = [self.seq_motif[y - 3:y].rjust(3, 'N') if y - 3 >= 0 else self.seq_motif[0:y].rjust(3, 'N')
               for y in range(motif_idx, 0, - 3)]
        l_w = []
        l_a = []
        aa_idx_lower_bound = None  # inclusive
        for i_aa, i_nt in [(i_aa, i_aa * 3)for i_aa in range(aa_idx - 1, aa_idx - len(l_m) - 1, -1)]:
            if i_nt < 0:
                raise ExtensionOutOfRangeException(i_nt)
            else:
                l_w.append(self.seq_wild_nt[i_nt: i_nt + 3])
                l_a.append(self.seq_wild_aa[i_aa])
                aa_idx_lower_bound = i_aa
        if aa_idx_lower_bound is None:
            raise ExtensionOutOfRangeException(l_m)
        return [(m, n, a) for m, n, a in zip(l_m, l_w, l_a)], aa_idx_lower_bound

    @staticmethod
    def align_aa_trimer(aa, trimer):
        """
        align the aa (codon) to a trimer. Only exact match is allowed. If the aa and the trimer cannot be aligned,
        we return 'none'. Otherwise a regular expression of the alignment is returned. Alignment is based on the
        reverse translation, so a nucleotide (iupac) can be aligned to a set of nucleotides - depending on he codon.
        This alignment is, in other words, the intersection of the codon of aa and the motif trimer
        :param aa: amino acid
        :param trimer: should be the motif trimer - the iupac annotation is accepted
        :return: regular expression of the trimer resultant from the intersection of the codon and motif
        """
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

    def update_results(self, align_optimized, align_wild, aa_idx_lower_bound, aa_idx_upper_bound):
        """
        if the alignment corresponds to exact match - align_optimized, this method updates both dictionaries:
        self.existent_motifs includes the existing motifs in the wild sequence, and self.motifs_to_integrate include
        the motifs that can be integrated in the wild sequence. These two dictionaries are populated by results of the
        following form: {'seq_motif':[(optimized alignment which is the intersection of the motif and aa, corresponding
        wild alignment, (aa index lower bound - inclusive, aa index upper bound -inclusive)), ]}
        :param align_optimized: the alignment of the intersection of the motif and aa
        :param align_wild: corresponding alignment on the wild sequence
        :param aa_idx_lower_bound: aa start index of the alignment - inclusive
        :param aa_idx_upper_bound: aa end index of the alignment - inclusive
        :return:
        """
        res = (align_optimized, align_wild, (aa_idx_lower_bound, aa_idx_upper_bound),
               (aa_idx_lower_bound * 3, aa_idx_upper_bound * 3 + 3))
        if re.match(align_optimized, align_wild):
            if self.seq_motif in self.existent_motifs.keys():
                if res not in self.existent_motifs[self.seq_motif]:
                    self.existent_motifs[self.seq_motif].append(res)
            else:
                self.existent_motifs[self.seq_motif] = [res]
        else:
            if self.seq_motif in self.motifs_to_integrate.keys():
                if res not in self.motifs_to_integrate[self.seq_motif]:
                    self.motifs_to_integrate[self.seq_motif].append(res)
            else:
                self.motifs_to_integrate[self.seq_motif] = [res]


def run_with_one_motif():
    seq_motif = 'WNNVYTAATTARYYNNN'
    with open('data/f8.orig.fasta', 'r') as f_seq_nt:
        seq_nt = str(SeqIO.read(f_seq_nt, 'fasta').seq)
    b = BlastIUPAC(seq_nt)
    b.run(seq_motif)
    print b.motifs_to_integrate


def run_with_multiple_motif():
    with open('data/f8.orig.fasta', 'r') as f_seq_nt,\
            open('data/liver_motifs.fasta') as f_motifs:
        seq_nt = str(SeqIO.read(f_seq_nt, 'fasta').seq)
        l_motif = [str(record.seq) for record in SeqIO.parse(f_motifs, 'fasta')]
    b = BlastIUPAC(seq_nt)
    b.run_batch(l_motif)
    print b.motifs_to_integrate


if __name__ == '__main__':
    # run_with_one_motif()
    run_with_multiple_motif()
