
__version__ = '1.1.0'
__author__ = 'Dhoha Abid'


import helper_functions as helper
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from random import randint
from random import sample
import traceback
from blast_iupac import BlastIUPAC
from helper_functions import translate_nt_sequence
from json import dump


class OptimizedSeqException(Exception):
    def __init__(self, msg):
        self.msg = msg


class OptimizeSeq:
    def __init__(self, seq_wild):
        self.seq_wild = seq_wild
        self.seq_opt = ['N' for i in range(len(seq_wild))]
        self.included_motifs = {}

    def run(self, d_motifs, f_name):
        self.optimize_wild_seq(d_motifs)
        self.check_optimized_seq()
        self.write_optimized_seq(f_name)

    def optimize_wild_seq(self, d_motifs):
        """
        create a string of same length as input. Assign all its nt the N character. update these characters using
        the motifs, then complete the missing values  from the wild sequence. Sometimes the same motif can be integrated
        using different nucleotide sequences, we pickup these sequences randomly. In a later step, we can be specific
        :param d_motifs: motif sequence
        :return:
        """
        for seq_motif, res in d_motifs.items():
            self.included_motifs[seq_motif] = []
            for opt, wild, aa_bound, nt_bound in res:
                start = nt_bound[0]
                end = nt_bound[1]
                if len(set(self.seq_opt[start:end])) == 1:  # check if this part of sequence has changed by other motifs
                    for i, nt in enumerate(helper.transform_reg_exp_to_list(opt), start):
                        try:
                            if type(nt).__name__ == 'str':
                                self.seq_opt[i] = nt
                            elif type(nt).__name__ == 'list':
                                self.seq_opt[i] = nt[randint(0, len(nt) - 1)]
                        except Exception:
                            traceback.print_exc()
                            sys.exit(1)
                    self.included_motifs[seq_motif].append((opt, wild, aa_bound, nt_bound))
        for i, opt, nt in zip([i for i in range(len(self.seq_wild))], self.seq_opt, self.seq_wild):
            if opt == 'N':
                self.seq_opt[i] = nt
        self.seq_opt = reduce(lambda x, y: x + y, self.seq_opt)

    def check_optimized_seq(self):
        """
        check if the sequence (can be string or Seq) translates to the same aa of the input (original) sequence
        :param seq: nucleotide sequence (can be string or Seq)
        :return: raise an error if the sequence does not translate to same aa of the input
        """
        if type(self.seq_opt).__name__ == 'str':
            seq = Seq(self.seq_opt)
        try:
            seq_translated = seq.translate()
            if seq_translated != translate_nt_sequence(self.seq_wild):
                raise OptimizedSeqException(seq_translated)
        except OptimizedSeqException as e:
            print e
            sys.exit(1)

    def write_optimized_seq(self, f_name):
        with open('res/' + f_name + '.fasta', 'w') as f_seq_opt\
                , open('res/' + f_name + '.json', 'w') as f_meta_seq_opt:
            SeqIO.write(SeqRecord(Seq(self.seq_opt), id='liver.motif', description=f_name), f_seq_opt, 'fasta')
            dump(self.included_motifs, f_meta_seq_opt)


def run_with_random_motifs():
    with open('data/f8.orig.fasta', 'r') as f_seq_nt, open('data/liver_motifs.fasta') as f_motifs:
        seq_nt = str(SeqIO.read(f_seq_nt, 'fasta').seq)
        l_motif = [str(record.seq) for record in SeqIO.parse(f_motifs, 'fasta')]
    b = BlastIUPAC(seq_nt)
    b.run_batch(l_motif)
    o = OptimizeSeq(seq_nt)
    o.run(b.motifs_to_integrate, 'test')
    # print len(b.res_motif_optimize)
    # sample(b.res_motif_optimize, 15)
    # TODO send the name of the file


if __name__ == '__main__':
    run_with_random_motifs()
