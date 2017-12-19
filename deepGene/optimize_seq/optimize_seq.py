
__version__ = '1.1.0'
__author__ = 'Dhoha Abid'


import helper_functions as helper
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from random import randint
import traceback


class OptimizedSeqException(Exception):
    def __init__(self, msg):
        self.msg = msg


class OptimizeSeq:
    def __init__(self, blast):
        self.blast = blast

    def optimize_seq_random_nt_choice(self, blast):
        """
        create a string of same length as input. Assign all its nt the N character. update these characters using
        the motifs, then complete the missing values  from the wild sequence. Sometimes the same motif can be integrated
        using different nucleotide sequences, we pickup these sequences randomly. In a later step, we can be specific
        :return: different fasta files of the enhanced sequence with metadata
        """
        self.blast = blast
        len_seq_nt = len(self.blast.seq_nt)
        seq_opt = ['N' for i in range(len_seq_nt)]
        included_motifs = {}
        for seq_motif, res in self.blast.res_motif_optimize.items():
            included_motifs[seq_motif] = []
            for opt, wild, aa_bound, nt_bound in res:
                start = nt_bound[0]
                end = nt_bound[1]
                if len(set(seq_opt[start:end])) == 1:  # check if this portion of the sequence has changed by other motifs
                    for i, nt in enumerate(helper.transform_reg_exp_to_list(opt), start):
                        try:
                            if type(nt).__name__ == 'str':
                                seq_opt[i] = nt
                            elif type(nt).__name__ == 'list':
                                seq_opt[i] = nt[randint(0, len(nt) - 1)]
                        except Exception:
                            traceback.print_exc()
                            sys.exit(1)
                    included_motifs[seq_motif].append(res)

        for i, opt, nt in zip([i for i in range(len_seq_nt)], seq_opt, self.blast.seq_nt):
            if opt == 'N':
                seq_opt[i] = nt
        seq_opt = reduce(lambda x, y: x + y, seq_opt)
        self.check_seq_opt(seq_opt)
        return seq_opt

    def check_seq_opt(self, seq):
        """
        check if the sequence (can be string or Seq) translates to the same aa of the input (original) sequence
        :param seq: nucleotide sequence (can be string or Seq)
        :return: raise an error if the sequence does not translate to same aa of the input
        """
        if type(seq).__name__ == 'str':
            seq = Seq(seq)
        try:
            seq_translated = seq.translate()
            if seq_translated != self.blast.seq_aa:
                raise OptimizedSeqException(seq_translated)
        except OptimizedSeqException as e:
            print e
            sys.exit(1)
        with open('res/res.liver.motif.v2.fasta', 'w') as f:
            SeqIO.write(SeqRecord(seq, id='liver.motif.v1'), f, 'fasta')

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
