import helper_functions as helper
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from random import randint
import traceback


def process_res(blast):
    """
    create a string of same length as input. Assign all its nt the N character. update these characters using
    the motifs, then complete the missing values  from the wild sequence

    :param self:
    :param blast:
    :return:
    """
    seq_opt = ['N' for i in range(len(blast.seq_nt))]
    for seq_motif, res in blast.res_motif_optimize.items():
        for opt, wild, aa_bound, nt_bound in res:
            start = nt_bound[0]
            end = nt_bound[1]
            if len(set(seq_opt[start:end])) == 1: # check if this portion of the sequence has changed by other motifs
                for i, nt in enumerate(helper.transform_reg_exp_to_list(opt[0]), start):
                    try:
                        if type(nt).__name__ == 'str':
                            seq_opt[i] = nt
                        elif type(nt).__name__ == 'list':
                            seq_opt[i] = nt[randint(0, len(nt) - 1)]
                    except:
                        traceback.print_exc()
                        sys.exit(1)
    j = 0
    l = [None] * len(blast.seq_nt)
    for s1, s2 in zip(seq_opt, blast.seq_nt):
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
