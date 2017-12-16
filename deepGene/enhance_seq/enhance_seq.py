import encode_iupac_to_nt as iupac
import reverse_translate as rv
import blast_iupac
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import math
from functools import reduce
import re


class ProcessFile:
    def __init__(self, seq_file='test/f8.orig.fasta', motif_file='test/enhancer_all.fasta'):
        self.seq_file = seq_file
        self.motif_file = motif_file
        with open(self.seq_file, 'r') as seq_f \
                , open(self.motif_file, 'r') as motif_f:
            self.nt_seq = SeqIO.read(seq_f, 'fasta').seq
            self.aa_seq = str(self.nt_seq.translate())
            self.l_motif = [str(enh.seq) for enh in SeqIO.parse(motif_f, 'fasta')]
        self.res = None

    def wrapper_blast_iupac(self):
        blast_iupac.run()

    def print_res(self):
        print self.res

    def dump_json(self):
        pass

    def run(self):
        pass


if __name__ == '__main__':

    v = ProcessFile()
    print(v.aa_seq)