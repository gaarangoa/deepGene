import unittest
import sys
from ..blast_iupac import BlastIUPAC as Blast
from Bio import SeqIO


class TestBlastIUPAC(unittest.TestCase):

    def __init__(self):
        with open('f8.orig.fasta', 'r') as f_orig:
            self.seq_nt = SeqIO.read(f_orig, 'fasta').seq
            self.seq_aa = str(self.seq_nt.translate())
            self.l_motif = [seq_motif for seq_motif in SeqIO.parse('enhancer_all.fasta')]

    def test_match_aa_trimer(self):
        self.assertEqual('T[TG]', Blast.match_aa_trimer(aa='V', trimer='GTK'))

    def test_get_trimers_right(self):
        pass

    def test_get_trimers_left(self):
        pass


if __name__ == '__main__':
    print('m here')
    unittest.main()
