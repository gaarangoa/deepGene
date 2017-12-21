
__version__ = '1.1.0'
__author__ = 'Dhoha Abid'

from unittest import TestCase
import unittest
from deepGene.optimize_seq.blast_iupac import BlastIUPAC as Blast
from deepGene.optimize_seq.blast_iupac import ExtensionOutOfRangeException


class TestBlastIUPAC(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBlastIUPAC, self).__init__(*args, **kwargs)
        self.blast = Blast('ATGCAAATAGAGCTCTCCACC')

    def test_match_aa_trimer(self):
        self.assertEqual('GT[TG]', Blast.align_aa_trimer(aa='V', trimer='GTK'))
        self.assertEqual('none', Blast.align_aa_trimer(aa='A', trimer='CCT'))

    def test_get_trimers_right(self):
        self.blast.seq_motif = 'CAGATC'
        self.assertEqual(([('CAG', 'CAA', 'Q'), ('ATC', 'ATA', 'I')], 2)
                         , self.blast.get_trimers_right(motif_idx=0, aa_idx=1))
        self.blast.seq_motif = 'TCTACCR'
        self.assertRaises(ExtensionOutOfRangeException, self.blast.get_trimers_right, 0, 5)
        with self.assertRaises(ExtensionOutOfRangeException):
            self.blast.get_trimers_right(0, 5)

    def test_extend_right(self):
        self.blast.seq_motif = 'TGTGGTKK'
        self.assertEqual('none', self.blast.extend_right(motif_idx=5, aa_idx=3))
        self.blast.seq_motif = 'TGTGGGARTTRT'
        self.assertEqual(('GA[AG]TT[AG]TC[ACTG]', 'GAGCTCTCC', 5), self.blast.extend_right(motif_idx=5, aa_idx=3))

    def test_get_trimers_left(self):
        self.blast.seq_motif = 'CAGATC'
        self.assertEqual(([('CAG', 'CAA', 'Q')], 1), self.blast.get_trimers_left(motif_idx=3, aa_idx=2))
        self.blast.seq_motif = 'RATGCAG'
        self.assertRaises(ExtensionOutOfRangeException, self.blast.get_trimers_left, 4, 1)
        self.blast.seq_motif = 'NNAGATC'
        self.assertEqual(([('NAG', 'CAA', 'Q'), ('NNN', 'ATG', 'M')], 0),
                         self.blast.get_trimers_left(motif_idx=4, aa_idx=2))

    def test_extend_left(self):
        self.blast.seq_motif = 'TGTGGGARTTRT'
        self.assertEqual('none', self.blast.extend_left(motif_idx=5, aa_idx=2))
        self.blast.seq_motif = 'AAGATCGARTTRT'
        self.assertEqual('none', self.blast.extend_left(motif_idx=5, aa_idx=2))
        self.blast.seq_motif = 'ARTTRTCN'
        self.assertEqual(('GA[AG]TT[AG]', 'GAGCTC', 3), self.blast.extend_left(motif_idx=5, aa_idx=5))


if __name__ == '__main__':
    unittest.main()

