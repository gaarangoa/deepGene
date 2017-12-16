from unittest import TestCase
from blast_iupac import BlastIUPAC as Blast


class TestBlastIUPAC(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBlastIUPAC, self).__init__(*args, **kwargs)
        self.blast = Blast(seq_nt='ATGCAAATAGAGCTCTCCACC'
                      , seq_aa='MQIELST')

    def test_match_aa_trimer(self):
        self.assertEqual('GT[TG]', Blast.match_aa_trimer(aa='V', trimer='GTK'))
        self.assertEqual('none', Blast.match_aa_trimer(aa='A', trimer='CCT'))

    def test_get_trimers_right(self):
        self.assertEqual([('M', 'TGT'), ('Q', 'GGT'), ('I', 'KNN')]
                         , self.blast.get_trimers_right(seq_motif='TGTGGTK', idx_motif=0, loc_aa=0))
        self.assertEqual([('M', 'TGT'), ('Q', 'GGT'), ('I', 'KKK')]
                         , self.blast.get_trimers_right(seq_motif='TGTGGTKKK', idx_motif=0, loc_aa=0))
        self.assertEqual([('M', 'TGT'), ('Q', 'GGT'), ('I', 'KKN')]
                         , self.blast.get_trimers_right(seq_motif='TGTGGTKK', idx_motif=0, loc_aa=0))
        self.assertEqual([('Q', 'TGG'), ('I', 'TKK')]
                         , self.blast.get_trimers_right(seq_motif='TGTGGTKK', idx_motif=2, loc_aa=1))

    def test_get_trimers_left(self):
        self.assertEqual([('E', 'TGG'), ('I', 'NTG')]
                         , self.blast.get_trimers_left(seq_motif='TGTGGTKK', idx_motif=5, loc_aa=3))
        self.assertEqual([('E', 'GGT'), ('I', 'TGT')]
                         , self.blast.get_trimers_left(seq_motif='TGTGGTKK', idx_motif=6, loc_aa=3))
        self.assertEqual([('I', 'ATC'), ('Q', 'NAG')]
                         , self.blast.get_trimers_left(seq_motif='AGATCGARTTRT', idx_motif=5, loc_aa=2))

    def test_extend_right(self):
        self.assertEqual('none'
                         , self.blast.extend_right(seq_motif='TGTGGTKK', idx_motif=5, loc_aa=3))
        self.assertEqual('GA[AG]TT[AG]TC[ACTG]'
                         , self.blast.extend_right(seq_motif='TGTGGGARTTRT', idx_motif=5, loc_aa=3))

    def test_extend_left(self):
        self.assertEqual('none'
                         , self.blast.extend_left(seq_motif='TGTGGGARTTRT', idx_motif=5, loc_aa=2))
        self.assertEqual('none'
                         , self.blast.extend_left(seq_motif='AAGATCGARTTRT', idx_motif=5, loc_aa=2))
