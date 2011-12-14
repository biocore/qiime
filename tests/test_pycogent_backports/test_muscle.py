#!/usr/bin/env python

from os import getcwd, remove, rmdir, mkdir, path
import tempfile, shutil
from cogent.core.moltype import RNA, DNA
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import flatten
from qiime.pycogent_backports.muscle import Muscle, muscle_seqs, aln_tree_seqs, \
        align_unaligned_seqs, build_tree_from_alignment, \
        align_and_build_tree, add_seqs_to_alignment, align_two_alignments

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Catherine Lozupone", "Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

class GeneralSetUp(TestCase):

    def setUp(self):
        """Muscle general setUp method for all tests"""
        self.seqs1 = ['ACUGCUAGCUAGUAGCGUACGUA','GCUACGUAGCUAC',
            'GCGGCUAUUAGAUCGUA']
        
        self.labels1 = ['>1','>2','>3']
        self.lines1 = flatten(zip(self.labels1,self.seqs1))

        self.seqs2=['UAGGCUCUGAUAUAAUAGCUCUC','UAUCGCUUCGACGAUUCUCUGAUAGAGA',
            'UGACUACGCAU']
        self.labels2=['>a','>b','>c']
        self.lines2 = flatten(zip(self.labels2,self.seqs2))
        
        self.temp_dir = tempfile.mkdtemp()
        self.temp_dir_spaces = '/tmp/test for muscle/'
        try:
            mkdir(self.temp_dir_spaces)
        except OSError:
            pass
        try:
            #create sequence files
            f = open(path.join(self.temp_dir, 'seq1.txt'),'w')
            f.write('\n'.join(self.lines1))
            f.close()
            g = open(path.join(self.temp_dir, 'seq2.txt'),'w')
            g.write('\n'.join(self.lines2))
            g.close()
        except OSError:
            pass

    def tearDown(self): 
        """cleans up all files initially created"""
        # remove the tempdir and contents
        shutil.rmtree(self.temp_dir)
        shutil.rmtree(self.temp_dir_spaces)

class MuscleTests(GeneralSetUp):
    """Tests for the Muscle application controller"""

    def test_base_command(self):
        """Muscle BaseCommand should return the correct BaseCommand"""
        c = Muscle()
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle']))
        c.Parameters['-in'].on('seq.txt')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle -in "seq.txt"']))
        c.Parameters['-cluster2'].on('neighborjoining')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "',getcwd(),'/"; ','muscle -cluster2 neighborjoining' +
            ' -in "seq.txt"']))

    def test_maxmb(self):
        """maxmb option should not break Muscle"""
        app = Muscle()
        app.Parameters['-maxmb'].on('250')
        outfile = tempfile.NamedTemporaryFile()
        app.Parameters['-out'].on(outfile.name)
        
        infile = tempfile.NamedTemporaryFile()
        infile.write(
            ">Seq1\nAAAGGGTTTCCCCT\n"
            ">Seq2\nAAAGGGGGTTTCCACT\n")
        infile.flush()
        result = app(infile.name)

        observed = result['MuscleOut'].read()
        expected = (
            ">Seq1\nAAA--GGGTTTCCCCT\n"
            ">Seq2\nAAAGGGGGTTTCCACT\n"
            )
        self.assertEqual(observed, expected)

    def test_changing_working_dir(self):
        """Muscle BaseCommand should change according to WorkingDir"""
        c = Muscle(WorkingDir='/tmp/muscle_test')
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/muscle_test','/"; ','muscle']))
        c = Muscle()
        c.WorkingDir = '/tmp/muscle_test2'
        self.assertEqual(c.BaseCommand,\
            ''.join(['cd "','/tmp/muscle_test2','/"; ','muscle']))
        
        #removing the dirs is proof that they were created at the same time
        #if the dirs are not there, an OSError will be raised
        rmdir('/tmp/muscle_test')
        rmdir('/tmp/muscle_test2')
    
    def test_aln_tree_seqs(self):
        "aln_tree_seqs returns the muscle alignment and tree from iteration2"
        tree, aln = aln_tree_seqs(path.join(self.temp_dir, 'seq1.txt'), 
                                   tree_type="neighborjoining",
                                   WorkingDir=self.temp_dir,
                                   clean_up=True)
        self.assertEqual(str(tree), '((1:1.125,2:1.125):0.375,3:1.5);')
        self.assertEqual(len(aln), 6)
        self.assertEqual(aln[-2], '>3\n')
        self.assertEqual(aln[-1], 'GCGGCUAUUAGAUCGUA------\n')

    def test_aln_tree_seqs_spaces(self):
        "aln_tree_seqs should work on filename with spaces"
        try:
            #create sequence files
            f = open(path.join(self.temp_dir_spaces, 'muscle_test_seq1.txt'),'w')
            f.write('\n'.join(self.lines1))
            f.close()
        except OSError:
            pass
        tree, aln = aln_tree_seqs(path.join(self.temp_dir_spaces,\
                                    'muscle_test_seq1.txt'), 
                                    tree_type="neighborjoining",
                                    WorkingDir=getcwd(),
                                    clean_up=True)
        self.assertEqual(str(tree), '((1:1.125,2:1.125):0.375,3:1.5);')
        self.assertEqual(len(aln), 6)
        self.assertEqual(aln[-2], '>3\n')
        self.assertEqual(aln[-1], 'GCGGCUAUUAGAUCGUA------\n')
        remove(self.temp_dir_spaces+'/muscle_test_seq1.txt')

    def test_align_unaligned_seqs(self):
        """align_unaligned_seqs should work as expected"""
        res = align_unaligned_seqs(self.seqs1, RNA)
        self.assertEqual(res.toFasta(), align1)

    def test_build_tree_from_alignment(self):
        """Muscle should return a tree built from the passed alignment"""
        tree_short = build_tree_from_alignment(build_tree_seqs_short, DNA)
        num_seqs = flatten(build_tree_seqs_short).count('>')
        self.assertEqual(len(tree_short.tips()), num_seqs)

        tree_long = build_tree_from_alignment(build_tree_seqs_long, DNA)
        seq_names = []
        for line in build_tree_seqs_long.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree_long.tips():
            if node.Name not in seq_names:
                self.fail()

    def test_align_and_build_tree(self):
        """Should align and build a tree from a set of sequences"""
        res = align_and_build_tree(self.seqs1, RNA)
        self.assertEqual(res['Align'].toFasta(), align1)

        tree = res['Tree']
        seq_names = []
        for line in align1.split('\n'):
            if line.startswith('>'):
                seq_names.append(line[1:])

        for node in tree.tips():
            if node.Name not in seq_names:
                self.fail()

    def test_add_seqs_to_alignment(self):
        """Should add sequences to an alignment"""
        res = add_seqs_to_alignment(seqs_to_add, align1)
        self.assertEqual(res.toFasta(), added_align_result)

    def test_align_two_alignments(self):
        """Should align to multiple sequence alignments"""
        res = align_two_alignments(align1, aln_to_merge)
        self.assertEqual(res.toFasta(), merged_align_result)

align1 = ">seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n---GCUACGUAGCUAC-------\n>seq_2\nGCGGCUAUUAGAUCGUA------"

# for use in test_add_seqs_to_alignment()
seqs_to_add = ">foo\nGCUACGUAGCU\n>bar\nGCUACGUAGCC"
added_align_result = ">bar\n---GCUACGUAGCC---------\n>foo\n---GCUACGUAGCU---------\n>seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n---GCUACGUAGCUAC-------\n>seq_2\nGCGGCUAUUAGAUCGUA------"

# for use in test_align_two_alignments()
aln_to_merge = ">foo\nGCUACGUAGCU\n>bar\n--UACGUAGCC"
merged_align_result = ">bar\n-----UACGUAGCC---------\n>foo\n---GCUACGUAGCU---------\n>seq_0\nACUGCUAGCUAGUAGCGUACGUA\n>seq_1\n---GCUACGUAGCUAC-------\n>seq_2\nGCGGCUAUUAGAUCGUA------"

build_tree_seqs_short = """>muscle_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>muscle_test_seqs_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>muscle_test_seqs_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>muscle_test_seqs_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>muscle_test_seqs_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>muscle_test_seqs_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>muscle_test_seqs_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>muscle_test_seqs_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>muscle_test_seqs_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>muscle_test_seqs_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""

build_tree_seqs_long = """>muscle_test_seqs_0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
AGCTTTAAATCATGCCAGTG
>muscle_test_seqsaaaaaaaa_1
GACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
TGCTTTCAATAATGCCAGTG
>muscle_test_seqsaaaaaaaa_2
AACCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
TGCTTTGAATCATGCCAGTA
>muscle_test_seqsaaaaaaaa_3
AAACCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
TGCTTTACATCATGCAAGTG
>muscle_test_seqsaaaaaaaa_4
AACCGCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
TGCTTTAAATCATGCCAGTG
>muscle_test_seqsaaaaaaaa_5
AACCCCCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
TGCTTTAAATCATGCCAGTT
>muscle_test_seqsaaaaaaaa_6
GACCCCCGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
TACTTTAGATCATGCCGGTG
>muscle_test_seqsaaaaaaaa_7
AACCCCCACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
TGCTTTAAATCATGCCAGTG
>muscle_test_seqsaaaaaaaa_8
AACCCCCACGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
TGCATTAAATCATGCCAGTG
>muscle_test_seqsaaaaaaaa_9
AAGCCCCACGGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
TGCTTTAAATCCTGACAGCG
"""


if __name__ == '__main__':
    main()
