#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import get_qiime_temp_dir, get_tmp_filename
from qiime.functional_assignment import usearch_function_assigner
from qiime.parse import parse_otu_map
from qiime.test import initiate_timeout, disable_timeout

class FunctionalAssignmentTests(TestCase):
    """
    """
    
    def setUp(self):
        """
        """
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_parallel_tests_',
                                         suffix='',
                                         result_constructor=str)
        create_dir(self.test_out)
        self.dirs_to_remove = [self.test_out]
        
        self.output_fp = join(self.test_out,'fmap.txt')
        self.failure_fp = join(self.test_out,'fail.txt')
        self.usearch_fp = join(self.test_out,'out.uc')
        self.bl6_fp = join(self.test_out,'out.bl6')
        self.log_fp = join(self.test_out,'fmap.log')
        self.files_to_remove = [self.output_fp, self.failure_fp, 
                                self.usearch_fp, self.log_fp, self.bl6_fp]
        
        self.refseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_refseqs',
                                            suffix='.fasta')
        refseqs1_f = open(self.refseqs1_fp,'w')
        refseqs1_f.write(refseqs1)
        refseqs1_f.close()
        self.files_to_remove.append(self.refseqs1_fp)
        
        self.inseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_inseqs',
                                            suffix='.fasta')
        inseqs1_f = open(self.inseqs1_fp,'w')
        inseqs1_f.write(inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)
        initiate_timeout(60)
    
    def tearDown(self):
        """ """
        disable_timeout()
        remove_files(self.files_to_remove,error_on_missing=False)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

class UsearchFunctionalAssignmentTests(FunctionalAssignmentTests):
    
    def test_usearch_function_assigner(self):
        """usearch_function_assigner functions as expected """
        usearch_function_assigner(query_fp=self.inseqs1_fp,
                              refseqs_fp=self.refseqs1_fp,
                              output_fp=self.output_fp,
                              failure_fp=self.failure_fp,
                              usearch_fp=self.usearch_fp,
                              blast6_fp=self.bl6_fp,
                              log_fp=self.log_fp,
                              evalue=1e-10,
                              min_id=0.75,
                              queryalnfract=0.35,
                              targetalnfract=0.0,
                              maxaccepts=1,
                              maxrejects=8,
                              temp_dir=get_qiime_temp_dir(),
                              HALT_EXEC=False)
        fmap = parse_otu_map(open(self.output_fp,'U'))
        self.assertEqual(len(fmap[0]),2)
        self.assertEqualItems(fmap[1],['eco:b0015', 'eco:b0122'])
        self.assertEqualItems(fmap[2],['eco:b0015-pr', 'eco:b0122-pr'])
        

refseqs1 = """>eco:b0001-pr
MKRISTTITTTITITTGNGAG
>eco:b0015-pr dnaJ
MAKQDYYEILGVSKTAEEREIRKAYKRLAMKYHPDRNQGDKEAEAKFKEIKEAYEVLTDS
QKRAAYDQYGHAAFEQGGMGGGGFGGGADFSDIFGDVFGDIFGGGRGRQRAARGADLRYN
MELTLEEAVRGVTKEIRIPTLEECDVCHGSGAKPGTQPQTCPTCHGSGQVQMRQGFFAVQ
QTCPHCQGRGTLIKDPCNKCHGHGRVERSKTLSVKIPAGVDTGDRIRLAGEGEAGEHGAP
AGDLYVQVQVKQHPIFEREGNNLYCEVPINFAMAALGGEIEVPTLDGRVKLKVPGETQTG
KLFRMRGKGVKSVRGGAQGDLLCRVVVETPVGLNERQKQLLQELQESFGGPTGEHNSPRS
KSFFDGVKKFFDDLTR
>eco:b0122-pr
MKTFFRTVLFGSLMAVCANSYALSESEAEDMADLTAVFVFLKNDCGYQNLPNGQIRRALV
FFAQQNQWDLSNYDTFDMKALGEDSYRDLSGIGIPVAKKCKALARDSLSLLAYVK
"""

inseqs1 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (N)
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (N)
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>eco:b0122
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa
"""


if __name__ == "__main__":
    main()