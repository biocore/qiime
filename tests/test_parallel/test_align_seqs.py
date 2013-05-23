#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from shutil import rmtree
from glob import glob
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.parallel.align_seqs import ParallelAlignSeqsPyNast
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename,
                        count_seqs_in_filepaths)
from qiime.test import initiate_timeout, disable_timeout

class ParallelAlignSeqsTests(TestCase):
    
    def setUp(self):
        """ """
        self.files_to_remove = []
        self.dirs_to_remove = []
        
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_parallel_tests_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)
        
        self.template_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='qiime_template',
                                            suffix='.fasta')
        template_f = open(self.template_fp,'w')
        template_f.write(pynast_test1_template_fasta)
        template_f.close()
        self.files_to_remove.append(self.template_fp)
        
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
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

class ParallelAlignSeqsPyNastTests(ParallelAlignSeqsTests):

    def test_parallel_align_seqs_pynast(self):
        """ parallel_align_seqs_pynast functions as expected """

        params = {
         'min_percent_id':0.75,
         'min_length':15,
         'template_fp':self.template_fp,
         'pairwise_alignment_method':'uclust',
         'blast_db':None
        }
        
        app = ParallelAlignSeqsPyNast()
        r = app(self.inseqs1_fp,
                self.test_out,
                params,
                job_prefix='PTEST',
                poll_directly=True,
                suppress_submit_jobs=False)
        # confirm that the total number of output sequences equals the total
        # number of input sequences
        num_input_seqs = count_seqs_in_filepaths([self.inseqs1_fp])[1]
        num_template_seqs = count_seqs_in_filepaths([self.template_fp])[1]
        num_output_seqs = \
         count_seqs_in_filepaths(glob(join(self.test_out,'*fasta')))[1] \
         - num_input_seqs - num_template_seqs
        self.assertEqual(num_input_seqs,num_output_seqs)


pynast_test1_template_fasta = """>1
ACGT--ACGTAC-ATA-C-----CC-T-G-GTA-G-T---
>2
AGGTTTACGTAG-ATA-C-----CC-T-G-GTA-G-T---
>3
AGGTACT-CCAC-ATA-C-----CC-T-G-GTA-G-T---
>4
TCGTTCGT-----ATA-C-----CC-T-G-GTA-G-T---
>5
ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---
"""

inseqs1 = """>1 description field
ACCTACGTTAATACCCTGGTAGT
>2
ACCTACGTTAATACCCTGGTAGT
>3
AA
"""


if __name__ == "__main__":
    main()