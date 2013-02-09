#!/usr/bin/env python
# File created on 09 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout

import os
from shutil import rmtree
from glob import glob
from tarfile import open as open_tarfile
from os.path import join, exists, getsize, split, splitext, dirname
from os import makedirs, system, chdir, getcwd
from numpy import array, absolute
from cogent import LoadTree, LoadSeqs
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import ApplicationNotFoundError
from qiime.util import get_tmp_filename
from cogent.parse.binary_sff import parse_binary_sff
from qiime.util import load_qiime_config, count_seqs
from qiime.parse import (parse_qiime_parameters,
                         parse_distmat_to_dict,
                         parse_distmat,
                         parse_taxa_summary_table)
from biom.parse import parse_biom_table
from qiime.test import initiate_timeout, disable_timeout, get_test_data_fps
from qiime.workflow import (no_status_updates,
                            WorkflowError,
                            get_params_str)
from qiime.core_qiime_analyses import run_core_qiime_analyses

class CoreQiimeAnalysesTests(TestCase):
    
    def setUp(self):
        
        self.files_to_remove = []
        self.dirs_to_remove = []
        
        # Create example output directory
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='core_qiime_analyses_test_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)
        
        # Get input data
        self.test_data = get_test_data_fps()
        
        self.qiime_config = load_qiime_config()
        
        # Define number of seconds a test can run for before timing out 
        # and failing
        initiate_timeout(120)

    
    def tearDown(self):
        
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_run_core_qiime_analyses_error_on_invalid_category(self):
        """run_core_qiime_analyses: error raised on invalid categories
        """
        # too few values in 'month' category
        self.assertRaises(ValueError,
                          run_core_qiime_analyses,
                          self.test_data['biom'][0],
                          self.test_data['map'][0],
                          20,
                          output_dir=self.test_out,
                          params=parse_qiime_parameters({}),
                          qiime_config=self.qiime_config,
                          categories=['SampleType','month'],
                          tree_fp=self.test_data['tree'][0],
                          parallel=False,
                          status_update_callback=no_status_updates)
        
        # invalid category name
        self.assertRaises(ValueError,
                          run_core_qiime_analyses,
                          self.test_data['biom'][0],
                          self.test_data['map'][0],
                          20,
                          output_dir=self.test_out,
                          params=parse_qiime_parameters({}),
                          qiime_config=self.qiime_config,
                          categories=['not-a-real-category'],
                          tree_fp=self.test_data['tree'][0],
                          parallel=False,
                          status_update_callback=no_status_updates)

    def test_run_core_qiime_analyses_serial(self):
        """run_core_qiime_analyses: functions (serially) using default qiime params
        """
        # this takes a long time, so use a longer sigalrm
        run_core_qiime_analyses(
                        self.test_data['biom'][0],
                        self.test_data['map'][0],
                        20,
                        output_dir=self.test_out,
                        params=parse_qiime_parameters({}),
                        qiime_config=self.qiime_config,
                        categories=['SampleType','days_since_epoch'],
                        tree_fp=self.test_data['tree'][0],
                        parallel=False,
                        status_update_callback=no_status_updates)
        
        # Basic sanity test of OTU table as details are tested 
        # in the pick_otus_through_otu_table tests
        self.assertTrue(exists('%s/bdiv_even20' % self.test_out))


    # def test_run_core_qiime_analyses_serial(self):
    #     """run_core_qiime_analyses: functions (serially) using default qiime params
    #     """
    #     # this takes a long time, so use a longer sigalrm
    #     restart_timeout(600)
    #     run_core_qiime_analyses(
    #         fna_fps=self.fasting_fna_fp,
    #         qual_fps=self.fasting_qual_fp,
    #         mapping_fp=self.fasting_mapping_fp,
    #         output_dir=self.wf_out,
    #         command_handler=call_commands_serially,
    #         params=self.run_core_qiime_analyses_params1,
    #         qiime_config=self.qiime_config,
    #         categories='BarcodeSequence',
    #         sampling_depth=100,
    #         arare_min_rare_depth=10,
    #         arare_num_steps=10,
    #         reference_tree_fp=None,
    #         parallel=False,
    #         status_update_callback=no_status_updates)
    #     
    #     # Basic sanity test of OTU table as details are tested 
    #     # in the pick_otus_through_otu_table tests
    #     otu_table_fp = join(self.wf_out,'otus','otu_table.biom')
    #     otu_table = parse_biom_table(open(otu_table_fp,'U'))
    #     sample_ids = list(otu_table.SampleIds)
    #     expected_sample_ids = ['PC.354','PC.355','PC.356','PC.481',
    #                            'PC.593','PC.607','PC.634','PC.635','PC.636']
    #     sample_ids.sort()
    #     expected_sample_ids.sort()
    #     self.assertEqual(sample_ids,expected_sample_ids)
    #     # even sampling directory exists
    #     self.assertTrue(exists('%s/bdiv_even100' % self.wf_out))
    # 
    # def test_run_core_qiime_analyses_serial_alt_params(self):
    #     """run_core_qiime_analyses: functions as expected (serially, alt params)
    #     """
    #     # Single category and sampling_depth=None
    #     # this takes a long time, so use a longer sigalrm
    #     restart_timeout(600)
    #     run_core_qiime_analyses(
    #         fna_fps=self.fasting_fna_fp,
    #         qual_fps=self.fasting_qual_fp,
    #         mapping_fp=self.fasting_mapping_fp,
    #         output_dir=self.wf_out,
    #         command_handler=call_commands_serially,
    #         params=self.run_core_qiime_analyses_params1,
    #         qiime_config=self.qiime_config,
    #         categories='BarcodeSequence',
    #         sampling_depth=None,
    #         arare_min_rare_depth=10,
    #         arare_num_steps=10,
    #         reference_tree_fp=None,
    #         parallel=False,
    #         status_update_callback=no_status_updates)
    #     
    #     # Basic sanity test of OTU table as details are tested 
    #     # in the pick_otus_through_otu_table tests
    #     otu_table_fp = join(self.wf_out,'otus','otu_table.biom')
    #     otu_table = parse_biom_table(open(otu_table_fp,'U'))
    #     sample_ids = list(otu_table.SampleIds)
    #     expected_sample_ids = ['PC.354','PC.355','PC.356','PC.481',
    #                            'PC.593','PC.607','PC.634','PC.635','PC.636']
    #     sample_ids.sort()
    #     expected_sample_ids.sort()
    #     self.assertEqual(sample_ids,expected_sample_ids)
    #     # even sampling directory exists (different depth may be chosen on 
    #     # different systems due to rounding)
    #     self.assertEqual(len(glob('%s/bdiv_even*' % self.wf_out)),1)
    # 
    # def test_run_core_qiime_analyses_parallel(self):
    #     """run_core_qiime_analyses: functions as expected in parallel
    #     """
    #     # this takes a long time, so use a longer sigalrm
    #     restart_timeout(600)
    #     run_core_qiime_analyses(
    #         fna_fps=self.fasting_fna_fp,
    #         qual_fps=self.fasting_qual_fp,
    #         mapping_fp=self.fasting_mapping_fp,
    #         output_dir=self.wf_out,
    #         command_handler=call_commands_serially,
    #         params=self.run_core_qiime_analyses_params1,
    #         qiime_config=self.qiime_config,
    #         categories='BarcodeSequence',
    #         sampling_depth=100,
    #         arare_min_rare_depth=10,
    #         arare_num_steps=10,
    #         reference_tree_fp=None,
    #         parallel=True,
    #         status_update_callback=no_status_updates)
    #     
    #     # Basic sanity test of OTU table as details are tested 
    #     # in the pick_otus_through_otu_table tests
    #     otu_table_fp = join(self.wf_out,'otus','otu_table.biom')
    #     otu_table = parse_biom_table(open(otu_table_fp,'U'))
    #     sample_ids= list(otu_table.SampleIds)
    #     expected_sample_ids = ['PC.354','PC.355','PC.356','PC.481',
    #                            'PC.593','PC.607','PC.634','PC.635','PC.636']
    #     sample_ids.sort()
    #     expected_sample_ids.sort()
    #     self.assertEqual(sample_ids,expected_sample_ids)


if __name__ == "__main__":
    main()