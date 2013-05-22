#!/usr/bin/env python
# File created on 09 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import sys
from StringIO import StringIO
from shutil import rmtree
from os.path import exists
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir,
                        get_tmp_filename,
                        load_qiime_config)
from qiime.parse import (parse_qiime_parameters)
from qiime.test import (initiate_timeout,
                        disable_timeout,
                        get_test_data_fps)
from qiime.workflow.util import (no_status_updates)
from qiime.workflow.core_diversity_analyses import run_core_diversity_analyses

class CoreDiversityAnalysesTests(TestCase):
    
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
        self.qiime_config['jobs_to_start'] = 2
        self.qiime_config['seconds_to_sleep'] = 1
        
        # suppress stderr during tests (one of the systems calls in the 
        # workflow prints a warning, and we can't suppress that warning with 
        # warnings.filterwarnings) here because it comes from within the code 
        # executed through the system call. Found this trick here:
        # http://stackoverflow.com/questions/9949633/suppressing-print-as-stdout-python
        self.saved_stderr = sys.stderr
        sys.stderr = StringIO()
        
        # Define number of seconds a test can run for before timing out 
        # and failing
        initiate_timeout(420)

    
    def tearDown(self):
        
        disable_timeout()
        
        # reset sys.stderr
        sys.stderr = self.saved_stderr
        
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_run_core_diversity_analyses_error_on_invalid_category(self):
        """run_core_diversity_analyses: error raised on invalid categories
        """
        # too few values in 'month' category
        self.assertRaises(ValueError,
                          run_core_diversity_analyses,
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
                          run_core_diversity_analyses,
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

    def test_run_core_diversity_analyses(self):
        """run_core_diversity_analyses functions with categories
        """
        run_core_diversity_analyses(
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
        
        # Basic sanity test that output directories and files exist
        fps = [
         '%s/bdiv_even20' % self.test_out,
         '%s/arare_max20' % self.test_out,
         '%s/taxa_plots' % self.test_out,
         '%s/bdiv_even20/unweighted_unifrac_dm.txt' % self.test_out,
         '%s/bdiv_even20/weighted_unifrac_pc.txt' % self.test_out,
         '%s/arare_max20/days_since_epoch_chao1.txt' % self.test_out,
         '%s/arare_max20/SampleType_PD_whole_tree.txt' % self.test_out,
         '%s/index.html' % self.test_out,
         '%s/table_mc%d.biom.gz' % (self.test_out,20)
        ]
        for fp in fps:
            self.assertTrue(exists(fp))

    def test_run_core_diversity_analyses_parallel(self):
        """run_core_diversity_analyses functions with categories in parallel
        """
        run_core_diversity_analyses(
                        self.test_data['biom'][0],
                        self.test_data['map'][0],
                        20,
                        output_dir=self.test_out,
                        params=parse_qiime_parameters({}),
                        qiime_config=self.qiime_config,
                        categories=['SampleType','days_since_epoch'],
                        tree_fp=self.test_data['tree'][0],
                        parallel=True,
                        status_update_callback=no_status_updates)
        
        # Basic sanity test that output directories and files exist
        fps = [
         '%s/bdiv_even20' % self.test_out,
         '%s/arare_max20' % self.test_out,
         '%s/taxa_plots' % self.test_out,
         '%s/bdiv_even20/unweighted_unifrac_dm.txt' % self.test_out,
         '%s/bdiv_even20/weighted_unifrac_pc.txt' % self.test_out,
         '%s/arare_max20/days_since_epoch_chao1.txt' % self.test_out,
         '%s/arare_max20/SampleType_PD_whole_tree.txt' % self.test_out,
         '%s/index.html' % self.test_out,
         '%s/table_mc%d.biom.gz' % (self.test_out,20)
        ]
        for fp in fps:
            self.assertTrue(exists(fp))

    def test_run_core_diversity_analyses_no_categories(self):
        """run_core_diversity_analyses functions without categories
        """
        # this takes a long time, so use a longer sigalrm
        run_core_diversity_analyses(
                        self.test_data['biom'][0],
                        self.test_data['map'][0],
                        20,
                        output_dir=self.test_out,
                        params=parse_qiime_parameters({}),
                        qiime_config=self.qiime_config,
                        categories=None,
                        tree_fp=self.test_data['tree'][0],
                        parallel=False,
                        status_update_callback=no_status_updates)
        
        # Basic sanity test that output directories and files exist
        fps = [
         '%s/bdiv_even20' % self.test_out,
         '%s/arare_max20' % self.test_out,
         '%s/taxa_plots' % self.test_out,
         '%s/bdiv_even20/unweighted_unifrac_dm.txt' % self.test_out,
         '%s/bdiv_even20/weighted_unifrac_pc.txt' % self.test_out,
         '%s/index.html' % self.test_out,
         '%s/table_mc%d.biom.gz' % (self.test_out,20)
        ]
        for fp in fps:
            self.assertTrue(exists(fp))
            
        # categorical output files don't exist
        self.assertFalse(exists(
         '%s/arare_max20/days_since_epoch_chao1.txt' % self.test_out))
        self.assertFalse(exists(
         '%s/arare_max20/SampleType_PD_whole_tree.txt' % self.test_out))

    def test_run_core_diversity_analyses_no_tree(self):
        """run_core_diversity_analyses functions without tree
        """
        # this takes a long time, so use a longer sigalrm
        run_core_diversity_analyses(
         self.test_data['biom'][0],
         self.test_data['map'][0],
         20,
         output_dir=self.test_out,
         params=parse_qiime_parameters(
          ['beta_diversity:metrics bray_curtis',
           'alpha_diversity:metrics observed_species,chao1']),
         qiime_config=self.qiime_config,
         categories=['SampleType'],
         tree_fp=None,
         parallel=False,
         status_update_callback=no_status_updates)
        
        # Basic sanity test that output directories and files exist
        fps = [
         '%s/bdiv_even20' % self.test_out,
         '%s/arare_max20' % self.test_out,
         '%s/taxa_plots' % self.test_out,
         '%s/bdiv_even20/bray_curtis_dm.txt' % self.test_out,
         '%s/arare_max20/SampleType_observed_species.txt' % self.test_out,
         '%s/index.html' % self.test_out,
         '%s/table_mc%d.biom.gz' % (self.test_out,20)
        ]
        for fp in fps:
            self.assertTrue(exists(fp))
        
        # phylogenetic diversity output files do not exist
        self.assertFalse(exists(
         '%s/bdiv_even20/unweighted_unifrac_dm.txt' % self.test_out))

if __name__ == "__main__":
    main()
