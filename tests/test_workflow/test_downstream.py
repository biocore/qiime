#!/usr/bin/env python
# File created on 20 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import sys
from StringIO import StringIO
from shutil import rmtree
from glob import glob
from os.path import join, exists, getsize, split, splitext
from tempfile import mkdtemp

from unittest import TestCase, main
from skbio.util import remove_files
from qiime.compare_alpha_diversity import compare_alpha_diversities
from qiime.util import load_qiime_config, get_qiime_temp_dir
from qiime.parse import (parse_qiime_parameters,
                         parse_distmat_to_dict)
from qiime.test import (initiate_timeout,
                        disable_timeout,
                        get_test_data_fps)
from qiime.workflow.util import (call_commands_serially,
                                 no_status_updates)
from qiime.workflow.downstream import (run_beta_diversity_through_plots,
                                       run_alpha_rarefaction,
                                       run_jackknifed_beta_diversity,
                                       run_summarize_taxa_through_plots)


class DownstreamWorkflowTests(TestCase):

    def setUp(self):
        """ """
        self.test_data = get_test_data_fps()
        self.files_to_remove = []
        self.dirs_to_remove = []

        # Create example output directory
        tmp_dir = get_qiime_temp_dir()
        self.test_out = mkdtemp(dir=tmp_dir,
                                prefix='core_qiime_analyses_test_',
                                suffix='')
        self.dirs_to_remove.append(self.test_out)

        self.qiime_config = load_qiime_config()
        self.params = parse_qiime_parameters(params_f1)

        # suppress stderr during tests (one of the systems calls in the
        # workflow prints a warning, and we can't suppress that warning with
        # warnings.filterwarnings) here because it comes from within the code
        # executed through the system call. Found this trick here:
        # http://stackoverflow.com/questions/9949633/suppressing-print-as-stdout-python
        self.saved_stderr = sys.stderr
        sys.stderr = StringIO()

        initiate_timeout(180)

    def tearDown(self):
        """ """
        disable_timeout()

        # reset sys.stderr
        sys.stderr = self.saved_stderr

        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_run_beta_diversity_through_plots(self):
        """ run_beta_diversity_through_plots generates expected results
        """
        run_beta_diversity_through_plots(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            tree_fp=self.test_data['tree'][0],
            parallel=False,
            status_update_callback=no_status_updates)

        unweighted_unifrac_dm_fp = join(
            self.test_out,
            'unweighted_unifrac_dm.txt')
        weighted_unifrac_dm_fp = join(self.test_out, 'weighted_unifrac_dm.txt')
        unweighted_unifrac_pc_fp = join(
            self.test_out,
            'unweighted_unifrac_pc.txt')
        weighted_unifrac_pc_fp = join(self.test_out, 'weighted_unifrac_pc.txt')
        weighted_unifrac_html_fp = join(self.test_out,
                                        'weighted_unifrac_emperor_pcoa_plot', 'index.html')

        # check for expected relations between values in the unweighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(unweighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)
        # check for expected relations between values in the weighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(weighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)

        # check that final output files have non-zero size
        self.assertTrue(getsize(unweighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_html_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_beta_diversity_through_plots_even_sampling(self):
        """ run_beta_diversity_through_plots functions with even sampling
        """

        run_beta_diversity_through_plots(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            sampling_depth=20,
            tree_fp=self.test_data['tree'][0],
            parallel=False,
            status_update_callback=no_status_updates)

        unweighted_unifrac_dm_fp = join(
            self.test_out,
            'unweighted_unifrac_dm.txt')
        weighted_unifrac_dm_fp = join(self.test_out, 'weighted_unifrac_dm.txt')
        unweighted_unifrac_pc_fp = join(
            self.test_out,
            'unweighted_unifrac_pc.txt')
        weighted_unifrac_pc_fp = join(self.test_out, 'weighted_unifrac_pc.txt')
        weighted_unifrac_html_fp = join(self.test_out,
                                        'weighted_unifrac_emperor_pcoa_plot', 'index.html')

        # check for expected relations between values in the unweighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(unweighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)
        # check for expected relations between values in the weighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(weighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)

        # check that final output files have non-zero size
        self.assertTrue(getsize(unweighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_html_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_beta_diversity_through_plots_parallel(self):
        """ run_beta_diversity_through_plots generates expected results in parallel
        """
        run_beta_diversity_through_plots(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            tree_fp=self.test_data['tree'][0],
            parallel=True,
            status_update_callback=no_status_updates)

        unweighted_unifrac_dm_fp = join(
            self.test_out,
            'unweighted_unifrac_dm.txt')
        weighted_unifrac_dm_fp = join(self.test_out, 'weighted_unifrac_dm.txt')
        unweighted_unifrac_pc_fp = join(
            self.test_out,
            'unweighted_unifrac_pc.txt')
        weighted_unifrac_pc_fp = join(self.test_out, 'weighted_unifrac_pc.txt')
        weighted_unifrac_html_fp = join(self.test_out,
                                        'weighted_unifrac_emperor_pcoa_plot', 'index.html')

        # check for expected relations between values in the unweighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(unweighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)
        # check for expected relations between values in the weighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(weighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)

        # check that final output files have non-zero size
        self.assertTrue(getsize(unweighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_html_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_alpha_rarefaction(self):
        """ run_alpha_rarefaction generates expected results """

        run_alpha_rarefaction(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            tree_fp=self.test_data['tree'][0],
            num_steps=5,
            parallel=False,
            min_rare_depth=3,
            max_rare_depth=18,
            status_update_callback=no_status_updates)

        html_fp = join(self.test_out, 'alpha_rarefaction_plots',
                       'rarefaction_plots.html')
        pd_averages_fp = join(self.test_out, 'alpha_rarefaction_plots',
                              'average_tables', 'PD_whole_treeSampleType.txt')
        pd_collated_fp = join(self.test_out, 'alpha_div_collated',
                              'PD_whole_tree.txt')

        # Confirm that palm and gut alpha diversities are different,
        # and suggestive of statistical significance (we only have a
        # few sequences, so we don't get significant results)
        ttest_res, alpha_avg = compare_alpha_diversities(open(pd_collated_fp),
                                                         open(
                                                             self.test_data[
                                                                 'map'][0]),
                                                         'SampleType',
                                                         18,
                                                         test_type='parametric')
        feces_palm_t = ttest_res[('feces', 'L_palm')][0]
        self.assertTrue(feces_palm_t < 0,
                        "t-statistic too high: %1.3f, but should be less than 0"
                        % feces_palm_t)

        # check that final output files have non-zero size
        self.assertTrue(getsize(html_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_alpha_rarefaction_stderr_and_stddev(self):
        """ run_alpha_rarefaction generates expected results """

        run_alpha_rarefaction(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            tree_fp=self.test_data['tree'][0],
            num_steps=5,
            parallel=False,
            min_rare_depth=3,
            max_rare_depth=18,
            status_update_callback=no_status_updates,
            plot_stderr_and_stddev=True)

        html_fp_stderr = join(self.test_out, 'alpha_rarefaction_plots_stderr',
                              'rarefaction_plots.html')
        pd_averages_fp_stderr = join(
            self.test_out, 'alpha_rarefaction_plots_stderr',
            'average_tables', 'PD_whole_treeSampleType.txt')
        html_fp_stddev = join(self.test_out, 'alpha_rarefaction_plots_stddev',
                              'rarefaction_plots.html')
        pd_averages_fp_stddev = join(
            self.test_out, 'alpha_rarefaction_plots_stddev',
            'average_tables', 'PD_whole_treeSampleType.txt')
        pd_collated_fp = join(self.test_out, 'alpha_div_collated',
                              'PD_whole_tree.txt')

        # Confirm that palm and gut alpha diversities are different,
        # and suggestive of statistical significance (we only have a
        # few sequences, so we don't get significant results)
        ttest_res, alpha_avg = compare_alpha_diversities(open(pd_collated_fp),
                                                         open(
                                                             self.test_data[
                                                                 'map'][0]),
                                                         'SampleType',
                                                         18,
                                                         test_type='parametric')
        feces_palm_t = ttest_res[('feces', 'L_palm')][0]
        self.assertTrue(feces_palm_t < 0,
                        "t-statistic too high: %1.3f, but should be less than 0"
                        % feces_palm_t)

        # check that final output files have non-zero size
        self.assertTrue(getsize(html_fp_stderr) > 0)
        self.assertTrue(getsize(html_fp_stddev) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_alpha_rarefaction_parallel(self):
        """ run_alpha_rarefaction generates expected results when run in parallel
        """

        run_alpha_rarefaction(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            tree_fp=self.test_data['tree'][0],
            num_steps=5,
            parallel=True,
            min_rare_depth=3,
            max_rare_depth=18,
            status_update_callback=no_status_updates)

        html_fp = join(self.test_out, 'alpha_rarefaction_plots',
                       'rarefaction_plots.html')
        pd_averages_fp = join(self.test_out, 'alpha_rarefaction_plots',
                              'average_tables', 'PD_whole_treeSampleType.txt')
        pd_collated_fp = join(self.test_out, 'alpha_div_collated',
                              'PD_whole_tree.txt')

        # Confirm that palm and gut alpha diversities are different,
        # and suggestive of statistical significance (we only have a
        # few sequences, so we don't get significant results)
        ttest_res, alpha_avg = compare_alpha_diversities(open(pd_collated_fp),
                                                         open(
                                                             self.test_data[
                                                                 'map'][0]),
                                                         'SampleType',
                                                         18,
                                                         test_type='parametric')
        feces_palm_t = ttest_res[('feces', 'L_palm')][0]
        self.assertTrue(feces_palm_t < 0,
                        "t-statistic too high: %1.3f, but should be less than 0"
                        % feces_palm_t)

        # check that final output files have non-zero size
        self.assertTrue(getsize(html_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_jackknifed_beta_diversity(self):
        """ run_jackknifed_beta_diversity generates expected results """

        run_jackknifed_beta_diversity(
            self.test_data['biom'][0],
            self.test_data['tree'][0],
            20,
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            self.test_data['map'][0],
            parallel=False,
            status_update_callback=no_status_updates)

        weighted_unifrac_upgma_tree_fp = join(self.test_out,
                                              'weighted_unifrac',
                                              'upgma_cmp', 'jackknife_named_nodes.tre')
        unweighted_unifrac_upgma_tree_fp = join(
            self.test_out, 'unweighted_unifrac', 'upgma_cmp',
            'jackknife_named_nodes.tre')
        weighted_unifrac_emperor_index_fp = join(
            self.test_out, 'weighted_unifrac', 'emperor_pcoa_plots',
            'index.html')
        unweighted_unifrac_emperor_index_fp = join(
            self.test_out, 'unweighted_unifrac', 'emperor_pcoa_plots',
            'index.html')

        input_file_basename = splitext(split(self.test_data['biom'][0])[1])[0]
        unweighted_unifrac_dm_fp = join(self.test_out, 'unrarefied_bdiv',
                                        'unweighted_unifrac_%s.txt' % input_file_basename)
        weighted_unifrac_dm_fp = join(self.test_out, 'unrarefied_bdiv',
                                      'weighted_unifrac_%s.txt' % input_file_basename)

       # check for expected relations between values in the unweighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(unweighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)
        # check for expected relations between values in the weighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(weighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)

        # check that final output files have non-zero size
        self.assertTrue(getsize(weighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_emperor_index_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_emperor_index_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_jackknifed_beta_diversity_parallel(self):
        """ run_jackknifed_beta_diversity generates expected results """

        run_jackknifed_beta_diversity(
            self.test_data['biom'][0],
            self.test_data['tree'][0],
            20,
            self.test_out,
            call_commands_serially,
            self.params,
            self.qiime_config,
            self.test_data['map'][0],
            parallel=True,
            status_update_callback=no_status_updates)

        weighted_unifrac_upgma_tree_fp = join(self.test_out,
                                              'weighted_unifrac',
                                              'upgma_cmp', 'jackknife_named_nodes.tre')
        unweighted_unifrac_upgma_tree_fp = join(
            self.test_out, 'unweighted_unifrac', 'upgma_cmp',
            'jackknife_named_nodes.tre')
        weighted_unifrac_emperor_index_fp = join(
            self.test_out, 'weighted_unifrac', 'emperor_pcoa_plots',
            'index.html')
        unweighted_unifrac_emperor_index_fp = join(
            self.test_out, 'unweighted_unifrac', 'emperor_pcoa_plots',
            'index.html')

        input_file_basename = splitext(split(self.test_data['biom'][0])[1])[0]
        unweighted_unifrac_dm_fp = join(self.test_out, 'unrarefied_bdiv',
                                        'unweighted_unifrac_%s.txt' % input_file_basename)
        weighted_unifrac_dm_fp = join(self.test_out, 'unrarefied_bdiv',
                                      'weighted_unifrac_%s.txt' % input_file_basename)

       # check for expected relations between values in the unweighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(unweighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)
        # check for expected relations between values in the weighted unifrac
        # distance matrix
        dm = parse_distmat_to_dict(open(weighted_unifrac_dm_fp))
        self.assertTrue(dm['f1']['f2'] < dm['f1']['p1'],
                        "Distance between pair of fecal samples is larger than distance"
                        " between fecal and palm sample (unweighted unifrac).")
        self.assertEqual(dm['f1']['f1'], 0)

        # check that final output files have non-zero size
        self.assertTrue(getsize(weighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_emperor_index_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_emperor_index_fp) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_summarize_taxa_through_plots(self):
        """ run_summarize_taxa_through_plots generates expected results
        """
        run_summarize_taxa_through_plots(
            self.test_data['biom'][0],
            self.test_data['map'][0],
            self.test_out,
            mapping_cat=None,
            sort=False,
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            status_update_callback=no_status_updates)

        # Check that summarized taxonomy files have non-zero size
        input_file_basename = splitext(split(self.test_data['biom'][0])[1])[0]
        for i in [2, 3, 4, 5, 6]:
            sum_taxa_file = join(self.test_out, input_file_basename + '_L%s.txt'
                                 % (str(i)))
            self.assertTrue(getsize(sum_taxa_file) > 0)

        # Check the html files are generated
        self.assertTrue(getsize(join(self.test_out, 'taxa_summary_plots',
                                     'area_charts.html')) > 0)

        self.assertTrue(getsize(join(self.test_out, 'taxa_summary_plots',
                                     'area_charts.html')) > 0)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out, 'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

params_f1 = """
multiple_rarefactions:num_reps	1
multiple_rarefactions_even_depth:num_reps	5
""".split('\n')

if __name__ == "__main__":
    main()
