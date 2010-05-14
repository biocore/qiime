#!/usr/bin/env python
# File created on 30 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

import signal
import os
from shutil import rmtree
from glob import glob
from os.path import join, exists, getsize, split, splitext
from os import makedirs
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename, ApplicationNotFoundError
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow import (run_qiime_data_preparation,
    run_beta_diversity_through_3d_plot,
    run_qiime_alpha_rarefaction,
    run_jackknifed_upgma_clustering,
    run_process_sra_submission,
    call_commands_serially,
    no_status_updates,WorkflowError,print_commands)

## The test case timing code included in this file is adapted from
## recipes provided at:
##  http://code.activestate.com/recipes/534115-function-timeout/
##  http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call
class TimeExceededError(Exception):
    pass


allowed_seconds_per_test = 240

def timeout(signum, frame):
    raise TimeExceededError,\
     "Test failed to run in allowed time (%d seconds)."\
      % allowed_seconds_per_test
    
class WorkflowTests(TestCase):
    
    def setUp(self):
        """ """
        self.qiime_config = load_qiime_config()
        self.dirs_to_remove = []
        self.files_to_remove = []
        
        tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'
        if not exists(tmp_dir):
            makedirs(tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(tmp_dir)
        
        self.wf_out = get_tmp_filename(tmp_dir=tmp_dir,
         prefix='qiime_wf_out',suffix='',result_constructor=str)
        self.dirs_to_remove.append(self.wf_out)
        
        self.fasting_mapping_fp = get_tmp_filename(tmp_dir=tmp_dir,
         prefix='qiime_wf_mapping',suffix='.txt')
        fasting_mapping_f = open(self.fasting_mapping_fp,'w')
        fasting_mapping_f.write(fasting_map)
        fasting_mapping_f.close()
        self.files_to_remove.append(self.fasting_mapping_fp)
        
        self.fasting_seqs_fp = get_tmp_filename(tmp_dir=tmp_dir,
            prefix='qiime_wf_seqs',suffix='.fasta')
        fasting_seqs_f = open(self.fasting_seqs_fp,'w')
        fasting_seqs_f.write(fasting_seqs_subset)
        fasting_seqs_f.close()
        self.files_to_remove.append(self.fasting_seqs_fp)
        
        self.fasting_seqs_denoiser_fp = get_tmp_filename(tmp_dir=tmp_dir,
            prefix='qiime_wf_seqs',suffix='.fasta')
        fasting_seqs_f = open(self.fasting_seqs_denoiser_fp,'w')
        fasting_seqs_f.write('\n'.join(fasting_seqs_subset.split('\n')[:44]))
        fasting_seqs_f.close()
        self.files_to_remove.append(self.fasting_seqs_denoiser_fp)
        
        self.fasting_otu_table_fp = get_tmp_filename(tmp_dir=tmp_dir,
            prefix='qiime_wf_otu_table',suffix='.txt')
        fasting_otu_table_f = open(self.fasting_otu_table_fp,'w')
        fasting_otu_table_f.write(fasting_subset_otu_table)
        fasting_otu_table_f.close()
        self.files_to_remove.append(self.fasting_otu_table_fp)
        
        self.fasting_tree_fp = get_tmp_filename(tmp_dir=tmp_dir,
            prefix='qiime_wf_tree',suffix='.tre')
        fasting_tree_f = open(self.fasting_tree_fp,'w')
        fasting_tree_f.write(fasting_subset_tree)
        fasting_tree_f.close()
        self.files_to_remove.append(self.fasting_tree_fp)
        
        self.template_aln_fp = get_tmp_filename(tmp_dir=tmp_dir,
         prefix='wf_template',suffix='.fasta')
        template_aln_f = open(self.template_aln_fp,'w')
        template_aln_f.write(template_alignment_subset)
        template_aln_f.close()
        self.files_to_remove.append(self.template_aln_fp)
        
        self.lanemask_fp = get_tmp_filename(tmp_dir=tmp_dir,
         prefix='wf_lanemask',suffix='.txt')
        lanemask_f = open(self.lanemask_fp,'w')
        lanemask_f.write(lanemask)
        lanemask_f.close()
        self.files_to_remove.append(self.lanemask_fp)
        
        self.sff_fp = get_tmp_filename(tmp_dir=tmp_dir,
         prefix='wf_sff',suffix='.txt')
        sff_f = open(self.sff_fp,'w')
        sff_f.write(fasting_subset_sff)
        sff_f.close()
        self.files_to_remove.append(self.sff_fp)
        
        working_dir = self.qiime_config['working_dir'] or './'
        jobs_dir = join(working_dir,'jobs')
        if not exists(jobs_dir):
            # only clean up the jobs dir if it doesn't already exist
            self.dirs_to_remove.append(jobs_dir)
        self.params = parse_qiime_parameters(qiime_parameters_f)
        self.params['align_seqs']['template_fp'] = self.template_aln_fp
        self.params['filter_alignment']['lane_mask_fp'] = self.lanemask_fp

        self.experiment_fp = get_tmp_filename(
            tmp_dir=tmp_dir, prefix='SRA_wf_exp', suffix='.txt')
        f = open(self.experiment_fp, 'w')
        f.write(sra_experiment_txt)
        f.close()
        self.files_to_remove.append(self.experiment_fp)

        self.submission_fp = get_tmp_filename(tmp_dir=tmp_dir,
            prefix='SRA_wf_sub',suffix='.txt')
        f = open(self.submission_fp, 'w')
        f.write(sra_submission_txt)
        f.close()
        self.files_to_remove.append(self.submission_fp)
        
        self.sra_params = parse_qiime_parameters(sra_submission_params_f)
        
        signal.signal(signal.SIGALRM, timeout)
        # set the 'alarm' to go off in allowed_seconds seconds
        signal.alarm(allowed_seconds_per_test)
        
    
    def tearDown(self):
        """ """
        # turn off the alarm
        signal.alarm(0)
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)
        
    def test_unsupported_options_handled_nicely(self):
        """WorkflowError raised on unsupported option """
        self.params['beta_diversity']['blah'] = self.fasting_otu_table_fp
        self.assertRaises(WorkflowError,run_beta_diversity_through_3d_plot,
         self.fasting_otu_table_fp, 
         self.fasting_mapping_fp,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         tree_fp=self.fasting_tree_fp,
         parallel=False, 
         status_update_callback=no_status_updates)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_qiime_data_preparation(self):
        """run_qiime_data_preparation runs without error"""
        run_qiime_data_preparation(
         self.fasting_seqs_fp, 
         self.wf_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.fasting_seqs_fp)[1])[0]
        otu_table_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'rdp_assigned_taxonomy','otu_table','%s_otu_table.txt' % 
         input_file_basename)
        tree_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'pynast_aligned_seqs','fasttree_phylogeny','%s_rep_set.tre' % 
         input_file_basename)
         
        # check that the two final output files have non-zero size
        self.assertTrue(getsize(tree_fp) > 0)
        self.assertTrue(getsize(otu_table_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)



    def test_run_qiime_data_preparation_denoise(self):
        """run_qiime_data_preparation denoises"""
        try:
            run_qiime_data_preparation(
             self.fasting_seqs_denoiser_fp, 
             self.wf_out, 
             call_commands_serially,
             self.params, 
             self.qiime_config, 
             self.sff_fp,
             self.fasting_mapping_fp,
             parallel=False,
             status_update_callback=no_status_updates)
        except WorkflowError:
            raise ApplicationNotFoundError,\
            "Denoiser (or other dependency) cannot be found."
         
        input_file_basename = 'denoised_seqs'
        
        otu_table_fp = join(self.wf_out,'uclust_picked_otus','denoised_otus',
         'rep_set','rdp_assigned_taxonomy','otu_table','%s_otu_table.txt' % 
         input_file_basename)
        tree_fp = join(self.wf_out,'uclust_picked_otus','denoised_otus',
         'rep_set','pynast_aligned_seqs','fasttree_phylogeny',
         '%s_rep_set.tre' % input_file_basename)
         
        # check that the two final output files have non-zero size
        self.assertTrue(getsize(tree_fp) > 0)
        self.assertTrue(getsize(otu_table_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_qiime_data_preparation_muscle(self):
        """run_qiime_data_preparation runs without error with muscle aligner"""
        self.params['align_seqs']['alignment_method'] = 'muscle'
        run_qiime_data_preparation(
         self.fasting_seqs_fp, 
         self.wf_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.fasting_seqs_fp)[1])[0]
        otu_table_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'rdp_assigned_taxonomy','otu_table','%s_otu_table.txt' % 
         input_file_basename)
        tree_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'muscle_aligned_seqs','fasttree_phylogeny','%s_rep_set.tre' % 
         input_file_basename)
         
        # check that the two final output files have non-zero size
        self.assertTrue(getsize(tree_fp) > 0)
        self.assertTrue(getsize(otu_table_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
    
    def test_run_qiime_data_preparation_parallel(self):
        """run_qiime_data_preparation runs in parallel without error"""
        run_qiime_data_preparation(
         self.fasting_seqs_fp, 
         self.wf_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=True,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.fasting_seqs_fp)[1])[0]
        otu_table_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'rdp_assigned_taxonomy','otu_table','%s_otu_table.txt' % 
         input_file_basename)
        tree_fp = join(self.wf_out,'uclust_picked_otus','rep_set',
         'pynast_aligned_seqs','fasttree_phylogeny','%s_rep_set.tre' % 
         input_file_basename)
         
        # check that the two final output files have non-zero size
        self.assertTrue(getsize(tree_fp) > 0)
        self.assertTrue(getsize(otu_table_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
         
         
    def test_run_beta_diversity_through_3d_plot(self):
        """ run_beta_diversity_through_3d_plot runs without error """
        run_beta_diversity_through_3d_plot(
         self.fasting_otu_table_fp, 
         self.fasting_mapping_fp,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         tree_fp=self.fasting_tree_fp,
         parallel=False, 
         status_update_callback=no_status_updates)
         
        unweighted_unifrac_pc_fp = join(self.wf_out,'unweighted_unifrac_pc.txt')
        weighted_unifrac_pc_fp = join(self.wf_out,'weighted_unifrac_pc.txt')
        weighted_unifrac_html_fp = join(self.wf_out,
        'weighted_unifrac_3d_continuous','weighted_unifrac_pc.txt_3D.html')
        
        # check that final output files have non-zero size
        self.assertTrue(getsize(unweighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_html_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
      
    def test_run_beta_diversity_through_3d_plot_parallel(self):
        """ run_beta_diversity_through_3d_plot runs in parallel without error """
        run_beta_diversity_through_3d_plot(
         self.fasting_otu_table_fp, 
         self.fasting_mapping_fp,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         tree_fp=self.fasting_tree_fp,
         parallel=True, 
         status_update_callback=no_status_updates)
         
        unweighted_unifrac_pc_fp = join(self.wf_out,'unweighted_unifrac_pc.txt')
        weighted_unifrac_pc_fp = join(self.wf_out,'weighted_unifrac_pc.txt')
        weighted_unifrac_html_fp = join(self.wf_out,
        'weighted_unifrac_3d_continuous','weighted_unifrac_pc.txt_3D.html')
        
        # check that final output files have non-zero size
        self.assertTrue(getsize(unweighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_pc_fp) > 0)
        self.assertTrue(getsize(weighted_unifrac_html_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_qiime_alpha_rarefaction(self):
        """ run_qiime_alpha_rarefaction runs without error """
    
        run_qiime_alpha_rarefaction(
         self.fasting_otu_table_fp, 
         self.fasting_mapping_fp,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         tree_fp=self.fasting_tree_fp,
         num_steps=10, 
         parallel=False, 
         min_seqs_per_sample=10,\
         status_update_callback=no_status_updates)
         
        pd_control_plot_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'html_plots','PD_whole_treeTreatmentControl_ave.png')
        pd_treatment_plot_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'average_plots','PD_whole_treeTreatment.png')
        pd_averages_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'average_tables','PD_whole_treeTreatment.txt')
        
        # check that final output files have non-zero size
        self.assertTrue(getsize(pd_control_plot_fp) > 0)
        self.assertTrue(getsize(pd_treatment_plot_fp) > 0)
        self.assertTrue(getsize(pd_averages_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_qiime_alpha_rarefaction_parallel(self):
        """ run_qiime_alpha_rarefaction runs in parallel without error """
    
        run_qiime_alpha_rarefaction(
         self.fasting_otu_table_fp, 
         self.fasting_mapping_fp,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         tree_fp=self.fasting_tree_fp,
         num_steps=10, 
         parallel=True, 
         min_seqs_per_sample=10,\
         status_update_callback=no_status_updates)
         
        pd_control_plot_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'html_plots','PD_whole_treeTreatmentControl_ave.png')
        pd_treatment_plot_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'average_plots','PD_whole_treeTreatment.png')
        pd_averages_fp = join(self.wf_out,'alpha_rarefaction_plots',
         'average_tables','PD_whole_treeTreatment.txt')
        
        # check that final output files have non-zero size
        self.assertTrue(getsize(pd_control_plot_fp) > 0)
        self.assertTrue(getsize(pd_treatment_plot_fp) > 0)
        self.assertTrue(getsize(pd_averages_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
         
    def test_run_jackknifed_upgma_clustering(self):
        """ run_jackknifed_upgma_clustering runs without error """
    
        run_jackknifed_upgma_clustering(
         self.fasting_otu_table_fp,
         self.fasting_tree_fp,
         100,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         parallel=False,
         status_update_callback=no_status_updates)
         
        weighted_unifrac_upgma_tree_fp = join(self.wf_out,
         'weighted_unifrac',
         'upgma_cmp','jackknife_named_nodes.tre')
        unweighted_unifrac_upgma_tree_fp = join(
         self.wf_out,'unweighted_unifrac','upgma_cmp',
         'jackknife_named_nodes.tre')
         
        # check that final output files have non-zero size
        self.assertTrue(getsize(weighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_upgma_tree_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_jackknifed_upgma_clustering_parallel(self):
        """ run_jackknifed_upgma_clustering runs in parallel without error """
    
        run_jackknifed_upgma_clustering(
         self.fasting_otu_table_fp,
         self.fasting_tree_fp,
         100,
         self.wf_out, 
         call_commands_serially,
         self.params,
         self.qiime_config,
         parallel=True,
         status_update_callback=no_status_updates)
         
        weighted_unifrac_upgma_tree_fp = join(self.wf_out,
         'weighted_unifrac',
         'upgma_cmp','jackknife_named_nodes.tre')
        unweighted_unifrac_upgma_tree_fp = join(
         self.wf_out,'unweighted_unifrac','upgma_cmp',
         'jackknife_named_nodes.tre')
         
        # check that final output files have non-zero size
        self.assertTrue(getsize(weighted_unifrac_upgma_tree_fp) > 0)
        self.assertTrue(getsize(unweighted_unifrac_upgma_tree_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.wf_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_process_sra_submission(self):
        """run_process_sra_submission produces a valid SRA submission package"""
        # TODO: remove dependence on external files, if possible
        test_dir = os.path.dirname(os.path.abspath(__file__))
        sff_dir = os.path.join(test_dir, 'sra_test_files', 'F6AVWTA')

        run_process_sra_submission(
                input_experiment_fp=self.experiment_fp,
                input_submission_fp=self.submission_fp,
                sff_dir=sff_dir,
                refseqs_fp=self.template_aln_fp,
                output_dir=self.wf_out,
                params=self.sra_params,
                qiime_config=self.qiime_config,
                command_handler=call_commands_serially,
                status_update_callback=no_status_updates,
                )

        tar_fp = os.path.join(self.wf_out, 'my_sffs.tgz')
        self.assertTrue(getsize(tar_fp) > 0)
        
        experiment_basename = os.path.splitext(
            os.path.basename(self.experiment_fp))[0]
        experiment_xml_fp = os.path.join(
            self.wf_out, experiment_basename + '.xml')
        self.assertTrue(getsize(experiment_xml_fp) > 0)
        run_xml_fp = os.path.join(
            self.wf_out, experiment_basename + '.xml')
        self.assertTrue(getsize(run_xml_fp) > 0)
        
sra_submission_params_f = """# split_libraries parameters
split_libraries:min-qual-score	5
split_libraries:min-seq-length	30
split_libraries:max-seq-length	1000
split_libraries:barcode-type	10
# by default, barcode-type should be 12 -- in this example it is 10
#split_libraries:barcode-type	12
split_libraries:max-homopolymer	1000
split_libraries:max-primer-mismatch	100
split_libraries:max-ambig	1000

# pick_otus parameters
pick_otus:otu_picking_method	cdhit
pick_otus:max_cdhit_memory	4000
pick_otus:prefix_prefilter_length	100
pick_otus:similarity	0.95

exclude_seqs_by_blast:word_size	10
exclude_seqs_by_blast:percent_aligned	0.25
exclude_seqs_by_blast:e_value	1e-20""".split('\n')

sra_experiment_txt = '''#EXPERIMENT_ALIAS	EXPERIMENT_CENTER	EXPERIMENT_TITLE	STUDY_REF	STUDY_CENTER	EXPERIMENT_DESIGN_DESCRIPTION	LIBRARY_CONSTRUCTION_PROTOCOL	SAMPLE_ALIAS	SAMPLE_CENTER	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	POOL_PROPORTION	BARCODE_READ_GROUP_TAG	BARCODE	LINKER	PRIMER_READ_GROUP_TAG	KEY_SEQ	PRIMER	RUN_PREFIX	RUN_ALIAS	REGION	PLATFORM	RUN_CENTER	RUN_DATE	INSTRUMENT_NAME
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015438	NCBI	F6AVWTA01_2878_700015438_V1-V3	B-2004-03-S1.sff	0.014492754	F6AVWTA01_ATGTTCGATG	AGACTCTGCT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	F6AVWTA01_2878	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015438	NCBI	F6AVWTA02_2878_700015438_V1-V3	B-2008-05-S1.sff	0.014492754	F6AVWTA02_ATGTTCTAGT	ATGTTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	F6AVWTA02_2878	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015470	NCBI	F6AVWTA01_2866_700015470_V1-V3	B-2004-04-S1.sff	0.014492754	F6AVWTA01_GCTCTACGTC	GCTCTACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	F6AVWTA01_2866	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015470	NCBI	F6AVWTA02_2866_700015470_V1-V3	B-2008-08-S1.sff	0.014492754	F6AVWTA02_GCTCTGTACT	GCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	F6AVWTA02_2866	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015766	NCBI	F6AVWTA01_2898_700015766_V1-V3	B-2004-08-S1.sff	0.014492754	F6AVWTA01_CATGAGCGTC	CATGAGCGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	F6AVWTA01_2898	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015766	NCBI	F6AVWTA02_2898_700015766_V1-V3	B-2009-06-S1.sff	0.014492754	F6AVWTA02_CATGAGCGTG	CATGAGCGTG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	F6AVWTA02_2898	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015468	NCBI	F6AVWTA01_2865_700015468_V1-V3	B-2005-06-S1.sff	0.014492754	F6AVWTA01_AGTACGTACT	AGTACGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	F6AVWTA01_2865	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015468	NCBI	F6AVWTA02_2865_700015468_V1-V3	B-2011-01-S1.sff	0.014492754	F6AVWTA02_AGTACACGTC	AGTACACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	F6AVWTA02_2865	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700016371	NCBI	F6AVWTA01_2907_700016371_V1-V3	B-2006-03-S1.sff	0.014492754	F6AVWTA01_TCTCTCTAGT	TCTCTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	F6AVWTA01_2907	0	FLX	JCVI 	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700016371	NCBI	F6AVWTA02_2907_700016371_V1-V3	B-2011-02-S1.sff	0.014492754	F6AVWTA02_TCTCTGTACT	TCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	F6AVWTA02_2907	0	FLX	JCVI 	NULL	NULL
'''

sra_submission_txt = '''#Field	Value	Example	Comments
accession	SRA003492	SRA003492	"leave blank if not assigned yet, e.g. if new submission"
submission_id	fierer_hand_study	fierer_hand_study	internally unique id for the submission
center_name	CCME	CCME	name of the center preparing the submission
submission_comment	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Free-text comments regarding submission
lab_name	Knight	Knight	"name of lab preparing submission, can differ from center (usually refers to the PI\'s info, not the sequencing center\'s)"
submission_date	2009-10-22T01:23:00-05:00	2009-10-22T01:23:00-05:00	timestamp of submission
CONTACT	Rob Knight;Rob.Knight@Colorado.edu	Rob Knight;Rob.Knight@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
CONTACT	Noah Fierer;Noah.Fierer@Colorado.edu	Noah Fierer;Noah.Fierer@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
study	study.xml	fierer_hand_study.study.xml	"leave blank if not submitting study, put in filename otherwise"
sample	sample.xml	fierer_hand_study.sample.xml	"leave blank if not submitting sample, put in filename otherwise"
experiment	experiment.xml	fierer_hand_study.experiment.xml	"leave blank if not submitting experiment, put in filename otherwise"
run	run.xml	fierer_hand_study.run.xml	"leave blank if not submitting run, put in filename otherwise"
file	my_sffs.tgz	fierer_hand_study.seqs.tgz	"leave blank if not submitting sequence data, put in filename otherwise"'''



qiime_parameters_f = """# qiime_parameters.txt
# WARNING: DO NOT EDIT OR DELETE Qiime/qiime_parameters.txt. Users should copy this file and edit copies of it.

# OTU picker parameters
pick_otus:otu_picking_method	uclust
pick_otus:clustering_algorithm	furthest
pick_otus:max_cdhit_memory	400
pick_otus:refseqs_fp
pick_otus:blast_db
pick_otus:similarity	0.97
pick_otus:max_e_value	1e-10
pick_otus:prefix_prefilter_length
pick_otus:trie_prefilter
pick_otus:prefix_length
pick_otus:suffix_length
pick_otus:optimal_uclust
pick_otus:exact_uclust
pick_otus:user_sort
pick_otus:suppress_presort_by_abundance_uclust
pick_otus:suppress_new_clusters

# Parallel options
parallel:jobs_to_start	2
parallel:retain_temp_files	False
parallel:seconds_to_sleep	1

# Representative set picker parameters
pick_rep_set:rep_set_picking_method	most_abundant
pick_rep_set:sort_by	otu

# Multiple sequence alignment parameters
align_seqs:template_fp
align_seqs:alignment_method	pynast
align_seqs:pairwise_alignment_method	uclust
align_seqs:blast_db
align_seqs:min_length	150
align_seqs:min_percent_id	75.0

# Alignment filtering (prior to tree-building) parameters
filter_alignment:lane_mask_fp
filter_alignment:allowed_gap_frac	 0.999999
filter_alignment:remove_outliers	False
filter_alignment:threshold	3.0

# Taxonomy assignment parameters
assign_taxonomy:id_to_taxonomy_fp
assign_taxonomy:reference_seqs_fp
assign_taxonomy:assignment_method	rdp
assign_taxonomy:confidence	0.8

# Phylogenetic tree building parameters
make_phylogeny:tree_method	fasttree
make_phylogeny:root_method	tree_method_default

# Beta diversity parameters
beta_diversity:metrics	weighted_unifrac,unweighted_unifrac

# Make 3D plot parameters
make_3d_plots:custom_axes

# Rarefaction parameters
multiple_rarefactions:num-reps	1
multiple_rarefactions:depth
multiple_rarefactions:lineages_included	False

# Even-depth rarefaction parameters
multiple_rarefactions_even_depth:num-reps	5

# Alpha diversity parameters
alpha_diversity:metrics	PD_whole_tree

# Make rarefaction plots parameters
make_rarefaction_plots:imagetype	png
make_rarefaction_plots:resolution	75

# Collate alpha
collate_alpha:example_path

""".split('\n')

fasting_map = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse__I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse__I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	Control_mouse__I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse__I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse__I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse__I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse__I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse__I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse__I.D._636
"""

fasting_seqs_subset = """>PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.636_6 FLP3FBN01A55WZ orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.634_7 FLP3FBN01D7O1S orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAGTGAAACCATGCGGAAATACTATACTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATC
>PC.634_8 FLP3FBN01CV5ZU orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_9 FLP3FBN01DQ783 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_10 FLP3FBN01DDPFF orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGGTACGCGTTACTC
>PC.634_11 FLP3FBN01CPD70 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.593_12 FLP3FBN01BBAE6 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCG
>PC.355_13 FLP3FBN01AWYZD orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.634_14 FLP3FBN01AM0P3 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCGTACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCG
>PC.635_15 FLP3FBN01BO4IE orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTCCCATTTATATCATGTGATATTCATGGCATATGCGGTATTAGCAGTCATTTCTAACTGTTGTTCCCCTCTGATAGG
>PC.635_16 FLP3FBN01BPX14 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATCGAGGGCAAGTTGCATACGCG
>PC.356_17 FLP3FBN01DB5I7 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTGCACCGCCGGAGCTTCCCCCGCCGCCCCATGCGGGGCTGCGGGCATATGCGGTATTAGCAGCCGTTTCCGGCTGTTGTCCCCCAGTGCAGGGCAGGTTGCCCACGCGTTACTCACCCGTCCG
>PC.634_18 FLP3FBN01AK9OO orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTCGAGCCGTTACCTCACCAACTAGCTAATGCGCCGCGGGCCCATCTCATAGCGGATTACTCCTTTAATTGCTACTTCATGCGAAGCTACAATCTTATGCGGTATTAATCTCCCTTTCGGAAGGCTATTCCCCTCTATGAGGTCAGGTTG
>PC.634_19 FLP3FBN01ANGF2 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.634_20 FLP3FBN01AF994 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTATCCCCTCCTGATAGGTAAGTTGCATACGCGTTACTTCACCCGTCGCGG
>PC.634_21 FLP3FBN01AHXU8 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCGGGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAGAGAG
>PC.593_22 FLP3FBN01APRWO orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
TTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTCTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>PC.634_23 FLP3FBN01C645E orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCAGTCG
>PC.634_24 FLP3FBN01A2K0Q orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCCCACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTATGGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCAT
>PC.634_25 FLP3FBN01DJY30 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACTAGCTAATCAGACGCGGGCCCATCTCATACCACCGGAGTTTTTGCCGCTGCACCATGCGGTGCTGTGGCCTTATGCGGTATTAGCAGCCATTTCTGGCTGTTATCCCCCTGTATGAGGCAGGTTGCCCACG
>PC.355_26 FLP3FBN01EB63O orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACCGGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTGTCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>PC.634_27 FLP3FBN01DUWHJ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCGGAATCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_28 FLP3FBN01ESJ5R orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.635_29 FLP3FBN01BJH6S orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGACCTCTCAGTCCGGCTACCGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCGGACGCGAGCCCATCTCCGAGCGATATAATCTTTGATACCAAAGGGATGTCCCTCCAGTATGTTATGCGGTATTAGCGACCGTTTCCAGCCGTTATTCCCCTCTCGAAGGTAGGTTGCTCACGTGTTACTCACCCGTCCG
>PC.634_30 FLP3FBN01BF988 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCATGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTC
>PC.634_31 FLP3FBN01DRASA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_32 FLP3FBN01EMWGJ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCG
>PC.634_33 FLP3FBN01BMERW orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGCGTTACTCACCCGTGCGCCGGTCG
>PC.634_34 FLP3FBN01A724Q orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCCGCGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.634_35 FLP3FBN01AKXJC orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCTTTTATCTCACCAACTAGCTAATCAGACGCGGATCCATCCCATACCACCGGAGTTTTTCACACAGGGCCATGCAGCCTCGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCGGTATGGGGCAGGTTATCCACGCGTT
>PC.634_36 FLP3FBN01EC04T orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.634_37 FLP3FBN01AX7NV orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.607_38 FLP3FBN01CS9UE orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGTCGGCTACCGATCGTCGCCTTGGTGGGCTGTTACCCCTCCAACTAGCTAATCGGACGCGGATCCATCTCTCACCGATTTCTCTTTCCCTCTCTCCTCATGCGACAAAAGAGGCTCATGCGGTATTAGCAACCGTTTCCAGCTGTTATCCCCCTGTGAAAGGCAGGTT
>PC.634_39 FLP3FBN01A38B0 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.634_40 FLP3FBN01CB1XF orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_41 FLP3FBN01ARU00 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCAGACGCGGGGCCATCTCATACCACCGGAGTTTTTCACACCGTGTCATGCGACACTGTGCGCTTATGCGGCATTAGCAGTCATTTCTAACTGTTATTCCCCTGTATGAGGCAGGTTCCCCACGCGTTACT
>PC.634_42 FLP3FBN01A24G3 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCGTTACG
>PC.593_43 FLP3FBN01AYGTR orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGGCCGTATTTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGATAAATCTTTGACATAAATGCCATGCGACACCTATGTGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTCCT
>PC.634_44 FLP3FBN01DPZMM orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTCTGGATGATCATCCTCTCAAACCATCTAACGATCGTCGACTTGGTGAGCCTTTACCTCACCAACTATCTAATCGTACGCAGGCCATTCCTAAAGCGCATAAATGCTTTAATCCGAAGATCATATGCGGTATTAGCCACTCTTTCGAGTAGTTATCCCTCACTTTAGGGTATGTTCCCACGCGTTACTCAGCCGTCCG
>PC.634_45 FLP3FBN01AJ7M7 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCCCTCAGGTCGGCTACTGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCTGGGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.634_46 FLP3FBN01DV3TK orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACAGAAACATGCGTCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGATTACCCACGTGTTACTCACCCG
>PC.634_47 FLP3FBN01AL42P orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTATCCCCCTCTGATAGGCAAGTTGCATACGCGTTACTC
>PC.634_48 FLP3FBN01BYPDD orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.481_49 FLP3FBN01CTMCP orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACCGGGGCATGCGCCCCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATATGGCAGGTTACCCACGCGTTACTCACCCG
>PC.607_50 FLP3FBN01BAJ58 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGACCTCTCAGTCCGGCTACCGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCGGACGCGAGCCCATCTTCGAGCGATAAAATCTTTGATACCAGCAGGATGTCCTCCCGGTATGTTATGCGGTATTAGCGACCGTTTCCAGCCGTTATTCCCCTCTCGAAGGCAGGTTGCTCACGTGTTACTCACCCGTCCG
>PC.634_51 FLP3FBN01A4K2C orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.634_52 FLP3FBN01DRBCR orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACAAACTAATGGAACGCATCCCCATCGATTATCGAAATTCTTTAATAACAAGAAGATGCCTTCTCGTTATGCTATCCAGTATTAATCTTTCTTTCGAAAGGCTATCCCGGAATAATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCAT
>PC.634_53 FLP3FBN01BWH0G orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGGTGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTC
>PC.634_54 FLP3FBN01D0JVY orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCGCCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCGCTGTTACCATGCGGTACTGCGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCACGGTTGCCCACGCGTTACTCACCCGTCCGCCGCT
>PC.634_55 FLP3FBN01EDR9T orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGGACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_56 FLP3FBN01CDTMI orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_57 FLP3FBN01AFERM orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGG
>PC.634_58 FLP3FBN01BFQFB orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCATCAACCACCGGAGTTTTCAACCCCATGCGATGCCGCATGAGATATTATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTTGATGGTAGGTTGCATACGCGTTACGCACCCG
>PC.635_59 FLP3FBN01D3JOZ orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGCTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGG
>PC.634_60 FLP3FBN01DVIVE orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATGCGCCGCATGCCCATCCATAGCCGCCTCAACTTTCCTCCCCAAGGGATGCCCCTCAGGGAGTGCACGCGGATTTAGACCGGATTTCTCCGGATTATTCCCCTGACAAGGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>PC.636_61 FLP3FBN01EQVOU orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGAGCCATTACCTCACCAACTAACTAATACGCCGCGGGATCATCTAAGAGCGACAGCAAAGCCGTCTTTCATTCATGAACCATGCGGTTCGCGAATCTATGCGGTATTAGCATCCGTTTCCGAATGTTATCCCCCACTCCTAGGCAGATTCCCC
>PC.634_62 FLP3FBN01DH4LG orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.354_63 FLP3FBN01EEU42 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.634_64 FLP3FBN01AE4G6 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_65 FLP3FBN01A559N orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCTTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.355_66 FLP3FBN01BYAA8 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGAACCATGCGGTTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.356_67 FLP3FBN01CSQSY orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGGCAGTCAGAGCCATGCGACCCAACTGCATTATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCACTCCAAGGCTAGG
>PC.634_68 FLP3FBN01DTXOJ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCGCCAACCAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAGCTTTCCCCACAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.481_69 FLP3FBN01EJD3A orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCACACCGCCCGAAGGCTTTTCACACCGCTCCATGCGAAGCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTGTGGGGCAGGTTGCCCACGCGTT
>PC.634_70 FLP3FBN01CIW8K orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTCTCTCAGTTCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACATGAGCCCATCTTTCAGCGGATTGCTCCTTTGATAACAGGATCATGCGATCCCGTCATTTCATGCGGTATTAGCACACCTTTCGATGTGTTATCCCCCTCTGAAAGG
>PC.481_71 FLP3FBN01BYVGX orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCATCGCATACCACCGGAGTTTTTCACACCAAGCCATGCGGCTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCGGTATGCGGCAGGTTGCCCACGCG
>PC.635_72 FLP3FBN01BKL5B orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCGATTAGACGCCGTTTCCTAGCGATTATTCCCCTCGATGAGGGCAAGTTGCTCACG
>PC.634_73 FLP3FBN01DX9GS orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCGG
>PC.593_74 FLP3FBN01BVHTK orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGAGCCGTTACCTCACCAACAACCTAATGGAACGCATCCCCATCCTTTACCGGAATCCTTTAATAATGAAACCATGCGGAATCATTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTAAAGGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCAG
>PC.607_75 FLP3FBN01A9PL2 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGTCACGCATGCCCATCCGTCTCCGGAGTCGCCTCCTTTGGCCGAAAGAAGATGCCTCCCCTCGGCCACATGCGGTATTAGGCCGGATTTCTCCGGGTTATCCCGCTGAAACGGGCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.634_76 FLP3FBN01C3PFP orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGGCAGTCAGAGCCATGCGACTCAACTGCATCATGCGGTATTAGCAGCTGTTTCCAACTGTTATTCCCCACTCCAAGGCAGGTTACT
>PC.634_77 FLP3FBN01A17GS orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCG
>PC.634_78 FLP3FBN01DM0XL orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACCCGTCATCGCCTTGGTGAGCCGTTACCTCACCAACTAACTAATAGGCCGCGAGCCCATCCCCAAGCGCATTGCTGCTTTAATCTTTCGATACTATGCGGTATTAATCCCAGTTTCCCAGGGCTATCCCGCTCTCGGGGGCAGGTTACTCACGTGTTACTCACCCGTGCGCC
>PC.634_79 FLP3FBN01DFQEI orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCATTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAACAAATATTCCCATGCGGGACCCCTGTTTTATGGGGTATTAATCCGACTTTCGCCGGGCTATCCCCTCTGATAGGTAAGTTGCATACGCGTTACTCACCCGTGCGCCGGTCCG
>PC.607_80 FLP3FBN01DN5MU orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGACCTCTCAGTCCGGCTACCCATCGTTGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGGGTCGCGAGCCCATCCTATACCGATAAATCTTTTACCCCTGTATCATGCGATACTGTGGTCTTATACCGTATTAGCACAAATTTCTCTGTGTTATCCCGTTGTATAGGGCAGGTTGCTCACGTGTTACTCACCCGTTCG
>PC.607_81 FLP3FBN01B1L15 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCGCTCTATCAGCGATAGCCTCTCTCGAGTCCGCCACCTTTCCTTCTGCCATCATGCGATGACTGAACCTTATGCGGTATTAGCACTGATTTCTCATTGTTATTCCCCT
>PC.634_82 FLP3FBN01AW01I orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGAACCATCGTATACCACCGGAGTTTTTCACACAGGGCCATGCGGCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCG
>PC.634_83 FLP3FBN01BI38I orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.593_84 FLP3FBN01DSC1C orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATGCGCCATAAGTCCATCCTCTACCAGTGCCTTGCAGCACTTTTAATACGGTCACCATGCAGTGTCCCTACCTATGCGGTCTTAGCTGCCGTTTCCAGCAGTTATCCCCCTGTAAAGGCCAGG
>PC.355_85 FLP3FBN01AIELL orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_86 FLP3FBN01D1T0J orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.607_87 FLP3FBN01EOHM5 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCTCGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCTCACAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCT
>PC.634_88 FLP3FBN01BMLLI orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAACAAATATTCCCATGCGGGACCCCTGTTTTATGGGGTATTAATCCGACTTTCGCCGGGCTATCCCCTCCTGATAGGTAAGTTGCATACGCGTTACTCACCCGTGCGCCGGTCG
>PC.634_89 FLP3FBN01D74QJ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
>PC.634_90 FLP3FBN01EAAV0 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_91 FLP3FBN01BC3GH orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACCTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAATGATTTAATTATTAAAAGATGCCTTCAAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.634_92 FLP3FBN01D9A6M orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.481_93 FLP3FBN01A8V4Q orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCGTTGCCGGAATCGCTCCCTTTGGCCGCAGGGCCATGCAGCCCCGCGGCATTACGCGGTATTAGACGGGATTTCTCCCGATTATCCCCCTGCAACGGGCAGGTCGCATACGCGTTACT
>PC.634_94 FLP3FBN01DCU7E orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_95 FLP3FBN01AM688 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.607_96 FLP3FBN01DMRPA orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTCCGGCTACCGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCAGGACATGCATCCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACT
>PC.607_97 FLP3FBN01BAZV6 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGAGGGCAGGTTGCATACGCGTTACTCACCCGTGCGCCGGTCGCCGGCGGAGCAAAG
>PC.634_98 FLP3FBN01DEBKQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGCAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGG
>PC.634_99 FLP3FBN01BAWOL orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGAAAGGAACATGCGTCCCCTTCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTTCGCCACTCGGGCAC
>PC.636_100 FLP3FBN01D7AVV orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAGTGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_101 FLP3FBN01B5T3T orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGGTGGGGGACCTTCCTCTCAGAACCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGTCG
>PC.635_102 FLP3FBN01B02WO orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.635_103 FLP3FBN01B2DRT orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCACTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGCCGCATGCCCATCCTTTACCGGATCGCTCCTTTGACATACCGATCATGCGGTCGGTATGTATTATGAGGTATTAGACGGAATTTCTCCCGATTATCCCTCTGTAAAGGGCAGGTCGCATACGTGTTACTCACCCG
>PC.634_104 FLP3FBN01A8KH9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_105 FLP3FBN01AU6ML orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_106 FLP3FBN01AJ4R9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACCAAAG
>PC.481_107 FLP3FBN01B32LN orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.634_108 FLP3FBN01COYKA orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTGGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.593_109 FLP3FBN01C4BE9 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCG
>PC.634_110 FLP3FBN01BU1CU orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGCCGTTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAACTAATCAGACGCAGGCCCATCCCTTACCGATAAATCTTTCACATCCCTGATATGCATCAGACTTGCATTATGAGGTTTTAATCCCGGTTTCCCGAGGCTATCCCTCTGTAAGGGGCAGGTTGCCTACGCGTTACTCACCCGTCCGCCG
>PC.634_111 FLP3FBN01AYTJU orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_112 FLP3FBN01DXJGO orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.636_113 FLP3FBN01BOOHE orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGCCGATCACCCTCTCAGGTCGGCTACGTATCGTTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATACGGCGCGGGTCCATCTATAAGTGACAGCCGAAACCGTCTTTCAACATTGAACCATGCGGTTCAATATATTATCCGGTATTAGCCCCGGTTTCCCGGAGTTATCCCAGTCTTATAGGTAGGTTACCCACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_114 FLP3FBN01BJL13 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTTCCAGTTCCAGTGTGACTGATCATCCTCTCAGACCAGTTAACCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACCAGCTAATGGTACGCGGACTCATCAAAAAGCAACAGCATATGCAGAGGCCATCTTTCCCACATAAGATAGATATGCAGCGTATCCGGTATTAGCAGCCGTTTCCAGCTGTTATCCCAAACTTCTTGG
>PC.634_115 FLP3FBN01CTRWB orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCCTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGGACCATGCAGTCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACT
>PC.636_116 FLP3FBN01BJCG1 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTC
>PC.481_117 FLP3FBN01CORWL orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCCCATACCACCGGAGTTTTTCACACAGGAGCATGCGCTCCCGTGCGCTTATGCGGTATTAGCACCTGTTTCCAGGTGGTATCCCCCGGTATGGGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.634_118 FLP3FBN01C68BF orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGAGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.635_119 FLP3FBN01B449W orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGAAAGAAACATGCGTCCCCTTCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCACTAACTCATACAT
>PC.607_120 FLP3FBN01EZ4SI orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCCCATACCACCGGAGTTTTTCACACAGGAGCATGCGCTTCCGTGCGCTTATGCGGTATTAGCACCTGTTTCCAGGTGGTATCCCCCCGGTGTGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCG
>PC.634_121 FLP3FBN01DD5NY orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCC
>PC.355_122 FLP3FBN01BHKSP orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.634_123 FLP3FBN01DHYS5 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGATAAATCTTTGACATAAATGCCATGCGACACCTATGTGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTCCT
>PC.634_124 FLP3FBN01D837F orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_125 FLP3FBN01DC9HR orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.634_126 FLP3FBN01DLDHU orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACCTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAACGATTTAATTATTATAAGATGCCTTACAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.634_127 FLP3FBN01D6Q98 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_128 FLP3FBN01CNKON orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTTACAGGTAGGTTGCCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_129 FLP3FBN01BCWX5 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCCTCCCTTTCGGAAGGTTCGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAG
>PC.635_130 FLP3FBN01EA91G orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.635_131 FLP3FBN01B06QS orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCTTCCACCACAAAAAGCTTTCAACCCAGAGAGATGCCTCTCCGAATTATATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTGAAAGGTAGGTTGCATACGCGTTACGCACCCGTCCGCCGGTCG
>PC.634_132 FLP3FBN01D0SFK orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGATAAATCTTTGACATAAATGCCATGCGACACCTATGTGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTCCT
>PC.634_133 FLP3FBN01DVTEG orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTATGGATCGTTGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCCAACGCGAGCCCATCCTTCGGCACCTCAGTTTTGGTATTTCTATCATGCGGTAAAAATACGTTATGCGGTATTACTGTCCGTTTCCAGACACTATCCCCCTCCGAAGGGCAGGTTGCTCACGCGTTACTCACCCGTCCGCC
>PC.634_134 FLP3FBN01CECHO orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCGTCGTCGGCTAGGTGAGCTGCTACCTCACCTACTACCTGATGGGCCGCGACCCCATCCCAGACCGCAAGAGCTTTCATCCCTCCGCCATGCGGTGGAAGGATAGTATCCGGTATTAGCTGCCGTTTCCGGCAGTTATCCCGATGTCTGGGGCAGGTTGGTCACGTGTT
>PC.634_135 FLP3FBN01ES2BT orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTAGTCCCCTCCTGATAGGTAAGTTGCATACGCGTTACTCACCCGTCGCG
>PC.636_136 FLP3FBN01BJO3F orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTATTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTTCATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.635_137 FLP3FBN01AMOUT orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCATCAACCACCGGAGTTTTCAACCCCATGCGATGCCGCATGAGATATTATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTTGATGGTAGGTTGCATACGCGTTACGCACCCG
>PC.356_138 FLP3FBN01CNT1Q orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATTATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>PC.634_139 FLP3FBN01DR739 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTCTCTCAGTTCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACATGAGCCCATCTTTCAGCGGATTGCTCCTTTGATAACAGGATCATGCGATCCCGTTCATTTCATTGCGGTATTAGCACACCTTTCGATGTGTTATCCCCCTCTGAAAGG
>PC.634_140 FLP3FBN01BBDL9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCTTTTACTCCACCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCTCAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATAAGGCAGGTTACCCACGCGTTACT
>PC.634_141 FLP3FBN01CLP3T orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
ATGGACCGTGTCTCAGTTCCATTGTGGCCGTTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATGCGCCGCAGGCTCATCCATCAGTGATGCCAGGAGCATCTTTAAACTTTCGTCCTATCCGGTATTAGCGATCGTTTCCAATCGTTGTCCCCGTCTGATGGGCAGATTACCTACGTGTTACTCACCCGTTCG
>PC.355_142 FLP3FBN01C1O95 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCGCCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGAACCATGCGGTTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGTAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.634_143 FLP3FBN01DHLQN orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCAGGTT
>PC.634_144 FLP3FBN01BA44D orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATGCGCCGCAGGCTCATCCATCAGTGATGCCAGGAGCATCTTTAAAACTTTCGTCCTATCCGGTATTAGCGATCGTTTCCAATCGTTGTCCCCGTCTGATGGGCAGATCACCTACGTGTTACTCACCCG
>PC.634_145 FLP3FBN01EHG6O orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTG
>PC.356_146 FLP3FBN01EE1C3 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACCGGGCCATGCGGCCCTGTGCGCTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCTGTGTAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.634_147 FLP3FBN01AER1M orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATGCGCCGCGGGCCCATCTCATAGCGGATTGCTCCTTTAATGCTACTTCATGCGAAGCTACAATCTTATGCGGTATTAATCTCCCTTTCGGAAGGCTATTCCCCTCTATGAGGCAGGTTCG
>PC.636_148 FLP3FBN01A1RXA orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>PC.634_149 FLP3FBN01B5Z5E orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATGCGCCGCAGGCTCATCCATCAGTGATGCCAGGAGCATCTTTAAACTTTCGTCCTATCCGGTATTAGCGATCGTTTCCAATCGTTGTCCCCGTCTGATGGGCAGATCACCTACGTGTTACTCACCCG
>PC.634_150 FLP3FBN01EB7G2 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTTCCAGTTCCAGTGTGACTGATCATCCTCTCAGACCAGTTAACCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACCAGCTAATGGTACGCGGACTCATCAAAAAGCAACAGCATATGCAGAGGCCATCTTTCCCACATAAGATAGATATGCAGCGTATTCGGTATTAGCAGCCGTTTCCAGCTGTTATCCCAAACTT
>PC.607_151 FLP3FBN01AYGA0 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCCTGTCCGGCCGAAGCCTTTCCTGCCTCCGGGATGCCCCGGTGGCATGTACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>PC.481_152 FLP3FBN01AL8JX orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGGGGCATGCGCCCCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATATGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.355_153 FLP3FBN01B7E8J orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.634_154 FLP3FBN01D6BHF orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACG
>PC.634_155 FLP3FBN01BLTKB orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACCTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAATGATTTAATTATTAAAAGATGCCTTCAAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.636_156 FLP3FBN01EK2JX orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCGGGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAGAGAG
>PC.634_157 FLP3FBN01EF15V orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTATCCCCCTCTGATAGGCAAGTTGCATACGCGTTACTC
>PC.634_158 FLP3FBN01BB8KH orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCCCACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACTAA
>PC.481_159 FLP3FBN01AMFGP orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGGAACATGTGTCCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.481_160 FLP3FBN01A6LEJ orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.634_161 FLP3FBN01BUCNS orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.634_162 FLP3FBN01ET0T4 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_163 FLP3FBN01B0AX7 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATGCGCCGCATGCCCATCCTTTGCCGGAATTGCTTCCTTTGACTCCCAAACCATGTGGTAAGGGAGTGTTATGCAGTATTAGTCGGAATTTCTCCCGGTTATCCTCCTGCAAAGGGCAGGTTGCATACGTGTTACTCACCCGTGCG
>PC.355_164 FLP3FBN01BM4QE orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAACTAATGCGCCGCATGCCCATCCATGACCGGATCGCTCCTTTGACTCCCGAGAGATGTCTCCCGGGGGTGTTATGCGGTATTAGACGGAATTTCTCCCGCTTATCCCCCTGTCATGGGCAAGTTGCATACGTGTTACTC
>PC.593_165 FLP3FBN01EPSX1 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACCAGCTAATCAGACGCGGGCTCATCTTATACTACCGGAGTTTTTCACACAGAAACATGCGTCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGATTACCCACGTGTTACT
>PC.634_166 FLP3FBN01DZ97Y orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTAGTCTTGGTGGGCCGTTACCCCCCAACAAACTAATGGAACGCATCCCCATCTCATACCGAATTTCTTTAATATAAAACCCATGCGGGAAATATATGCTATCGGATATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTATGAGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATCGG
>PC.634_167 FLP3FBN01BFDT5 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_168 FLP3FBN01BKTRL orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATGTCGCGCGTGCCCGTCCCGTACCACCGGAATTTTAAATCGAGAGCCATGCGGCTCTCGAGTATCATGGGATGTTAGTCCACGTTTCCGCGGGTTATCTCCCGGTACGGGGTTGGTTGCACACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAGG
>PC.634_169 FLP3FBN01C9NT9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCGTTACG
>PC.634_170 FLP3FBN01BE65F orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAGCTTTCCCCACAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.607_171 FLP3FBN01AY2QI orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAACTTTCCCCAGAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.636_172 FLP3FBN01BYYFZ orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCACGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.634_173 FLP3FBN01AFUF1 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACCAAAG
>PC.634_174 FLP3FBN01B7I8E orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACCAAAG
>PC.356_175 FLP3FBN01EIDVP orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.634_176 FLP3FBN01DOR0M orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTTGACTTGGTGGGCCGTTACCCCGCCAACTATCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAGTGAAACCATGCGGAAATACTATACTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATC
>PC.634_177 FLP3FBN01DO8MU orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_178 FLP3FBN01ATR31 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_179 FLP3FBN01CO4O9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCGTGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.356_180 FLP3FBN01CZUXO orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTCTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCACGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>PC.634_181 FLP3FBN01A9EX0 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.634_182 FLP3FBN01BYSJT orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_183 FLP3FBN01AV5YV orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.634_184 FLP3FBN01DU3OR orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCTGGGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTACGGGCAGGTTGCCA
>PC.636_185 FLP3FBN01EWW35 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_186 FLP3FBN01A3LZ0 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTATCCCCTCTGATAGGTAAGTTGCATACGCGTTACTCACCCGTGCG
>PC.634_187 FLP3FBN01D55K5 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGAGCTCATCTCAGAGCGATAAGTCTTTGACAGCCAACCCCATGCGGGATAGCTGTATCATGCGGTATTAGCGGCTGTTTCCAGCCGTTATTCCCCACTCCAAGGCAGATTG
>PC.636_188 FLP3FBN01CHNXI orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTATGGATCGTTGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCCAACGCGAGCCCATCCTTCGGCACCTCAGTTTTGGTATTTCTATCATGCGGTAAAAATACGTTATGCGGTATTACTGTCCGTTTCCAGACACTATCCCCCTCCGAAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCG
>PC.634_189 FLP3FBN01EETLW orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATGCGCCGCATGCCCATCCATAGCCGCCTCAACTTTCCTCCCCAAGGGATGCCCCTCAGGGAGTGTACGCGGATTTAGACCGGATTTCTCCGGATTATTCCCCTGCTATGGGCAGGTTGCATACGTGTTACTCACCCGTGCG
>PC.634_190 FLP3FBN01B4E5A orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCGATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACG
>PC.356_191 FLP3FBN01ET99E orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACCACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCACCGGACCATGCGGTCCTGCGGTCTTATGCGGTATTAGCAGTCATTTCTGACTGTTGTCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.634_192 FLP3FBN01BS7JZ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCC
>PC.634_193 FLP3FBN01DW7MZ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.634_194 FLP3FBN01BKWR9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAACTAATGCGCCGCGGGCCCATCTCATAGCGGATTGCTCCTTTAATTGCTACTTCATGCGAAGCTGCAATCTTATGCGGTATTAATCTCCCTTTCGGAAGGCTATTCCCCTCTATGAGGCAGGTT
>PC.636_195 FLP3FBN01DLUYP orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.634_196 FLP3FBN01DBVV5 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCAGGTT
>PC.607_197 FLP3FBN01B5DKJ orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTATACCACCGGAGTTTTTCACACCGGAGCATGCGCTCCTGTGCGCTTATGCGGTATTAACAGTCGTTTCCAACTGTTATCCCCCTGTATAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>PC.634_198 FLP3FBN01ASQES orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.634_199 FLP3FBN01DPNWR orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAGGAGATGCCTCTCGGAGTTATTACGCGATATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_200 FLP3FBN01BLIY9 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAACGATTTAATTATTATAAGATGCCTTACAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.634_201 FLP3FBN01DQGLV orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCGTCGTCGGCTAGGTGGGCCACTGCCCCGCCTACTACCTGATGGGCCGCGACCCCATCCTGCACCGCCGGAGCTTTCATCCGCTCCCCATGCGGGGTGCGGATAGTATCCGGTATTAGCTGCCGTTTCCGGCAGTTATCCCGATGTGCAGGGCAGGTTGGTCACGTGTTACT
>PC.354_202 FLP3FBN01CUQNU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTT
>PC.634_203 FLP3FBN01B0QGX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGGTGGCCGATCGCCCTCCTCAGGCCGGCTACCCATCGCAGGCTAGGTGGGCCGTTGCCCCGCCTACTACCTAATGGGCCGCGACCCCATCCCGCACCGTCGGAGCTTTCCTCCGTGGCGCATGCGCGCCTCGGAGAGTATCCGGTATTAGCCGCCGTTTCCGGCGGTTGTCCCGGGGTGCGGGGCAGGTTGGTCACGTGTTACTCACCCGTTCGCC
>PC.634_204 FLP3FBN01A3PGP orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>PC.634_205 FLP3FBN01ALY2B orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACG
>PC.634_206 FLP3FBN01BT8QN orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTCTGGATGATCATCCTCTCAAACCATCTAACGATCGTCGACTTGGTGAGCCTTTACCTCACCAACTATCTAATCGTACGCAGGCCATTCCTAAAGCGCATAAATGCTTTAATCCGAAGATCATATGCGGTATTAGCCACTCTTTCGAGTAGTTATCCCTCACTTTAGGGTATGTTCCCACGCGTTACTCAGCCGTCCG
>PC.634_207 FLP3FBN01DATE7 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGCGGGGGGTTAACCTCTCAGTCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACCAAAG
>PC.635_208 FLP3FBN01A5XYK orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTAGGCTTCTACCCCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCTCAGTTTTTACCACTGTACCATGCAGTACCGTGGTCTTATGCGGTATTAGCAATCATTTCTAATTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.607_209 FLP3FBN01BAT5U orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAACTTTCCCCAGAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.634_210 FLP3FBN01ECTYY orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGCGCCGCATGCCCATCCTCCCCCGATGAATCTTTCCTCCTCCGGAGATGCCTCCGGAGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGGGAAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_211 FLP3FBN01A2T1M orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCAGGTT
>PC.634_212 FLP3FBN01CBDAE orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.634_213 FLP3FBN01D8Y20 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGACCTTCCTCTCAGAACCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCGCCAACCAGCTAATGCGCCGCATGCCCATCCTTTACCGGATCGCTCCTTTGACATACCGATCATGCGGTCGGTATGTGTTATGAGGTATTAGACGGAATTTTCTCCCGATTATCCCTCTGTAAAGGGCAGGTCGCATACGTGTTACTCACCCG
>PC.634_214 FLP3FBN01A6MXL orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGATGTCTTGGTGGGCCGTTACCCCGCCAACAAACTAATGGAACGCATCCCCATCGATTATCGAAATTCTTTAATAACAAGAAGATGCCTTCTCGTTATGCTATCCAGTATTAATCTTTCTTTCGAAAGGCTATCCCGGAATAATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.635_215 FLP3FBN01ED1E9 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATTGGAACCAGATGCCTGATTCCATCGCATGGGGCATTAGACGCCGTTTCCAGCCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.356_216 FLP3FBN01AYIUY orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAGACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTGAGGCAGGTTATCCACG
>PC.634_217 FLP3FBN01A0UCW orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCC
>PC.636_218 FLP3FBN01CDJ70 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCTTTGGTGGGCCTCTACCCCGCCAACTGGCTAATCAGACGCGGGCCCCTCCCATACCACTATTGTTTTTCACACAGAACCATGCGGTCCCGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCG
>PC.481_219 FLP3FBN01EFY7W orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACAAGCTAATCAGACGCGGATCCATCGTATACCACCAAAAGCTTTAGCTTTTTGTTTTCCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCT
>PC.634_220 FLP3FBN01AGI0O orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.356_221 FLP3FBN01CUN6D orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTATTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.636_222 FLP3FBN01EFPI1 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.634_223 FLP3FBN01A8V8L orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCAGGCCCATCCTTTGCCGATAAATCTTTGACCAAACAGCCATGTGACTATTTAGTATTATGAGGTTTTAATCATCGTTTCCAATGGCTATCCCTCTGCAAAGGGCAGGTTGCCTACGCGTTACGT
>PC.354_224 FLP3FBN01BZOEE orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCG
>PC.593_225 FLP3FBN01EGX4W orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTATTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.635_226 FLP3FBN01D2AYU orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCACCGGACCATGCGGTCCTGTGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCGCT
>PC.634_227 FLP3FBN01BHVDZ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTCACGCGGTATTAGACGTAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGCGTTACTCACCCGTGCGCCGGTCG
>PC.481_228 FLP3FBN01ED5UR orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.634_229 FLP3FBN01CU37V orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACCGATCGTCGACTTGGTGGGCCGTTACCCCGCCAACTATCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATCCCCAAGGGATGTCCCTCAAGGATGTCATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGTGTTACTCACCCGTCCG
>PC.354_230 FLP3FBN01B86XL orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCACCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.355_231 FLP3FBN01DJ1VW orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.634_232 FLP3FBN01AKD0L orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACTCTCTCAAGCCGGCTACTGATCGTTGTCTTGGTAGGCCGTTACCCTGCCAACTAACTAATCAGACGCGGGCCCATCCTGTACCACCGTGGTTTTCCCTGCTGTGCCATGCGGCACAGCAGGCTTATGCGGTATTAGCAGCCATTTCTGGCTGTTGTCCCCCGGTACAGGGCAGGTTGCCCACGCG
>PC.356_233 FLP3FBN01D6H8C orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACT
>PC.635_234 FLP3FBN01DQFJG orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGGGCCTCTACCCCGCCAACTATCTAATCAGCCGCGAGTCCATCTCCCATCGGATTTCTCCTTTGACTATCCAAGGATGCCCTCAAATAGTATTATGCGGTATTAGCGTCCGTTTCCAGACGTTATCCCCCTTTGAAAGGTAGGTTACTCACGTGTTACTCACCCGTCCG
>PC.635_235 FLP3FBN01AGC0G orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.634_236 FLP3FBN01D0FVD orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCCACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.634_237 FLP3FBN01BG1KW orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.356_238 FLP3FBN01BW4TT orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGACCGTATCTCAGTTCCAATGTGGCCGATCAGCCTCTCAGCTCGGCTATGCATCGTTGCCTTGGTAGGCCATTGCCCCACCAACTAGCTAATACACCGCAAGCTCATCCTAAGGTGAAGCAAACGCTTCTTTTAACATATCTAGTTATCTAGTTATGTATCATCCGGTATTAGCGTTCGTTTCCAAACGTTATCCCAGTCCCTAGGGTAGATTACCC
>PC.635_239 FLP3FBN01BOFGD orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCCCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTAGGCTTCTACCCCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCTCAGTTTTTACCTCTGTACCATGCGGTACTGGGGTCTTATGCGGTATTAGCAATCATTTCTAATTGTTATCCCCCTGTGTAAGGCAGGTTGCTCACGCGTTACTCACCCGTCCGCC
>PC.355_240 FLP3FBN01CWLZU orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.593_241 FLP3FBN01AMXMK orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGGTCCATCTTATACCACCGGAGTTTTTTCACACTGTACCATGCAGTACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGGTTACCCACGCGTTACTCACCCG
>PC.354_242 FLP3FBN01BHNXZ orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.356_243 FLP3FBN01DRT2O orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>PC.607_244 FLP3FBN01ATU7A orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCACCAACTACCTAATGCGCCGCATGCCCATCCTTTACCGGATCGCTCCTTTGACATACCGATCATGCGGTCGGTATGTGTTATGAGGTATTAGACGGAATTTCTCCCGATTATCCCTCTGTAAAGGGCAGGTCGCATACGTGTTACTC
>PC.607_245 FLP3FBN01DQMNE orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACAAGCTAATCAGACGCGGGTCCATCTCGCACCACCGGAGTTTTCAGGGCAGGGGCATGCGCCCCCTCCCGTTATGCGGTGTTAGCACCTATTTCTGGGTGTTATCCCCCAGTGTAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.635_246 FLP3FBN01BGO68 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.355_247 FLP3FBN01EC56I orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGAGCTCTTACCTCACCAACTAGCTAATCAGACGCGGGCCCATCTTACACCACCTCAGTTTTTTCCACAAGGTCATGCGACCCTGTGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACT
>PC.636_248 FLP3FBN01D9EK7 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTCACTCACGTGTTACT
>PC.356_249 FLP3FBN01BXX60 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGCCGGCTGCTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACT
>PC.636_250 FLP3FBN01A208L orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.635_251 FLP3FBN01EBSVN orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTCCGCCGCT
>PC.636_252 FLP3FBN01BDPN5 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCATGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTC
>PC.635_253 FLP3FBN01DF6NE orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTTAGACCGGCAGGTTGGATACGTGTTACTCACCCG
>PC.355_254 FLP3FBN01BSM84 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.354_255 FLP3FBN01BJIH7 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCGTTGCCGGAATCGCTTCCTTTGGCCGCAGGGCCATGCAGCCCCGCGGCATTACGCGGTATTAGACGGGATTTCTCCCGATTATCCCCCTGCAACGGGCAGGTCGCATACGCGTTACT
>PC.354_256 FLP3FBN01D8FO0 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.354_257 FLP3FBN01AGYWU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTT
>PC.481_258 FLP3FBN01D6UF4 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGCGCCGCATGACCATCCGCAGCCGGATCGCTCCTTTGAACCAACAGGGATGCCCCCGTCGGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTTCATACGTGTTACTCACCCGTGCG
>PC.635_259 FLP3FBN01AEB9M orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.355_260 FLP3FBN01CIT14 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGCCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAGGAGATGCCTCTCGGAGTTATTACGCGATATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.636_261 FLP3FBN01DZG2M orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.635_262 FLP3FBN01CN7JX orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTCGGTGGGCCGTTACCCCGCCGACTGCCTAATGCGCCGCATGCCCATCCTCCACCACCGGAGTTTTCCTCCCACGGAGATGCCTCCGCGGGATTTACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCACTAACTCAAGG
>PC.636_263 FLP3FBN01B74TP orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTCCGCCACTAGGGCG
>PC.636_264 FLP3FBN01BXZ1T orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGCCGATCACCCTCTCAGGTCGGCTACGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATACGGCGCGGGTCCATCTATAAGTGATAGCAAAACCATCTTTCACTTTAGAACCATGCGGTTCTAAATGTTATCCGGTATTAGCTCCGGTTTCCCGAAGTTATCCCAGTCTTATAGGTAGGTTACCCACGTGTTACTCACCCGTCCGCCGCTAAG
>PC.355_265 FLP3FBN01C6LRW orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.355_266 FLP3FBN01AOR4O orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGGATCATGCAATCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTACCCACGCGTTACTCACCCG
>PC.607_267 FLP3FBN01CSZZD orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATGCGCCGCATGCCCATCCTCCGCCGGAATCGCTTCCTTTTAACACCCGTGGGATGCCCCACAGATGTACCACGCGGTATTAGTCCGGATTTCTCCGGGTTATCCCCCTGCGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.636_268 FLP3FBN01DNS9D orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.356_269 FLP3FBN01EGW4M orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTCTACCCGCCAACAAGCTAATCAGACGCGGGTCCATCGTATACCACCGGAGTTTTTCACACCGGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTCGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.636_270 FLP3FBN01C9TUA orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTATGGATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCCAACGCGGACCCATCCTATGCCGCCTCAGCTTTTCACACCGTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCACCCGTTTCCGGATGTTATCCCCCTGCATAGGGCAGGTTGTCCACGCGTTACTCACCCGTCCGCCG
>PC.635_271 FLP3FBN01DZ0UC orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.607_272 FLP3FBN01D9E65 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCTTCACCGAAATTCTTTCCCCCTCGGAAGATGCCTCCCAAGGAGTATATGCGGTATTAGGCGAAATTTCTTCCGGTTATCCCGCTGTAAAGGGTGCTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.636_273 FLP3FBN01ECWT9 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCTCACCAACTACCTAATCAGACGCGGGTCCCTCCTATACCACTATCGTTTTTCACACAGGGCCATGCGGCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCGCT
>PC.607_274 FLP3FBN01AHH1K orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACTATACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAACTGTTTCCAGTTGGTATCCCCCTGCATACGGCAGGTTATCCACGCGTTACTCACCCG
>PC.636_275 FLP3FBN01DQLWU orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTC
>PC.636_276 FLP3FBN01ALCMH orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCGTTGCCGGAATCGCTTCCTTTGGCCGCAGGGCCATGCAGCCCCGCGGCATTACGCGGTATTAGACGGGATTTCTCCCGATTATCCCCCTGCAACGGGCAGGTCGCATACGCGTTACTCACCCG
>PC.635_277 FLP3FBN01C4N17 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
>PC.607_278 FLP3FBN01B8F4K orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCAGAGTTTTTCACACCGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGTGTTACTCACCCG
>PC.481_279 FLP3FBN01CATW7 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
>PC.593_280 FLP3FBN01B97H8 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGATTGCTCCTTTGACCGCCCCGGGATGTCCCGGAATGGTGTTACGCGGAATTAGTCGGAATTTCTTCCGGTTATTCCCCTGTGACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.355_281 FLP3FBN01BA934 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCGCCACCGGAATCAACCTTTGGCACCAACAGGATGTCCCATAGGTGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTT
>PC.636_282 FLP3FBN01A3YNV orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGACCCATCCCGCACCGCATGCGCTTTCCGCGCGGCCCCATGCGGGGCCGTGCGCATATGCGGTATTAGCAGCCGTTTCCAGCTGTTGTCCCCCAGTGCGGGGTAGGTTGTCCACGCGTTACTCACCCG
>PC.481_283 FLP3FBN01BJOS4 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGAACCATGCGGTTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.607_284 FLP3FBN01EP59G orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACAAACTAATGCACCGCAGGTCCATCCGCACCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>PC.354_285 FLP3FBN01CTU4A orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCGGGACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.607_286 FLP3FBN01DD28G orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACAAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>PC.635_287 FLP3FBN01BW0YB orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGCGCCGCATGCCCAGCCGCCACCGGATTACTCCTTTCGCCCGGACCGGATGCCCGGTCCGGGCGGCATGGGGTATTAGGCCGGGTTTCCCCGGGTTATCCCCCTGTGGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAGGTG
>PC.355_288 FLP3FBN01AZT44 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTTCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.355_289 FLP3FBN01AGMG5 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.635_290 FLP3FBN01ARUVC orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACTAATGTTTTTCACTCTGTCCCATGCGGGACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTGTAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.355_291 FLP3FBN01D86I8 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAACTGTTATCCCCTGTATGAGGCAGGTTACCCACGCGTT
>PC.481_292 FLP3FBN01COOOS orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAG
>PC.355_293 FLP3FBN01D8J4C orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCACCAACAGGATGTCCCATAGGTGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTTACT
>PC.636_294 FLP3FBN01BR61D orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGACGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTTTAATC
>PC.356_295 FLP3FBN01BF0UO orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGCGCCGCATGCCCAGCCGCCACCGGATTCCTCCTTTCGCCCGGTCCGGATGCCCGGTCCGGGCGGCATGGGGTATTAGGCCGGGTTTCCCCGGGTTATCCCCCTCGTGGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAGGTG
>PC.354_296 FLP3FBN01DB7BE orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCTCCAACTAGCTAATCAGACGCGGGTCCATCTCATACCGTCTCGGCTTTTCACCCCGAACCATGCGGTTCTGTGTGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.607_297 FLP3FBN01AZSCT orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>PC.636_298 FLP3FBN01D3EWI orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.354_299 FLP3FBN01A55LW orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCCGCACCGGATCGCTCCTTTGACCGCTCCCCCATGCAGAGGAACGGTGTCATGCCGTATTAGTCCGGATTTCTCCGGGTTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTT
>PC.355_300 FLP3FBN01DZ6P8 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTAAGGAAAAGACATGCATCTTCTCCTGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTGTAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.355_301 FLP3FBN01AEB2N orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCCGCGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.607_302 FLP3FBN01B9U1E orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCTTCACCGAAATTCTTTCCCCCTCGGAAGATGCCTCCCAAGGAGTATATGCGGTATTAGGCGAAATTTCTTCCGGTTATCCCGCTGTAAAGGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAG
>PC.354_303 FLP3FBN01A77ZB orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.481_304 FLP3FBN01BCQ7B orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATGCGCCGCAGGCTCATCCATCAGTGATGCCAGGAGCATCTTTAAACTTTCGTCCTATCCGGTATTAGCGATCGTTTCCAATCGTTGTCCCCGTCTGATGGGCAGATCACCTACGTGTTACTCACCCG
>PC.635_305 FLP3FBN01D5TAZ orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTAGGCTTCTACCCCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCTCAGTTTTTACCTCTGTACCATGCGGTACTGGGGTCTTATGCGGTATTAGCAATCATTTCTAATTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.635_306 FLP3FBN01BX26A orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCTGTGCGCCAGTCGCCGG
>PC.481_307 FLP3FBN01DTX5C orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.354_308 FLP3FBN01BCTWA orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCACCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.636_309 FLP3FBN01BA22V orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.636_310 FLP3FBN01A4UET orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.636_311 FLP3FBN01EII4M orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTCCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
>PC.607_312 FLP3FBN01DPG83 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTGGTCTTGGTGAGCCGTTACCCCACCAACTAACTAATACAACGCATGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTATGGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCGACAAT
>PC.593_313 FLP3FBN01AYGIQ orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACAAACTAATGCACCGCAGGTCCATCCGCACCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTT
>PC.593_314 FLP3FBN01AGF7L orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGCCGATCACCCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAAGCCGTTACCCTACCAACTAGCTAATCAGACGCGGGTCCATCCTGTACCGCAAAAGCTTTGATACTTCTACCATGCGATAAAGTATTTTATCTCGTATTAGCATACCTTTCGGTATGTTATCCGTGTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.635_315 FLP3FBN01A4NV1 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGGAGTCCCCATGCAGGGTCCCCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTAAAACAGTCCGGGG
>PC.635_316 FLP3FBN01AZ963 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAACTTTCCCTAGAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACT
>PC.481_317 FLP3FBN01DMBQY orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGAATCATGCGATTCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTT
>PC.636_318 FLP3FBN01A6SLO orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCATGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCG
>PC.635_319 FLP3FBN01CQHPF orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTT
>PC.635_320 FLP3FBN01CSQG7 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTTGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCGTATACCTCCGGAAATTTTCACACTCTGTCATGCGACAGTGTGCGCTTATGCGGTATTAGCAGTTGTTTCCAACTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCG
>PC.481_321 FLP3FBN01AMENY orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGAAAGGAACATGCGTCCCCTTCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCACTAAGGCCCGCGCCTTCCGGGT
>PC.481_322 FLP3FBN01D1PPR orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTATCATGCGATACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCG
>PC.607_323 FLP3FBN01CPKVQ orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.355_324 FLP3FBN01BZK6L orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGTAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTCAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.356_325 FLP3FBN01DF6XB orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACCGATCGTCGCCTTGGTGGGCCTCTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCGTATACCACCGGAGTTTTTCACACCGGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.481_326 FLP3FBN01CF2XK orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCG
>PC.355_327 FLP3FBN01CXD3D orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTAAGGAAAAGACATGCATCTTTCCCTGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAGAT
>PC.355_328 FLP3FBN01CQGY3 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTAGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTCA
>PC.481_329 FLP3FBN01B52I9 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTTCGGCTACGCATCATTGCCTTGGTAAGCCTTTACCCCACCAACTAGCTAATGCGCCGCGGGCCCATCCAAAAGCGGTAGCATAGCCACCTTTTACATAGTTACCATGCGGTAACTATGGTTATGCGGTATTAGCACCTGTTTCCAAGTGTTATCCCCCTCTTTTGGGCAGGTTGCCCACGTGTTACTCACCCG
>PC.607_330 FLP3FBN01AS2O8 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATATGGCAGGTTACCCACGCGTTACTCACCCG
>PC.636_331 FLP3FBN01D7MLX orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>PC.607_332 FLP3FBN01EUKQQ orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAATCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>PC.635_333 FLP3FBN01AYCPW orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTAGGTGGGCCGTTACCCCGCCTACCAGCTAATCAGACGCGAGGCCATCTTCCAGCGATAAATCTTTGGTGTCTCGATGATGCCATCGAAACACATCATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCTCACG
>PC.354_334 FLP3FBN01BDSSU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.636_335 FLP3FBN01ELLAQ orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGACCATCCGCAGCCGGATCGCTCCTTTGAATCACCAGAGATGCCTCCGGTGATTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTTCATACGTGTTACTCACCCG
>PC.356_336 FLP3FBN01CEUYF orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGCCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTCTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCGTATACCACCGGAGTTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCG
>PC.636_337 FLP3FBN01CI8Z2 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.481_338 FLP3FBN01ENO6X orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCACCAAACCATGCGGTTTTGTGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>PC.636_339 FLP3FBN01B33P1 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>PC.636_340 FLP3FBN01ALS1C orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCATCGGCTTGGTGGGCCGTTACCTCACCAACTACCTAATCAGACGCGGGTCCCTCCTATACCACTATCGTTTTTCACACAGGGCCATGCGGCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCC
>PC.356_341 FLP3FBN01AY1NE orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATAAATCTTTTCCGTCAGTACCATGCGGTACCAGCGGGTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.593_342 FLP3FBN01D9HWD orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATGCGCCATAAGTCCATCCTCTACCAGTGCCTTGCAGCACTTTTAATACGGTCACCATGCAGTGTCCCTACCTATGCGGTCTTAGCTGCCGTTTCCAGCAGTTATCCCCCTGTAAAGGCCAGGTT
>PC.593_343 FLP3FBN01B0EIF orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCTTTTATCTCACCAACTAGCTAATCAGACGCAGATCCATCCCATACCACCGGAGTTTTTCACACAGGGCCATGCAGCCTCGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCGGTATGGGGCAGGTTATCCACG
>PC.354_344 FLP3FBN01BPI1O orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>PC.635_345 FLP3FBN01BBDRA orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAACTTTCCCCAGAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.354_346 FLP3FBN01C3FYG orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.356_347 FLP3FBN01DCLBL orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGACCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCTTGGGATGTCCCTAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>PC.607_348 FLP3FBN01D7ZOI orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>PC.636_349 FLP3FBN01C0173 orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCATGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTC
>PC.635_350 FLP3FBN01DPEUG orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGACCATCCGCAGCCGGATCGCTCCTTTGAATCTCCGGAGATGCCTCCGGAGATTATTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTTCATACGTGTT
>PC.636_351 FLP3FBN01B7KTL orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTC
>PC.607_352 FLP3FBN01CXEBD orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>PC.355_353 FLP3FBN01BVDVK orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.636_354 FLP3FBN01CMYHR orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGCGAGAATGCTCACTCATGATACTATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCG
>PC.356_355 FLP3FBN01D6N8J orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.607_356 FLP3FBN01COUSC orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTAGGCCATTACCCCACCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCAGAGTTTTTCACACCGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTTTCCACGCGTTACTCACCCGTCCG
>PC.481_357 FLP3FBN01AZPBJ orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCGCTGTACCATGCGGTACTGCGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCGCT
>PC.354_358 FLP3FBN01A8PO2 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.354_359 FLP3FBN01CVZXE orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAG
>PC.607_360 FLP3FBN01D1IC0 orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCACCAACCAGCTAATCAGACGCGGGCCCATCTTGCACCACCGGAATCAACCTTTGGCACCAACAGGATGTCCCGTCGATGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCCTCTGTAAAGGGCAGGTCGCATACGTGTTACT
>PC.355_361 FLP3FBN01AF0HW orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACTCTCTCAAGCCGGCTACTGATCGTTGCTTTGGTAGGCCATTACCCTGCCAACTGGCTAATCAGACGCGGGGCCATCGTATGCCGATAACTCTTTTCACACCATGCCATGCAGCATTGTGTGCTTATGCGGTATTAGCACCTATTTCTAACTGTTATCCCCCTGTGTAAGGCAGG
>PC.593_362 FLP3FBN01DS5AU orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTCACCGGATTGCTCCTTTGACCGCTCCGGGATGCCCCGGAATGGTGTTACGCGGAATTAGTCGGAATTTCTTCCGGTTATTCCCCTGTGACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>PC.356_363 FLP3FBN01DVSSP orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCGCCACCGGAATCAACCTTTGGCACCAACAGGATGTCCCATAGGTGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTT
>PC.635_364 FLP3FBN01AYYJN orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTGGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGACCCATCGTGTACCGTACTAGATAAGATCTAGGCTTTCCGCCCTGTGCCATGCGGCACTGTGCGCATATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTACACGGCAGGTTG
>PC.607_365 FLP3FBN01BQL2S orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
TTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATGCGCCGCATGCCCATCCTTTGCCGGAATTGCTTCCTTTGACTCCCAAACCATGTGGTAAGGGAGTGTTATGCAGTATTAGTCGGAATTTCTCCCGGTTATCCTCCTGCAAAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>PC.636_366 FLP3FBN01C6OIE orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACTCTCTCAAGCCGGCTACTGATCGTTGCTTTGGTAGGCCATTACCCTGCCAACTGGCTAATCAGACGCGGGGCCATCGTATGCCGATAACTCTTTTCACACCATGCCATGCAGCATTGTGTGCTTATGCGGTATTAGCAGTCATTTCTGACTGTTGTCCCCCTGCATACGG
>PC.635_367 FLP3FBN01D8PY9 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATCAGACGCGAGCTCATCTCAGAGCGATAAATCTTTGGCGTCCAGAGAGATGCCTCCCAGACGCATCATGCGGTATTAGCGGCTGTTTCCAACCGTTATTCCCCACTCCAAGG
>PC.635_368 FLP3FBN01BRKZB orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACCTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAATGATTTAATTATTAAAAGATGCCTTCAAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>PC.593_369 FLP3FBN01EHLMD orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTCTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCG
>PC.354_370 FLP3FBN01BOOJJ orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCTCCAACTAGCTAATCGGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>PC.481_371 FLP3FBN01CO1SB orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAACTAATGCGCCGCATGCCCATCCATGACCGGATCGCTCCTTTGACTCCCGAGAGATGTCTCCCGGGGGTGTTATGCGGTATTAGACGGAATTTCTCCCGCTTATCCCCCTGTCATGGGCAGGTTGCATACGTGTTACTC
>PC.355_372 FLP3FBN01D9QTP orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCGTCCG
>PC.636_373 FLP3FBN01CRC8T orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>PC.355_374 FLP3FBN01DMQWE orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.355_375 FLP3FBN01EU78N orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
>PC.636_376 FLP3FBN01AN2IX orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>PC.481_377 FLP3FBN01AMBQN orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCTGTACCACCGTGGTTTTCCCTGCTGTGCCATGCGGCACCGCAGGCTTATGCGGTATTAGCAGCCATTTCTGGCTGTTGTCCCCCGGTACAGGGCAGGTTGCCCACGCGTTACTCACCCG
>PC.354_378 FLP3FBN01B8NID orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAATCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.636_379 FLP3FBN01A2TMU orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGGCCATCCGTAGCCGGTGTTACCCTTTAAACCCCAAGAGATGCCTCTCGGAGTTATTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTACGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAG
>PC.607_380 FLP3FBN01AELDH orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGGCCTTTACCCCGCCAACCAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCG
>PC.635_381 FLP3FBN01ED2F6 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTCGGTGGGCCGTTACCCCGCCGACTGCCTAATGCGCCGCATGCCCATCCTCCACCACCGGAGTTTTCCTCCCAAGGAGATGCCTCCATGGGATTTACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCCCCCGTGCGCCAGTCGCCG
>PC.355_382 FLP3FBN01EABDR orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGGATCATGCAATCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACT
>PC.635_383 FLP3FBN01DIYZ3 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACAAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>PC.355_384 FLP3FBN01APLDI orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACCAGCTAATCAGACGCGGGCTCATCTTATACTACCGGAGTTTTTCACACAGAAACATGCGTCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGATTACCCACGTGTTACTCACCCGTCCG
>PC.355_385 FLP3FBN01D0DDK orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.635_386 FLP3FBN01CDWQ4 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGCCCATCCTTGTCCGGATAAATCCTTTGATCGAATTCTCATGCGAGAACCCGATGTCACGCGGTATTAGACCGGATTTCTCCGGCTTATCCCCCTGACAAGGGTAGGTTGCATACGTGTTACTC
>PC.354_387 FLP3FBN01AGMY0 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT
>PC.635_388 FLP3FBN01DQSWF orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTAGGCTTCTACCCCACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCTCAGTTTTTACCTCTGTACCATGCGGTACTGGGGTCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTGTAAGGCAGGTTCTCCACGCGTTACTCACCCGTCCGCCACTAAG
>PC.636_389 FLP3FBN01CXUVC orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTCAAGTAAAAGACATGCGTCTCCTACTGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTAAAGGCCAGGTTACTTATGTATTACTCACCCGTTCGCCACTCGGGC
>PC.635_390 FLP3FBN01EFNNG orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCTACCTCTCAGTACGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGAGCCCATCTTTAAGCGATAAATCTTTGATACACAAACCATGCGATTCATGTATATTATGCGGTATTAGCGGTCGTTTCCGACCGTTATCCCACTCTTAAAGGCAGGTTGCTC
>PC.607_391 FLP3FBN01EAC1O orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGACCTCTCAGTCCGGCTACCGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCGGACGCGAGCCCATCTCCGAGCGATAAATCTTTGATACCAAAGGGATGTCCCTCCAGTATGTTATGCGGTATTAGCGACCGTTTCCAGCCGTTATTCCCCTCTCGAAGGCAGGTTGCTCACGTGTTACTCACCCGTCCG
>PC.636_392 FLP3FBN01DO7JM orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGATAAATCTTTGACATAAATGCCATGCGACACCTATGTGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTCCTCACG"""

fasting_subset_otu_table = """#Full OTU Counts
#OTU ID	PC.354	PC.355	PC.356	PC.481	PC.593	PC.607	PC.634	PC.635	PC.636	Consensus Lineage
0	0	0	0	0	0	0	0	0	1	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
1	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
2	0	0	0	0	2	2	1	5	1	Root;Bacteria
3	0	0	0	0	0	0	0	0	1	Root;Bacteria
4	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
5	0	0	2	0	0	0	0	1	0	Root
6	0	0	0	0	0	0	0	1	0	Root;Bacteria
7	0	1	2	0	9	1	1	1	3	Root;Bacteria;Bacteroidetes
8	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
9	0	0	0	0	0	0	2	7	22	Root;Bacteria;Bacteroidetes
10	1	2	0	2	1	6	0	2	4	Root;Bacteria;Bacteroidetes
11	0	0	0	0	0	0	0	0	2	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
12	0	1	0	0	0	3	1	1	1	Root;Bacteria;Bacteroidetes
13	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
14	0	0	1	0	0	0	1	1	0	Root;Bacteria;Bacteroidetes
15	0	0	0	0	0	0	1	3	0	Root;Bacteria
16	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
17	0	0	0	0	0	5	0	0	0	Root;Bacteria;Bacteroidetes
18	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
19	0	2	2	4	0	5	1	5	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
20	0	0	0	0	0	0	1	1	1	Root;Bacteria;Bacteroidetes
21	0	0	0	0	0	0	0	1	0	Root;Bacteria
22	0	0	0	0	0	0	1	2	6	Root;Bacteria;Bacteroidetes
23	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";Bacillales;"Staphylococcaceae";Staphylococcus
24	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
25	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
26	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
27	0	0	0	1	0	0	1	1	9	Root;Bacteria;Bacteroidetes
28	0	0	0	0	0	1	0	0	0	Root;Bacteria
29	0	1	0	0	0	0	0	0	0	Root;Bacteria
30	0	1	1	1	1	0	3	0	4	Root;Bacteria;Bacteroidetes
31	0	0	0	0	9	1	0	0	0	Root;Bacteria
32	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
33	0	0	0	0	0	2	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
34	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
35	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
36	0	0	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes
37	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
38	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
39	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Prevotellaceae
40	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
41	6	0	5	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
42	0	0	1	0	0	0	0	0	0	Root;Bacteria
43	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
44	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes
45	0	0	0	0	2	0	0	0	0	Root;Bacteria;Bacteroidetes
46	0	0	0	0	0	5	9	5	3	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
47	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
48	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia"
49	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
50	0	0	0	0	3	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
51	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
52	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
53	0	0	0	0	0	0	13	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
54	0	0	0	0	0	2	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae;Olsenella
55	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
56	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
57	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes
58	0	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae
59	0	0	0	0	1	0	0	0	0	Root;Bacteria
60	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
61	0	0	0	0	2	0	0	4	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
62	0	2	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
63	0	0	0	0	0	0	2	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
64	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
65	5	9	0	3	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
66	0	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
67	0	2	0	0	0	1	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
68	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
69	0	1	2	0	3	0	1	0	0	Root;Bacteria;Bacteroidetes
70	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
71	0	0	0	0	0	1	1	0	0	Root;Bacteria;Bacteroidetes
72	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
73	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
74	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
75	0	11	0	0	1	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
76	2	1	10	2	24	0	0	1	1	Root;Bacteria
77	0	0	1	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
78	0	0	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
79	0	1	0	0	0	0	0	0	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
80	0	4	3	0	1	2	0	2	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
81	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
82	0	0	0	0	0	1	0	0	0	Root;Bacteria
83	0	0	0	1	0	0	1	2	19	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
84	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
85	0	1	1	1	0	0	1	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
86	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
87	0	0	1	0	0	0	0	0	0	Root;Bacteria
88	0	0	0	0	1	0	0	0	0	Root;Bacteria
89	0	0	0	0	0	0	0	0	1	Root;Bacteria
90	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
91	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
92	0	7	1	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
93	0	0	0	0	1	0	0	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
94	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
95	0	0	0	1	0	0	1	0	0	Root;Bacteria;Bacteroidetes
96	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
97	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
98	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
99	0	0	0	0	0	0	0	2	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
100	0	0	0	0	0	0	1	0	0	Root;Bacteria
101	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
102	0	0	0	0	0	1	0	0	0	Root;Bacteria
103	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
104	0	1	1	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
105	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
106	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
107	0	1	1	4	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
108	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
109	0	0	0	0	0	0	5	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Porphyromonadaceae;Parabacteroides
110	0	0	1	0	0	0	0	0	0	Root;Bacteria
111	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
112	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
113	0	0	1	0	0	0	0	0	0	Root;Bacteria
114	0	0	0	0	0	6	0	3	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
115	0	0	0	0	0	0	2	1	0	Root;Bacteria;Bacteroidetes
116	0	0	0	0	0	1	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales
117	1	0	0	5	17	20	0	0	0	Root;Bacteria
118	0	0	0	0	0	3	5	2	5	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
119	0	1	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
120	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
121	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
122	1	3	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
123	0	1	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
124	2	1	0	5	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
125	0	1	0	3	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
126	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
127	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Bryantella
128	0	1	0	1	1	0	0	0	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
129	1	0	3	0	0	4	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
130	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
131	1	0	0	0	0	0	0	0	0	Root;Bacteria
132	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
133	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
134	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
135	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
136	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
137	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII;Anaerovorax
138	0	0	0	0	0	1	0	0	0	Root;Bacteria
139	0	0	0	0	0	0	1	0	0	Root;Bacteria;Bacteroidetes
140	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
141	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
142	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
143	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
144	0	0	0	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
145	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
146	2	3	8	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
147	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
148	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
149	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
150	1	0	0	2	4	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
151	0	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
152	1	0	0	0	0	4	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
153	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
154	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
155	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
156	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
157	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
158	0	0	0	0	0	0	0	0	1	Root;Bacteria
159	1	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia"
160	0	1	1	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
161	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
162	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
163	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
164	0	0	0	0	0	0	1	0	0	Root;Bacteria
165	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
166	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
167	8	0	0	0	3	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
168	0	0	0	0	0	1	0	0	0	Root;Bacteria
169	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
170	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
171	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
172	1	4	2	6	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
173	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
174	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
175	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
176	0	0	0	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
177	3	0	0	0	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
178	29	1	10	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
179	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
180	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
181	0	0	0	0	0	0	0	1	1	Root;Bacteria;Bacteroidetes
182	0	0	0	1	0	0	0	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
183	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
184	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
185	2	0	2	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
186	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
187	0	0	0	4	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
188	1	0	0	0	10	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Peptostreptococcaceae";"Peptostreptococcaceae Incertae Sedis"
189	0	0	0	0	1	0	0	0	0	Root;Bacteria
190	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
191	0	0	0	1	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
192	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes
193	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
194	0	0	0	1	0	0	0	0	0	Root;Bacteria
195	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
196	9	0	0	0	5	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
197	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
198	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
199	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
200	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
201	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
202	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
203	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
204	0	1	0	3	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
205	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
206	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
207	0	0	0	0	0	0	0	1	0	Root;Bacteria
208	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
209	1	0	0	0	0	0	0	0	0	Root;Bacteria
210	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
211	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
212	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
213	0	0	0	0	0	1	0	0	0	Root;Bacteria
214	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
215	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
216	1	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
217	0	0	0	1	0	0	0	0	0	Root;Bacteria
218	0	0	0	0	0	1	0	0	0	Root;Bacteria;Proteobacteria;Deltaproteobacteria
219	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
220	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
221	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
222	0	0	0	0	0	7	0	2	2	Root;Bacteria;Bacteroidetes
223	0	0	0	0	0	0	0	2	0	Root;Bacteria
224	0	0	0	0	0	0	0	0	1	Root;Bacteria
225	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
226	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
227	2	18	0	1	0	0	21	4	4	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
228	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
229	1	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
230	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
231	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
232	0	0	0	0	1	0	0	0	0	Root;Bacteria
233	0	0	0	1	0	2	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
234	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
235	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
236	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
237	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
238	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
239	0	0	0	0	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
240	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
241	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
242	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
243	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
244	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
245	0	0	0	4	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
246	0	0	0	0	0	0	1	0	0	Root;Bacteria
247	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
248	0	0	0	0	0	0	0	1	0	Root;Bacteria
249	0	0	0	0	1	0	0	0	0	Root;Bacteria
250	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
251	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
252	0	2	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
253	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
254	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
255	1	3	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
256	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
257	0	0	0	0	0	1	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII
258	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
259	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
260	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
261	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
262	0	1	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
263	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
264	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
265	1	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
266	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
267	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
268	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
269	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
270	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
271	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
272	0	1	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
273	1	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
274	0	0	0	1	0	0	1	5	2	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
275	0	0	0	0	0	0	0	1	0	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Rikenellaceae;Alistipes
276	1	4	3	2	0	0	0	0	0	Root;Bacteria;Bacteroidetes
277	0	0	1	0	0	0	0	0	0	Root;Bacteria
278	0	0	0	0	0	1	0	0	0	Root;Bacteria
279	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
280	2	2	0	1	0	0	0	0	0	Root;Bacteria;Bacteroidetes
281	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
282	0	0	0	1	4	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
283	0	0	0	1	0	2	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
284	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
285	1	0	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
286	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes
287	0	0	0	0	5	2	0	0	0	Root;Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter
288	0	0	0	0	0	1	0	0	0	Root;Bacteria
289	0	0	1	0	0	0	0	0	0	Root;Bacteria
290	0	1	0	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
291	0	0	0	2	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
292	9	12	5	13	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
293	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
294	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
295	4	2	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
296	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
297	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
298	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
299	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
300	0	0	0	0	0	0	0	2	0	Root;Bacteria
301	0	0	0	0	0	0	8	8	3	Root;Bacteria
302	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
303	0	0	0	3	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
304	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
305	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
306	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
307	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
308	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
309	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
310	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
311	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
312	0	0	0	8	0	0	1	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Turicibacter
313	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
314	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
315	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
316	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
317	1	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
318	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
319	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
320	1	0	0	0	0	0	0	1	0	Root;Bacteria
321	0	0	0	0	2	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
322	0	0	0	0	2	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
323	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia"
324	0	2	1	0	0	0	0	0	0	Root;Bacteria;Bacteroidetes
325	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
326	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
327	0	0	0	1	0	0	0	0	0	Root;Bacteria
328	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
329	0	0	2	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Acetanaerobacterium
330	0	0	0	0	0	1	0	0	0	Root;Bacteria
331	0	0	0	0	4	0	0	0	2	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
332	0	0	0	0	0	1	0	0	0	Root;Bacteria
333	0	0	0	0	0	0	0	0	1	Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
334	0	0	0	0	0	0	1	0	0	Root;Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Verrucomicrobiaceae;Akkermansia
335	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
336	0	0	0	0	1	0	0	0	0	Root;Bacteria
337	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
338	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes
339	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
340	0	0	0	0	0	0	1	0	0	Root;Bacteria;Actinobacteria;Actinobacteria;Coriobacteridae;Coriobacteriales;Coriobacterineae;Coriobacteriaceae
341	0	0	0	0	0	0	0	0	1	Root;Bacteria;Deferribacteres;Deferribacteres;Deferribacterales;Deferribacteraceae;Mucispirillum
342	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
343	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
344	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
345	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
346	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
347	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
348	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
349	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
350	0	1	0	0	0	0	0	0	0	Root;Bacteria
351	0	0	0	0	0	0	0	1	0	Root;Bacteria
352	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
353	0	0	0	0	0	2	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
354	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";"Carnobacteriaceae";"Carnobacteriaceae 1";Atopostipes
355	0	0	0	0	0	0	0	2	0	Root;Bacteria;Bacteroidetes
356	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Coprobacillus
357	0	0	0	0	0	0	2	3	1	Root;Bacteria;Bacteroidetes
358	0	0	0	0	0	0	1	0	0	Root;Bacteria
359	0	0	0	0	0	0	2	0	0	Root;Bacteria;TM7;TM7_genera_incertae_sedis
360	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
361	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
362	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
363	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
364	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";"Lachnospiraceae Incertae Sedis"
365	1	0	0	0	2	0	1	0	2	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
366	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
367	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
368	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
369	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Butyrivibrio
370	0	0	1	0	0	0	0	0	0	Root;Bacteria
371	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
372	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
373	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes
374	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
375	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
376	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
377	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
378	0	0	0	0	0	0	4	0	0	Root;Bacteria;Firmicutes;"Erysipelotrichi";"Erysipelotrichales";Erysipelotrichaceae;Erysipelotrichaceae Incertae Sedis
379	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
380	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
381	0	0	1	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
382	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
383	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
384	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
385	0	0	0	1	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
386	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
387	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
388	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
389	0	0	0	0	0	0	1	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
390	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
391	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
392	0	1	0	0	0	0	3	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Clostridiaceae;"Clostridiaceae 1";Clostridium
393	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
394	0	0	1	0	0	0	0	0	0	Root;Bacteria
395	0	0	0	0	0	0	0	1	0	Root;Bacteria
396	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
397	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
398	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
399	0	0	0	0	0	0	0	1	0	Root;Bacteria;Actinobacteria;Actinobacteria
400	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
401	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
402	0	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
403	0	0	0	0	0	0	2	0	1	Root;Bacteria;Proteobacteria;Deltaproteobacteria
404	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
405	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Streptococcaceae;Streptococcus
406	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
407	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
408	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
409	0	0	0	0	0	0	1	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
410	1	0	0	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
411	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
412	0	0	0	0	0	0	0	0	1	Root;Bacteria;Firmicutes
413	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
414	0	0	3	0	0	0	0	1	3	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
415	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes
416	0	0	0	0	0	1	0	0	0	Root;Bacteria
417	14	1	14	1	0	0	0	0	0	Root;Bacteria;Firmicutes;"Bacilli";"Lactobacillales";Lactobacillaceae;Lactobacillus
418	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
419	1	0	0	0	0	0	0	2	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";Ruminococcus
420	0	0	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
421	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
422	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
423	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae";"Ruminococcaceae Incertae Sedis"
424	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
425	0	0	0	0	0	0	2	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Ruminococcaceae"
426	0	0	0	0	0	0	0	1	0	Root;Bacteria
427	0	0	0	0	0	1	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
428	0	1	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
429	0	0	1	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
430	1	1	1	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
431	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
432	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
433	1	0	0	0	0	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales
434	0	0	0	0	0	0	0	1	0	Root;Bacteria;Firmicutes
435	1	0	1	0	0	0	0	1	1	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae"
436	0	0	0	0	1	0	0	0	0	Root;Bacteria;Firmicutes;"Clostridia";Clostridiales"""

fasting_subset_tree = """(91:0.03138,((56:0.01554,(171:0.00015,41:0.00014)0.868:0.01019)0.786:0.00551,(80:0.01059,(((260:0.00581,244:0.02632)0.871:0.00571,((319:0.01626,70:0.00524)0.764:0.00509,(298:0.03941,369:0.01625)0.858:0.01096)0.925:0.00016)0.907:0.0156,(272:0.00015,219:0.02761)0.371:0.00527)0.886:0.01577)0.904:0.01538)0.213:0.00506,(((((124:0.02331,193:0.01557)0.550:0.02651,390:0.02213)0.922:0.02481,((376:0.02606,(420:0.01061,(375:0.01758,(398:0.05303,(((429:0.00014,(427:0.00569,409:0.00793)0.175:0.00015)0.715:0.00016,(434:0.00016,426:0.02959)0.962:0.01738)0.942:0.02633,(414:0.01138,(413:0.01803,424:0.02362)0.839:0.01221)0.715:0.00567)0.706:0.00547)0.860:0.01803)0.748:0.00612)0.961:0.04689)0.973:0.0476,(((152:0.06102,400:0.06529)0.948:0.04668,((((188:0.14301,207:0.09317)0.630:0.00562,((((318:0.00503,320:0.00632)0.963:0.06151,((421:0.02097,(430:0.00014,(((432:0.03044,366:0.01302)0.756:0.01196,361:0.02203)0.807:0.01147,431:0.02457)0.835:0.0118)0.928:0.02756)0.899:0.03039,(401:0.00901,433:0.09069)0.881:0.02632)0.941:0.06451)0.909:0.04762,(170:0.04503,187:0.02247)0.806:0.01706)0.835:0.0181,347:0.03498)0.900:0.02275)0.843:0.01899,(((((21:0.06701,(137:0.00015,257:0.07002)0.343:0.01903)0.900:0.04272,((248:0.00016,(139:0.03037,(115:0.00014,95:0.03129)0.822:0.02203)0.998:0.07935)0.999:0.11876,(82:0.00014,(59:0.0964,(117:0.00015,332:0.02889)0.836:0.01403)0.926:0.01929)1.000:0.23864)0.050:0.00109)0.759:0.03556,((87:0.01453,(57:0.00252,(58:0.04076,(((((((116:0.04298,46:0.00016)0.942:0.02413,84:0.02412)0.751:0.00986,((274:0.01098,275:0.0279)0.767:0.00508,74:0.01398)0.758:0.00703)0.746:0.00561,51:0.07406)0.864:0.01536,16:0.00015)0.955:0.03768,35:0.00015)0.997:0.03175,19:0.00015)0.861:0.02236)0.996:0.07959)0.974:0.06082)0.901:0.04353,(278:0.06061,(((28:0.02168,(88:0.00612,(20:0.00544,((((10:0.0,12:0.0,14:0.0):0.00014,(6:0.04332,36:0.0156)0.911:0.00016)0.947:0.01021,26:0.00014)0.722:0.00016,24:0.00508)0.995:0.06124)0.842:0.0489)0.997:0.0834)0.999:0.0992,((131:0.02646,((30:0.00014,(100:0.01017,69:0.00014)0.928:0.00015)0.808:0.0051,(((((((((45:0.00506,((224:0.00014,144:0.01568)0.399:0.01048,76:0.00016)0.787:0.00519)0.815:0.00508,44:0.00016)0.860:0.01021,7:0.00526)0.754:0.00768,113:0.04271)0.692:0.01687,((2:0.01311,4:0.03731)0.748:0.01016,(104:0.00015,276:0.03287)0.973:0.03782)0.892:0.02397)0.854:0.01508,((((((285:0.03762,(15:0.10332,((358:0.04675,(((158:0.06378,354:0.1522)0.803:0.03759,(194:0.00973,((189:0.01519,249:0.02483)0.954:0.03416,68:0.01853)0.789:0.01036)0.934:0.05047)0.819:0.02898,(164:0.03164,277:0.02089)0.136:0.01345)0.829:0.0189)0.942:0.06408,((301:0.01074,(355:0.01673,300:0.02274)0.867:0.01124)0.785:0.00582,357:0.00539)1.000:0.17151)0.830:0.03922)0.939:0.05382)0.701:0.01345,67:0.02688)0.823:0.017,(((109:0.04428,63:0.03926)0.979:0.05805,(37:0.08014,79:0.0163)0.220:0.02124)0.923:0.05077,(39:0.10195,(93:0.01059,((((53:0.00352,(((227:0.00016,(83:0.00014,286:0.06319)0.991:0.03771)0.788:0.00014,333:0.06987)0.998:0.06107,86:0.00015)0.994:0.04517)0.349:0.00596,50:0.04512)0.503:0.0083,(34:0.0426,(32:0.02608,(49:0.06727,90:0.00496)0.974:0.02959)0.523:0.00358)0.863:0.02028)0.936:0.03486,136:0.07971)0.734:0.00702)0.961:0.08516)0.910:0.05397)0.866:0.0375)0.894:0.02444,((8:0.02694,(3:0.03192,5:0.00016)0.973:0.04398)0.978:0.04594,(33:0.02159,(106:0.03443,96:0.03596)0.864:0.0251)0.831:0.02373)0.304:0.00014)0.978:0.02591,(71:0.051,(222:0.02838,246:0.00436)0.937:0.0244)0.284:0.01734)0.429:0.00645,(((((27:0.0,9:0.0):0.00014,(25:0.01023,38:0.02084)0.885:0.01026)0.133:0.00014,22:0.01547)0.827:0.0208,(324:0.03204,17:0.0441)0.123:0.01)0.947:0.03497,280:0.02253)0.718:0.00493)0.842:0.00611)0.000:0.00016,181:0.02181)0.925:0.01592,31:0.02127)0.477:0.00997,(289:0.01603,290:0.01098)0.959:0.02745)0.924:0.00014)0.998:0.07977)0.749:0.0055,(((209:0.01801,325:0.03269)0.853:0.01642,(42:0.05299,29:0.01425)0.827:0.014)0.417:0.00015,232:0.02092)0.982:0.03398)0.649:0.00014)0.953:0.02633,(89:0.0342,(138:0.0459,327:0.04257)0.940:0.04586)0.883:0.02242)0.907:0.0455)0.948:0.06336)0.574:0.04246)0.813:0.02689,(((102:0.01072,54:0.05569)0.961:0.07824,((0:0.0468,111:0.02902)0.786:0.01542,((340:0.03302,(302:0.02503,399:0.03369)0.781:0.01109)0.945:0.03337,182:0.0403)0.832:0.01488)0.992:0.10706)0.904:0.05948,(((211:0.02581,213:0.02321)0.969:0.06007,(142:0.04144,223:0.05483)0.800:0.02186)0.955:0.05472,((((267:0.03308,(134:0.07285,150:0.00015)0.813:0.00015)0.366:0.00872,((391:0.08217,417:0.00818)0.524:0.03815,(330:0.07239,370:0.05715)0.980:0.09661)0.925:0.06443)0.974:0.0939,((11:0.0269,23:0.00015)1.000:0.15357,((378:0.00014,(168:0.05313,97:0.00066)0.999:0.08303)0.999:0.11217,((331:0.04314,356:0.08848)0.965:0.05433,(405:0.1024,(312:0.01087,(308:0.03192,307:0.0212)0.910:0.00015)0.995:0.08955)0.688:0.01126)0.228:0.00016)0.806:0.02968)0.833:0.03579)0.642:0.00917,(334:0.12276,392:0.05097)0.847:0.02572)0.877:0.02757)0.426:0.02497)0.933:0.04732)0.871:0.02613,(((287:0.26862,(217:0.08165,218:0.00929)0.884:0.05245)0.705:0.05092,((403:0.09144,359:0.22127)0.852:0.04039,351:0.07049)0.676:0.02494)0.959:0.07101,((350:0.10412,(402:0.00014,410:0.00551)0.945:0.06389)0.892:0.04421,(412:0.00448,415:0.03109)0.992:0.06861)0.462:0.0102)0.914:0.03206)0.981:0.00016)0.831:0.01613,(((((((((((395:0.0607,((((228:0.06268,419:0.02151)0.928:0.04763,((173:0.00015,176:0.01044)0.978:0.05847,(130:0.10381,(((165:0.00015,((225:0.03875,226:0.00536)0.937:0.01609,((281:0.01056,(386:0.02218,326:0.02751)0.987:0.00015)0.867:0.02806,329:0.0339)0.744:0.00533)0.338:0.01064)0.794:0.00618,284:0.021)0.594:0.01533,163:0.06412)0.670:0.00331)0.833:0.01554)0.827:0.01685)0.796:0.01647,425:0.06687)0.769:0.00951,374:0.00443)0.714:0.01154)0.714:0.00872,(363:0.01064,387:0.03462)0.889:0.01546)0.000:0.00502,396:0.06189)0.480:0.01186,(360:0.04307,(416:0.05007,237:0.02283)0.740:0.00734)0.835:0.01155)0.921:0.01801,((365:0.00014,389:0.00014)0.846:0.00538,(407:0.00015,(385:0.01085,(233:0.00728,404:0.00015)0.166:0.00014)0.963:0.02238)0.828:0.00729)0.942:0.00015)0.151:0.01071,((352:0.00013,380:0.01643)0.954:0.02126,(247:0.01614,243:0.04517)0.430:0.00014)0.987:0.03338)0.758:0.00496,(353:0.00014,382:0.00014)0.999:0.00014)0.871:0.0346,379:0.00364)0.970:0.03438,(337:0.04417,(422:0.00664,((341:0.0223,(118:0.03295,279:0.00014)0.888:0.02922)1.000:0.11353,(411:0.02454,423:0.0304)0.693:0.01914)0.931:0.05373)0.904:0.03125)0.719:0.00734)0.909:0.02779,(159:0.04393,221:0.00015)0.999:0.07613)0.313:0.02026,259:0.01429)0.933:0.02439)0.521:0.00282)0.846:0.00221,((((127:0.04567,(((191:0.00016,(236:0.01619,408:0.01671)0.757:0.00515)0.542:0.00016,(256:0.00014,342:0.01065)0.997:0.03284)0.721:0.02916,(394:0.05783,(151:0.04317,(103:0.01259,161:0.02538)0.818:0.01669)0.468:0.02681)0.415:0.01576)0.916:0.03996)0.922:0.04147,(323:0.01447,338:0.00321)0.853:0.02723)0.389:0.01609,((406:0.08533,(162:0.04584,(153:0.06999,322:0.02039)0.809:0.02113)0.098:0.01333)0.844:0.01478,((((123:0.02634,(169:0.01045,174:0.01039)0.949:0.00015)0.693:0.01633,269:0.03999)0.838:0.01016,(((201:0.02365,283:0.00356)0.852:0.01226,(114:0.01093,346:0.01663)0.691:0.00343)0.884:0.01235,(112:0.0158,(195:0.02712,212:0.02152)0.872:0.01073)0.750:0.0053)0.753:0.00522)0.753:0.00535,(180:0.04103,(155:0.03597,242:0.04333)0.767:0.00883)0.919:0.0276)0.323:0.00014)0.958:0.02654)0.123:0.00016,((192:0.03384,(((110:0.13051,(373:0.00534,(254:0.0107,264:0.0218)0.763:0.00528)0.953:0.00015)0.830:0.01036,156:0.05032)0.141:0.00015,206:0.01589)0.844:0.01058)0.738:0.00516,((335:0.01999,372:0.02628)0.912:0.02056,(266:0.00015,294:0.02745)0.853:0.01082)0.856:0.01148)0.902:0.01565)0.816:0.00518)0.736:0.00672)0.886:0.01398)0.396:0.00015,((13:0.00015,(1:0.0051,(48:0.03225,18:0.00014)0.610:0.00519)0.912:0.01029)0.908:0.01953,(179:0.0515,(52:0.01038,149:0.06621)0.812:0.0109)0.739:0.00921)0.783:0.00867)0.952:0.01558,((261:0.05074,262:0.02236)0.919:0.03024,(((238:0.03629,166:0.00767)0.097:0.00537,(303:0.01324,383:0.03403)0.887:0.02784)0.870:0.01529,(((72:0.00982,(47:0.00014,(146:0.0082,239:0.0327)0.763:0.00803)0.841:0.01299)0.138:0.00826,((62:0.00015,((((105:0.02651,186:0.00015)0.775:0.00513,(215:0.00874,(203:0.00016,305:0.03337)0.824:0.00879)0.766:0.00877)0.933:0.01583,(94:0.0389,(65:0.00014,(184:0.0267,185:0.02183)0.938:0.00015)0.883:0.01609)0.804:0.01015)0.799:0.00516,((((((((((((200:0.01046,(98:0.01556,(73:0.01564,(122:0.02128,175:0.00015)0.889:0.01398)0.762:0.00707)0.970:0.03203)0.938:0.00014,349:0.02174)0.532:0.00015,(((((((299:0.01635,306:0.01635)0.288:0.00542,135:0.03053)0.776:0.007,(345:0.01083,(119:0.00512,310:0.00014)0.789:0.00517)0.810:0.00514)0.458:0.00014,255:0.00015)0.916:0.00016,343:0.01618)0.804:0.0102,(388:0.06006,(190:0.03767,315:0.01441)0.751:0.00755)0.878:0.01614)0.867:0.00016,(((183:0.01056,258:0.0162)0.763:0.00521,(241:0.01533,(202:0.01143,362:0.0273)0.551:0.01038)0.782:0.00615)0.801:0.00519,(304:0.01611,364:0.03412)0.898:0.00015)0.763:0.00015)0.857:0.00514)0.782:0.00513,252:0.01635)0.754:0.00514,(((268:0.0139,(147:0.00998,265:0.03126)0.228:0.01262)0.793:0.00914,((296:0.00014,(245:0.00528,291:0.01087)0.884:0.01072)0.564:0.02251,((348:0.00837,384:0.00827)0.951:0.04585,(428:0.02172,435:0.02712)0.838:0.02582)0.743:0.03345)0.023:0.00016)0.926:0.016,(229:0.00015,((309:0.02156,(157:0.01646,(235:0.0739,288:0.02349)0.935:0.04288)0.888:0.02448)0.806:0.01075,((143:0.01222,99:0.02577)0.065:0.0049,(120:0.02212,148:0.03357)0.732:0.00443)0.943:0.02146)0.851:0.01089)0.855:0.01021)0.889:0.01046)1.000:0.00014,((101:0.01909,(141:0.00014,154:0.01582)0.763:0.00487)0.875:0.01291,((128:0.01036,253:0.00519)0.734:0.00539,(((240:0.02482,(126:0.01979,(140:0.0054,199:0.03843)0.738:0.01722)0.550:0.00726)0.066:0.00598,(40:0.01546,((107:0.00014,((263:0.02739,((230:0.00525,251:0.02704)0.668:0.00015,(((231:0.02647,297:0.01435)0.843:0.01391,(393:0.01089,(75:0.00014,160:0.00014)0.397:0.00015)0.843:0.00513)0.635:0.00014,92:0.00016)0.838:0.00511)0.957:0.01032)0.321:0.00014,85:0.0051)0.840:0.00508)0.877:0.01022,(61:0.00016,66:0.03726)0.983:0.03719)0.018:0.00015)0.797:0.00886)0.893:0.01188,132:0.00015)0.974:0.02665)0.931:0.01579)0.788:0.00514)0.871:0.01027,((133:0.00464,381:0.02942)0.932:0.01625,(121:0.01046,344:0.01099)0.921:0.01599)0.908:0.00015)0.964:0.00015,(313:0.03321,314:0.04605)0.895:0.01872)0.862:0.02644,(((177:0.00015,(((418:0.01756,((129:0.0,214:0.0):0.01191,(273:0.0172,316:0.02655)0.704:0.00543)0.747:0.00614)0.880:0.01298,282:0.01152)0.761:0.00625,(234:0.02239,436:0.03053)0.762:0.00988)0.898:0.00015)0.838:0.00517,((167:0.0,196:0.0):0.00016,321:0.01066)0.928:0.01041)0.975:0.00015,(295:0.00015,(292:0.00016,(((108:0.00824,(210:0.00509,336:0.06592)0.913:0.02095)0.752:0.00816,328:0.00552)0.249:0.01037,271:0.00015)0.959:0.02084)0.580:0.0103)0.920:0.01027)0.865:0.01052)0.751:0.00597,((293:0.01657,270:0.04174)0.741:0.00513,((311:0.05063,371:0.01364)0.418:0.01153,(205:0.00787,367:0.01549)0.887:0.01792)0.700:0.00769)0.935:0.0179)0.871:0.01345,(250:0.01588,339:0.01855)0.940:0.02734)0.837:0.0123,(((((178:0.01106,((125:0.0,204:0.0):0.00531,220:0.0385)0.740:0.00464)0.672:0.00524,(197:0.04008,81:0.00398)0.917:0.01677)0.527:0.00015,(((64:0.04233,(317:0.00014,(368:0.03503,377:0.01132)0.757:0.00501)0.725:0.00686)0.876:0.01641,172:0.00015)0.785:0.00554,145:0.00934)0.874:0.01052)0.674:0.00549,78:0.00469)0.829:0.00531,(198:0.03931,(216:0.05763,77:0.0047)0.807:0.0104)0.553:0.0107)0.802:0.00014)0.903:0.01049)0.741:0.00014)0.859:0.00015,(55:0.0211,(397:0.0598,60:0.03851)0.760:0.01044)0.752:0.00488)0.762:0.00647)0.821:0.01394,(208:0.04985,43:0.00598)0.937:0.02237)0.829:0.02204)0.606:0.00014)0.335:0.00016)0.788:0.00544);"""

template_alignment_subset = """>114239
............................................................................................................aGAGTTT-GA--T-CC-T-G-GCTC-AG-AT-TGAA-C-GC--TGG-C--G-GT-A-TG--C----T-T--AACACA-T-GC-A-AGT-CGA-A-CG----------G-TAA-CA-G-----------------------------CAG-A-AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTT-CT------------------GGCT--G--AC--G--AG-T-GG-C-GG-A--C-------------GGG-TGAGT-A--AC-GC-G-T-A-GG---A-A--T-CT-G--C-CTTA---CA-G------------------------------------------------------------------T-GG----GGG-AT-AG-CCC-------------------------G-G-T-----------------------GAA-A---ACC-GGA-TTAA-TA---CC-G--C-AT-A----------C--------------------G-------------------------------------CC-C-----------------------------------------------------------------------------------------------------------------------T-AC-G--------------------------------------------------------------------------------------------------------------------------------------G-G-G---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-CTG--G-----G--GA-T--C--------------------------------------------------------------------------------------------------------------------TTC-G----------------------------------------------------------------------------------------------------------------------G-A--CC-TG--G---C-A--------------C----T-G---T-TA-G---AT---G-A-----G-CCT-GCG--T-GAG--A------TT--A--G-CT-T----G---TTGG-T-G-GG-G-T----AAT-GG---C-CTACCA--A-GG-C-A--A-CG-A------------TCT-C-T------AG-CT-G-G-TCT-G-AG----A--GG-AC--G-AT-C-AG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-GCCCAG-A-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GC-GA-G----A-G-CC-T-GA-TG-CA-GCAA-TACC-G-CG-T---G-T-G--T--GA-A-G--A--A-G-G--C----CTG-AG---------G-G-T-T-G-T--A---AA-G-CAC--------TT-TC-A-A--T--TG---TGA-A--G---AAAAGCT---T-TT-GG----T--T--AA-T---A----------AC-C-TTGAGTC-TT-GA-CA-TTAA-C-A--A-TA-C---------AA-----------GAAGC-ACC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GGT-GCA-A-G-CG-TTAA-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GAGT-AC--G-TA-G-G-T-G------------G--T-TC-G-T-T-AA----GTC---A---G-ATG-TG-A-AA-GC--CC-CGG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-T-G-GG-AA-C----T-G-C-A-T-T--------T---GAA-A-C-T-G-GCA--A-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-A-----T-GG--TA-G-A------------G-GC-G-GG-T----AG--AA-T-TCA-TAGT--GT-A-GCG-GTGAAA-TG-CGT-AGAT-ATTA-T-G-A--GG-A-AT-A-CC-AG--T--G--GC-GAA-G--G-C---G----G--C-C-CGCTG------G-AC-CA--------------------------------------------------------------AT-A-C-C--GA--CA-----CT-GA-GG--T-A-CGA--AA-G-C--------------G-TGGG-GAG-C-G-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-TC--AA-CT---------A-GC--C--G-T-TG-G-GA-C-----------------------------------------------------------------------------------------CTT-GA--------------------------------------------------------------------------------------------------------------------------------------------------G-G-T-CT--C-A-G-T-GG-T------GC--A----GC-TAA--CG-C---G--T-GAA-GT--T----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGC-T--AAA-ACTC-AAA---------TGAA-TTG-ACGAG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-ATCCC-TT-G-AC-ATC-A---------------TC-A-G-------------A-A--CTT-G--TT--A-GA-G-A-T--A-A-C----T-CGG--T-G-----CC-------------------------------------T--TC-G------------------------------------------GG------------AA--CTGA-AT--GA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TC--C-TTAG--T-T-G-C-C---AG-C-A----CG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAT------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GG---T----G-G-----------------G---A-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---TAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-G-C-C-CTT----AT-G--GG-A-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-T-CA-GTA--C-AAA-GG-GC--------------------------------------------------------------------------------------------------A-G-C-C-A--A-CTCG-C--G---------------------------------------A-GA-G-T-----------G--C-G-CA---A----------A--TCC-C------A-T-AAAGCTGA---T-C-G-TAG-TCC--------GGA-T-TGGAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTCC-A-T-G-AA-G-TT-GGAAT-CG-C-TA--G-TA-AT-C-G-T----GAA-TC-A-G--A------AT--GTC-AC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGCA----CACACCG-CCC-GTC-----a---ca--cca-tg-gg-a--g---tgg-g-ct-gc-aaa--a-gaa------g--t-agg-ta-g-t-t-t-aa-c-c--------------------------------------------------------------ttc-g------------------------------------------------------------------------------------------------------gg-a--ga-a--c---gc-tta--cc--act-t----t-gtg-gt-tca------------------------tg--act-gggg-tg-aag-tcgtaacaa-ggtag-ccct-aggggaa-cctg-gggc-tggatcacctcctt.................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
>200959
............................................................................................................agagttt-ga--t-cc-t-g-gctc-ag-at-tgaa-c-gc--tgg-c--g-gc-a-gg--c----c-t--aacaca-t-gc-a-agt-cga-a-cg----------g-taa-ca-g----------------------------gaag-c-ag----------------------------------------------------ctt-g----------------------------------------------------------------------------------ctg-cttt----------------g-ct--g--ac--g--ag-t-gg-c-gg-a--c-------------ggg-tgagt-a--at-gt-c-t-g-GA---TAA--T-CT-A--C-C-CTT--GA-G------------------------------------------------------------------T-CA----AGG-AT-AA-GCA-------------------------C-G-G-----------------------GAA-A---CTG-TGG-CTAA-CA---CT-T--G-AT-A---------------------------------G--G-G-A--A--T-----------------TA-GGCC-T-----------------------------------------------------------------------------------------------------------------------G-GA-A--------------------------------------------------------------------------------------------------------------------------------------A-G-GC-A-T---------------T--T-T-C-C-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAAC-TCGT---------------------------------------------------------------------------------------------------------------------------------------GCT-A----------------------------------------------------------------------------------------------------------------------------------C-------------------------------A-AG-G---AT---G-A-----G-TCT-GCG--G-ATG--A------TT--A--G-CT-T----G---TTGG-T-G-GG-G-T----AAC-TG-C-C-T-ACCA--A-GG-C-G--A-TG-A------------TCA-T-T------AC-GG-G-C-TCT-T-AG----C--GG-AG--T-TG-C-CC-GGAG-A-AGGAA--T-C-TG-A-GA-C-AC-T-A-TTCCTA-GCTC-TAC-G--G-A-G-T-GC-A-GC-A-G-TC---GC-G-G-AAC-ATTTA-C-AA-T-GC--AC-GA-A----A-G-TG-T-GA-TA-AT-GCAA-GCCA-A-AG-T---G-C-T--T--A----------T-C-A--------TT-T--------------A-G-A--------T--A-AGC--------TT-TT-T-C--T-TAG----T------------------------------G--T--AA-A---A----A------------------------------T-C-TA-G-GA-G---------AA-----------TAAGGATCTGGG-A--AA---G-AC-TGGT--GCCA--G-C---C--GCCG---C-GG--TA-AT--AC---CA-GC-AGT-TCA-A-G-TT-GCAT-C-CAG-TT-TT-A--T-T--GGGT-CTA----AA-ACAT-CC--G-TA-G-C-T-T------------G--T-TT-A-C-T-AA----G-T-C-T---C-CTG-TG-A-AA-TC--CT-ACG-T--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------GC-G-T-GG-GG-CT---T-G-C-A-G-G--------A--GA-T-A-C-T-G-GTA--A-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-A-T-----C-GG--TA-G-A------------C-GT-A-AG-G----AG--TACG-CCG-G-GG--GT-A-ACG-GTTAAA-TG-TGT-TAAT-C-CT-C-GGT--GG-A-CT-A-AC-AA--T--G-AGC-GAA-G--G-C---A----C--C-T-TACGA------G-GA-C---------------------------------------------------------------GAAT-C-T--GA--CA-----GT-GA-GG--G-A-TGA--AG-G-C--------------C-AGGG-GCG-C-A-AAAC--GG-ATTA-G-ATA-C-----CC-G-T-GTA-G-T----C-CT--G-G-CAG-T-AAA--C-GCTG-CT--CA-CT---------A-AA--C--A-T-CG-G-GC-C--C--TC----------------------------------------------------------------------------------TTC-GG---------------------------------------------------------------------------------------------------------------------------------------------GA---G-G-A-TT--C-G-G-T-GC-T------GA--A----GC-GAA--GG-C-G-A--T--AA-GT--G----A-GCT-ACC-T-G-GG-AAG-TA---TAG-----T-C--G-C-A-A-GAT-C--GAA-ACTT-AAA---------GGAA-TTG-GC-GG-G-G-AGAC----ATA--C-A-A-CGA-GT-G--A--CG-CG-T--GC-GGT-TC-RATT-A-G-ATT-TTAC-A-CC-G-TG-A-A-CC-TC-A-CC-AGGAG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CGACAGC--AGGATGAAGGTCAGTCTGAAGGGCTTACC-TGA-CACGCTGAG-C-G-G-A-GATG-CA-TGG-CC--GAC-GTC-A-GC-CT---G-TG-CT-G--TGA-AG-TGA-C-CC-G-TT-AA-AT-CAGGT-AA--------C-AGG-CGA-G-ACC-C-G-TG--TC--G-TTAA--T-T-A-C-T---AC-A-G--A--AA-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-G------------G----G---T-A--CC---------------T-T-A-A-C-GA-G--AC-T-G-CC---T-T------------------------------------G-G---TAA----------------------------------C-A-A-G--G-A-GG-A--AGG-T--GCGG-G-CTAC-GTT--AGGT-G---AGT-A-C-G-T-C-C-CTA----AT-C--CC-C-T-GG-GC-TA-CAC-GCGCG-C--GG--CAATG---G-TATG-T-T--C-AAA-TG-GA--------------------------------------------------------------------------------------------------T-G-C-A-A--C-TCCG-A--A---------------------------------------A-GG-A-G-----------A--A-G-CT---A----------A--TCC-C------C-T-AACGC-AT-A-T-C-T-CAG-TCC--------AGA-T-TGAGG-GT--T-GCAA-CC-C-------------------------------------------------------------------------------------------------A-CCCTC-A-T-G-AC-G-AT-GGAAT-TC-G-TA--G-TA-AT-C-G-G----GTA-TT-A-G--C------AC--TGC-CC-G-GT-G-AAT-AAGT-C-CCTGTCTCT-TGCA----CACACCG-CCC-GTC-----A---AA--CCA-AC-CG-A--G---TTA-T-GT-AG-GGG--C-GA-------A--G-CCC-T--T-A-T-T-aa-c-c-------------------------------------------------------------t-tc-g------------------------------------------------------------------------------------------------------gg-a--gg-g--c---gc-tta--cc--act-t----t-gtg-at-tca------------------------tg--act-gggg-tg-aag-tcgtaacaa-ggtaa-ccgt-aggggaa-cctg-cggc-tggatcacctcctt.................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
>83155
............................................................................................................agagttt-ga--t-cc-t-g-gctc-ag-aa-cgaa-c-gc--tgg-c--g-gc-g-cg--t----c-t--taagca-t-gc-a-agt-cga-a-cgg---------c-aa-gag-g------------------------------aa---ag----------------------------------------------------ctt-g----------------------------------------------------------------------------------ct--tt------------------tctC--C--TA--G--AG-T-GG-C-GG-A--C-------------TGG-TGAGT-A--AC-AC-G-T-G-GG---TGA--T-AT-A--C-C-TTT--TG-G------------------------------------------------------------------A-CG----GGG-AT-AG-CCC-------------------------C-T-A-----------------------GAA-A---TAG-GGG-GTAA-TA---CC-G--G-AT-A----------A--------------------G-G--C-C-G--T--G----------------CG----GG-T-----------------------------------------------------------------------------------------------------------------------T-GG-A--------------------------------------------------------------------------------------------------------------------------------------G-C-C----G---------------C--A-C-G-G-G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-G-C-GC-----------------------------------------------------------------------------------------------------------------------------------TAT-G----------------------------------------------------------------------------------------------------------------------------G-CG--C---C-G--------------C----C-G---A-AA-G---AA---T-A-----G-CCT-GCG--A-CCT--A------TC--A--G-CT-T----G---TTGG-T-G-AG-G-T----AAA-GG-C-C-C-ACCA--A-GG-C-G--A-TG-A------------CAG-G-T------AT-CC-G-G-CCT-G-AG----A--GG-GT--G-AA-C-GG-ACAC-A-TTGGG--A-C-TG-A-GA-T-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-CT---AA-G-A-ATA-TTCCG-C-AA-T-GG--AG-GA-A----A-C-TC-T-GA-CG-GA-GCAA-CGCC-G-CG-T---G-G-A--C--GA-T-G--A--A-G-G-CC-----GG-AA---------G-G-T-T-G-T--A---AA-G-TCC--------TT-TT-A-T--A-ATT----GA-G--G---AATAAGC---GGGA-CA----G--G--GA-A--------T-----GG-T-TCC-GCG-GT-GA-CT-GT-A-G-AT-TAT--G---------AA-----------TAAGC-ACC-GG-C-TAA---T--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AC--AC---GT-AA-GGT-GCG-A-G-CG-TTAT-T-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GGGC-AC--G-CA-G-G-C-G------------G--C-TT-T-G-C-AA----G-C-T-T---G-GTG-TG-A-AA-TC--TC-AGG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-C-T-GA-AA-C----T-G-C-A-T-T--------G--AG-A-A-C-T-G-CAT--G-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-T-T-----A-GT--GA-G-G------------G-GA-A-AC-C----GG--AATT-CCA-G-GT--GT-A-GGG-GTGAAA-TC-TGT-AGAT-A-TC-T-GGA--AG-A-AC-A-CC-AA--T--G--GC-GAA-G--G-C---A----G--G-T-TTCCA------G-CA-CA--------------------------------------------------------------TA-A-C-T--GA--CG-----CT-CA-GG--T-G-CGA--AG-G-T--------------G-CGGG-GAT-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CG--C-A-CAG-T-AAA--C-GATG-TG--CA-CT---------A-GG--T--G-T-CT-G-GG-C------------------------------------------------------------------------------------------AT-AA----------------------------------------------------------------------------------------------------------------------------------------------------G-C-CC--G-G-G-T-GC-C------AA--A----GC-AAA--CG-C-G-A--T--AA-GT--G----C-ACC-GCC-T-G-GG-GAG-TA---TGC-----C-C--G-C-A-A-GGG-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-GTAC-G-CG-A-GG-A-A-CC-TT-A-CC-TGGGT-TT-G-AC-A-T-A------------C-ATA-C-C-------------G-AT-T-A-T--TT--A-GA-G-A-T--A-A-G--T-A-A-G--C-G-----TA-------------------------------------G--CA-A------------------------------------------TA----C----G---GG-TAT---G--AA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-CC-G--TGA-GG-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-C-TA--CT--G-CCAG--T-T-A-C-T---AA-C-A--G--G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TGA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AG---C----T-G------------A----G---G-A--CT---------------C-T-G-G-C-GG-A--AC-T-G-CCG--G-C------------------------------------G-A---CAA----------------------------------G-C-C-G--G-A-GG-A--AGG-C--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-G-C-C-CTT----AT-G--TC-C-A-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-TCGG-T-A--C-AGA-GT-GA--------------------------------------------------------------------------------------------------T-G-C-G-A--G-GCAG-T--G---------------------------------------A-TG-C-G-----------G--A-G-CA---A----------A--ACG-C------A-G-AAAAC-CG-A-T-C-G-TAG-TCC--------GGA-T-TGGAG-TC--T-GAAA-CC-C-------------------------------------------------------------------------------------------------G-ACTCC-A-T-G-AA-G-TT-GGAAT-CG-C-TA--G-TA-AT-C-G-C----ACA-TC-A-G-C-------AC--GGT-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TC-CG-A--G---TCG-A-GG-AT-ACC--C-AAA------G--T-CGG-TA-G-G-C-T-AA-C-C---------------------------------------------------------------GC-A---------------------------------------------------------------------------------------------------A--GG-G--GG-C--C---GC-TGC--CC--AAG-G----T-ATG-CT-TGG------------------------TA--AGG-GGGG-TG-AAG-TCGTAACAA-GGTAG-CCGT-ACTGGAA-AGTG-Cggc-tggatcacctcctt.................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
>105109
............................................................................................................AGAGTTT-GA--T-CC-T-G-GCTC-AG-GA-TAAA-C-CT--TGG-A--A-GGAG-CA--C----A-T--AAGACA-T-GC-A-AGT-CGA-A-CGA---------A-C---AG-A-------------A--------------GCTA-T-TC----------------------------------------------------TTA-T----------------------------------------------------------------------------------AGA-TTAGT-----------------GG--A-AGT--T--AG-T-GG-C-GG-A--C-------------TGG-TGAGT-A--AT-GT-A-T-A-AG---TAA--C-CT-G--C-C-TAT--CA-G------------------------------------------------------------------A-GG----GGA-AC-AA-CAG-------------------------T-T-G-----------------------GAA-A---CGA-CTG-CTAA-TA---CC-G--C-AT-A----------T--------------------G-C--C-A-T-AA--G-----------------GT---TC-G-------------------------------------------------------------------------------------------------------------------------CA-T--------------------------------------------------------------------------------------------------------------------------------------G-G-A--C-C---------------A--A-GT--G-G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-G-A--------------------------------------------------------------------------------------------------------------------------------------GCA-A-------------------------------------------------------------------------------------------------------------------------------T--C---C-G--------------C----T-G---A-TA-G---AT---G-G-----A-CTT-ATA--T-CTG--A------TT--A--G-CT-A----G---TAGG-T-AGAG-G-T----AAT-GG-C-T-C-ACCT--A-GG-C-G--A-CG-A------------TCA-G-T------AG-CC-G-G-ACT-G-A-G---A--GG-TT--A-AA-C-GG-CCAC-A-TTGGT--A-C-TG-A-GA-T-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TC---GG-G-A-ATA-TTGCG-C-AA-T-GG--AG-GA-A----A-C-TC-T-GA-CG-CA-GTGA-CGCC-G-CG-T---G-C-A--G--GA-A-G--A--A-G-G-TT-----TT-CG---------G-A-T-T-G-T--A---AA-C-TGC--------TTTAG-A-C--A-GGG----AA-G--A---A-----------------------G--AA------------------------------GTGA-CA-GT-A-C-CT-G-TA-G---------AA-----------TAAGC-TCC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GGA-GCG-A-G-CG-TTAT-C-CGG-AT-TT-A--T-T--GGGT-GTA----AA-GGGT-GC--G-TA-G-A-C-G------------G--G-AA-G-T-T-AA----G-T-T-A---G-TTG-TG-A-AA-TC--CC-TCG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CT-G-A-GG-AA-C----T-G-C-A-A-C--------T--AA-A-A-C-T-G-ATT--T-T-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-A-----C-TG--GA-G-A------------G-GA-A-AG-T----GG--AATT-CCT-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAT-A-TT-A-GGA--GG-A-AC-A-CC-AG--T--G--GCG--AAG--G-C---G----A--C-T-TTCTG------G-AC-AG--------------------------------------------------------------AA-A-C-T--GA--CG-----TT-GA-GG--C-A-CGA--AA-G-T--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-C-T-GGT-A-G----TCCA--C-A-CCG-T-AAA--C-GATG-GA--TA-CT---------A-GG--T--G-T-AG-G-GC-A----------------------------------------------------------------------------------------T-GA-AT-AA-------------------------------------------------------------------------------------------------------------------------------------------------T-G-CT--C-T-G-T-GC-C------GT--C----GC-AAA--CG-C-A-A--T--AA-GT--A----T-CCC-ACC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGT-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCA-GT-G--G--AG-TA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGGC-TT-G-AC-A-T-A-------------TAAG-T-G-------------A-AT-A-A-C--TA--A-GA-G-A-T----T-A--G-T-T-A--G-------TT-------------------------------------C--TT-C------------------------------------------GG---------A---AC-ACT-T-A--TA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-AG-G-TT-AA-GT-CCTGC-AA--------C-GAG-CGC-A-ACC-C-C-TG--TT--C-TTAG--T-T-G-C-C---AG-C-ATGT----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AAAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------A----T-G------------G----G---A-A--CT---------------C-T-A-G-G-AA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-T-C-G--G-A-GG-A--AGG-T--GGGG-A-CGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-TTT----AT-G--TC-C-T-GG-GC-TA-CAC-ACGTA-C--TA--CAATG---G-CTAT-A-A--C-AAA-GG-GA--------------------------------------------------------------------------------------------------A-G-C-G-A--A-GTAG-T--G---------------------------------------A-TA-T-G-----------G--A-G-CA---A----------A--ACC-C------A---AAAAAGTA-G-T-C-T-CAG-TTC--------GGA-T-TGAAG-GC--T-GAAA-TT-C-------------------------------------------------------------------------------------------------G-CCTTC-A-T-G-AA-G-CT-GGAAT-TG-C-TA--G-TA-AT-G-G-C----AGG-TC-A-G-C-------AT--ACT-GC-C-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-AG-A--G---TTG-G-GA-AT-ACC--C-GAA------G--C-CTG-TG-A-G-C-C-AA-C-CG------------------------------------------------------------T-AA-G------------------------------------------------------------------------------------------------------GG-G--GC-A------GC-AGT--CG--AAG-G----T-AGA-AT-CAA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGTAG-CCGT-ATCGGAA-GGTG-CGGC-TGGATCACCTCCTTTCT..............................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
"""

lanemask = """00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110111101011001110100101101011001000010100111111010110101110111010111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011010110101101001000000000000011101111101001101101010101100011100101101001010111001101000000000000000000000000000000000000000000000000000000000000000000101100001110110110111000000000000000000000000010101000000000000000000000001110100011101110111101100011010010110100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100001010001011010001100010100000101110111001011100100000011001001011010000100011110101011010100001110110101010111100101101010010110100000000000011101010000001101101010111010110000100110110010110101101111010111110010101101011010110101011111101111011101001010101011010110101011000110101011101111101011010110011011010000101011010110110110111101111010110100010101001001101010010010101011000001101100000000010101010101001000110101110000000011011010100101110000110100100000000000000000000000000000000000000000000000000000000000000000001101101010110101101000000000110000000000011111011101101011100010110111100111100101000100111100010110011011001100011011011101110101011011110101110110110100101001111011100001101111011001011010101010000000000001001011010101011000010101010001011101101011011001101110100000000000000000000000000000000000000000000000000000000000000000000110101100000000000000000000000000000000000000000000000000000000000000000000000001101010110110100001010101010100000000100110101010101011100101010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001010101010101000001011001101010000000000001011010110100001100111101110101100110101110111111011011101111010110101110011010110101101100100100110111010010100010000100101011111000000101101100000000000000000000000000000000000000000000000000000000000000110101010011001100000110110110010101110011010100000000000000101111011101010111100110111101011101000001101010111010100001011001010111010111001011110110011011000000000101100100101011010100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101010101101000000110010000110111001101010100100110110010000101110111010101101110110001110000010100101010101110100111011110111000000000111101110111110101011110000101001010101110110100100110110100110111011011110101011101111010110101101010110110101101111101101011010101000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001010101010111101101110110011101110101101100010110110100111011011101011010110110110111110110000000010111011101011101010110011001011110010101010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101001100000000000000010101010101101001101010111001010000000000000000000000000000000000001010001110000000000000000000000000000000000101010100101011010011101001111010111101110011110100011101010101010101110000110100110101011011011011101111101001100111110001011110101001011101101100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101010101001011110100100000000000000000000000000000000000000010110101000000000001001010110001000000000010011101000000101011111011010101010111011100000000111010111110110010111101101000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010111110101010110101101111101101011001011011010101000011101101010100000001100111011010110101110111101011111111101111000011111110111011100000100011001110110110100100011101011011011100101110000001001011100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111001100111010000101110110111000000000000000000000000110011101111011011101111111110111110111100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
"""
fasting_subset_sff = """Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  647488
  Index Length:  9251
  # of Reads:    392
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FLP3FBN01ELBSX
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1766_0111

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       254
  Clip Qual Left:   5
  Clip Qual Right:  254
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.13	1.07	0.19	0.21	0.86	0.20	0.97	0.20	1.01	1.06	0.13	0.11	0.95	0.15	1.05	0.15	1.03	0.15	1.04	1.08	0.06	1.01	1.89	0.14	0.05	1.04	0.15	1.00	0.10	0.98	0.15	0.16	1.02	0.18	0.15	1.00	0.13	0.19	0.95	0.18	0.14	1.03	0.16	1.04	0.11	0.17	1.00	0.14	0.11	1.81	0.14	1.05	0.13	2.79	0.94	1.01	1.02	0.19	1.87	0.13	0.99	0.18	1.01	1.03	0.13	1.06	0.10	1.01	0.10	0.08	2.89	0.12	0.09	1.95	0.87	0.94	0.10	0.10	1.05	1.01	0.08	1.04	0.13	0.98	0.10	1.06	0.15	0.08	1.03	0.11	1.03	0.97	0.11	2.84	0.13	0.08	1.95	0.13	0.13	1.01	0.01	0.18	0.83	1.02	0.10	0.09	1.74	0.11	0.09	1.80	1.07	2.75	1.11	2.72	0.13	0.94	0.12	0.84	0.15	1.05	0.13	1.02	0.16	0.12	1.03	0.06	1.88	0.20	0.06	1.85	1.88	0.10	0.08	0.87	0.05	1.09	1.07	0.91	1.06	0.10	0.07	0.99	0.11	0.06	1.00	0.11	0.11	1.07	0.06	1.05	0.09	0.03	1.09	0.07	0.03	0.93	0.06	1.04	0.95	0.02	0.11	2.09	0.00	1.91	0.09	0.03	1.94	0.92	0.06	0.00	2.78	0.08	0.13	2.14	1.04	2.08	1.16	1.95	0.14	1.04	0.06	0.81	0.19	0.11	0.93	2.02	0.13	0.12	1.83	1.08	0.04	1.12	0.83	0.08	1.04	0.08	0.00	1.06	0.16	1.01	1.90	0.06	0.10	1.06	0.00	0.09	0.98	0.13	0.11	0.97	0.96	0.17	0.09	1.91	1.02	0.10	0.12	0.95	0.00	0.18	1.07	0.06	2.20	1.08	0.09	2.03	0.10	0.09	0.96	0.00	0.17	1.05	0.08	2.01	0.16	0.01	1.04	0.08	0.09	1.00	0.03	0.11	1.06	1.85	0.03	1.10	0.00	0.08	0.92	0.78	1.12	0.01	0.04	2.01	0.01	1.93	0.00	0.04	1.04	0.10	1.04	0.16	0.00	1.10	0.00	0.04	3.00	0.12	0.10	1.13	0.82	0.13	0.05	0.99	0.09	2.88	1.87	0.16	0.13	0.76	0.61	0.12	0.17	0.92	0.99	1.18	0.06	1.04	0.00	0.09	1.03	0.10	1.03	0.07	1.02	0.14	0.08	1.10	0.10	0.04	1.03	0.13	0.01	1.02	0.03	0.09	1.08	0.06	0.07	0.96	0.97	0.07	0.09	1.00	0.10	0.98	0.11	1.05	0.11	0.93	0.10	0.15	1.21	1.06	0.99	0.03	0.06	1.04	0.92	2.32	0.06	1.04	0.99	0.01	0.15	0.78	0.00	2.28	2.11	4.18	1.11	0.00	0.90	0.13	0.13	1.14	0.04	0.85	1.05	1.87	1.15	2.85	0.00	1.79	0.00	0.11	1.01	0.02	1.07	0.09	0.12	1.14	0.00	0.11	0.98	0.08	1.20	2.03	0.97	0.16	0.00	1.02	0.00	3.20	2.20	0.14	0.89	0.72	0.07	0.17	0.93	1.20	0.08	0.02	0.98	0.17	0.12	1.21	0.12	0.09	3.20	0.05	0.00	1.14	0.00	0.64	1.33	0.35	2.36
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	72	75	75	76	77	80	81	83	85	87	90	92	93	95	95	95	98	98	101	104	105	108	108	111	111	112	113	113	113	114	115	115	115	117	119	121	123	126	128	128	131	131	132	132	135	137	138	139	140	143	146	149	151	154	157	159	160	163	163	165	165	168	168	169	172	172	172	175	175	176	177	177	178	179	179	181	183	186	187	187	190	190	191	193	194	196	199	201	202	202	205	208	211	212	215	215	216	219	222	224	224	225	227	227	230	233	235	235	238	241	244	245	245	247	250	251	252	255	255	257	257	260	262	265	268	268	268	271	272	275	277	277	277	278	278	281	282	285	286	287	289	292	294	296	299	302	305	308	311	312	315	317	319	321	324	325	326	329	330	331	331	333	334	337	339	339	340	340	341	341	341	341	342	344	347	349	350	351	351	352	353	353	353	355	355	358	360	363	366	368	369	369	370	373	375	375	375	376	376	378	379	382	383	386	389	392	392	392	395	397	398	400	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
Quality Scores:	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	33	33	33	36	37	37	37	37	37	37	40	40	40	39	39	38	40	40	40	40	40	40	40	37	37	37	37	37	35	35	35	37	37	37	37	37	35	35	35	31	31	23	23	23	31	21	21	21	35	35	37	37	37	36	36	36	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	28	28	28	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	35	32	32	32	32	35	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	32	32	32	36	37	35	32	32	32	32	32	32	32	32	36	37	37	37	37	36	36	31	31	32	32	36	36	36	36	36	36	36	36	36	36	36	28	27	27	27	26	26	26	30	29	30	29	24	24	24	21	15	15	13	13

>FLP3FBN01EG8AX
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1719_1463

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       280
  Clip Qual Left:   5
  Clip Qual Right:  280
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.05	0.13	1.03	0.18	0.21	0.93	0.20	0.94	0.17	1.04	1.04	0.10	0.13	1.05	0.18	1.01	0.14	1.04	0.14	1.03	1.08	0.06	1.01	2.02	0.14	0.04	0.95	0.12	1.04	0.08	0.98	0.15	0.14	1.08	0.19	0.12	1.04	0.12	0.15	0.96	0.13	0.11	1.01	0.12	0.93	0.10	0.16	1.02	0.13	0.12	2.15	0.13	1.03	0.14	2.57	1.08	1.01	1.04	0.18	1.83	0.20	0.91	0.17	0.96	2.68	0.11	0.19	1.92	0.14	1.07	1.95	0.99	1.05	0.14	0.14	0.98	1.07	0.11	1.04	0.11	0.97	0.15	0.95	0.15	0.14	0.96	0.12	1.06	2.02	0.08	1.93	0.13	0.14	1.96	0.12	0.12	1.03	0.04	0.17	1.26	0.93	0.08	0.19	5.41	0.12	0.10	2.00	0.00	1.93	0.08	1.83	0.13	1.01	0.08	1.00	0.15	1.04	0.09	0.83	0.19	0.13	1.07	0.18	1.09	0.11	2.03	3.90	0.11	1.08	0.82	0.16	0.10	1.02	0.11	1.88	0.09	0.07	1.04	0.13	0.13	1.04	0.05	1.09	1.02	0.11	1.93	0.14	1.88	0.11	0.03	1.01	0.04	1.85	0.09	0.13	2.12	0.93	0.10	0.11	2.96	0.10	0.13	2.16	0.83	2.13	1.17	3.97	1.01	0.07	0.02	1.72	0.08	0.13	2.15	0.92	0.12	0.13	1.89	1.88	0.09	1.00	1.81	0.06	0.16	1.06	0.14	0.12	2.05	0.11	1.89	1.05	1.01	0.13	0.05	0.96	0.10	0.07	0.88	0.21	0.11	1.10	0.14	3.70	0.17	0.06	0.97	0.13	0.06	1.09	0.13	1.05	1.04	0.06	1.03	0.15	0.00	1.03	0.03	0.06	1.04	0.13	0.91	1.96	0.94	0.09	2.04	0.16	1.09	2.00	0.00	1.09	0.13	2.80	2.05	0.11	0.06	0.87	0.91	0.11	1.11	2.00	0.14	1.02	0.07	0.90	1.07	1.91	0.16	0.03	0.95	0.12	0.04	1.06	0.05	0.11	1.04	0.09	0.01	1.10	1.98	0.17	1.96	0.09	0.96	0.08	1.79	1.07	0.10	1.00	0.79	0.09	0.06	1.03	0.07	0.15	1.07	0.09	0.14	2.06	0.11	0.04	1.09	0.01	0.11	1.01	0.00	1.11	2.95	0.90	1.08	0.16	0.00	1.89	2.04	0.08	0.17	0.96	0.01	1.03	0.13	2.77	0.06	0.92	0.04	2.77	0.08	1.00	0.95	0.04	2.95	0.11	1.95	0.06	0.07	1.02	0.08	1.08	1.10	0.12	0.05	0.92	0.12	3.68	0.93	0.02	0.98	0.15	1.11	1.03	0.14	1.04	0.09	0.01	1.02	0.04	0.11	1.07	0.03	0.98	1.94	0.11	0.04	1.11	0.00	0.09	0.94	0.11	1.95	2.10	0.00	0.16	1.87	0.11	0.88	0.03	0.00	1.07	1.11	1.16	1.05	0.88	0.11	0.00	0.97	1.84	1.04	1.10	0.14	0.98	0.15	0.85	0.00	0.00	1.05	2.79	1.07	1.17	0.02	0.00	0.74	0.08	0.10	1.01	1.02	0.10	0.10	1.78	2.02	1.12	0.10	0.83	0.71	0.15	0.12	1.96	0.16	0.25	1.09	0.16	0.19
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	146	146	148	148	151	153	153	156	156	157	160	160	160	163	163	164	165	165	166	167	167	167	167	168	171	171	174	174	175	178	178	179	179	181	182	182	185	188	188	190	190	191	192	195	198	201	203	203	203	203	206	209	211	212	214	217	220	222	223	223	224	226	226	228	229	229	231	233	233	233	234	234	237	238	240	241	241	243	245	246	247	247	250	253	256	259	260	260	262	262	264	266	266	267	269	270	273	276	279	279	282	285	287	288	288	288	289	290	293	293	294	294	297	299	301	301	301	303	305	305	305	307	308	310	310	310	312	312	315	317	318	321	323	323	323	323	324	326	328	329	331	334	337	339	340	340	343	346	348	348	349	349	352	352	354	357	358	359	360	361	364	365	365	366	367	369	371	374	375	375	375	376	377	380	383	384	387	387	388	388	389	391	392	395	395	398
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	37	33	33	21	21	21	26	33	37	36	36	40	33	24	24	29	33	33	39	39	39	40	39	39	39	40	37	37	37	37	37	37	37	37	37	37	37	32	32	20	20	20	20	20	35	35	37	37	37	37	37	37	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	36	37	37	37	37	37	36	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	33	28	28	28	28	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	28	28	28	37	28	28	28	37	37	37	37	37	36	36	36	36	36	28	26	26	26	26	28	36	36	36	36	36	36	36	37	38	38	38	38	38	37	37	37	37	37	31	31	31	31	31	31	31	31	31	31	31	31	30	22	22	22	25	25	31	31	31	31	31	31	31	25	25	25	25	25	28

>FLP3FBN01EEWKD
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1692_3531

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       249
  Clip Qual Left:   5
  Clip Qual Right:  248
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.16	1.06	0.23	0.23	0.88	0.23	0.95	0.22	0.88	0.14	1.01	0.16	0.07	1.12	0.15	0.14	1.04	1.07	1.06	0.15	1.01	0.15	1.03	0.16	0.09	1.96	0.15	1.05	1.07	0.94	0.19	0.18	0.96	0.20	0.17	0.98	0.16	0.17	0.96	0.15	0.15	1.09	0.15	1.03	0.16	0.13	1.01	0.14	0.15	1.71	0.17	1.03	0.14	2.65	1.00	1.11	0.95	0.17	1.92	0.10	0.94	0.20	0.92	2.80	0.13	0.13	2.83	0.11	0.11	1.98	0.92	1.07	0.13	0.10	1.03	0.92	0.13	1.04	0.11	0.92	0.14	0.89	0.18	0.13	1.06	0.10	1.04	1.07	0.12	2.78	0.13	0.12	1.92	0.08	0.17	1.02	0.01	0.16	0.94	0.92	0.10	0.14	1.72	0.10	0.12	1.98	0.85	0.10	1.13	0.05	0.11	1.04	0.08	0.85	0.15	0.12	1.06	0.08	1.10	1.01	0.11	1.02	0.14	1.02	0.10	0.90	0.15	2.01	2.06	1.00	0.09	0.83	0.00	0.89	2.06	0.11	0.09	1.00	0.05	0.97	1.09	0.08	0.10	0.97	0.04	0.12	1.02	0.08	0.11	1.01	0.09	0.12	1.17	0.13	0.09	1.02	0.10	1.01	0.06	0.09	1.04	0.03	0.11	1.98	0.05	0.03	1.06	0.10	0.10	1.94	0.02	2.05	0.11	0.00	2.18	1.05	2.29	0.06	0.87	0.01	0.00	2.08	1.14	2.23	1.00	2.03	0.16	1.97	0.96	1.88	0.17	0.10	1.81	0.92	0.11	1.12	0.75	0.13	1.06	0.10	0.09	1.05	0.13	1.08	2.08	0.06	0.10	1.00	0.00	0.06	1.01	0.13	0.10	0.91	0.06	0.15	1.12	1.79	1.01	0.14	0.01	0.90	0.00	0.13	1.05	0.03	2.13	1.10	0.02	1.88	0.06	0.09	1.02	0.13	0.12	1.07	0.05	2.09	0.02	0.08	2.07	0.04	1.06	0.08	0.68	0.17	0.89	1.04	0.08	0.12	1.02	0.07	0.96	0.20	0.04	1.01	0.99	0.06	1.10	0.03	0.05	1.15	0.12	0.12	0.86	0.00	1.12	0.20	1.97	2.18	0.09	0.13	0.81	0.07	0.12	1.02	0.07	1.12	0.14	3.15	0.00	1.02	0.11	0.03	3.24	0.85	0.25	0.97	0.00	1.04	0.05	1.08	0.93	0.14	1.17	0.06	0.82	1.14	0.07	0.04	0.97	0.11	0.02	1.03	0.03	0.07	1.04	0.11	0.12	1.07	0.90	1.02	0.12	1.10	0.07	0.85	1.28	0.08	0.89	0.85	0.00	0.13	0.96	1.95	0.04	0.15	1.10	1.85	1.12	0.06	0.00	1.02	0.10	2.27	2.18	0.85	0.93	0.00	0.10	2.05	1.01	0.03	1.07	0.00	0.07	1.16	0.05	0.01	1.11	0.02	0.16	1.06	0.00	1.07	0.04	0.73	0.10	0.00	1.10	2.94	0.16	2.27	0.00	0.05	1.10	0.00	2.15	1.09	0.05	0.12	0.73	2.08	1.13	0.02	0.07	0.93	0.07	3.33	0.12	0.15	1.00	0.02	1.14	1.05	0.04	0.92	0.14	1.01	0.10	1.05	0.03	2.00	0.11	0.14	3.02	0.01	0.48	1.03	0.13
Flow Indexes:	1	3	6	8	10	12	15	18	19	20	22	24	27	27	29	30	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	68	71	71	72	73	76	77	79	81	83	86	88	89	91	91	91	94	94	97	100	101	104	104	107	107	108	110	113	115	118	120	121	123	125	127	129	129	130	130	131	133	135	136	136	139	141	142	145	148	151	154	157	159	162	165	165	168	171	171	173	173	176	176	177	178	178	180	183	183	184	185	185	186	187	187	189	189	190	191	191	194	194	195	197	198	200	203	205	206	206	209	212	215	218	219	219	220	223	226	228	228	229	231	231	234	237	239	239	242	242	244	246	248	249	252	254	257	258	260	263	266	268	270	270	271	271	274	277	279	281	281	281	283	286	286	286	287	289	291	293	294	296	298	299	302	305	308	311	312	313	315	317	318	320	321	324	325	325	328	329	329	330	333	335	335	336	336	337	338	341	341	342	344	347	350	353	355	357	360	361	361	361	363	363	366	368	368	369	372	373	373	374	377	379	379	379	382	384	385	387	389	391	393	393	396	396	396	399
Bases:	tcagAGCACGAGCCTACATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGGc
Quality Scores:	36	35	35	35	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	33	33	26	21	21	31	36	37	37	37	36	36	33	34	34	36	37	36	37	37	37	37	37	37	37	37	37	37	37	36	28	28	28	36	36	37	37	35	35	35	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	32	32	32	35	35	37	35	32	32	32	37	37	37	37	37	37	36	36	36	36	36	36	36	36	37	37	37	37	37	37	37	35	35	35	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	35	35	35	37	37	37	37	37	37	37	37	37	37	36	36	36	37	36	35	35	35	37	28	28	28	32	35	37	37	37	36	36	36	37	37	37	37	37	37	35	35	35	35	35	37	37	37	37	36	36	36	37	28	28	28	28	35	36	37	37	37	37	37	37	37	37	37	37	36	33	33	32	31	36	36	33	33	27	27	27	36	31	25	25	25	32	36	36	36	36	36	36	36	36	36	36	36

>FLP3FBN01DEHK3
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1278_0245

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       256
  Clip Qual Left:   5
  Clip Qual Right:  256
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.01	0.02	0.93	0.03	0.02	1.07	0.05	1.06	0.04	1.05	1.74	0.11	0.07	1.02	0.10	0.92	0.06	0.13	0.83	0.90	0.08	1.02	1.01	0.11	0.97	1.12	0.14	0.90	0.12	0.09	1.00	0.07	0.05	1.05	0.07	0.06	1.02	0.08	0.05	1.01	0.12	0.01	0.99	0.12	0.96	0.02	0.16	1.01	0.11	0.02	2.10	0.12	1.05	0.05	2.93	1.07	1.00	1.02	0.13	1.92	0.15	0.89	0.08	1.08	1.14	0.07	1.00	0.14	0.97	0.03	0.11	2.97	0.15	0.02	1.88	1.01	1.06	0.01	0.12	1.03	0.93	0.03	1.00	0.09	0.98	0.08	1.01	0.10	0.16	1.03	0.10	1.03	1.06	0.11	2.82	0.11	0.14	1.82	0.06	0.14	0.98	0.13	0.08	1.17	0.95	0.08	0.11	2.00	0.09	0.04	2.00	1.01	1.92	0.12	1.02	0.09	0.12	1.84	1.97	0.07	0.99	0.05	1.01	0.05	1.01	0.13	1.07	0.14	0.13	1.01	0.11	1.03	0.98	0.10	1.90	1.87	0.12	0.09	1.01	0.08	0.98	1.07	1.01	0.14	1.03	0.09	0.14	1.04	0.13	0.88	0.12	0.10	1.05	0.08	1.01	1.04	1.02	0.07	1.00	1.08	0.06	0.98	1.07	0.12	1.98	0.09	0.13	2.14	1.03	0.11	0.21	0.92	0.09	1.07	0.14	1.02	0.05	0.14	2.15	0.96	2.11	1.12	2.04	0.15	0.91	0.16	0.84	0.11	0.08	1.14	1.95	0.09	0.17	1.87	0.89	0.09	1.03	1.02	0.14	0.09	1.03	0.15	1.14	0.09	0.91	1.90	0.11	0.09	0.97	0.09	1.11	0.13	0.16	0.89	0.17	1.08	0.12	0.95	1.10	0.93	0.11	0.13	1.08	0.96	0.11	1.06	0.08	0.78	0.13	0.13	2.95	0.07	0.13	1.07	0.11	0.12	1.01	0.10	1.02	0.11	3.02	0.10	0.93	0.15	0.05	1.01	0.17	1.08	0.12	0.13	1.03	1.88	0.08	1.01	0.09	0.12	2.13	0.07	0.11	1.06	0.12	0.09	0.95	0.08	1.00	0.10	2.11	0.08	3.00	0.15	0.06	2.19	0.82	1.21	0.05	0.15	2.03	0.06	1.96	2.07	0.08	0.16	0.93	0.90	0.09	0.97	0.09	0.15	1.14	0.04	0.00	1.06	0.09	0.13	1.96	1.00	0.18	0.00	1.80	0.05	0.06	3.75	0.08	0.14	1.23	0.49	1.05	0.12	0.00	0.71	0.08	0.12	2.07	0.91	0.14	0.15	1.02	0.07	0.12	1.11	0.18	0.06	1.10	1.95	1.07	0.92	0.11	0.13	2.01	0.89	0.07	0.99	0.08	0.12	1.08	0.04	0.05	1.01	0.14	1.02	1.05	0.06	0.96	1.07	3.08	0.00	1.98	0.07	0.12	2.15	1.02	0.04	0.97	0.00	0.13	1.08	2.10	0.15	0.00	1.10	0.81	0.12	5.20	0.15	0.95	0.14	0.51	0.06	1.07	0.04	0.18	1.06	0.07	3.03	0.09	2.04	0.16	0.02	1.10	0.14	0.00	0.97	0.18	1.85	2.07	0.03	0.06	0.78	0.08	0.11	1.00	0.05	1.10	0.14	0.83	0.04	0.16	1.10	1.07	1.13
Flow Indexes:	1	3	6	8	10	11	11	14	16	19	20	22	23	25	26	28	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	72	75	75	76	77	80	81	83	85	87	90	92	93	95	95	95	98	98	101	104	105	108	108	111	111	112	113	113	115	118	118	119	119	121	123	125	127	130	132	133	135	135	136	136	139	141	142	143	145	148	150	153	155	156	157	159	160	162	163	165	165	168	168	169	172	174	176	179	179	180	181	181	182	183	183	185	187	190	191	191	194	194	195	197	198	201	203	205	206	206	209	211	214	216	218	219	220	223	224	226	228	231	231	231	234	237	239	241	241	241	243	246	248	251	252	252	254	257	257	260	263	265	267	267	269	269	269	272	272	273	274	277	277	279	279	280	280	283	284	286	289	292	295	295	296	299	299	302	302	302	302	305	307	310	313	313	314	317	320	323	324	324	325	326	329	329	330	332	335	338	340	341	343	344	345	345	345	347	347	350	350	351	353	356	357	357	360	361	363	363	363	363	363	365	367	369	372	374	374	374	376	376	379	382	384	384	385	385	388	391	393	395	398	399	400
Bases:	tcagACCAGCGACTAGCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
Quality Scores:	37	37	37	37	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	38	38	39	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	33	33	33	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	33	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	28	28	28	28	33	33	33	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	36	31	31	27	27	28	28	28	27	21	31	31	36	36	36	36	36	36	36	36	36	36	36	31	31	31	31	31	31	31

>FLP3FBN01DGFYQ
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1300_1344

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       273
  Clip Qual Left:   5
  Clip Qual Right:  273
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.01	0.00	1.01	0.00	0.00	1.02	0.00	1.09	0.00	1.12	0.79	0.15	0.14	1.04	0.03	1.11	0.00	0.95	0.14	0.84	0.70	0.19	1.05	1.87	0.03	0.33	1.00	0.28	1.01	0.13	0.92	0.10	0.08	1.10	0.14	0.07	1.02	0.13	0.06	1.00	0.17	0.00	0.98	0.15	0.94	0.00	0.21	0.82	0.15	0.00	2.13	0.14	1.03	0.05	3.29	1.11	1.17	1.04	0.31	2.11	0.25	1.01	0.01	1.22	2.73	0.12	0.08	2.73	0.19	0.03	1.74	1.05	0.99	0.00	0.11	1.15	1.00	0.00	1.05	0.12	1.04	0.03	0.83	0.13	0.17	0.86	0.14	0.81	1.11	0.08	2.81	0.06	0.10	2.00	0.13	0.11	1.07	0.08	0.17	1.02	0.90	0.12	0.13	2.15	0.03	0.07	1.70	1.07	2.08	0.01	1.01	0.14	0.11	0.88	2.91	0.07	1.09	0.07	1.10	0.01	0.85	0.14	1.03	0.10	0.12	1.01	0.29	1.72	0.14	0.11	1.74	1.82	0.10	0.08	1.04	0.12	0.98	1.02	1.09	0.13	1.01	0.11	0.13	0.94	0.09	1.15	0.09	0.14	0.85	0.11	0.83	0.83	1.19	0.00	0.89	1.16	0.07	0.00	1.01	0.11	2.98	0.07	0.11	2.12	1.04	0.93	0.10	2.10	0.12	0.05	2.07	1.16	2.20	1.02	3.12	0.08	1.04	0.15	0.32	0.77	0.07	0.23	2.08	0.11	0.08	2.09	1.12	0.07	0.98	0.18	0.13	1.91	0.11	0.15	0.92	0.09	1.06	2.09	0.11	0.15	0.82	0.08	1.01	0.06	0.02	1.01	0.07	0.95	0.08	0.84	1.17	1.05	0.10	0.03	1.13	2.70	1.06	0.10	1.82	0.00	0.16	1.01	0.15	0.13	0.97	0.07	1.06	0.07	1.00	0.08	1.04	0.10	0.18	1.08	0.94	0.10	0.11	1.15	1.85	1.01	0.17	0.76	0.04	0.14	2.24	2.13	0.20	0.11	0.73	0.20	1.08	0.06	4.64	0.19	1.98	0.22	0.00	2.06	2.09	0.15	0.21	1.10	0.16	0.90	0.24	0.95	0.17	1.00	0.17	0.16	1.09	0.09	0.15	1.05	0.15	0.07	1.33	0.08	0.08	0.96	0.10	0.23	0.84	0.98	0.12	0.11	3.50	0.09	1.01	0.11	0.01	1.07	2.19	0.20	0.12	1.88	0.02	0.05	1.05	0.01	2.09	1.01	0.11	0.12	1.18	0.06	0.01	0.88	0.39	0.05	1.11	1.93	0.92	0.98	0.07	0.13	2.04	0.98	0.08	0.92	0.06	0.09	0.96	1.62	1.20	0.15	1.01	1.05	2.95	0.01	1.68	0.15	0.15	2.03	0.93	0.15	1.04	0.09	0.21	1.04	2.03	1.06	0.00	0.19	1.08	0.02	4.70	0.14	1.00	0.10	0.08	1.03	0.88	0.10	0.19	0.79	1.11	0.14	0.12	1.12	0.17	1.10	0.13	1.93	0.14	0.13	1.21	0.00	0.07	1.06	0.10	2.14	2.23	0.76	2.99	0.03	0.01	0.70	0.78	1.00	0.09	0.13	1.22	1.18	2.00	0.96	1.12	0.07	0.88	0.10	1.09	0.05	0.08	0.98	2.85	1.05	1.06	0.07	2.26	1.04
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	68	71	71	72	73	76	77	79	81	83	86	88	89	91	91	91	94	94	97	100	101	104	104	107	107	108	109	109	111	114	115	115	115	117	119	121	123	126	128	128	131	131	132	132	135	137	138	139	141	144	146	149	151	152	153	155	156	159	161	161	161	164	164	165	166	168	168	171	171	172	173	173	174	175	175	175	177	180	183	183	186	186	187	189	192	192	195	197	198	198	201	203	206	208	210	211	212	215	216	216	216	217	219	219	222	225	227	229	231	234	235	238	239	239	240	242	245	245	246	246	249	251	253	253	253	253	253	255	255	258	258	259	259	262	264	266	268	271	274	277	280	283	284	287	287	287	287	289	292	293	293	296	296	299	301	301	302	305	308	311	312	312	313	314	317	317	318	320	323	324	324	325	327	328	329	329	329	331	331	334	334	335	337	340	341	341	342	345	347	347	347	347	347	349	352	353	356	357	360	362	364	364	367	370	372	372	373	373	374	375	375	375	378	379	380	383	384	385	385	386	387	389	391	394	395	395	395	396	397	399	399	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
Quality Scores:	35	35	35	35	35	35	35	35	33	31	31	31	33	35	35	35	35	35	35	35	35	35	35	35	35	35	23	20	20	31	31	33	33	33	35	23	17	17	21	20	20	20	31	31	33	35	35	35	35	35	33	33	33	35	31	31	31	35	35	35	35	35	35	35	31	31	31	33	35	35	35	35	35	35	35	35	35	35	31	31	31	26	26	26	26	35	35	35	35	35	35	35	33	31	31	31	35	35	35	35	35	35	35	35	35	35	35	35	35	35	31	31	31	35	35	35	33	33	33	33	33	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	30	26	26	26	30	33	35	35	35	35	35	35	35	35	33	33	33	35	33	27	27	25	25	25	27	14	14	14	14	14	25	25	34	34	35	35	35	32	33	33	32	35	35	32	25	25	15	20	20	20	28	35	33	33	33	33	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	29	24	24	24	29	35	35	35	35	33	33	31	31	34	34	34	34	34	34	31	20	20	20	20	20	31	34	31	31	31	31	32	31	31	33	34	25	25	20	20	18	25	28	28	22	20	22	28	28	28	30	30	29	29	29	30	25	25	25	29	29	26	26	25

>FLP3FBN01A55WZ
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0363_0481

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       276
  Clip Qual Left:   5
  Clip Qual Right:  276
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.00	0.00	1.01	0.02	0.01	1.04	0.04	1.05	0.01	1.16	0.79	1.82	0.90	0.11	0.09	0.98	0.11	1.02	0.09	0.97	0.97	0.07	0.13	0.96	1.05	0.14	1.92	0.05	0.12	1.12	0.17	0.03	1.00	0.16	0.06	0.93	0.15	0.03	0.92	0.11	1.00	0.01	0.16	0.94	0.16	0.00	1.97	0.14	0.94	0.02	3.19	0.97	1.07	1.01	0.36	2.12	0.15	1.10	0.12	1.14	2.76	0.12	0.08	2.04	0.22	0.78	1.74	1.14	0.93	0.06	0.16	1.09	0.99	0.01	1.01	0.15	1.05	0.03	1.01	0.11	0.14	0.91	0.17	0.90	1.94	0.12	1.97	0.13	0.12	2.04	0.22	0.07	0.97	0.16	0.16	1.01	1.06	0.14	0.11	4.82	0.10	1.01	1.77	0.25	1.95	0.09	1.85	0.16	0.95	0.02	1.05	0.12	1.03	0.07	1.15	0.08	0.08	1.01	0.20	0.89	0.15	1.92	3.80	0.12	1.05	0.97	0.12	0.10	0.92	0.09	2.04	0.18	0.08	1.02	0.14	0.10	1.07	0.12	1.09	1.04	0.10	2.00	0.14	0.99	0.14	0.93	0.95	0.12	1.06	0.92	0.11	2.07	1.07	0.12	0.10	3.06	0.12	0.08	2.07	0.93	2.17	1.04	3.75	1.12	0.12	0.06	1.83	0.18	1.00	1.02	1.08	0.14	1.01	1.03	0.20	0.14	0.94	0.16	1.33	0.06	1.05	2.20	0.21	0.09	0.86	0.08	0.16	2.10	0.06	2.06	1.14	0.78	0.11	0.13	0.97	0.02	0.12	1.12	0.04	0.07	1.07	0.10	3.77	0.12	0.07	1.08	0.08	0.10	0.98	0.12	1.17	0.11	2.06	1.04	0.14	0.16	0.76	1.04	2.14	2.18	0.95	2.92	2.05	0.06	2.92	1.94	0.12	0.13	0.77	0.09	1.11	0.09	0.14	1.11	0.13	0.10	1.08	0.08	0.18	1.08	0.14	1.04	0.16	0.80	0.12	3.60	0.21	0.06	1.06	0.00	0.12	1.00	0.09	0.14	1.04	0.03	1.09	0.09	1.00	0.12	0.16	0.99	0.88	0.12	1.03	0.16	0.85	0.21	0.15	1.07	0.08	0.11	1.12	0.16	0.12	1.21	0.13	0.90	0.13	0.02	1.08	0.93	1.96	0.13	0.06	0.87	0.13	0.15	1.12	0.02	1.03	0.11	2.06	0.11	0.12	1.08	0.87	1.01	0.13	0.05	2.01	1.99	0.14	0.12	1.12	0.00	1.05	0.17	1.02	0.11	2.97	0.13	3.08	0.15	0.63	0.17	0.08	0.89	0.06	1.03	0.17	1.80	0.15	1.79	0.15	0.09	1.07	0.00	1.07	1.00	0.00	0.16	1.09	0.02	2.04	0.10	0.09	1.86	0.05	0.99	0.20	0.60	0.14	0.99	1.12	0.75	0.17	0.00	0.99	1.92	0.18	1.98	0.14	0.00	1.09	0.00	0.15	1.11	0.13	1.89	2.10	0.13	0.18	1.79	0.11	1.02	0.19	0.13	1.11	0.75	1.14	0.91	0.04	0.12	0.82	1.01	2.03	1.01	1.12	0.06	0.87	0.14	1.07	0.15	0.11	1.06	2.79	1.08	1.11	0.03	0.14	0.94	0.13	0.11	1.04	0.79	0.23	0.14	2.21	2.05
Flow Indexes:	1	3	6	8	10	11	12	12	13	16	18	20	21	24	25	27	27	30	33	36	39	41	44	47	47	49	51	51	51	52	53	54	56	56	58	60	61	61	61	64	64	66	67	67	68	69	72	73	75	77	79	82	84	85	85	87	87	90	90	93	96	97	100	100	100	100	100	102	103	103	105	105	107	107	109	111	113	115	118	120	122	122	123	123	123	123	125	126	129	131	131	134	137	139	140	142	142	144	146	147	149	150	152	152	153	156	156	156	159	159	160	161	161	162	163	163	163	163	164	167	167	169	170	171	173	174	177	179	181	182	182	185	188	188	190	190	191	192	195	198	201	203	203	203	203	206	209	211	213	213	214	217	218	219	219	220	220	221	222	222	222	223	223	225	225	225	226	226	229	231	234	237	240	242	244	246	246	246	246	249	252	255	257	259	262	263	265	267	270	273	276	278	281	282	283	283	286	289	291	293	293	296	297	298	301	301	302	302	305	307	309	311	311	311	313	313	313	315	318	320	322	322	324	324	327	329	330	333	335	335	338	338	340	342	344	345	346	349	350	350	352	352	355	358	360	360	361	361	364	364	366	369	370	371	372	375	376	377	377	378	379	381	383	386	387	387	387	388	389	392	395	396	399	399	400	400
Bases:	tcagACGGTGAGTGTCCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
Quality Scores:	35	35	35	35	32	32	32	32	35	35	35	35	35	35	35	35	35	35	35	35	35	35	38	38	39	39	32	32	32	35	35	35	35	35	34	31	21	21	25	35	32	25	25	25	32	35	35	37	39	35	35	35	35	35	35	35	35	35	35	35	35	35	32	32	32	32	32	35	32	32	32	32	35	35	35	35	35	35	35	35	35	35	32	32	32	32	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	34	34	26	26	26	26	32	35	35	35	35	35	35	34	34	34	34	34	34	34	35	35	35	35	35	35	35	35	35	26	26	26	26	35	35	35	35	35	35	35	35	34	34	34	34	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	30	22	24	21	21	21	30	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	35	29	29	26	34	35	27	27	27	27	35	35	35	35	35	35	35	32	32	32	32	32	35	35	35	35	35	35	35	35	35	35	31	32	32	25	28	25	25	25	25	30	30	30	30	30	30	30	30	28	22	22	22	28	28	30	25	22	22	22	30	30

>FLP3FBN01D7O1S
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1610_3070

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       287
  Clip Qual Left:   5
  Clip Qual Right:  283
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.06	0.11	1.06	0.15	0.19	0.84	0.20	1.04	0.16	0.89	1.02	0.13	0.13	1.01	0.17	1.05	0.13	0.97	0.16	1.02	1.07	0.05	1.05	2.03	0.14	0.03	1.01	0.13	1.10	0.08	0.88	0.13	0.14	0.94	0.17	0.14	0.94	0.12	0.17	1.01	0.14	0.12	1.02	0.13	1.00	0.10	0.18	1.04	0.13	0.10	1.83	0.12	1.02	0.10	2.62	1.01	1.06	1.07	0.18	2.06	0.08	0.96	0.19	0.91	3.07	0.14	0.17	2.07	0.10	1.08	2.03	0.86	0.97	0.14	0.15	1.05	1.03	0.10	0.99	0.08	1.02	0.10	1.03	0.12	0.14	1.05	0.17	1.05	2.09	0.12	1.95	0.09	0.08	1.92	0.12	0.11	1.00	0.05	0.16	0.97	1.03	0.12	0.15	4.98	0.08	1.09	2.04	0.00	2.12	0.14	1.78	0.14	0.82	0.11	0.92	0.13	1.02	0.10	0.85	0.09	0.06	1.09	0.06	1.19	0.10	1.97	3.89	0.10	1.04	0.94	0.11	0.06	0.96	0.12	1.90	0.08	0.08	1.06	0.13	0.14	0.96	0.12	1.08	1.05	0.10	2.11	0.11	1.03	0.07	0.90	1.07	0.02	0.90	1.07	0.13	2.01	1.04	0.08	0.11	2.96	0.05	0.13	1.89	0.81	2.07	1.19	3.63	1.09	0.10	0.09	1.74	0.06	1.06	1.03	0.98	0.14	0.95	0.96	0.11	0.12	0.85	0.09	1.08	0.11	1.01	2.00	0.12	0.16	1.09	0.09	0.16	2.03	0.11	1.91	0.92	0.98	0.14	0.04	1.05	0.16	0.06	0.98	0.01	0.14	1.08	0.14	3.78	0.15	0.08	0.99	0.10	0.08	1.07	0.05	0.99	1.03	0.11	0.96	0.04	0.08	1.05	1.85	1.89	1.04	0.04	2.95	0.12	0.13	1.84	0.00	0.98	0.14	2.72	1.87	0.16	0.06	0.98	1.02	0.17	1.01	1.06	0.00	0.08	1.06	0.07	2.77	1.87	0.14	0.08	1.03	0.23	0.03	1.02	0.10	0.10	0.99	0.14	0.00	1.02	2.09	0.07	3.01	0.13	0.11	1.00	1.04	1.08	0.16	1.04	1.02	0.16	0.10	1.09	0.97	1.03	0.12	0.85	1.06	0.09	0.15	1.01	0.09	1.15	3.06	1.06	0.88	0.16	0.12	1.85	2.20	0.14	0.14	0.99	0.09	0.91	0.15	2.73	0.13	1.06	0.08	2.94	0.00	0.97	0.93	0.19	2.82	0.20	1.96	0.05	0.06	1.07	0.17	0.89	1.08	0.13	0.06	1.14	0.04	3.51	1.01	0.04	0.85	0.26	1.01	1.76	1.05	0.07	0.03	1.10	0.11	0.98	2.06	0.15	0.14	1.10	0.05	0.00	1.07	0.17	1.80	2.00	0.13	0.16	2.24	0.11	0.95	0.12	0.04	1.03	1.01	1.14	0.92	0.93	0.04	0.16	1.19	2.10	1.05	0.87	0.17	1.05	0.17	1.02	0.09	0.07	0.89	2.79	0.89	1.11	0.15	0.22	1.05	0.17	0.12	0.86	1.03	0.07	0.08	1.83	1.88	1.05	0.08	0.92	1.08	0.15	0.17	1.85	0.17	0.07	0.89	0.07	0.08	1.09	0.06	0.90	0.16	0.29	3.32	1.12	0.22
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	146	146	148	150	151	153	154	156	156	157	160	160	160	163	163	164	165	165	166	167	167	167	167	168	171	171	173	174	175	177	178	181	183	185	186	186	189	192	192	194	194	195	196	199	202	205	207	207	207	207	210	213	215	216	218	221	222	222	223	223	224	226	226	226	229	229	231	233	233	233	234	234	237	238	240	241	244	246	246	246	247	247	250	253	256	259	260	260	262	262	262	265	266	267	269	270	273	274	275	277	278	281	283	284	284	284	285	286	289	289	290	290	293	295	297	297	297	299	301	301	301	303	304	306	306	306	308	308	311	313	314	317	319	319	319	319	320	322	324	325	325	326	329	331	332	332	335	338	340	340	341	341	344	344	346	349	350	351	352	353	356	357	357	358	359	361	363	366	367	367	367	368	369	372	375	376	379	379	380	380	381	383	384	387	387	390	393	395	398	398	398	399
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAGTGAAACCATGCGGAAATACTATACTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATCaaac
Quality Scores:	37	36	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	37	36	40	28	28	28	33	40	38	37	37	39	37	40	37	38	38	38	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	33	33	33	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	35	35	37	31	24	21	21	21	21	23	32	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	28	28	28	28	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	28	28	28	33	33	37	37	37	37	37	28	28	28	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	33	33	33	33	37	36	28	28	26	36	37	37	37	37	36	31	31	26	31	31	36	36	30	17	16	16	19	19	19	26	30	30	32	33	36	36	36	36	36	33	33	33	36	36	33	33	28	31	31	31	31	31	28	28	28	31	31	31	31	25	22	22	22	25	28	28	28	31	28	28	28	28	31	29	28	28	28	29	25	16	16	16	16	15

>FLP3FBN01CV5ZU
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1069_1480

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       278
  Clip Qual Left:   5
  Clip Qual Right:  278
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.00	0.04	1.04	0.07	0.08	1.01	0.10	1.05	0.07	1.07	1.03	0.06	0.04	1.04	0.08	1.00	0.06	1.07	0.07	1.00	1.02	0.02	0.99	2.06	0.07	0.02	1.01	0.10	1.00	0.03	1.01	0.08	0.09	1.03	0.10	0.07	0.94	0.06	0.08	1.01	0.09	0.04	1.04	0.10	1.01	0.04	0.14	0.99	0.12	0.03	2.02	0.12	1.00	0.06	3.03	1.02	1.04	1.03	0.16	2.04	0.15	0.99	0.09	1.01	3.03	0.10	0.10	2.02	1.01	0.03	2.07	1.00	1.02	0.04	0.13	1.02	0.99	0.03	1.01	0.10	0.99	0.06	1.00	0.13	0.11	1.02	0.14	1.01	1.03	1.05	2.02	0.10	0.13	1.99	0.14	0.13	1.01	0.10	0.10	1.07	0.99	0.14	0.15	6.14	2.04	2.11	2.02	0.07	0.98	0.13	0.95	0.09	0.90	0.09	0.87	0.13	0.12	1.00	0.10	1.08	1.15	0.15	4.79	0.15	1.00	1.04	0.07	0.10	1.03	0.12	0.08	0.79	0.97	1.03	0.08	0.12	1.01	0.08	1.01	1.06	0.96	0.03	1.04	1.04	0.17	0.04	2.07	0.07	2.01	0.15	0.08	2.06	0.91	2.16	0.09	0.96	0.09	0.11	2.14	1.01	2.19	1.12	2.07	0.12	1.98	1.10	2.03	0.05	0.06	2.07	2.00	0.12	0.09	0.94	0.08	1.06	0.12	0.17	0.93	0.13	1.05	2.02	0.15	0.10	1.03	0.98	1.02	0.13	0.15	1.94	1.03	1.08	0.09	0.08	1.05	0.08	0.11	0.95	0.11	0.10	0.98	0.07	0.06	1.00	0.16	0.09	2.91	0.13	0.14	1.01	0.16	0.07	1.07	0.08	1.02	0.11	0.93	0.02	0.07	1.04	1.02	2.00	1.89	0.98	0.08	0.01	1.91	0.97	0.06	2.97	1.02	0.09	2.92	0.05	1.05	0.10	0.04	2.03	2.00	0.08	0.10	0.97	1.02	0.13	0.06	2.01	0.13	1.06	0.06	0.97	0.16	1.07	0.12	0.99	0.14	0.15	1.03	0.07	0.10	0.92	0.15	0.07	2.06	0.10	1.06	0.05	1.07	0.09	0.97	0.12	1.02	0.09	0.14	1.02	0.11	0.14	0.96	1.04	0.14	1.09	1.04	0.15	0.10	0.94	2.01	0.95	0.10	0.17	1.05	0.06	0.17	1.03	0.19	0.18	1.12	2.05	1.11	1.00	0.11	0.17	1.92	0.95	0.06	0.96	1.03	1.04	2.04	1.08	0.06	0.87	0.11	0.03	3.05	0.10	1.02	0.11	1.06	0.15	0.99	0.12	0.20	1.05	0.13	0.97	1.08	0.10	0.17	1.01	1.97	1.02	0.21	0.09	1.03	0.06	4.99	0.12	1.03	0.18	0.14	1.05	2.05	0.97	1.10	0.05	0.06	0.80	0.12	2.06	1.14	1.00	0.14	1.90	2.00	0.07	0.16	1.02	0.09	0.17	1.12	0.10	0.19	0.85	0.11	0.10	1.15	0.97	0.93	1.00	0.04	0.15	1.05	1.01	2.10	1.06	1.00	1.01	0.10	0.13	1.01	0.05	0.10	0.95	3.02	1.08	1.02	0.03	0.15	0.97	0.10	0.07	1.01	0.99	0.10	0.05	1.99	1.92	0.99	0.05	1.21	0.93
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	69	71	71	72	73	76	77	79	81	83	86	88	89	90	91	91	94	94	97	100	101	104	104	104	104	104	104	105	105	106	106	107	107	109	111	113	115	118	120	121	123	123	123	123	123	125	126	129	132	133	134	137	139	140	141	143	144	147	147	149	149	152	152	153	154	154	156	159	159	160	161	161	162	163	163	165	165	166	167	167	170	170	171	171	174	176	179	181	182	182	185	186	187	190	190	191	192	195	198	201	204	207	207	207	210	213	215	217	220	221	222	222	223	223	224	227	227	228	230	230	230	231	233	233	233	235	238	238	239	239	242	243	246	246	248	250	252	254	257	260	263	263	265	267	269	271	274	277	278	280	281	284	285	285	286	289	292	295	296	296	297	298	301	301	302	304	305	306	307	307	308	310	313	313	313	315	317	319	322	324	325	328	329	329	330	333	335	335	335	335	335	337	340	341	341	342	343	346	348	348	349	350	352	352	353	353	356	359	362	365	366	367	368	371	372	373	373	374	375	376	379	382	383	383	383	384	385	388	391	392	395	395	396	396	397	399	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACCAGCTAATACAACGCATGCCCATCTGTAACCGCCGAAACTTTCAACCACAAGAGATGCCTCTCATAGTGTTATGCGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTACAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	31	31	31	31	31	31	37	37	37	37	37	37	37	37	37	37	37	38	38	33	33	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	40	40	40	40	40	40	40	40	40	40	38	38	38	38	40	38	38	38	38	38	38	38	38	38	38	38	38	38	30	30	31	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	32	28	28

>FLP3FBN01DQ783
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1423_0757

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       268
  Clip Qual Left:   5
  Clip Qual Right:  268
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.00	0.03	1.01	0.08	0.06	1.04	0.07	0.93	0.06	1.07	1.01	0.07	0.03	1.02	0.06	0.97	0.06	1.05	0.06	1.02	1.00	0.02	0.93	1.96	0.08	0.02	1.03	0.09	1.01	0.03	1.00	0.07	0.07	1.09	0.09	0.07	1.00	0.04	0.06	1.01	0.08	0.02	1.02	0.11	0.97	0.03	0.12	1.01	0.08	0.02	1.96	0.11	1.03	0.04	2.97	0.98	1.02	1.07	0.14	1.97	0.15	1.01	0.07	1.03	3.01	0.09	0.09	1.94	0.10	1.03	2.00	0.99	0.95	0.06	0.16	0.98	1.02	0.02	1.05	0.11	1.01	0.04	1.03	0.09	0.14	1.03	0.14	1.01	2.03	0.12	1.97	0.11	0.15	2.05	0.13	0.10	1.03	0.12	0.07	1.00	1.01	0.05	0.12	4.89	0.13	1.04	2.07	0.06	2.05	0.08	1.94	0.09	1.02	0.02	0.98	0.08	1.00	0.07	1.02	0.11	0.10	1.01	0.09	1.04	0.10	2.02	4.02	0.11	1.04	1.01	0.16	0.06	1.00	0.07	2.01	0.09	0.13	1.03	0.11	0.10	1.04	0.11	1.03	1.01	0.12	1.96	0.12	0.97	0.11	0.97	1.00	0.14	1.08	1.04	0.15	1.99	1.02	0.12	0.16	2.92	0.14	0.07	2.13	1.01	2.05	1.06	3.84	1.00	0.09	0.07	1.92	0.11	0.98	1.05	0.99	0.10	0.99	1.03	0.13	0.09	1.03	0.13	1.04	0.08	1.02	2.06	0.09	0.15	1.03	0.06	0.12	1.99	0.07	1.85	1.02	1.00	0.07	0.09	1.02	0.09	0.10	0.97	0.06	0.09	1.08	0.11	3.94	0.09	0.11	1.02	0.12	0.11	1.04	0.03	0.97	1.05	0.98	0.11	1.04	0.07	1.07	1.01	2.04	1.99	0.14	2.00	0.07	0.09	1.02	0.92	2.03	0.13	2.99	2.05	0.19	0.08	1.01	0.11	1.02	0.08	0.17	1.03	0.11	0.13	0.99	0.07	0.10	0.98	1.05	0.02	0.13	1.00	0.10	1.91	1.10	0.12	0.14	0.93	0.10	0.15	1.03	0.06	0.16	1.05	0.07	0.09	0.97	1.94	0.14	1.06	0.94	0.10	1.08	0.11	1.04	0.12	0.14	0.99	0.09	0.14	0.95	0.09	0.12	1.12	0.14	0.95	0.10	0.17	1.03	0.06	0.06	1.01	0.09	0.13	2.05	0.09	0.17	1.06	0.17	0.12	0.94	0.05	1.02	0.14	2.14	0.14	0.11	1.16	1.01	0.92	0.14	0.06	2.17	2.13	0.12	0.13	0.97	0.06	1.01	0.10	2.03	0.11	1.98	0.09	3.08	0.08	1.07	0.16	0.10	0.84	0.09	1.05	0.07	1.99	0.17	1.98	0.07	0.07	0.99	0.09	1.09	0.13	0.06	1.04	1.01	0.03	2.06	0.16	0.10	1.98	0.12	1.04	0.05	1.01	0.16	0.94	1.10	0.92	0.18	1.00	0.06	0.92	1.17	1.97	0.14	0.15	0.98	0.04	0.04	0.98	0.07	1.95	2.19	0.16	0.14	2.02	0.15	0.87	0.17	0.07	1.04	0.98	1.25	0.98	1.04	0.11	0.14	1.04	1.87	1.05	1.09	0.10	1.03	0.15	0.81	0.07	0.09	0.97	3.06	1.25
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	146	146	148	150	151	153	154	156	156	157	160	160	160	163	163	164	165	165	166	167	167	167	167	168	171	171	173	174	175	177	178	181	183	185	186	186	189	192	192	194	194	195	196	199	202	205	207	207	207	207	210	213	215	216	217	219	221	222	223	223	224	224	226	226	229	230	231	231	233	233	233	234	234	237	239	242	245	248	249	252	254	254	255	258	261	264	267	268	268	270	271	273	275	278	281	284	286	289	292	295	295	298	301	303	305	305	308	309	310	313	313	314	314	317	319	321	321	323	323	325	325	325	327	330	332	334	334	336	336	339	341	344	345	347	347	350	350	352	354	356	357	358	360	362	363	364	364	367	370	372	372	373	373	376	376	378	381	382	383	384	385	388	389	389	390	391	393	395	398	399	399	399	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	37	37	38	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	38	38	38	38	40	38	38	38	38	38	38	38	38	38	36	38	37	37	36	36	36	30	31	31	31	31	31	31	31	31	31	32	31	31	31	31	31	25	25	25

>FLP3FBN01DDPFF
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1269_0617

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       270
  Clip Qual Left:   5
  Clip Qual Right:  263
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.02	0.02	1.03	0.08	0.05	1.00	0.06	1.03	0.05	1.05	1.00	0.05	0.02	1.03	0.06	1.00	0.04	1.07	0.03	1.01	1.00	0.00	0.92	2.04	0.04	0.01	0.95	0.08	1.06	0.00	0.97	0.08	0.05	1.06	0.08	0.07	0.98	0.03	0.07	0.99	0.07	0.01	1.04	0.11	0.99	0.01	0.12	0.99	0.10	0.02	1.98	0.13	1.03	0.04	2.96	1.01	1.02	1.05	0.16	2.02	0.15	1.02	0.06	1.03	3.10	0.10	0.11	2.04	0.09	1.01	2.03	0.95	1.03	0.04	0.11	1.01	1.04	0.01	1.07	0.10	1.02	0.05	0.96	0.11	0.13	0.99	0.13	0.99	2.03	0.09	1.97	0.10	0.13	1.97	0.10	0.12	1.06	0.09	0.08	1.06	0.98	0.07	0.14	4.80	0.14	1.00	2.05	0.05	2.03	0.05	1.86	0.10	1.02	0.02	1.01	0.07	1.02	0.07	0.99	0.14	0.10	1.06	0.10	1.03	0.09	2.05	4.03	0.12	1.06	0.92	0.12	0.06	0.97	0.10	2.06	0.04	0.12	1.08	0.12	0.12	1.01	0.10	1.06	1.04	0.08	2.10	0.14	0.94	0.10	1.01	1.03	0.16	1.01	1.04	0.14	2.03	0.92	0.11	0.11	2.88	0.17	0.09	2.10	1.06	1.97	1.05	3.84	1.04	0.08	0.13	1.99	0.15	1.06	1.03	1.07	0.09	0.98	1.00	0.16	0.07	0.86	0.13	1.00	0.08	1.07	2.07	0.15	0.12	0.99	0.13	0.15	2.04	0.17	2.02	1.04	0.95	0.09	0.12	1.03	0.19	0.08	0.98	0.07	0.07	1.02	0.09	4.00	0.12	0.12	1.02	0.16	0.13	1.03	0.15	1.02	1.05	1.05	0.10	1.01	0.20	1.00	1.03	2.12	2.02	0.17	2.06	0.07	0.04	1.02	0.99	2.06	0.13	3.02	2.02	0.18	0.07	1.02	0.10	1.04	0.09	0.13	1.02	0.07	0.12	1.05	0.04	0.11	0.99	1.02	0.15	0.12	0.97	0.14	1.95	1.09	0.07	0.13	0.89	0.08	0.15	1.04	0.05	0.14	1.04	0.06	0.12	0.99	1.98	0.12	1.08	1.00	0.11	1.01	0.10	1.06	0.13	0.13	0.96	0.05	0.15	0.99	0.11	0.11	1.08	0.15	0.93	0.08	0.17	1.00	0.16	0.07	1.05	0.06	0.12	2.08	0.07	0.14	1.05	0.14	0.13	1.11	0.05	0.95	0.16	2.05	0.17	0.19	1.03	0.94	0.94	0.12	0.09	1.94	2.07	0.17	0.07	1.06	0.05	1.09	0.05	1.01	0.09	2.90	0.10	2.76	0.14	0.94	0.15	0.08	0.99	0.04	1.03	0.15	1.91	0.12	1.84	0.16	0.14	1.08	0.07	1.14	0.97	0.03	0.15	0.93	0.02	2.03	0.14	0.05	1.99	0.08	0.91	0.16	1.12	0.11	1.04	1.15	0.93	0.09	0.04	0.92	2.11	0.15	2.00	0.10	0.14	1.00	0.05	0.03	0.95	0.09	1.89	1.88	0.14	0.17	2.74	1.11	1.03	1.11	1.01	0.02	0.10	0.92	1.05	1.96	1.10	1.01	0.11	0.98	0.13	0.86	0.22	0.16	0.95	2.88	0.92	1.02	0.00	0.23	1.26
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	146	146	148	150	151	153	154	156	156	157	160	160	160	163	163	164	165	165	166	167	167	167	167	168	171	171	173	174	175	177	178	181	183	185	186	186	189	192	192	194	194	195	196	199	202	205	207	207	207	207	210	213	215	216	217	219	221	222	223	223	224	224	226	226	229	230	231	231	233	233	233	234	234	237	239	242	245	248	249	252	254	254	255	258	261	264	267	268	268	270	271	273	275	278	281	284	286	289	292	295	295	298	301	303	305	305	308	309	310	313	313	314	314	317	319	321	323	323	323	325	325	325	327	330	332	334	334	336	336	339	341	342	345	347	347	350	350	352	354	356	357	358	361	362	362	364	364	367	370	372	372	373	373	376	376	376	377	378	379	380	383	384	385	385	386	387	389	391	394	395	395	395	396	397	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGGTACGCGTTACTCacccgtg
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	34	33	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	31	31	31	38	37	37	36	36	33	33	33	36	37	36	38	38	38	37	38	38	38	38	38	38	36	36	36	36	38	36	36	36	28	28	29	21	21	22	28	32	31	31	31	31	31	32	31	31	28	28	28	31	31	31	30	30

>FLP3FBN01CPD70
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0992_0942

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       280
  Clip Qual Left:   5
  Clip Qual Right:  278
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.02	0.01	1.01	0.07	0.04	1.03	0.09	1.04	0.05	1.10	0.92	0.08	0.06	1.02	0.08	0.99	0.04	1.05	0.07	0.93	1.01	0.04	1.02	1.98	0.09	0.05	1.02	0.07	1.01	0.04	0.99	0.06	0.07	1.07	0.09	0.07	0.99	0.06	0.09	0.99	0.09	0.02	1.00	0.10	0.91	0.02	0.16	1.00	0.10	0.02	2.06	0.12	1.09	0.07	2.96	1.06	1.01	1.08	0.14	2.10	0.12	1.01	0.10	0.98	0.91	0.08	1.01	0.14	1.01	0.05	0.16	2.03	0.13	1.02	2.04	1.06	1.03	0.06	0.17	0.96	1.01	0.03	1.03	0.12	1.03	0.03	1.00	0.09	0.14	1.04	0.15	1.01	2.08	0.10	1.95	0.13	0.14	1.93	0.12	0.10	1.01	0.05	0.10	1.07	1.01	0.04	0.14	5.81	0.08	0.02	1.98	0.07	2.02	0.03	2.02	0.12	0.94	0.04	0.97	0.08	1.05	0.07	1.02	0.10	0.12	1.08	0.15	1.02	0.09	1.96	3.98	0.14	1.06	0.91	0.12	0.08	0.97	0.10	2.04	0.09	0.14	0.98	0.11	0.13	1.06	0.12	1.03	1.05	0.12	1.93	0.12	2.04	0.13	0.02	1.02	0.12	2.07	0.04	0.17	2.11	1.03	0.04	0.10	3.00	0.13	0.06	2.08	0.90	2.13	1.02	3.61	1.04	0.03	0.13	1.93	0.10	0.07	2.01	1.06	0.13	0.10	1.89	2.06	0.11	1.03	2.00	0.14	0.13	0.97	0.05	0.14	1.99	0.06	1.94	1.03	0.97	0.08	0.14	0.96	0.11	0.10	1.01	0.09	0.10	1.08	0.14	3.87	0.11	0.10	1.00	0.14	0.14	1.03	0.12	1.01	1.13	0.14	1.01	0.14	0.11	1.02	0.11	0.11	0.94	0.07	1.00	1.99	1.03	0.08	2.10	0.15	0.98	2.12	0.02	1.04	0.08	3.01	2.07	0.15	0.15	0.88	1.02	0.19	1.13	2.06	0.12	1.00	0.09	1.08	1.03	1.98	0.11	0.13	1.00	0.11	0.14	1.00	0.08	0.16	0.96	0.06	0.05	1.03	2.03	0.12	2.02	0.13	1.02	0.14	1.84	1.06	0.05	1.11	0.88	0.13	0.17	0.92	0.06	0.06	1.08	0.08	0.14	2.12	0.09	0.10	1.01	0.12	0.09	0.97	0.07	1.06	3.06	1.07	0.97	0.10	0.05	2.06	2.04	0.15	0.09	0.95	0.05	1.00	0.16	2.99	0.14	1.06	0.15	2.91	0.12	0.99	1.00	0.11	2.96	0.19	1.93	0.11	0.13	1.06	0.06	0.98	1.00	0.03	0.07	1.04	0.11	4.04	1.05	0.15	0.95	0.08	0.90	1.11	0.11	1.01	0.11	0.14	1.01	0.20	0.04	1.09	0.05	1.07	2.00	0.11	0.14	0.99	0.01	0.05	1.03	0.10	1.99	2.10	0.13	0.17	1.87	0.12	0.95	0.17	0.13	1.06	1.02	1.08	1.09	1.01	0.09	0.10	0.86	2.06	1.06	1.17	0.08	0.88	0.08	0.80	0.07	0.12	1.07	3.08	1.10	1.05	0.02	0.08	1.01	0.06	0.09	1.02	0.91	0.05	0.06	2.13	2.01	1.17	0.07	0.95	0.86	0.04	0.12	1.93	0.06
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	74	75	75	76	77	80	81	83	85	87	90	92	93	93	95	95	98	98	101	104	105	108	108	108	108	108	108	111	111	113	113	115	115	117	119	121	123	126	128	130	130	131	131	131	131	133	134	137	139	139	142	145	147	148	150	150	152	152	155	157	157	160	160	161	164	164	164	167	167	168	169	169	170	171	171	171	171	172	175	175	178	178	179	182	182	183	183	185	186	186	189	192	192	194	194	195	196	199	202	205	207	207	207	207	210	213	215	216	218	221	224	226	227	227	228	230	230	232	233	233	235	237	237	237	238	238	241	242	244	245	245	247	249	250	251	251	254	257	260	263	264	264	266	266	268	270	270	271	273	274	277	280	283	283	286	289	291	292	292	292	293	294	297	297	298	298	301	303	305	305	305	307	309	309	309	311	312	314	314	314	316	316	319	321	322	325	327	327	327	327	328	330	332	333	335	338	341	343	344	344	347	350	352	352	353	353	356	356	358	361	362	363	364	365	368	369	369	370	371	373	375	378	379	379	379	380	381	384	387	388	391	391	392	392	393	395	396	399	399
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTCTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGcc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	32	32	32	32	32	32	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	34	34	28	28	28	26	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	40	40	40	40	40	40	40	38	38	38	38	38	38	38	40	40	40	40	40	40	40	38	38	38	38	31	31	32	32	32	31	31	31	31	31	30	30	30	31	31	32	32	32	32	32	32	32	32	32	32	32	31	31	31	31	31

>FLP3FBN01BBAE6
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0421_2032

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       273
  Clip Qual Left:   5
  Clip Qual Right:  273
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.02	0.03	1.02	0.05	0.04	1.00	0.08	1.03	0.05	1.13	0.20	1.29	0.17	0.02	0.86	0.15	0.05	1.15	0.06	1.09	0.16	0.04	0.73	0.19	0.18	0.94	1.33	0.04	1.82	0.19	0.04	0.81	1.17	0.00	0.89	0.19	0.05	0.90	0.22	0.05	1.08	0.04	0.08	1.09	0.06	0.05	1.35	0.07	1.18	0.04	0.46	0.88	0.16	0.01	2.11	0.18	0.91	0.00	3.19	1.00	1.05	1.04	0.14	2.12	0.04	1.03	0.11	0.98	1.05	0.06	0.92	0.07	0.92	0.05	0.14	2.11	0.99	0.06	2.02	1.05	0.97	0.04	0.12	1.02	1.06	0.01	1.02	0.08	0.98	0.08	0.90	0.10	0.11	0.99	0.18	1.01	1.04	1.06	2.08	0.13	0.08	1.07	0.16	0.96	1.04	0.12	0.07	0.97	1.05	0.10	0.10	5.03	0.15	1.03	1.95	0.08	1.92	0.08	1.83	0.11	1.02	0.06	1.03	0.07	1.03	0.08	1.05	0.10	0.13	1.04	0.09	0.94	0.12	2.02	4.08	0.13	0.93	0.96	1.03	1.01	0.10	0.16	1.07	0.13	0.10	1.03	0.18	0.11	1.02	0.11	0.95	1.04	0.93	0.05	1.04	0.87	0.10	0.06	1.97	0.06	0.99	0.06	1.07	1.98	0.94	0.06	0.15	2.73	0.08	0.09	1.84	1.03	1.92	1.01	3.62	1.01	0.08	0.13	1.83	1.07	0.13	1.00	1.04	0.10	1.04	0.88	0.19	1.10	0.16	0.14	1.02	0.09	0.88	1.86	0.13	0.15	0.95	0.04	0.13	1.02	0.09	0.05	1.02	1.00	0.11	0.05	1.93	0.98	0.12	0.12	0.91	0.12	0.15	0.99	0.15	0.10	1.07	0.13	0.09	1.03	0.14	0.09	2.91	0.15	0.15	1.02	0.06	0.07	1.05	0.12	1.86	1.12	0.11	0.11	1.93	0.04	0.15	0.91	1.97	2.01	1.11	1.83	0.09	0.11	0.77	0.04	3.01	0.11	3.04	0.10	0.20	1.86	0.08	0.09	1.12	2.08	0.17	0.14	1.08	0.10	0.09	1.03	2.03	2.99	0.13	0.85	0.20	0.10	1.02	0.08	0.13	1.03	0.12	0.12	4.02	1.11	0.17	1.02	0.16	0.17	0.93	0.14	0.15	1.03	0.09	0.08	2.07	0.97	0.09	0.06	0.92	0.96	1.03	0.14	1.05	0.05	0.06	1.02	1.04	1.02	0.06	0.12	1.14	1.93	0.92	0.97	0.18	0.17	1.81	0.92	0.14	0.94	0.10	0.71	1.15	1.91	0.06	1.85	0.09	0.07	2.72	0.04	0.96	0.09	1.75	0.12	1.95	1.10	0.00	0.15	1.05	0.08	2.08	1.11	0.15	0.13	0.75	0.17	4.94	0.16	0.94	0.09	0.13	1.02	0.97	0.08	0.17	1.86	0.09	0.07	1.05	2.84	0.12	0.14	0.98	0.01	0.07	1.08	0.23	1.90	2.02	0.07	0.17	0.90	0.04	0.13	1.16	0.08	0.17	1.01	0.12	0.11	1.07	0.92	1.12	0.99	0.83	0.12	0.19	1.05	1.96	1.04	1.09	0.08	0.97	0.10	1.04	0.09	0.14	0.98	2.90	0.99	1.04	0.15	0.14	0.89	0.04	0.08	1.39	1.37
Flow Indexes:	1	3	6	8	10	12	15	18	20	23	26	27	29	29	32	33	35	38	41	44	47	49	52	55	55	57	59	59	59	60	61	62	64	64	66	68	69	71	73	76	76	77	79	79	80	81	84	85	87	89	91	94	96	97	98	99	99	102	104	105	108	109	112	112	112	112	112	114	115	115	117	117	119	119	121	123	125	127	130	132	134	134	135	135	135	135	137	138	139	140	143	146	149	151	152	153	155	156	159	159	161	163	164	164	165	168	168	168	171	171	172	173	173	174	175	175	175	175	176	179	179	180	182	183	185	186	188	191	193	194	194	197	200	203	204	207	207	208	211	214	217	220	223	223	223	226	229	231	231	232	235	235	238	239	239	240	240	241	242	242	245	247	247	247	249	249	249	252	252	255	256	256	259	262	263	263	264	264	264	266	269	272	275	275	275	275	276	278	281	284	287	287	288	291	292	293	295	298	299	300	303	304	304	305	306	309	309	310	312	314	315	316	316	318	318	321	321	321	323	325	325	327	327	328	331	333	333	334	337	339	339	339	339	339	341	344	345	348	348	351	352	352	352	355	358	360	360	361	361	364	367	370	373	374	375	376	377	380	381	381	382	383	385	387	390	391	391	391	392	393	396	399	400
Bases:	tcagAGCAGCACTTGTCATGCTGCCTCCCGTAGGAGTCTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCG
Quality Scores:	37	37	37	36	35	35	35	36	36	32	32	32	34	34	35	35	37	37	37	36	34	32	34	36	37	37	33	33	33	37	37	38	38	38	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	28	26	26	35	35	37	37	35	31	21	21	21	26	32	35	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	35	35	35	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	35	32	32	32	37	36	36	26	26	26	36	36	33	33	33	36	36	36	33	33	33	33	33	36	36	36	36	36	36	36	36	36	33	33	33	36	33	33	33	36	36	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	31	28	22	20	19

>FLP3FBN01AWYZD
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0258_1671

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       248
  Clip Qual Left:   5
  Clip Qual Right:  248
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.05	1.03	0.09	0.08	0.96	0.09	0.99	0.08	2.02	1.05	0.07	1.03	0.03	1.02	1.01	1.02	0.01	1.03	0.97	0.08	1.03	0.07	0.08	0.98	0.02	0.07	1.02	0.07	0.03	0.95	0.08	0.08	1.05	0.07	0.08	1.01	0.04	0.09	0.99	0.12	0.05	1.01	0.10	1.01	0.05	0.11	1.00	0.10	0.04	1.99	0.11	1.02	0.07	2.93	1.01	1.03	1.08	0.13	2.06	0.15	0.99	0.09	1.01	3.06	0.07	0.10	3.05	0.12	0.02	2.01	1.04	1.02	0.04	0.13	0.94	1.04	0.06	0.98	0.12	1.00	0.07	0.96	0.08	0.10	1.03	0.17	1.02	1.05	0.12	2.97	0.13	0.11	2.01	0.09	0.09	0.97	0.10	0.05	0.95	1.03	0.01	0.14	2.05	0.10	0.06	2.08	1.00	0.08	1.07	0.10	0.10	1.03	0.12	0.88	0.11	0.13	1.01	0.10	1.06	1.03	0.07	0.98	0.13	1.01	0.03	1.00	0.07	1.01	0.10	1.01	0.06	0.12	2.01	1.03	0.10	1.04	0.05	1.01	2.01	0.11	0.07	0.99	0.09	1.02	1.01	0.14	0.11	0.94	0.10	0.12	0.97	0.12	0.07	1.05	0.10	0.08	1.07	0.16	0.06	1.04	0.11	1.02	0.04	0.06	1.01	0.12	0.09	1.97	0.11	0.08	0.96	0.08	0.05	1.93	0.06	2.01	0.05	0.06	2.02	1.09	2.00	0.07	1.03	0.05	0.10	2.00	1.04	2.01	1.04	1.87	0.11	1.97	0.93	1.93	0.11	0.13	2.02	1.10	0.12	1.03	1.03	0.14	0.98	0.11	0.15	0.97	0.16	0.97	2.06	0.13	0.13	1.04	0.09	0.14	1.01	0.10	0.10	1.03	0.14	0.09	0.95	2.02	1.06	0.08	0.12	1.03	0.08	0.12	1.02	0.07	2.07	1.06	0.11	2.05	0.12	0.13	1.01	0.10	0.11	1.00	0.04	1.96	0.16	0.07	2.02	0.12	1.09	0.15	1.02	0.11	1.01	1.03	0.13	0.07	1.04	0.13	0.97	0.08	0.14	1.01	1.01	0.14	1.06	0.08	0.10	1.09	0.17	0.06	0.95	0.07	1.04	0.07	2.01	2.12	0.07	0.13	1.04	0.12	0.10	0.97	0.11	1.01	0.06	3.01	0.03	1.00	0.09	0.14	3.08	1.09	0.10	1.02	0.08	1.03	0.13	1.02	0.99	0.20	1.03	0.06	0.92	1.04	0.07	0.13	1.01	0.11	0.15	1.02	0.09	0.17	1.03	0.18	0.08	1.06	1.02	1.01	0.14	0.99	0.17	1.01	1.01	0.10	0.97	0.95	0.05	0.13	1.04	1.94	0.14	0.12	1.09	2.05	1.13	0.11	0.07	0.94	0.15	2.02	1.96	1.03	0.99	0.14	0.17	2.00	0.99	0.18	1.05	0.09	0.13	1.00	0.08	0.17	1.02	0.20	0.14	1.08	0.05	1.01	0.15	0.98	0.08	0.05	1.06	2.91	0.10	1.99	0.09	0.07	0.90	0.09	1.97	1.15	0.14	0.17	0.98	1.89	1.06	0.15	0.10	1.06	0.13	2.99	0.11	0.08	0.91	0.10	1.04	1.03	0.06	0.87	0.16	1.08	0.09	1.00	0.06	1.89	0.08	0.16	3.11
Flow Indexes:	1	3	6	8	10	10	11	13	15	16	17	19	20	22	25	28	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	68	71	71	72	73	76	77	79	81	83	86	88	89	91	91	91	94	94	97	100	101	104	104	107	107	108	110	113	115	118	120	121	123	125	127	129	131	134	134	135	137	139	140	140	143	145	146	149	152	155	158	161	163	166	169	169	172	175	175	177	177	180	180	181	182	182	184	187	187	188	189	189	190	191	191	193	193	194	195	195	198	198	199	201	202	204	207	209	210	210	213	216	219	222	223	223	224	227	230	232	232	233	235	235	238	241	243	243	246	246	248	250	252	253	256	258	261	262	264	267	270	272	274	274	275	275	278	281	283	285	285	285	287	290	290	290	291	293	295	297	298	300	302	303	306	309	312	315	316	317	319	321	322	324	325	328	329	329	332	333	333	334	337	339	339	340	340	341	342	345	345	346	348	351	354	357	359	361	364	365	365	365	367	367	370	372	372	373	376	377	377	378	381	383	383	383	386	388	389	391	393	395	397	397	400	400	400
Bases:	tcagAACTCGTCGATGCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	38	40	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	38	37

>FLP3FBN01AM0P3
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0145_0217

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       268
  Clip Qual Left:   5
  Clip Qual Right:  268
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.02	0.01	1.03	0.00	0.00	1.00	0.04	1.02	0.00	1.14	0.86	0.08	0.07	0.87	0.05	1.01	0.09	1.06	0.05	1.02	0.89	0.04	1.02	1.76	0.13	0.11	0.95	0.13	1.04	0.04	1.01	0.04	0.09	1.03	0.08	0.02	1.01	0.07	0.05	1.01	0.07	0.00	0.99	0.10	0.97	0.00	0.12	0.89	0.12	0.01	2.03	0.13	1.02	0.05	3.03	1.09	1.02	1.03	0.17	2.13	0.11	0.99	0.07	1.13	1.05	0.07	0.74	0.19	1.04	0.02	0.11	2.93	0.12	0.00	1.77	1.18	0.97	0.00	0.14	0.97	1.00	0.02	1.00	0.08	0.93	0.09	1.01	0.08	0.11	1.04	0.19	0.85	1.04	0.12	2.71	0.12	0.11	1.78	0.11	0.09	1.00	0.10	0.06	1.23	1.00	0.04	0.12	2.20	0.16	0.04	2.07	1.73	0.21	0.04	1.75	0.72	0.11	0.05	2.78	0.09	0.85	0.08	1.08	0.09	0.94	0.12	1.21	0.06	0.11	1.03	0.24	1.82	0.07	0.15	1.95	2.02	0.10	0.10	1.02	0.13	0.98	0.95	1.08	0.08	0.98	0.15	0.11	1.03	0.13	0.79	0.10	0.10	1.06	0.08	1.00	1.02	1.06	0.02	0.92	1.16	0.15	0.01	1.95	0.21	1.90	0.05	0.20	1.89	1.02	0.07	0.20	2.79	0.04	0.14	2.02	0.96	1.83	1.13	3.79	0.97	0.11	0.11	1.97	0.08	0.10	1.91	1.92	0.12	0.07	1.10	0.09	0.96	0.10	0.05	1.08	0.08	0.93	2.08	0.11	0.10	0.95	0.12	1.09	0.10	0.05	0.92	0.16	0.95	0.05	1.04	1.03	1.02	0.04	0.10	1.06	2.97	0.07	0.14	2.82	0.17	0.13	1.05	0.22	0.10	0.96	0.16	2.90	1.06	1.06	1.05	2.01	0.03	0.17	0.96	1.95	1.98	0.15	1.04	0.11	0.96	4.77	0.15	0.99	0.06	0.05	1.03	1.05	0.11	0.18	1.01	0.97	0.12	1.13	0.11	0.10	0.99	0.09	0.11	1.03	0.12	1.97	0.17	0.97	0.12	0.10	1.01	0.13	0.11	1.01	0.13	0.13	1.01	0.06	0.06	1.01	0.95	0.12	1.87	0.08	1.00	0.07	0.08	1.01	0.12	1.11	0.16	0.10	1.01	1.00	0.15	0.08	1.01	0.19	0.13	0.91	0.87	0.17	0.15	1.06	0.15	2.10	1.06	0.14	0.11	1.08	0.11	0.09	1.04	0.16	0.08	1.04	2.00	1.07	0.84	0.15	0.11	2.15	1.01	0.16	1.06	0.16	0.17	0.99	0.06	0.15	1.01	2.16	0.21	1.10	0.92	0.08	0.18	3.00	0.06	1.07	0.09	0.85	1.96	0.18	1.09	1.09	0.14	0.20	0.94	2.03	1.08	0.16	0.16	1.11	0.12	4.70	0.11	0.05	0.93	0.04	1.06	1.12	0.86	0.17	0.02	1.02	0.98	1.09	2.21	0.12	0.16	1.02	0.00	0.04	1.04	0.11	2.07	2.11	0.08	1.07	0.01	0.73	0.08	2.01	0.18	0.16	0.92	1.06	1.07	0.14	0.15	1.00	1.07	1.86	0.89	1.08	0.13	0.95	0.14	1.06	0.08	0.18	0.82	3.25	1.25
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	72	75	75	76	77	80	81	83	85	87	90	92	93	95	95	95	98	98	101	104	105	108	108	111	111	112	112	115	115	116	119	119	119	121	123	125	127	130	132	132	135	135	136	136	139	141	142	143	145	148	150	153	155	156	157	159	160	163	163	165	165	168	168	169	172	172	172	175	175	176	177	177	178	179	179	179	179	180	183	183	186	186	187	187	190	192	195	197	198	198	201	203	206	208	210	211	212	215	216	216	216	219	219	219	222	225	227	227	227	228	229	230	231	231	234	235	235	236	236	238	240	241	241	241	241	241	243	246	247	250	251	253	256	259	261	261	263	266	269	272	275	276	278	278	280	283	285	288	289	292	295	296	299	301	301	302	305	308	311	312	312	313	314	317	317	318	320	323	326	327	327	329	330	333	333	333	335	337	338	338	340	341	344	345	345	346	349	351	351	351	351	351	354	356	357	358	361	362	363	364	364	367	370	372	372	373	373	375	377	379	379	382	383	384	387	388	389	389	390	391	393	395	398	399	399	399	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCGTACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	39	38	38	39	39	38	38	38	33	33	33	36	39	39	39	37	37	37	36	36	29	29	29	33	33	37	37	37	37	37	35	35	31	31	25	25	25	21	21	21	36	37	37	36	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	33	33	33	37	37	37	36	36	36	33	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	33	33	33	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	28	28	28	28	28	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	37	37	37	37	36	36	36	36	37	37	37	37	37	37	36	36	36	36	36	36	30	26	22	22	22	22	22	33	36	36	36	36	36	31	31	32	36	36	36	36	36	36	36	30	30	30	31	31	31	31	31	30	30	30	31	31	25	22	16	16	16	16

>FLP3FBN01BO4IE
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0579_0388

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       245
  Clip Qual Left:   5
  Clip Qual Right:  244
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.03	0.97	0.07	0.07	1.02	0.06	1.01	0.06	1.05	2.05	0.97	0.06	0.01	0.99	0.09	0.01	1.07	0.04	0.91	0.06	1.06	0.04	0.96	1.01	0.03	1.00	0.08	0.08	1.07	0.95	0.06	0.08	1.05	0.07	0.07	1.00	0.05	0.06	1.00	0.08	0.03	1.03	0.09	0.98	0.02	0.11	1.00	0.08	0.03	2.04	0.11	1.04	0.05	3.09	0.99	1.02	1.06	0.14	2.03	0.11	1.01	0.07	0.96	1.04	0.06	1.02	0.09	1.05	0.03	0.14	3.18	0.10	0.05	2.20	0.98	0.96	1.03	0.13	0.10	0.94	0.11	0.93	0.08	1.05	0.07	1.06	0.11	0.13	1.03	0.08	0.99	1.01	0.08	3.01	0.14	0.13	1.99	0.21	0.11	0.91	0.14	0.08	1.01	1.03	0.04	0.12	1.97	0.08	0.03	1.91	1.09	1.94	0.11	0.95	0.11	0.14	1.92	1.96	0.10	0.81	0.13	1.09	0.08	1.04	0.08	1.01	0.10	0.13	1.02	0.08	0.93	0.96	0.07	1.90	1.97	0.14	0.08	1.07	0.15	0.98	1.05	1.04	0.12	1.01	0.07	0.11	1.01	0.14	1.02	0.09	0.14	0.97	0.05	0.98	0.97	1.05	0.08	0.91	1.07	0.14	0.83	1.23	0.14	1.91	0.09	0.05	1.91	1.12	0.91	0.13	1.93	0.14	0.08	2.05	1.00	2.01	1.03	3.90	0.10	0.11	0.95	1.97	0.14	0.12	1.99	1.10	0.09	1.01	1.10	0.10	0.14	0.91	0.07	1.01	0.10	0.96	2.00	0.13	0.08	1.05	0.16	0.96	0.14	0.07	0.95	0.14	1.06	0.11	1.03	0.96	1.02	0.07	0.05	1.07	0.13	0.06	2.01	0.11	0.94	0.07	0.14	3.02	0.13	0.11	1.02	0.09	0.12	1.09	0.16	0.90	0.09	1.09	1.03	0.08	0.07	1.01	0.10	0.94	0.12	0.10	1.01	0.06	1.08	0.11	0.07	1.06	1.92	0.08	0.96	0.05	0.11	2.11	0.05	0.04	1.00	0.14	0.06	0.99	0.11	1.03	0.05	2.05	0.14	3.33	0.12	2.65	0.11	0.10	1.13	0.10	0.16	2.68	0.85	0.22	0.08	1.14	1.07	0.15	0.06	1.07	0.13	1.02	0.07	0.13	0.90	0.21	0.14	0.89	0.18	0.12	1.05	1.10	0.16	0.08	1.10	0.07	0.74	0.26	0.18	0.89	1.17	0.16	0.13	2.06	0.08	1.06	0.10	0.15	1.05	0.24	0.09	1.04	0.15	0.10	1.84	0.14	0.12	1.07	0.12	0.32	1.03	0.12	0.15	1.10	0.85	0.11	0.21	0.95	0.03	0.12	1.20	0.20	0.13	0.86	1.84	1.20	0.95	0.11	0.11	1.80	1.09	0.08	0.94	0.12	0.18	0.95	0.15	0.17	0.86	0.17	1.16	1.00	0.14	0.97	0.05	0.44	0.98	0.14	0.15	2.87	0.00	1.01	0.14	0.86	2.02	1.22	0.14	1.04	0.02	0.09	1.08	2.15	0.19	0.11	0.79	2.07	0.11	4.32	0.13	0.90	0.08	0.83	0.08	1.04	0.12	0.00	0.98	0.02	1.31	0.10	0.00	1.16	0.98	0.04	2.06	0.09	0.15	1.06	0.30
Flow Indexes:	1	3	6	8	10	11	11	12	15	18	20	22	24	25	27	30	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	72	75	75	76	77	78	81	83	85	87	90	92	93	95	95	95	98	98	101	104	105	108	108	111	111	112	113	113	115	118	118	119	119	121	123	125	127	130	132	133	135	135	136	136	139	141	142	143	145	148	150	153	155	156	157	159	160	162	163	165	165	168	168	169	170	172	172	175	175	176	177	177	178	179	179	179	179	182	183	183	186	186	187	189	190	193	195	197	198	198	201	203	206	208	210	211	212	215	218	218	220	223	223	223	226	229	231	233	234	237	239	242	244	247	248	248	250	253	253	256	259	261	263	263	265	265	265	267	267	267	270	273	273	273	274	277	278	281	283	286	289	292	293	296	298	301	302	305	305	307	310	313	316	316	319	322	325	326	329	332	335	336	336	337	338	341	341	342	344	347	350	352	353	355	358	361	361	361	363	365	366	366	367	369	372	373	373	376	377	377	379	379	379	379	381	383	385	388	390	393	394	396	396	399
Bases:	tcagACCGCAGAGTCACATGCTGCCTCCCGTAGGAGTCTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTCCCATTTATATCATGTGATATTCATGGCATATGCGGTATTAGCAGTCATTTCTAACTGTTGTTCCCCTCTGATAGGc
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	33	33	33	38	38	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	36	36	36	36	36	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	33	33	23	20	20	20	20	20	30	20	20	20	22	28	37	37	37	37	37	37	37	36	33	33	33	37	37	37	37	37	37	36	36	36	37	37	37	36	36	33	33	33	36	36	36	36	37	36	36	36	37	36	36	36	37	37	37	36	36	36	36	36	36	30	30	30	30	30	22	22	22	22	26	30	36	36	36	36	36	36	36	36

>FLP3FBN01BPX14
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0588_1818

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       258
  Clip Qual Left:   5
  Clip Qual Right:  258
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.01	0.00	1.00	0.03	0.03	1.03	0.04	1.01	0.04	1.12	1.85	0.63	0.18	0.02	0.60	0.22	0.12	1.02	0.14	1.01	0.11	1.20	0.15	1.02	1.00	0.10	0.91	0.12	0.24	0.77	1.05	0.02	0.02	1.14	0.05	0.05	1.04	0.03	0.02	1.00	0.05	0.00	1.01	0.06	0.97	0.01	0.07	0.99	0.06	0.02	2.09	0.06	1.02	0.03	3.04	1.01	1.06	1.05	0.11	2.04	0.10	1.07	0.03	1.08	1.29	0.14	0.67	0.24	0.84	0.03	0.14	2.12	0.87	0.00	1.91	1.08	1.07	0.00	0.13	0.95	1.06	0.02	1.12	0.09	0.98	0.10	0.92	0.11	0.13	0.96	0.24	0.89	1.02	1.07	2.14	0.13	0.10	1.02	0.12	1.01	1.19	0.06	0.05	1.25	0.91	0.01	0.17	5.57	0.15	0.01	2.05	0.11	1.84	0.05	2.11	0.07	1.00	0.01	0.97	0.05	0.96	0.05	1.00	0.10	0.10	1.06	0.11	1.09	0.11	1.88	3.90	0.16	1.01	1.09	1.21	0.96	0.12	0.16	0.99	0.12	0.16	0.84	0.14	0.20	0.71	0.11	1.06	0.93	0.99	0.13	1.10	1.10	0.17	0.05	2.33	0.10	0.27	0.89	0.90	2.26	1.31	0.07	0.08	3.31	0.14	0.09	2.35	1.35	2.24	1.14	4.47	1.08	0.13	0.16	2.49	0.94	0.81	0.25	1.32	0.11	0.08	2.28	0.11	1.02	0.15	0.13	1.00	0.16	1.00	2.13	0.14	0.07	0.91	0.18	0.12	1.15	0.87	0.15	1.07	0.12	0.12	0.93	1.01	1.06	0.04	0.14	1.05	1.02	0.04	0.86	0.24	1.03	0.12	0.22	1.93	0.06	1.18	0.91	0.31	0.15	0.94	0.09	1.96	0.18	1.00	0.11	1.29	0.18	0.06	1.00	0.33	0.11	1.19	0.28	1.79	0.72	0.49	0.93	1.21	1.63	0.14	1.22	0.17	0.19	1.01	0.15	1.02	0.12	3.26	0.80	0.23	0.97	0.13	1.02	0.22	0.09	1.07	0.08	0.26	1.89	0.13	2.04	1.86	0.08	0.13	0.98	0.12	1.14	0.15	0.87	0.28	0.22	0.97	0.15	0.10	0.95	0.08	0.23	2.07	0.09	1.05	0.22	0.12	1.00	0.07	1.10	0.02	0.08	2.14	0.12	2.25	0.15	0.10	1.03	0.11	0.04	1.21	0.09	0.87	1.43	0.18	0.09	1.77	0.11	0.24	1.03	0.16	0.13	1.07	0.07	0.16	3.95	0.25	0.04	1.08	0.05	0.05	0.98	0.12	0.13	2.39	0.73	0.35	1.00	0.14	1.04	1.02	1.10	0.29	0.14	1.67	1.16	2.74	0.14	2.11	0.07	0.14	0.88	0.15	1.25	0.14	0.15	1.14	0.92	0.00	1.06	0.09	0.30	2.11	0.94	0.12	0.09	1.64	0.23	3.53	0.06	1.21	0.13	0.13	1.08	0.20	0.86	0.22	0.10	1.15	0.13	0.58	0.64	0.22	0.88	0.15	2.64	0.17	0.09	1.12	0.12	0.02	1.62	0.20	1.03	1.99	0.19	0.16	1.06	0.08	0.13	1.01	0.15	0.01	0.87	0.31	0.11	1.31	0.83	1.10	0.88	0.21	0.10	1.40	1.18
Flow Indexes:	1	3	6	8	10	11	11	12	15	18	20	22	24	25	27	30	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	73	75	75	76	77	80	81	83	85	87	90	92	93	94	95	95	98	100	101	104	105	108	108	108	108	108	108	111	111	113	113	115	115	117	119	121	123	126	128	130	130	131	131	131	131	133	134	135	136	139	142	145	147	148	149	151	152	155	155	158	159	160	160	161	164	164	164	167	167	168	169	169	170	171	171	171	171	172	175	175	176	177	179	182	182	184	187	189	190	190	193	196	197	199	202	203	204	207	208	210	212	215	215	217	218	221	223	223	225	227	230	233	235	235	236	238	239	240	240	242	245	247	249	249	249	250	252	254	257	260	260	262	262	263	263	266	268	270	273	276	279	279	281	284	286	289	289	291	291	294	297	299	300	303	303	306	309	312	312	312	312	315	318	321	321	322	324	326	327	328	331	331	332	333	333	333	335	335	338	340	343	344	346	349	349	350	353	353	355	355	355	355	357	360	362	365	367	368	370	372	372	372	375	378	378	380	381	381	384	387	390	393	394	395	396	399	400
Bases:	tcagACCGCAGAGTCACATGCTGCCTCCCGTAGGAGTCTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATCGAGGGCAAGTTGCATACGCG
Quality Scores:	34	34	34	34	32	22	22	22	22	27	31	34	34	34	32	32	32	34	34	34	34	34	40	40	40	40	40	40	40	40	39	40	39	39	36	28	23	23	22	33	36	36	36	36	36	36	36	36	37	34	34	34	34	34	34	33	32	32	32	32	27	27	14	14	14	14	14	14	29	29	32	32	32	32	33	34	34	34	34	34	34	34	34	34	34	34	34	34	29	34	30	31	31	31	34	33	30	30	24	24	25	24	18	18	18	15	15	15	18	18	14	17	17	17	14	14	14	14	16	12	12	16	17	22	25	25	29	33	34	34	34	32	32	32	34	34	34	34	34	34	32	32	32	34	31	31	30	32	32	32	25	22	16	16	16	16	16	15	16	16	16	22	29	14	14	14	18	24	32	34	34	34	34	34	34	34	34	30	30	30	32	34	34	34	34	34	29	29	29	29	32	25	25	18	18	18	29	32	29	34	29	29	32	29	17	17	17	17	29	24	21	14	14	14	14	14	17	29	24	24	24	24	24	24	24	24	24	14	14	14	14	14	14	14	19	19	14	14	14	14	14	14	14	19	19	19	19	24	26	22	22	19	20	20	20	20	16	16

>FLP3FBN01DB5I7
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  1251_1905

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       272
  Clip Qual Left:   5
  Clip Qual Right:  272
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.00	0.00	1.03	0.03	0.02	1.03	0.03	1.05	0.01	1.15	0.94	0.11	0.07	0.98	0.09	0.94	0.08	1.04	1.83	0.11	0.08	1.07	0.92	0.09	1.06	0.16	0.96	0.12	0.10	1.00	1.03	0.03	0.04	1.07	0.07	0.05	1.02	0.09	0.06	0.98	0.10	0.00	0.98	0.12	0.91	0.01	0.15	0.91	0.12	0.01	2.10	0.12	1.03	0.05	2.94	1.07	1.04	1.07	0.15	1.99	0.10	1.00	0.06	1.09	2.74	0.09	0.12	2.83	0.13	0.04	1.72	1.14	0.92	0.02	0.17	1.00	1.04	0.01	0.91	0.10	1.01	0.04	1.02	0.07	0.10	1.00	0.14	0.90	1.00	0.14	2.77	0.13	0.09	2.05	0.19	0.07	0.93	0.10	0.08	1.02	1.02	0.06	0.12	2.15	0.10	0.03	1.87	1.14	1.02	0.06	1.86	0.11	0.11	0.90	2.73	0.06	0.91	0.12	1.10	0.05	0.95	0.10	1.02	0.10	0.10	1.03	0.23	1.87	0.07	0.09	1.84	1.81	0.10	0.09	1.09	0.16	0.93	1.03	1.09	0.11	0.93	0.10	0.14	1.00	0.10	1.00	0.14	0.12	0.95	0.11	0.92	1.01	1.06	0.06	0.97	1.08	0.10	0.01	1.92	0.19	1.99	0.04	0.19	2.01	1.04	0.06	0.16	2.82	0.10	0.07	1.91	0.14	2.84	1.07	3.75	1.00	0.13	0.09	1.99	0.09	0.12	1.97	2.08	0.13	0.09	1.08	0.12	0.95	0.11	0.10	1.02	0.11	0.96	2.06	0.15	0.11	1.01	0.11	1.03	0.12	0.04	1.02	0.20	1.00	0.06	0.97	1.03	0.96	0.10	0.08	1.03	2.78	1.11	0.13	1.96	0.14	0.09	1.05	0.18	0.11	0.96	0.12	1.95	0.11	1.06	0.14	0.17	0.95	0.11	0.17	1.09	0.07	0.17	1.03	1.95	0.92	0.14	0.12	1.90	1.87	0.17	1.02	0.15	1.02	0.14	0.10	1.02	0.14	2.15	0.13	4.55	1.08	0.04	0.14	1.71	0.94	0.17	0.13	3.71	0.13	0.13	0.93	0.09	0.09	1.00	0.09	0.16	0.95	0.11	0.12	1.04	3.68	0.12	0.07	1.05	0.07	1.03	0.09	0.09	1.01	0.13	0.09	1.00	2.85	0.14	0.10	0.99	0.04	0.10	1.02	0.14	0.17	1.02	0.92	0.17	0.09	1.10	0.05	0.13	1.02	0.09	0.14	1.03	1.89	1.15	0.95	0.13	0.09	2.04	0.91	0.12	1.09	0.09	0.10	1.05	0.04	0.17	0.96	0.05	1.07	0.13	0.04	1.96	0.89	2.95	0.14	1.84	1.89	0.13	0.11	0.92	0.13	1.08	0.21	0.15	1.08	2.00	0.10	0.18	0.88	1.04	0.10	4.53	0.13	0.09	0.90	0.03	1.04	1.14	0.06	0.20	0.91	0.10	0.07	1.04	0.08	0.18	1.01	0.12	2.80	0.20	0.09	1.16	0.00	0.06	1.01	0.12	2.00	2.17	0.06	0.18	0.89	0.00	0.08	2.84	0.08	0.18	0.97	0.89	1.15	0.10	0.11	1.10	1.01	2.02	0.98	1.02	0.11	0.83	0.14	1.04	0.11	0.03	0.94	2.86	1.07	1.15	0.05	2.06	0.98
Flow Indexes:	1	3	6	8	10	11	14	16	18	19	19	22	23	25	27	30	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	68	71	71	72	73	76	77	79	81	83	86	88	89	91	91	91	94	94	97	100	101	104	104	107	107	108	109	111	111	114	115	115	115	117	119	121	123	126	128	128	131	131	132	132	135	137	138	139	141	144	146	149	151	152	153	155	156	159	159	161	161	164	164	165	168	168	168	171	171	173	173	173	174	175	175	175	175	176	179	179	182	182	183	183	186	188	191	193	194	194	197	199	202	204	206	207	208	211	212	212	212	213	215	215	218	221	223	223	225	228	231	234	235	235	236	239	239	240	240	242	244	247	249	249	251	251	251	251	251	252	255	255	256	259	259	259	259	262	265	268	271	272	272	272	272	275	277	280	283	284	284	284	287	290	293	294	297	300	303	304	304	305	306	309	309	310	312	315	318	320	323	323	324	325	325	325	327	327	328	328	331	333	336	337	337	340	341	343	343	343	343	343	346	348	349	352	355	358	360	360	360	363	366	368	368	369	369	372	375	375	375	378	379	380	383	384	385	385	386	387	389	391	394	395	395	395	396	397	399	399	400
Bases:	tcagACAGACCACTCACATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTGCACCGCCGGAGCTTCCCCCGCCGCCCCATGCGGGGCTGCGGGCATATGCGGTATTAGCAGCCGTTTCCGGCTGTTGTCCCCCAGTGCAGGGCAGGTTGCCCACGCGTTACTCACCCGTCCG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	39	39	39	39	39	37	36	30	30	30	38	38	38	36	36	33	39	39	39	39	37	37	37	37	37	31	31	31	37	37	37	37	37	37	37	38	38	38	37	37	38	34	30	28	28	38	37	37	37	37	37	37	34	34	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	33	33	33	37	37	38	38	38	38	28	28	28	28	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	31	31	31	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	22	22	15	15	15	15	15	22	23	23	28	17	17	17	17	32	37	37	34	21	21	21	21	38	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	33	36	33	33	32	32	28	28	15	15	15	15	15	28	32	35	33	35	35	26	26	26	31	36	36	36	36	36	36	25	25	25	31	31	31	31	31	31	31	32	31	31	31	32	28	28	28	30	30	32	31	31

>FLP3FBN01AK9OO
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0125_0438

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       257
  Clip Qual Left:   5
  Clip Qual Right:  253
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	0.94	0.02	1.09	0.05	0.02	0.96	0.03	1.06	0.04	1.14	0.57	0.21	0.06	1.08	0.04	1.08	0.12	1.08	0.06	0.98	0.83	0.07	0.94	1.71	0.19	0.26	0.93	0.21	0.67	0.13	1.03	0.06	0.09	1.09	0.07	0.08	1.05	0.06	0.07	1.06	0.08	0.02	1.04	0.11	0.94	0.01	0.15	1.01	0.08	0.03	2.06	0.09	1.01	0.07	3.04	1.07	1.02	1.07	0.14	2.10	0.12	0.92	0.05	1.05	0.94	0.08	1.08	0.13	1.06	0.01	0.11	3.21	0.11	0.01	2.19	0.95	1.04	0.08	0.09	0.92	0.94	0.02	0.93	0.09	1.08	0.03	1.04	0.07	0.08	1.07	0.12	1.03	1.09	0.11	3.07	0.12	0.09	2.00	0.05	0.09	0.93	0.12	0.07	0.99	1.11	0.03	0.14	2.14	0.09	0.03	2.22	1.18	0.16	0.86	0.24	0.11	0.77	0.28	1.16	0.06	0.14	0.97	2.95	0.06	1.08	0.10	1.05	0.06	1.01	0.12	0.77	0.13	0.09	0.71	0.36	2.19	1.01	0.10	1.14	1.89	0.15	0.17	1.20	0.08	0.96	1.11	0.73	1.01	0.14	0.11	1.00	0.15	0.15	0.81	0.14	0.17	0.86	0.08	1.12	1.01	1.36	0.03	0.94	1.20	0.11	0.03	1.98	0.40	1.91	0.05	0.20	2.30	1.23	0.13	0.60	0.95	0.22	1.09	0.11	0.75	0.03	0.11	2.21	1.03	2.08	1.17	2.01	0.16	0.93	0.11	1.01	0.14	0.05	1.08	1.84	0.10	0.13	1.97	0.99	0.14	1.07	0.84	0.13	1.07	0.04	0.10	1.13	0.16	0.93	1.87	0.19	0.09	1.02	0.10	0.12	1.00	0.08	0.14	1.09	0.97	0.13	0.14	1.79	0.86	0.24	0.25	1.06	2.68	0.16	0.12	2.61	0.11	0.05	1.11	0.14	0.22	0.93	0.16	0.99	0.10	1.10	0.14	0.98	0.12	0.14	1.08	0.22	0.13	1.08	0.89	0.18	0.96	0.40	0.06	1.20	1.99	0.33	0.74	0.05	0.00	1.99	1.09	1.02	0.10	0.51	0.01	2.03	0.14	2.82	2.00	0.20	0.13	1.92	0.07	0.07	1.17	0.02	0.16	0.98	0.19	0.77	1.00	0.80	0.22	2.00	0.07	1.11	0.21	0.00	1.18	0.10	0.12	1.11	0.13	0.09	0.89	0.27	0.00	0.94	0.96	0.24	1.86	0.12	1.12	0.04	0.13	1.02	0.23	1.13	0.70	1.03	0.15	0.09	1.93	0.12	0.13	1.01	0.13	1.06	0.24	1.70	0.88	0.12	0.12	1.19	0.11	0.03	0.82	0.33	0.00	1.13	1.84	0.65	0.98	0.19	0.03	1.91	1.79	0.08	0.12	1.17	0.07	1.03	0.13	1.05	0.14	3.08	0.00	3.12	0.25	0.81	1.96	0.13	2.00	0.04	2.03	0.04	0.04	1.13	0.21	0.88	1.00	0.12	0.19	1.93	0.15	3.69	0.27	0.98	0.03	1.15	0.16	1.22	1.14	0.21	0.03	0.81	0.17	0.26	0.90	0.11	1.02	0.07	2.26	0.77	0.03	1.49	0.20	0.00	1.13	0.12	2.21	2.27	0.14	0.26	1.19	0.29	0.19	3.94	0.26
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	72	75	75	76	77	80	81	83	85	87	90	92	93	95	95	95	98	98	101	104	105	108	108	111	111	112	114	117	119	122	123	123	123	125	127	129	131	134	136	136	137	139	140	140	143	145	146	147	148	151	154	157	159	160	161	163	164	167	167	169	169	172	172	173	175	176	178	180	183	183	184	185	185	186	187	187	189	191	194	195	195	198	198	199	201	202	204	207	209	210	210	213	216	219	220	223	223	224	227	228	228	228	231	231	231	234	237	239	241	243	246	249	250	252	255	256	256	258	261	261	262	263	265	267	267	269	269	269	270	270	273	273	276	279	281	282	283	285	285	287	290	293	296	299	300	302	302	304	307	309	310	311	314	314	317	319	321	321	322	325	328	331	332	332	333	334	337	337	338	338	341	343	345	347	347	347	349	349	349	351	352	352	354	354	356	356	359	361	362	365	365	367	367	367	367	369	371	373	374	377	380	382	384	384	385	387	390	392	392	393	393	396	399	399	399	399
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTCGAGCCGTTACCTCACCAACTAGCTAATGCGCCGCGGGCCCATCTCATAGCGGATTACTCCTTTAATTGCTACTTCATGCGAAGCTACAATCTTATGCGGTATTAATCTCCCTTTCGGAAGGCTATTCCCCTCTATGAGGTCAGGTTGcccc
Quality Scores:	34	34	34	31	27	27	27	31	34	34	32	27	27	22	22	23	30	34	34	34	34	34	37	40	40	40	40	40	40	40	40	39	38	38	38	39	37	36	36	30	30	30	36	36	36	36	37	37	37	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	27	27	26	27	27	27	31	34	34	34	34	32	32	27	21	21	22	34	34	30	30	30	32	32	32	34	32	32	32	31	31	30	30	30	30	30	31	31	19	19	18	19	24	26	30	31	31	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	34	33	32	26	22	24	24	18	18	18	18	14	14	27	27	32	30	32	32	32	30	31	31	25	25	25	34	34	25	25	25	25	25	25	25	25	34	34	34	32	34	34	32	32	32	32	32	34	34	34	34	34	32	32	32	34	30	29	29	29	34	32	32	24	24	24	24	29	29	29	17	17	17	24	29	29	29	29	32	34	34	34	34	34	34	24	24	18	28	34	34	34	34	34	34	34	34	23	23	18	18	16	16	19	12	24	24	24	24	24	14	14	14	14	19	14	12	12	12	15	14	14	20	20

>FLP3FBN01ANGF2
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0150_0112

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       257
  Clip Qual Left:   5
  Clip Qual Right:  257
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	0.93	0.06	1.05	0.06	0.07	1.03	0.08	1.02	0.06	1.09	1.00	0.08	0.08	1.06	0.07	0.90	0.07	1.01	0.05	1.02	0.96	0.09	0.89	1.96	0.09	0.01	1.08	0.11	1.03	0.03	1.05	0.07	0.09	1.07	0.08	0.08	0.98	0.07	0.05	0.93	0.11	0.04	0.98	0.10	0.97	0.03	0.13	0.97	0.13	0.04	2.05	0.14	1.03	0.05	3.00	1.01	1.01	1.05	0.12	1.95	0.13	1.02	0.10	1.02	1.02	0.09	0.97	0.14	1.01	0.07	0.15	1.87	0.13	1.01	1.95	1.05	1.03	0.07	0.15	0.98	1.04	0.03	1.03	0.09	1.01	0.07	0.92	0.11	0.14	1.01	0.12	1.01	1.94	0.14	1.90	0.12	0.13	0.98	0.14	1.04	1.00	0.08	0.07	1.12	0.90	0.08	0.13	4.70	0.16	1.00	2.11	0.09	1.94	0.08	1.94	0.08	1.01	0.02	1.11	0.07	0.95	0.11	0.95	0.12	0.10	1.07	0.18	1.04	0.11	2.06	3.98	0.12	1.04	1.01	0.15	0.89	0.13	1.02	0.96	0.10	0.13	1.00	0.14	0.14	1.07	0.11	0.99	1.11	0.95	0.04	1.02	0.99	0.09	0.07	2.00	0.12	0.11	1.02	1.07	2.09	1.01	0.12	0.10	2.96	0.12	0.09	2.15	0.91	2.02	1.03	4.13	1.04	0.10	0.11	2.07	0.94	0.96	0.10	1.00	0.13	0.16	1.97	0.16	1.00	0.12	0.06	1.05	0.07	1.06	2.00	0.06	0.15	1.00	0.00	0.11	1.02	0.96	0.13	1.09	0.14	0.12	0.98	1.11	0.94	0.09	0.10	1.02	0.95	0.11	0.96	0.07	0.94	0.10	0.09	2.15	0.09	1.10	1.07	0.16	0.08	0.98	0.14	2.17	0.07	0.99	0.12	0.93	0.08	0.09	1.07	0.15	0.14	1.06	0.09	1.96	1.06	0.06	1.02	1.02	1.99	0.15	1.01	0.15	0.01	1.09	0.10	0.98	0.13	2.87	1.01	0.07	1.05	0.12	0.98	0.19	0.17	1.05	0.06	0.03	2.01	0.07	1.97	2.09	0.18	0.11	0.92	0.19	1.08	0.14	0.91	0.16	0.16	1.05	0.09	0.12	1.03	0.16	0.11	2.04	0.16	1.04	0.12	0.11	1.00	0.08	1.03	0.05	0.12	1.97	0.15	2.04	0.14	0.09	1.01	0.18	0.15	0.98	0.04	0.95	1.11	0.05	0.14	1.97	0.10	0.08	1.03	0.14	0.00	1.07	0.05	0.10	3.99	0.12	0.15	1.04	0.01	0.05	1.04	0.20	0.07	2.09	0.91	0.18	1.02	0.06	1.00	1.04	1.04	0.05	0.04	2.10	1.08	3.00	0.17	2.02	0.07	0.05	1.02	0.13	0.92	0.09	0.11	1.08	1.00	0.01	1.06	0.17	0.05	2.13	0.96	0.07	0.12	1.94	0.09	4.20	0.16	0.99	0.16	0.10	0.94	0.18	1.10	0.05	0.05	0.99	0.14	0.15	1.07	0.12	0.99	0.13	2.82	0.16	0.17	1.01	0.13	0.05	2.02	0.16	1.05	2.04	0.07	0.06	1.01	0.14	0.13	0.95	0.12	1.01	0.12	0.94	0.08	0.11	0.96	0.91	1.00	0.28	0.15	1.11	1.13
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	67	69	72	72	74	75	75	76	77	80	81	83	85	87	90	92	93	93	95	95	98	100	101	104	105	108	108	108	108	108	110	111	111	113	113	115	115	117	119	121	123	126	128	130	130	131	131	131	131	133	134	136	138	139	142	145	147	148	149	151	152	155	155	158	159	160	160	161	164	164	164	167	167	168	169	169	170	171	171	171	171	172	175	175	176	177	179	182	182	184	187	189	190	190	193	196	197	199	202	203	204	207	208	210	212	215	215	217	218	221	223	223	225	227	230	233	235	235	236	238	239	240	240	242	245	247	249	249	249	250	252	254	257	260	260	262	262	263	263	266	268	270	273	276	279	279	281	284	286	289	289	291	291	294	297	299	300	303	303	306	309	312	312	312	312	315	318	321	321	322	324	326	327	328	331	331	332	333	333	333	335	335	338	340	343	344	346	349	349	350	353	353	355	355	355	355	357	360	362	365	368	370	372	372	372	375	378	378	380	381	381	384	387	389	391	394	395	396	399	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTCTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGGATCTTTAGATGGAACCAGATGCCTGATTCCATCGCCATGGGGCATTAGACGCCGTTTCCAGCGATTATTCCCCTGATGAGGGCAAGTTGCTCACGCG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	38	38	28	28	28	28	28	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	37	37	36	36	36	29	28	28	28	36	36	36	36	36	36	28	28	29	36	38	38	38	38	38	36	32	32	32	32	30	28	28	28

>FLP3FBN01AF994
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0068_1402

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       276
  Clip Qual Left:   5
  Clip Qual Right:  276
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.02	0.06	1.01	0.06	0.07	1.02	0.08	1.01	0.07	1.06	1.00	0.05	0.03	1.02	0.07	0.99	0.04	1.06	0.03	1.01	1.03	0.01	1.03	2.09	0.08	0.04	0.99	0.10	1.04	0.03	1.00	0.05	0.06	1.04	0.08	0.06	0.99	0.05	0.07	0.93	0.09	0.05	1.08	0.10	0.99	0.05	0.09	0.98	0.07	0.03	1.97	0.11	1.01	0.05	3.09	1.01	1.06	1.03	0.14	2.11	0.14	1.02	0.07	0.98	3.05	0.07	0.08	2.07	1.00	0.05	2.09	1.01	1.00	0.05	0.12	0.98	1.05	0.04	1.00	0.11	1.02	0.05	1.03	0.09	0.09	0.98	0.14	1.01	1.02	1.07	2.06	0.10	0.10	2.00	0.12	0.09	1.01	0.12	0.07	0.99	0.98	0.06	0.11	4.99	0.14	1.05	2.06	0.12	2.02	0.10	2.04	0.11	1.06	0.01	0.92	0.07	1.00	0.07	1.01	0.11	0.10	1.01	0.07	1.07	0.10	2.02	4.17	0.10	1.05	1.05	0.13	0.03	0.99	0.16	2.04	0.05	0.13	1.01	0.10	0.11	0.96	0.09	0.96	1.03	1.05	0.03	1.04	2.02	3.00	0.06	0.11	1.96	1.05	0.05	0.14	3.03	0.09	0.08	2.08	0.99	2.09	1.08	4.05	0.94	0.04	0.16	2.05	0.13	0.09	2.09	0.94	0.09	0.95	0.06	0.13	1.08	0.16	0.09	1.96	0.15	0.97	1.98	0.07	0.16	1.06	0.10	0.10	2.10	0.14	1.90	1.06	1.03	0.07	0.08	0.95	0.17	0.09	0.95	0.14	0.11	0.99	0.12	0.14	1.06	0.11	0.07	2.21	0.16	1.03	1.02	0.13	0.01	0.99	0.12	1.02	0.00	1.05	0.97	0.11	0.04	1.07	0.17	1.01	0.11	0.15	1.02	0.18	1.04	0.05	0.10	1.07	1.05	0.09	0.91	0.19	0.15	1.06	0.07	0.07	0.96	0.15	2.00	0.14	0.15	1.01	0.10	1.04	0.15	3.05	1.05	0.13	1.06	0.10	0.02	0.92	0.09	0.16	3.09	0.12	0.10	1.03	1.03	0.08	0.04	0.93	0.17	4.07	0.08	0.11	1.04	0.17	0.16	1.03	0.16	0.20	1.13	0.10	0.04	1.01	3.81	0.11	0.06	3.83	0.13	0.96	0.05	0.18	0.94	0.06	0.07	1.05	0.13	1.96	0.09	0.93	0.11	0.10	1.03	0.13	0.19	0.97	0.17	0.18	1.09	0.19	0.12	1.17	2.00	1.04	0.98	0.14	0.10	1.94	0.94	0.06	1.00	1.05	0.15	1.99	1.07	0.08	0.93	1.05	0.14	2.79	0.08	1.01	0.94	0.23	0.14	2.05	2.87	1.67	0.96	0.33	0.32	1.12	0.13	4.48	0.11	1.01	0.10	1.57	0.20	1.10	0.13	0.23	1.04	0.18	0.99	0.18	0.14	1.03	0.91	0.10	1.77	0.87	1.75	0.16	1.15	1.78	0.36	0.24	1.08	0.44	0.22	0.74	0.30	0.19	0.72	0.23	0.20	0.96	1.16	0.86	0.97	0.27	0.24	1.23	1.17	1.62	1.00	1.10	0.29	1.57	0.21	1.35	0.23	0.34	0.91	2.89	0.83	0.91	0.19	0.76	1.15	0.33	0.14	1.32	1.58
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	69	71	71	72	73	76	77	79	81	83	86	88	89	90	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	133	135	135	138	141	143	144	145	147	148	148	149	149	149	152	152	153	156	156	156	159	159	160	161	161	162	163	163	163	163	164	167	167	170	170	171	173	176	179	179	181	182	182	185	188	188	190	190	191	192	195	198	201	204	207	207	209	210	213	215	217	218	221	223	226	228	231	232	234	237	240	242	242	245	247	249	249	249	250	252	255	258	258	258	261	262	265	267	267	267	267	270	273	276	279	280	280	280	280	283	283	283	283	285	288	291	293	293	295	298	301	304	307	308	308	309	310	313	313	314	316	317	319	319	320	322	323	325	325	325	327	328	331	331	332	332	332	333	333	334	337	339	339	339	339	341	343	343	345	348	350	353	354	356	356	357	358	358	360	361	361	364	367	370	373	374	375	376	379	380	381	381	382	383	385	385	387	390	391	391	391	392	393	395	396	399	400	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGTCCGTGTCTCAGTACCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATGGAACGCATGCCTATCTATCAGCGATGAATCTTTAGCAAATATCCCCATGCGGGGCCCCTGCTTCATGCGGTATTAGTCCGACTTTCGCCGGGTTATCCCCTCCTGATAGGTAAGTTGCATACGCGTTACTTCACCCGTCGCGG
Quality Scores:	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	37	37	37	37	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	31	31	31	31	35	35	35	35	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	36	31	31	31	36	36	35	30	28	16	16	16	16	16	16	16	16	16	16	16	14	14	16	24	35	35	30	24	24	24	24	24	24	16	16	14	14	14	21	21	21	21	14	14	13	13	21	14	13	13	14	14	13	24	24	21	21	13	14	14	13	13

>FLP3FBN01AHXU8
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0087_0802

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       289
  Clip Qual Left:   5
  Clip Qual Right:  289
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.03	0.03	0.99	0.06	0.05	1.03	0.07	1.02	0.06	1.07	1.00	0.08	0.06	1.00	0.06	0.95	0.07	1.04	0.10	0.98	1.01	0.07	0.88	1.98	0.10	0.04	1.04	0.09	1.03	0.05	1.00	0.06	0.08	1.06	0.08	0.06	1.00	0.05	0.07	1.00	0.09	0.03	1.01	0.08	0.98	0.02	0.14	1.01	0.11	0.01	2.02	0.10	1.02	0.06	3.00	1.01	1.05	1.05	0.12	2.06	0.15	1.01	0.08	1.04	3.07	0.08	0.10	2.02	0.10	0.97	1.96	1.12	0.94	0.06	0.16	0.92	1.01	0.04	0.99	0.11	1.05	0.03	1.01	0.09	0.12	1.03	0.15	1.02	2.05	0.13	2.04	0.10	0.08	2.03	0.11	0.09	1.05	0.14	0.07	1.01	0.99	0.06	0.11	4.99	0.15	1.00	2.03	0.14	2.04	0.08	2.03	0.12	1.04	0.06	1.04	0.07	1.01	0.07	1.01	0.11	0.15	0.98	0.08	1.00	0.09	1.98	3.88	0.14	1.03	1.04	1.04	1.00	0.12	0.09	1.01	0.10	0.09	0.96	0.15	0.11	1.04	0.12	0.99	1.02	1.00	0.10	0.93	1.96	2.04	1.07	0.10	1.93	1.02	0.11	0.15	2.90	0.08	0.07	2.00	1.03	2.03	1.02	3.76	1.04	0.15	0.13	1.89	0.12	0.92	1.01	1.06	0.15	1.02	0.14	0.11	1.03	0.12	0.14	2.00	0.08	1.02	2.11	0.15	0.11	1.00	0.08	0.16	1.02	0.07	0.10	1.00	0.96	0.09	0.06	1.88	1.00	0.09	0.09	0.95	0.15	0.11	1.04	0.11	0.12	1.06	0.13	0.02	1.05	0.14	0.10	3.00	0.10	0.12	1.08	0.07	0.11	0.99	0.05	1.89	0.09	1.02	0.12	1.90	0.09	0.13	1.04	2.04	2.03	1.08	2.09	0.04	0.09	0.95	0.02	2.08	0.14	2.96	0.11	1.78	0.05	1.00	0.12	4.56	3.87	0.09	1.02	0.15	0.15	1.08	0.05	0.04	1.05	0.12	0.04	4.97	0.06	0.12	2.08	0.12	2.95	0.11	0.70	0.10	0.06	0.98	1.01	0.15	0.07	1.10	1.03	1.18	1.03	0.09	0.10	1.05	2.98	0.05	1.03	0.13	0.12	2.13	0.91	0.07	1.01	0.13	0.05	2.02	0.15	0.98	0.09	2.80	0.09	2.96	0.09	1.03	1.98	0.03	1.94	0.18	1.86	2.03	0.04	0.07	1.01	1.02	0.12	4.66	0.08	1.01	0.09	0.03	1.04	1.01	0.13	0.10	1.81	0.08	1.05	0.09	2.82	0.18	0.08	1.00	0.01	0.07	1.01	0.13	1.94	2.06	0.07	0.19	1.01	0.05	0.15	1.05	0.09	0.09	1.07	0.19	0.17	1.04	1.01	1.07	0.96	0.87	0.11	0.13	1.07	2.05	1.03	1.09	0.09	0.95	0.09	0.77	0.18	0.13	1.06	3.05	1.13	1.09	0.02	0.12	0.93	0.06	0.08	0.89	1.05	0.17	0.11	2.01	0.15	0.17	0.96	0.17	1.09	1.08	0.01	1.03	0.84	0.08	0.09	1.97	2.06	0.13	0.12	0.96	0.10	0.09	1.12	0.08	1.09	0.12	0.70	0.16	0.94	0.14	0.96	0.26	1.14
Flow Indexes:	1	3	6	8	10	11	14	16	18	20	21	23	24	24	27	29	31	34	37	40	43	45	48	51	51	53	55	55	55	56	57	58	60	60	62	64	65	65	65	68	68	70	71	71	72	73	76	77	79	81	83	86	88	89	89	91	91	94	94	97	100	101	104	104	104	104	104	106	107	107	109	109	111	111	113	115	117	119	122	124	126	126	127	127	127	127	129	130	131	132	135	138	141	143	144	145	147	148	148	149	149	150	152	152	153	156	156	156	159	159	160	161	161	162	163	163	163	163	164	167	167	169	170	171	173	176	179	179	181	182	182	185	188	191	192	195	195	196	199	202	205	208	211	211	211	214	217	219	219	221	223	223	226	227	227	228	228	229	230	230	233	235	235	237	237	237	239	239	241	243	243	243	243	243	244	244	244	244	246	249	252	255	255	255	255	255	258	258	260	260	260	262	265	266	269	270	271	272	275	276	276	276	278	281	281	282	284	287	287	289	291	291	291	293	293	293	295	296	296	298	298	300	300	301	301	304	305	307	307	307	307	307	309	312	313	316	316	318	320	320	320	323	326	328	328	329	329	332	335	338	341	342	343	344	345	348	349	349	350	351	353	355	358	359	359	359	360	361	364	367	368	371	371	374	376	377	379	380	383	383	384	384	387	390	392	394	396	398	400
Bases:	tcagACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCGGGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAGAGAG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	32	32	32	32	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	32	28	15	15	15	15	15	22	22	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	34	34	34	37	37	37	37	37	37	37	36	36	36	36	36	36	31	24	22	22	22	22	33	36	36	36	36	36	26	27	28	37	38	38	38	37	38	38	38	32	32	31	32	32	31	31	31	32	32	31	30	30	30	30	32	32	32	32	32	32	32	31	31	31	31	31	32	31	31	31	31	31	31	28	28	28	28	28	28

>FLP3FBN01APRWO
  Run Prefix:   R_2008_12_09_13_51_01_
  Region #:     1
  XY Location:  0176_1782

  Run Name:       R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted
  Analysis Name:  D_2009_06_10_10_15_54_node176_fullProcessingAmplicons
  Full Path:      /srv/seq/454/2008_12_09/R_2008_12_09_13_51_01_FLX04080350_Administrator_PC_Calrestricted/D_2009_06_10_10_15_54_node176_fullProcessingAmplicons/

  Read Header Len:  32
  Name Length:      14
  # of Bases:       276
  Clip Qual Left:   5
  Clip Qual Right:  276
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.06	0.03	1.02	0.05	0.04	0.93	0.05	1.01	0.04	1.07	0.03	1.06	0.06	0.02	0.92	0.13	0.01	0.99	0.06	0.93	0.04	0.05	1.02	0.05	0.10	1.09	1.02	0.07	2.07	0.05	0.05	0.97	1.03	0.06	0.93	0.09	0.06	1.04	0.07	0.06	1.02	0.05	0.03	1.00	0.08	0.02	1.02	0.08	0.97	0.03	0.15	0.97	0.09	0.04	1.90	0.14	1.00	0.08	2.86	1.07	1.03	1.04	0.13	1.96	0.11	0.92	0.04	1.10	2.82	0.06	0.11	2.71	0.13	0.01	1.77	1.06	1.02	0.04	0.15	0.92	1.04	0.03	0.98	0.12	1.00	0.08	0.98	0.10	0.12	1.01	0.15	0.94	1.03	1.14	1.92	0.14	0.13	0.95	0.08	1.04	1.04	0.09	0.07	1.29	1.03	0.09	0.15	4.81	0.19	1.00	1.98	0.16	1.95	0.08	2.02	0.09	1.04	0.03	1.07	0.06	1.01	0.09	0.97	0.12	0.12	1.07	0.25	1.03	0.11	2.01	3.88	0.16	0.91	1.04	1.15	0.88	0.13	0.13	0.97	0.13	0.12	0.88	0.13	0.12	0.99	0.11	1.00	0.92	1.03	0.03	1.08	1.02	0.11	0.02	1.97	0.17	0.95	0.03	0.95	2.05	1.09	0.03	0.13	2.85	0.14	0.09	2.18	1.15	2.01	1.02	4.07	1.05	0.10	0.13	2.27	1.05	0.11	0.98	0.89	0.14	0.98	1.00	0.15	0.99	0.16	0.14	0.96	0.09	1.03	1.93	0.09	0.17	0.98	0.06	0.14	1.00	0.14	0.04	1.10	1.00	0.03	0.03	2.08	1.01	0.04	0.07	0.98	0.10	0.11	0.96	0.10	0.09	1.03	0.14	0.11	2.00	0.09	0.10	1.98	0.13	0.14	1.04	0.12	0.15	1.01	0.13	1.85	1.08	0.05	0.10	1.03	0.09	0.04	0.97	0.09	0.94	0.07	0.13	1.88	1.03	0.13	1.12	0.09	0.07	1.10	2.92	0.10	0.07	1.04	0.12	0.99	0.12	3.01	3.17	1.09	0.10	0.17	0.95	0.18	0.13	0.89	0.16	1.02	3.05	0.12	1.00	0.08	1.06	0.12	0.95	0.13	0.17	1.05	0.17	0.09	1.03	0.09	0.06	2.22	0.09	1.02	0.07	2.96	0.11	0.10	1.89	1.06	1.03	2.04	0.16	1.03	0.10	2.09	1.03	0.94	1.01	0.07	0.19	1.02	1.96	1.09	0.96	0.10	0.08	2.12	1.04	0.18	0.97	0.11	0.95	1.05	1.96	0.23	1.96	0.11	0.06	2.95	0.11	1.05	0.16	2.07	0.11	1.96	1.01	0.15	0.14	1.03	0.07	2.10	1.04	0.23	0.11	0.92	0.04	4.86	0.14	1.00	0.07	0.08	0.92	0.17	0.01	1.08	0.10	1.05	0.03	0.06	0.98	0.04	0.03	1.08	3.02	0.14	0.12	0.86	0.01	0.10	1.08	0.06	2.01	2.16	0.18	2.25	0.11	0.15	0.98	0.15	0.19	1.05	1.00	1.10	1.04	1.04	0.17	0.19	0.90	1.96	0.97	1.08	0.10	0.95	0.13	1.05	0.19	0.07	1.04	3.01	1.15	1.07	0.01	0.09	0.94	0.09	0.09	0.92	0.88	0.20	0.04	2.28	2.28
Flow Indexes:	1	3	6	8	10	12	15	18	20	23	26	27	29	29	32	33	35	38	41	44	47	49	52	55	55	57	59	59	59	60	61	62	64	64	66	68	69	69	69	72	72	72	75	75	76	77	80	81	83	85	87	90	92	93	94	95	95	98	100	101	104	105	108	108	108	108	108	110	111	111	113	113	115	115	117	119	121	123	126	128	130	130	131	131	131	131	133	134	135	136	139	142	145	147	148	149	151	152	155	155	157	159	160	160	161	164	164	164	167	167	168	169	169	170	171	171	171	171	172	175	175	176	178	179	181	182	184	187	189	190	190	193	196	199	200	203	203	204	207	210	213	216	216	219	219	222	225	227	227	228	231	234	236	239	239	240	242	245	246	246	246	249	251	253	253	253	254	254	254	255	258	261	263	264	264	264	266	268	270	273	276	279	279	281	283	283	283	286	286	287	288	289	289	291	293	293	294	295	296	299	300	300	301	302	305	305	306	308	310	311	312	312	314	314	317	317	317	319	321	321	323	323	324	327	329	329	330	333	335	335	335	335	335	337	340	343	345	348	351	352	352	352	355	358	360	360	361	361	363	363	366	369	370	371	372	373	376	377	377	378	379	381	383	386	387	387	387	388	389	392	395	396	399	399	400	400
Bases:	tcagAGCAGCACTTGTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTCTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	38	38	38	39	39	39	39	39	39	39	34	34	36	30	30	30	39	39	39	39	39	39	39	37	37	37	37	37	37	37	37	37	37	37	38	38	33	33	33	33	33	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	38	38	37	37	37	37	37	37	37	38	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	36	37	38	38	37	37	38	38	38	38	38	38	38	38	38	38	38	38	36	36	36	38	38	38	38	29	29	29	23	28	30	31	32	31	31	31	31	32	32	32	32	32	32	31	31	31	32	32	31	28	25	19	18	18	16
"""

if __name__ == "__main__":
    main()
