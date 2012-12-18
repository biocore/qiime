#!/usr/bin/env python
# File created on 02 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

import signal
from os import makedirs, system, chdir, getcwd
from os.path import exists
from shutil import rmtree
from numpy import array
from cogent import LoadTree, LoadSeqs
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.util import load_qiime_config, count_seqs, get_tmp_filename
from qiime.workflow import (call_commands_serially,no_status_updates,WorkflowError)
from qiime.parse import (parse_qiime_parameters,
                         fields_to_dict)
from qiime.pick_subsampled_reference_otus_through_otu_table import (
                                pick_subsampled_open_reference_otus,
                                iterative_pick_subsampled_open_reference_otus,
                                final_repset_from_iteration_repsets)
from biom.parse import parse_biom_table

## The test case timing code included in this file is adapted from
## recipes provided at:
##  http://code.activestate.com/recipes/534115-function-timeout/
##  http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call
class TimeExceededError(Exception):
    pass


allowed_seconds_per_test = 500

def timeout(signum, frame):
    raise TimeExceededError,\
     "Test failed to run in allowed time (%d seconds)."\
      % allowed_seconds_per_test


class PickSubsampledReferenceOtusThroughOtuTableTests(TestCase):
    """ """
    
    def setUp(self):
        """
        """
        self.start_dir = getcwd()
        self.qiime_config = load_qiime_config()
        self.dirs_to_remove = []
        self.files_to_remove = []
        
        # use 'fast' uclust settings for reduced runtime
        self.params = parse_qiime_parameters([
         "pick_otus:max_rejects	8",
         "pick_otus:word_length	8",
         "pick_otus:max_accepts	1",
         "pick_otus:stepwords	8"])
        
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'
        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(self.tmp_dir)
        
        self.wf_out = get_tmp_filename(tmp_dir=self.tmp_dir,
         prefix='qiime_wf_out',suffix='',result_constructor=str)
        self.dirs_to_remove.append(self.wf_out)
        
        self.fasting_seqs_fp1 = get_tmp_filename(tmp_dir=self.tmp_dir,
            prefix='qiime_wf_seqs',suffix='.fasta')
        fasting_seqs_f = open(self.fasting_seqs_fp1,'w')
        fasting_seqs_f.write(fasting_seqs_subset)
        fasting_seqs_f.close()
        self.files_to_remove.append(self.fasting_seqs_fp1)
        
        self.pick_ref_otus_refseqs2 = get_tmp_filename(
            tmp_dir=self.tmp_dir,prefix='ref_otus_wf',suffix='.fna')
        f = open(self.pick_ref_otus_refseqs2,'w')
        f.write(pick_ref_otus_refseqs2)
        f.close()
        self.files_to_remove.append(self.pick_ref_otus_refseqs2)
        
        self.fasting_seqs_fp2 = get_tmp_filename(tmp_dir=self.tmp_dir,
            prefix='qiime_wf_seqs',suffix='.fasta')
        fasting_seqs_f = open(self.fasting_seqs_fp2,'w')
        fasting_seqs_f.write(subsample_otus_seqs)
        fasting_seqs_f.close()
        self.files_to_remove.append(self.fasting_seqs_fp2)
        
        signal.signal(signal.SIGALRM, timeout)
        # set the 'alarm' to go off in allowed_seconds seconds
        signal.alarm(allowed_seconds_per_test)
    
    def tearDown(self):
        """ """
        # turn off the alarm
        signal.alarm(0)
        # change back to the start dir - some workflows change directory
        chdir(self.start_dir)
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)
    
    def test_final_repset_from_iteration_repsets(self):
        """ final_repset_from_iteration_repsets functions as expected """
        repset1 = """>o1
ACCGT
>o2
AGG
>o3
ACCGTT""".split('\n')

        repset2 = """>o4
TACCGT
>o5
TAGG
>o6
TACCGTT""".split('\n')

        repset3 = """>o4
CAT
>o7
AAAA
>o1
A""".split('\n')
        
        exp = [("o1","ACCGT"), ("o2","AGG"), ("o3","ACCGTT"),
                ("o4","TACCGT"), ("o5","TAGG"), ("o6","TACCGTT")]
        actual = list(final_repset_from_iteration_repsets([repset1,repset2]))
        self.assertEqual(actual,exp)

        exp = [("o1","ACCGT"), ("o2","AGG"), ("o3","ACCGTT"),
                ("o4","TACCGT"), ("o5","TAGG"), ("o6","TACCGTT"),('o7','AAAA')]
        actual = list(final_repset_from_iteration_repsets([repset1,repset2,repset3]))
        self.assertEqual(actual,exp)
    

    def test_iterative_pick_subsampled_open_reference_otus_no_prefilter(self):
        """pick_subsampled_open_reference_otus functions as expected without prefilter
        """
        iterative_pick_subsampled_open_reference_otus(
          input_fps=[self.fasting_seqs_fp1,self.fasting_seqs_fp2],
          refseqs_fp=self.pick_ref_otus_refseqs2,
          output_dir=self.wf_out,
          percent_subsample=0.5,
          new_ref_set_id='wf.test.otu',
          command_handler=call_commands_serially,
          params=self.params,
          qiime_config=self.qiime_config,
          prefilter_refseqs_fp=None,
          prefilter_percent_id=None,
          step1_otu_map_fp=None,
          step1_failures_fasta_fp=None,
          parallel=False,
          suppress_step4=False,
          logger=None,
          status_update_callback=no_status_updates)
        
        for i in (0,1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out,i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out,i)
            self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
            # No final failures file should be created
            self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
            
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out,i)
            self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
            
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out,i)
            self.assertTrue(exists(repset_fp),"Rep set doesn't exist.")
            
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out,i)
            self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist.")
            
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(set(otu_table.SampleIds),
                         set(['PC.355','PC.481','PC.636','PC.354','PC.635',
                              'PC.593','PC.607','PC.356','PC.634','sample1']))
        self.assertTrue('r102' in otu_table.ObservationIds,\
         "Reference OTU (r102) is not in the final OTU table.")
        # at least one reference and one cleanup otu from the first iteration are in the 
        # otu table - this may very rarely fail by chance
        self.assertTrue(array([o.startswith('wf.test.otu.0.Reference')\
          for o in otu_table.ObservationIds]).any())
        self.assertTrue(array([o.startswith('wf.test.otu.0.CleanUp.Reference')\
          for o in otu_table.ObservationIds]).any())
        
        # At least one otu from the second iteration exists
        self.assertTrue(array([o.startswith('wf.test.otu.1.Reference')\
                               for o in otu_table.ObservationIds]).any()
                        or 
                        array([o.startswith('wf.test.otu.1.CleanUp.Reference')\
                               for o in otu_table.ObservationIds]).any())
        
        # No singletons in the final OTU table
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        
        # number of otus in otu table equals sum of seqs in rep_set.fna
        self.assertEqual(len(otu_table.ObservationIds),
                         count_seqs('%s/rep_set.fna' % self.wf_out)[0])
        
        # confirm that tax/align/tree files exist
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        self.assertTrue(exists(tree_fp))
        self.assertTrue(exists(aln_fp))
        num_tree_tips = len(LoadTree(tree_fp).tips())
        num_align_seqs = LoadSeqs(aln_fp).getNumSeqs()
        self.assertEqual(num_tree_tips,num_align_seqs)
        self.assertTrue(num_tree_tips > 30)
        # OTU table minus pynast failures has same number of otus as aligned
        # sequences
        otu_table_fp = '%s/otu_table_mc2_w_tax_no_pynast_failures.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp))
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(len(otu_table.ObservationIds),num_align_seqs)

    def test_iterative_pick_subsampled_open_reference_otus(self):
        """pick_subsampled_open_reference_otus functions as expected with prefilter
        """
        iterative_pick_subsampled_open_reference_otus(
          input_fps=[self.fasting_seqs_fp1,self.fasting_seqs_fp2],
          refseqs_fp=self.pick_ref_otus_refseqs2,
          output_dir=self.wf_out,
          percent_subsample=0.5,
          new_ref_set_id='wf.test.otu',
          command_handler=call_commands_serially,
          params=self.params,
          qiime_config=self.qiime_config,
          prefilter_refseqs_fp=None,
          prefilter_percent_id=0.80,
          step1_otu_map_fp=None,
          step1_failures_fasta_fp=None,
          parallel=False,
          suppress_step4=False,
          logger=None,
          status_update_callback=no_status_updates)
        
        for i in (0,1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out,i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out,i)
            self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
            # No final failures file should be created
            self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
            
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out,i)
            self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
            
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out,i)
            self.assertTrue(exists(repset_fp),"Rep set doesn't exist.")
            
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out,i)
            self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist.")
            
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(set(otu_table.SampleIds),
                         set(['PC.355','PC.481','PC.636','PC.354','PC.635',
                              'PC.593','PC.607','PC.356','PC.634']))
        self.assertTrue('r102' in otu_table.ObservationIds,\
         "Reference OTU (r102) is not in the final OTU table.")
        # at least one reference and one cleanup otu from the first iteration are in the 
        # otu table - this may very rarely fail by chance
        self.assertTrue(array([o.startswith('wf.test.otu.0.Reference')\
          for o in otu_table.ObservationIds]).any())
        self.assertTrue(array([o.startswith('wf.test.otu.0.CleanUp.Reference')\
          for o in otu_table.ObservationIds]).any())
        
        # One archael sequence from second iteration formed a new OTU
        self.assertTrue(array([o.startswith('wf.test.otu.1.Reference')\
                               for o in otu_table.ObservationIds]).any()
                        ^ #xor 
                        array([o.startswith('wf.test.otu.1.CleanUp.Reference')\
                               for o in otu_table.ObservationIds]).any())
        
        # No singletons in the final OTU table
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        
        # number of otus in otu table equals sum of seqs in rep_set.fna
        self.assertEqual(len(otu_table.ObservationIds),
                         count_seqs('%s/rep_set.fna' % self.wf_out)[0])
        
        # confirm that tax/align/tree files exist
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        self.assertTrue(exists(tree_fp))
        self.assertTrue(exists(aln_fp))
        num_tree_tips = len(LoadTree(tree_fp).tips())
        num_align_seqs = LoadSeqs(aln_fp).getNumSeqs()
        self.assertEqual(num_tree_tips,num_align_seqs)
        self.assertTrue(num_tree_tips > 30)
        # OTU table minus pynast failures has same number of otus as aligned
        # sequences
        otu_table_fp = '%s/otu_table_mc2_w_tax_no_pynast_failures.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp))
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(len(otu_table.ObservationIds),num_align_seqs)


    def test_iterative_pick_subsampled_open_reference_otus_parallel(self):
        """pick_subsampled_open_reference_otus functions as expected in parallel
        """
        iterative_pick_subsampled_open_reference_otus(
          input_fps=[self.fasting_seqs_fp1,self.fasting_seqs_fp2],
          refseqs_fp=self.pick_ref_otus_refseqs2,
          output_dir=self.wf_out,
          percent_subsample=0.5,
          new_ref_set_id='wf.test.otu',
          command_handler=call_commands_serially,
          params=self.params,
          prefilter_refseqs_fp=None,
          qiime_config=self.qiime_config,
          step1_otu_map_fp=None,
          step1_failures_fasta_fp=None,
          parallel=True,
          suppress_step4=False,
          logger=None,
          status_update_callback=no_status_updates)
        
        for i in (0,1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out,i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out,i)
            self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
            # No final failures file should be created
            self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
            
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out,i)
            self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
            
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out,i)
            self.assertTrue(exists(repset_fp),"Rep set doesn't exist.")
            
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out,i)
            self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist.")
            
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(set(otu_table.SampleIds),
                         set(['PC.355','PC.481','PC.636','PC.354','PC.635',
                              'PC.593','PC.607','PC.356','PC.634']))
        self.assertTrue('r102' in otu_table.ObservationIds,\
         "Reference OTU (r102) is not in the final OTU table.")
        # at least one reference and one cleanup otu from the first iteration are in the 
        # otu table - this may very rarely fail by chance
        self.assertTrue(array([o.startswith('wf.test.otu.0.Reference')\
          for o in otu_table.ObservationIds]).any())
        self.assertTrue(array([o.startswith('wf.test.otu.0.CleanUp.Reference')\
          for o in otu_table.ObservationIds]).any())
        
        # At least one otu from the second iteration exists
        self.assertTrue(array([o.startswith('wf.test.otu.1.Reference')\
                               for o in otu_table.ObservationIds]).any()
                        or 
                        array([o.startswith('wf.test.otu.1.CleanUp.Reference')\
                               for o in otu_table.ObservationIds]).any())
        
        # No singletons in the final OTU table
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        
        # number of otus in otu table equals sum of seqs in rep_set.fna
        self.assertEqual(len(otu_table.ObservationIds),
                         count_seqs('%s/rep_set.fna' % self.wf_out)[0])
        
        # confirm that tax/align/tree files exist
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        self.assertTrue(exists(tree_fp))
        self.assertTrue(exists(aln_fp))
        num_tree_tips = len(LoadTree(tree_fp).tips())
        num_align_seqs = LoadSeqs(aln_fp).getNumSeqs()
        self.assertEqual(num_tree_tips,num_align_seqs)
        self.assertTrue(num_tree_tips > 30)
        # OTU table minus pynast failures has same number of otus as aligned
        # sequences
        otu_table_fp = '%s/otu_table_mc2_w_tax_no_pynast_failures.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp))
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(len(otu_table.ObservationIds),num_align_seqs)
        
    def test_pick_subsampled_open_reference_otus(self):
        """pick_subsampled_open_reference_otus functions as expected
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  prefilter_refseqs_fp=None,
                                  qiime_config=self.qiime_config,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates)
        final_otu_map_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        
        self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(exists(final_repset_fp),"Final representative set doesn't exist")
        self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist")
        
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        self.assertEqual(len(otu_table.ObservationIds),count_seqs(final_repset_fp)[0])
        
        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.pick_ref_otus_refseqs2)[0]+
                         len([o for o in otu_table.ObservationIds 
                              if not o.startswith('r')]))
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map. This is done on the OTU map singletons get filtered
        # before building the otu table
        final_otu_map = fields_to_dict(open(final_otu_map_fp))
        self.assertTrue('r102' in final_otu_map,\
         "Reference OTU (r102) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU5' in final_otu_map,\
         "Failure OTU (wf.test.otu.ReferenceOTU5) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU1' in final_otu_map,\
         "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU1) is not in the final OTU map.")

        # confirm that tax/align/tree files exist
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        self.assertTrue(exists(tree_fp))
        self.assertTrue(exists(aln_fp))
        num_tree_tips = len(LoadTree(tree_fp).tips())
        num_align_seqs = LoadSeqs(aln_fp).getNumSeqs()
        self.assertEqual(num_tree_tips,num_align_seqs)
        self.assertTrue(num_tree_tips > 30)
        # OTU table minus pynast failures has same number of otus as aligned
        # sequences
        otu_table_fp = '%s/otu_table_mc2_w_tax_no_pynast_failures.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp))
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(len(otu_table.ObservationIds),num_align_seqs)

    def test_pick_subsampled_open_reference_otus_suppress_assign_tax(self):
        """pick_subsampled_open_reference_otus functions without assign tax step
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  prefilter_refseqs_fp=None,
                                  qiime_config=self.qiime_config,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates,
                                  run_assign_tax=False)

        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table_w_tax_fp = '%s/otu_table_mc2_w_tax.biom' % self.wf_out
        self.assertFalse(exists(otu_table_w_tax_fp),
                         "OTU table w tax exists (it shouldn't).")

    def test_pick_subsampled_open_reference_otus_suppress_align_and_tree(self):
        """pick_subsampled_open_reference_otus functions without align/tree step
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  prefilter_refseqs_fp=None,
                                  qiime_config=self.qiime_config,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates,
                                  run_align_and_tree=False)

        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        tree_fp = '%s/rep_set.tre' % self.wf_out
        self.assertFalse(exists(tree_fp),
                         "Tree exists (it shouldn't).")


    def test_pick_subsampled_open_reference_otus_no_prefilter(self):
        """pick_subsampled_open_reference_otus functions as expected without prefilter
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  qiime_config=self.qiime_config,
                                  prefilter_refseqs_fp=None,
                                  prefilter_percent_id=None,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates)
        final_otu_map_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        
        self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(exists(final_repset_fp),"Final representative set doesn't exist")
        self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist")
        
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        self.assertEqual(len(otu_table.ObservationIds),count_seqs(final_repset_fp)[0])
        
        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.pick_ref_otus_refseqs2)[0]+
                         len([o for o in otu_table.ObservationIds 
                              if not o.startswith('r')]))
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map. This is done on the OTU map singletons get filtered
        # before building the otu table
        final_otu_map = fields_to_dict(open(final_otu_map_fp))
        self.assertTrue('r102' in final_otu_map,\
         "Reference OTU (r102) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU5' in final_otu_map,\
         "Failure OTU (wf.test.otu.ReferenceOTU5) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU1' in final_otu_map,\
         "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU1) is not in the final OTU map.")


    def test_pick_subsampled_open_reference_otus_parallel(self):
        """pick_subsampled_open_reference_otus functions as expected
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  qiime_config=self.qiime_config,
                                  prefilter_refseqs_fp=None,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=True,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates)
        final_otu_map_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        
        self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),\
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(exists(final_repset_fp),"Final representative set doesn't exist")
        self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist")
        
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        self.assertEqual(len(otu_table.ObservationIds),count_seqs(final_repset_fp)[0])
        
        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.pick_ref_otus_refseqs2)[0]+
                         len([o for o in otu_table.ObservationIds if not o.startswith('r')]))
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map. This is done on the OTU map singletons get filtered
        # before building the otu table
        final_otu_map = fields_to_dict(open(final_otu_map_fp))
        self.assertTrue('r102' in final_otu_map,\
         "Reference OTU (r102) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU5' in final_otu_map,\
         "Failure OTU (wf.test.otu.ReferenceOTU5) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU1' in final_otu_map,\
         "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU1) is not in the final OTU map.")

        # confirm that tax/align/tree files exist
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        self.assertTrue(exists(tree_fp))
        self.assertTrue(exists(aln_fp))
        num_tree_tips = len(LoadTree(tree_fp).tips())
        num_align_seqs = LoadSeqs(aln_fp).getNumSeqs()
        self.assertEqual(num_tree_tips,num_align_seqs)
        self.assertTrue(num_tree_tips > 30)
        # OTU table minus pynast failures has same number of otus as aligned
        # sequences
        otu_table_fp = '%s/otu_table_mc2_w_tax_no_pynast_failures.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp))
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        self.assertEqual(len(otu_table.ObservationIds),num_align_seqs)


    def test_pick_subsampled_open_reference_otus_suppress_step4(self):
        """pick_subsampled_open_reference_otus functions as expected wo step 4
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  prefilter_refseqs_fp=None,
                                  qiime_config=self.qiime_config,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=True,
                                  logger=None,
                                  status_update_callback=no_status_updates)
        final_otu_map_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        
        self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
        self.assertTrue(exists(final_failure_fp),"Final failures file doesn't exist")
        self.assertTrue(exists(final_repset_fp),"Final representative set doesn't exist")
        self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist")
        
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        self.assertEqual(len(otu_table.ObservationIds),count_seqs(final_repset_fp)[0])
        
        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.pick_ref_otus_refseqs2)[0]+
                         len([o for o in otu_table.ObservationIds if not o.startswith('r')]))
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map. This is done on the OTU map singletons get filtered
        # before building the otu table
        final_otu_map = fields_to_dict(open(final_otu_map_fp))
        self.assertTrue('r102' in final_otu_map,\
         "Reference OTU (r102) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU5' in final_otu_map,\
         "Failure OTU (wf.test.otu.ReferenceOTU5) is not in the final OTU map.")
        self.assertFalse('wf.test.otu.CleanUp.ReferenceOTU1' in final_otu_map,\
         "Step4 (clean-up) OTUs are present in the final OTU map but shouldn't be.")

    def test_pick_subsampled_open_reference_otus_parallel_suppress_step4(self):
        """pick_subsampled_open_reference_otus functions as expected in parallel wo step4
        """
        pick_subsampled_open_reference_otus(input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=self.params,
                                  qiime_config=self.qiime_config,
                                  prefilter_refseqs_fp=None,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=True,
                                  suppress_step4=True,
                                  logger=None,
                                  status_update_callback=no_status_updates)
        final_otu_map_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        
        self.assertTrue(exists(final_otu_map_fp),"Final OTU map doesn't exist")
        self.assertTrue(exists(final_failure_fp),"Final failures file doesn't exist")
        self.assertTrue(exists(final_repset_fp),"Final representative set doesn't exist")
        self.assertTrue(exists(new_refseqs_fp),"New refseqs file doesn't exist")
        
        otu_table_fp = '%s/otu_table_mc2.biom' % self.wf_out
        self.assertTrue(exists(otu_table_fp),"OTU table doesn't exist.")
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        for row in otu_table.iterObservationData():
            self.assertTrue(sum(row) >= 2,"Singleton OTU detected in OTU table.")
        self.assertEqual(len(otu_table.ObservationIds),count_seqs(final_repset_fp)[0])
        
        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.pick_ref_otus_refseqs2)[0]+
                         len([o for o in otu_table.ObservationIds if not o.startswith('r')]))
        
        # spot check a few of the otus to confirm that we're getting reference and new
        # otus in the final otu map. This is done on the OTU map singletons get filtered
        # before building the otu table
        final_otu_map = fields_to_dict(open(final_otu_map_fp))
        self.assertTrue('r102' in final_otu_map,\
         "Reference OTU (r102) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU5' in final_otu_map,\
         "Failure OTU (wf.test.otu.ReferenceOTU5) is not in the final OTU map.")
        self.assertFalse('wf.test.otu.CleanUp.ReferenceOTU1' in final_otu_map,\
         "Step4 (clean-up) OTUs are present in the final OTU map but shouldn't be.")

    def test_pick_subsampled_open_reference_otus_invalid_input(self):
        """pick_subsampled_open_reference_otus raises error on refseqs in params file
        """
        self.assertRaises(WorkflowError,
         pick_subsampled_open_reference_otus,input_fp=self.fasting_seqs_fp1, 
                                  refseqs_fp=self.pick_ref_otus_refseqs2,
                                  output_dir=self.wf_out,
                                  percent_subsample=0.5,
                                  new_ref_set_id='wf.test.otu',
                                  command_handler=call_commands_serially,
                                  params=parse_qiime_parameters(
                                   ['pick_otus:refseqs_fp %s' % self.pick_ref_otus_refseqs2]),
                                  prefilter_refseqs_fp=None,
                                  qiime_config=self.qiime_config,
                                  step1_otu_map_fp=None,
                                  step1_failures_fasta_fp=None,
                                  parallel=False,
                                  suppress_step4=False,
                                  logger=None,
                                  status_update_callback=no_status_updates)

subsample_otus_seqs = """>sample1_1 r0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTTGACTTGGTGGGCCGTTACCCCGCCAACTATCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAGTGAAACCATGCGGAAATACTATACTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATCAA
>PC.607_1 r1
CTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTCCGCCACTAGGGCG
>PC.607_2 hits new otu
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTGCCCCGCCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCTCACCATGCGGCAAAATAATGTCATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTT
>PC.607_1000 85% id to r100
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACCGGTCTGTACGGTAAGGTTACCTATACGCT
>PC.607_1001 85% id to r100
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACCGGTCTGTACGGTAAGGTTACCTATACGCT
>PC.607_1002 85% id to r100
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACCGGTCTGTACGGTAAGGTTACCTATACGCT
>PC.607_1003 85% id to r100
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACCGGTCTGTACGGTAAGGTTACCTATACGCT
>PC.607_1004 85% id to r100
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACCGGTCTGTACGGTAAGGTTACCTATACGCT
>sample1_100 random sequence
AGCGGATGACCGTGAAGCGTATGACGTATCCTAACAGACCCTCATAACAGGGTGAATTTGATACGTCCTCAGTCGCAAATAGCCCAAGGTCTCCTCACTGGCGAATCTGCAAGGCCGTTCCGGGCGGCTCCTTGGTTTAAGCCGCCGGAGGCCCCCACTAATAAGTGTGATAATTTGAACTATAAGCGGACGTGTTGTGCAGCGTTACATCAGAAGTGTTTGATCGGGTAGACGCTCCCGGAAGTCACGG
>sample1_101 random sequence
TATGGGCCACAAGGTCACATGGCTCCTGCAATGCTCGAAGATAACGTTGGGAAGCAATTGCCAGAAGTAGCTCAGTGTGGGCATTTCGTTAAGTTGCCCGTGGTTATGCAAGTTGTTAGCTCCAGAATACAGTTAGCGGGGAGAATCTGGAAGTAAACACGTTAGTTGGCCGAAATGGGGCTCTTGAGTTCGCAAACTAAGTAATATCTTAGGAATGGGCGGCTGACATCGAGATTTTATTTTCAGGATA
>sample1_102 random sequence
CATCGCTCAGTATCTCCCTTATCCCCAGGCAAGCGAGCATACGGGCCCCGACGTCCTGTTTCAGTTGATGTGTGTTCACTGCCATATGAGGCTTGTCCCACCGGAGCGTGCAAATTACGTTAGCCAAACGGAATCGCGTGACAGGAGGGGCAGATCTCTGGCCCATCGGCCTTATTCAACGGGGTTAAACCTCCTTTCAGCACACCGTCTGAGCTAGGCGCACGACCTATGTAAGCGGATTAGGGGATCG
>sample1_103 random sequence
TTGACCTACACGGACTTTACTCGCGCATGTCTGCACCCAAGGGACCAAAGCGTTGTAGCTAGACGCCTCCTCTACCTATTCGAATATTGCTACACCCAAAATCTTCCATTCGTTCGAGGCACATGCGAGTGGCAACTCGGAGGATAGAGGTGCCTTATGAGGGTCATTTTGCTTGCAGCTCTTGCCTGGATGCCACCCAGCCCGGCTATAACATGAAAACTGTACGCCCCCAAAGTGTCAATATATGTGG
>sample1_104 random sequence
TAGTAATTGAGTAGGGCACTGGCAAAACCGACTCCACGCCCCGGATCTCATCGCACTATGCACTGGCAATTGCTGGCTGTTTTCGATATGTTTACTTTGGAACTCCTGAGGCAGTCCTGGTACCGAAGTGTCGGACTCCCAAGATGTCTTTACATAGAGTCAGCAAGACGGAGCTACGCTACATGTAACGGGGCCCGCTGTAATCCCTCCCCCAGAATCGTCGGCGCGGCTGATTAACGTCCTCACTCGA
>sample1_105 random sequence
ATATAGAACCGCGCGTTCTAAGTGCGTGCACTTAAGTGACGATTCATGGCATGTCCTGAGCTTGCCACCGCGCTTCCAGGCGTAGTCTCAGTAGCTTGACTGTACTTTGCGGGCGCAGAGGTCCCCAGTGCCAAGCCTCAAGTGTTCTGAATTACCCCACCGATCTTTAGCGTAGGTCCCGCACTCGCATATCGTAGAGCGACTCGTAGTAAAAGGGGTAGGCTTTAATAGAACACCTTTAGCAATAGCG
>sample1_106 random sequence
CCCGCACCAGATCACGAGTATCCACTCGGTCACACTGCGAGGCAATCTGTGTCAAAGGTCCACTGTAAGTATGAAATGTCCCCGTTCTCGTTAGATCCGGTAGGCATCAGAAGATCACTGGATCTACCCGAACTGTATGTCACGGTCATATAAGCGACTAAATGATCCTGCGGACACGTGCAACGAGCACGGTGATGGATGGGGAATACAGCACCGAATCATGATATCTCTTGGCGTACGTAACATAGTC
>sample1_107 random sequence
CTTGGAGCCTATAAGGATAGGAACAGCCCTAATAAGAATGAGTGCCGAACACATTGTGCAGGGCTATTGTTAGTATCGCACCCGGAGGGGAGCCTGAAATGCAAACTGGCCTCGACATCTGGATATCGCACAGACCCCCGCAGTGATTTCACGGTAGACAGTAAAGCGTTAGACACAGCGTCCGGATTCGGACTAAGTTCCCCTCAGTTTAGTTGGATCTTATTCCTATGACGCGACCTCCACTGCTCAG
>sample1_108 random sequence
GAGTAGCTATAACTCCGCGCGTTGTCTGAGCAATAAAAAGCAACTTTATTTTCGTAGCACGGGCTATCCGAACCAGCGTAATGGACTCATGGCCTAGCAGATGGGACAAAAGTAGCTCCCTGCGCCCAGCTGTATGGGTAGGACGTACAGGAGGAATTTGATTTTATCCGGCTAACTTCCGAGAAGTCGTCTAGGCGCGTTCCATATAGGATACTCGCGGTTTGCTCACGAGGTACCAACCGTTACCAAG
>sample1_109 random sequence
AAGTGTGTTTCCCGTCGCATGGGATTTGGCTAGTCCTGGTTTGTTCGTCCGTGTAGAAGTCATAATTATTTAAAAGCAGCTTTTCCTAAGTGAGTTTTTGTCGTGCTGACACATATCCAGGGAGAGCGCCTGTATGCGTGATAGACTACCAAACAATGTCACCCAATCCATACACAATGGTGCTAGAGACTCAGGACAGGTCGTATGCCATGAGGCCATTCTGCCTTCAGGTCCAGGGATCCCTTATCGG
>sample1_110 random sequence
CCTTATTCTTACCTCCCACGGACGCCACAATTAGCCAAGTCACCCGTACCTAATGGTTCCGAAAGGAGAAGCCTGCCTCACGTGGACCCAACCGCCCTGACATACGGAACAAGCTTCAGGTAAGCCGGATGATACATAATGCCCCTAGGATATAGTAGACTTGCGTCGGCAGTAGCTGCCACGTCTCACTTCGTGATAAGAGAGCATCAAATCAATTTCACCGACCGATTTCGCCCGGGACTATGCGGTA
>sample1_111 random sequence
CCACACTGAGGCCAGGTCGGTTACACCTTAAATTCTAACGCCTGATAGGTTACTTGGGGTTGCCCTTTGCAGGCGCCTAAAAGATACGATGGAAGACGAGCGTAGGACGGCCGAACACAGGCCTTCGAAACTTAATACCAAATGAGGGACGGAGGGAGCAGCTCTTCGTCGAGTACAGCGCTACGACCGGCATAGAATAACGACTACACTTTTCCAACTGAAAAGATAGCTTAGAGAAGGTTTGCGAGGT
>sample1_112 random sequence
CCTTGACGTATCCCAAAGACCACTTGGTGTTGCCCTACACGGGCAAACGGTTGGTTAAATAGCACGCTCCAGTGAAGGGCTTTCGTATGACCGCGCGCGCCGACCTCGTACACGTTATCGGAGCCGTCTAGTTCTGTCCAGGTCCCATCGTCCTGCCGAGTTATGAGATAATGAGTAGTACACGAGGGCACCGATCTGAAGCGTTGCTTATATAATAACGGCGGGTGTTAATGAGTGAACACCTAGAACT
>sample1_113 random sequence
TAGTACAGGCCCTATGCCCAACATCACCCGCCAGTGTGCCGCTACTTTGTACTGGATCATTTATTTTCGATCGAGACCGGGTCTAATCGCTGGCGGGATGGTTAGAATTGGATTTGATTTGCCTACATATATGACCAGAGTGCCGTTGAGGAGGAGGATGGACACGCCGGACCTTCTCGAAGTGTTGTTTGCTCGACATACTAGGCGTTCTTGTTGTATGGTTCCCCCTGATAACACCATGGTCGCGATT
>sample1_114 random sequence
CCCTAAGATTAACCCTGGGAGGGTTTGAACTAGCGCCGGTCTACTACGGTGAATTCCTAAAGTCGAACTATATAGTATCGGTGACTTGACAGGAGAAGGCCCAGCTAGAGGTAGTGTTTGCTAGAGAGGGGCCGGGAAAACCAGCATAACACAATTAACGACAAGCATATTGATCGTATAGCTTGGAAGGAAATATCCTGGAACCGCATCGTGGTAAAACTCCGGTAGCGAACCGCAGAGAGCTTTTGGG
>sample1_115 random sequence
AATTCTTTTTATTTGGAGGCTCTGACCTCTTCGAGGCACCTCTTGAGCGCCTTCAGGTTTGCTAAAAGTTAGTCGTTAGTCCGGAATCCCACCATCAGAGGTAGTAATATCTTTTAGACTACGGCTTCACATAGTGTCGCTCCATTGGTGTGCCTCGGTGGCCACTGAGAACTCGGCCTGCAGGACTCACTGCTTCGGCAATTATAATCGACGTCTGAGCGCTATCAAAGCGGAGCTGACGAAAGTTGGG
>sample1_116 random sequence
GGCTCGTGATCCTATAATAGACGTATGCGGTAGGAAGTTGGTCCATAGAAGGTTCACGCGGGTGTGGCCCCGACGTAGGGATATACTTTTTGAGTCAGAGGAAGCGAATGTTCACTCGAGCTTGATCTACGACCATTCTACGTTGTTAGAGTCCTCCAGATCTCTAAATGAGTGCCTGGTCAACGTTTCACTCACTCACTTCTTCCATAAGAATCTGTCCAATCGGAATGGATCAACCCAGACTTGTGCT
>sample1_117 random sequence
AGGAGCCCGCTAGTGTGACCGGATAAAACGATTCATAGTCGGATTGGGAATTGCATTTGAGTGTGCCCGAACAAGATAGATAGCTACTTACCCAGTTCCCGCTGACACGAGGGTCCTGTCAGCTAACTCCGTTGGGAAGCCCCAGACATCCCTGTGGTCCCATACACGTATTAGGCTCCAAGTGCACGCCGTCAATGTTATTGATCAGTGATATTCCCCTTTGGGGCACCTGGCCGTTCAGAAGCATCGG
>sample1_118 random sequence
AAACGGTACTATGCAAATCAATGACCGGATTGAGCAGGTCAAGCTTAATATCCCCCGCATGAAACGCATACATTACACACACGGAGACATGTCTCCCTCTCAACCAGTCACGTCTGAAAGTTTGTCCCGATTGACCGGCGCCCTAAACCGCGCACATTAGTGAGACCGTCCCCTACCACCCACCGCGACTTAACAGATCTAGGCAGGGCCGAACTGACCTTCCGGGGTGTTTTCTGTAATCGGCAACCGC
>sample1_119 random sequence
CACGTCGAGTCTGATGCGCGGATAGGTTTCAAAGCTTAGGATAATATATAGTCGAACCTTAATGATGGAACTGCCGAGTTTTTAAGACTTGTAGACATAACGGAATTCGCCGGTATGACGGATGACGCGTCTGCTCCAGCTGTAGTATTTTGACATCGAAGCTTCCAGTGTGGACACAGGTGAGTGTACGTAGCCGTGCACCCTTTCGACTTTAGACCCCGATAATTGGGGACCCTGTAAACGAAGCCCC
>sample1_120 random sequence
CTCATGCGCACCTCGCATTTCAACAATACATCGAATAGGCTCCGGTCTAGGTGATCATACTGGGGGGCCAGTCGTCCGGAGTGGGCCAGGGCTTAATAACACTCTGTTGTCTCGACTAGATGTGGCTGGACAGTCGACCAAGACGTACTACACCCATGTTATCGGCCAAAGCTGTACTTTAACTAAAATGACCCCTTAATTCCTGTTCCGACACAGGCTTACACCAGGAGATTGTCTTTCTTGGACGCAA
>sample1_121 random sequence
AACGCGATCTGGGCGTTGCCGTATTTTATCTTGCCGTAAGGTCCTTCACGCCCGCACATGGCATTGGTCGATCGGGGTGAGTATTGACCAAACCTGTCGCTAGATATGGATGTAAGTTTGGGTCGTTCTCAGGCTAGATCGTATACAGGGAGGATAGAGTAGATGGGGCCGAACGCCGAGCTGGCTCGACTATAGGGGCGACGCGTGCGAAGGGGGAAAGCGCGTTAGTACACCCCGCCGATAGGGTGGC
>sample1_122 random sequence
TGTTGGTAATGCGCAAGCAAAGATACAAGCCGAGTAGCCCTTGTTACAAGGGGGACTGGCATTTCGGACTCATCCGAGAAGAGGGGCTGCAGACATTCTACGGGCCAGGGTCTACACGGTCCTTGGTGTATGTGCGCGGCGGAAAGTTGCTAACGACTCCAGTTAGCTGTAGGGTCTGTCTGTCATATGCATAGTCGCGTCTGATTGTGCTCTCATTTTCTAGAGACAAGCGTTTGGTTACTATGCCAGG
>sample1_123 random sequence
ACTCGCGCGATGACTGACCCTGACTGAAAAGAGTCCGATGGTAACCCCATGGCGATCATGGGTAGTACTACGAGTAAGGCCACACGAAAGAGACGCGGTGGCAACCTGCGTCCGCTTTCTGAATGCTAACTTCTTCAGATATTGTCTGAGTTAACTGGTGTGTGTAAAGTTGCCGAGGCTGGGCCCAAACAGCAGCGGGAGAGTTTGTCTACGAAGTTGGGTAGGTGGCACCGTGTGACGCTTACATAGT
>sample1_124 random sequence
TTCTAGCTGGCGAGATCCAAACGGACTCATGTAGGCATACGTGGGTCGTAGTATTCCTTTTGTGATGGCAGCGTGGCTCAACACTTGGTGGGCGCCGTAACTCTAGGTCGCCACGGCTTTATTTTATGAACGCCTAGCGATGTGAGTATCTCTGGGACAATGAACGACGTAGTTCCATTGGAATCCCAAAGCAAGTGGTAAGGGGAGCCTTCCGCCTCGACACTCTTACGTATAAAAGTCTTCATGTTCA
>sample1_125 random sequence
AGACACCTCACTGTTCGAAAGGGAGACGTGATGTCGCAAAGTTTAATTGAGTGCATCACCAATTGCTATAGCGATGTCTGTTCAGTATGGGGTGATCAGAATCTGTATAGAATTTTCGCTGCTCTCAGATACCTCCTACTCACATTATGAGCTCCTGGACCGCCGATAGACCTTAACAGGGGTCGAGGCATCCAAACGGCCCTTCAGGAGCCACGCCGGAGCATCGTCAGTGCAAGAATTGGTGACGTGC
>sample1_126 random sequence
GTTCCTAGCTAGCCAAGTGGCTTGGTGGGTAGAGGGCCACCCATACCACAATTCACTTATGAAGGCCAGCTATCCAGACGCGTTCCCGTAGTTGTATATTATCGTAGCTTGCTTTCCAAATGCTTATCTAATCAATTAGGGAGTGGGCTTAGGTCACGCATTTAAGTTAATATGGGACAGTGGTCCTTGTAATAGAACGAACGGACCATGGTCAGCTTTGCTTAGCCAACATAAGAACATAAACGCACGC
>sample1_127 random sequence
ATGGACAAATGACCAGCACCACATAGCAAGGGTCGTAGCAGTTTCGGGTGGGCGTGGATAGTCACGGCCAGTCAAAGGTGGGGGCCTTGCATAAATCTTCAGCCGTAGCGCGGGCTGTTGTAAGGTCTTGGCTACAAAAGCGTAAAAGTAAGAAACACATGACAATACCCCGCTGCGATACCTGCCCATCTATTCGAGTAGAGATTCAGCGGATTTCTAAACGGTAGGGAGTCGGAGTTGTAACAGGCAA
>sample1_128 random sequence
CAACTCTCTCTTATGGTTGGGGTGTTCCTTGGCGAAGTCAGCTCTTACGGTCCCGACAGCTAACGGTAACTGGCTGCCTATATGCAAGGCGTACTCCATTACCTTATCCCGAAGCACGCACGCATAAACTTCCGCTTCGTGCTTGGTAAAGAATGACGGCATCTTGGTCACATTTGTGCAACCGGTTCTTTAGTTATCGGTTGCGCGATATATATCTCACTGATTCTAGTGTTAATACTCTGTACCCTAT
>sample1_129 random sequence
ACCGAGCTAGCTCTAATTTAGAAAGACCGACAACAACACTCGGTGTTTCACAAAAATAAGTCATCAATTAAGAAACGCCGAAATCATCAGAAACAATGCTTGAATAGTTGTGCAGATGTGAAACCCGGGTCCTGGAAAGTGTCATGAGGAATGGTGGTGTGCCGAGGGACTGCTGGTTAATATATGGAGCCTCGCCGTCCTGGCGCTGCGTATCTAGCCAGGCGAGAGGGTCGCCTGGTTCAAAGATATC
>sample1_130 random sequence
AGTAATTGACCCACTCGGGCTGACGGAAGAAGGCTTAAATGGCATCCAACCCTGAGACAGGCCCCGTGGTATACTCACCGAGAGTGACGAAATCCGACGATCTGACGTTTAAACGTGCACATGAATCGTGGTAATACGGATTGCGGCAAGGCATACCCACTAAAGCAAAAGGTCGCAACTGGCTATGCCGATACTTTGTGGGTTCGCGGCGAAAGTGAGATATCACACTGTCATGTTTTCAACGTAGGAG
>sample1_131 random sequence
TACACTGAGCTGAAGCATTCCATCACTGCGAGTTTGCAGACTTTGATCGGGGCTGTAATACCAATAAAGAACAAGCACAGGCCTAAAAGTTACTCTGCTGTGTATCGGGTCGCGACAGTGGAACTCTTGTTAGGGATACCCCTGGGAGGATCGTTCCGCAGTATGTGCTCGGCTAGGAACTAACTGGTCCATGCTAGGTTGCCAATATTCGCCAAATGGTACGTAGCACCGGACAGCTGTCACACATATC
>sample1_132 random sequence
TAAAGGTTCGGCGCGGGAGGAGACTCGCCTGGGTTGCAGGTACCATGGTTGCCTCTTTTACAGCTTATCAGATAAGGTGGACAGTGTGTGTGGCACAACGTGACTAGAAGAACATCTGTACTACACAGCCTTCATCCTTTATGAGTGCACCAACAAAGCTCCCCGAACTGCTTGAAGCCATTGACGGTCAGACGACCAATCGAACGAGGTGATTCGAACGTACACCACTCCGACCCCACCTCACCTTACA
>sample1_133 random sequence
GACACTGCATGGGCGGCCCCATAGCCAGTGTCAAAATTCACGTGTTGATTAAACAAACGTAGGACCATTAAACTACAAATAACCACACCCTAGCAAGTTACAGCTTGACCCTACTCTAATCCGCGGCCATTTGAGAAGTCCAACTCGGTCTCCCAGAGATTCGGAGTGATCTTATACGGATGGCACATAGAAAAGACGGTGATTTGAAGAGTTAATGGCTGAATAGCAACTCCACGAAGTGCGCCTATGT
>sample1_134 random sequence
AAGCCAGACGGCTCAAATATAAATATAGCCTCTAGCACCTTTCTGCCGGCAGCGACACCGTCACGCCGGGAAGCGGAAACATAATTACTACGAGATCTCGAAGCCAATCCGTTATTCATCCCTGGCGATCGACGCATCTATTCACGGATCGATACAAAGGCTATTGAGAAACCACAGCCGCGAGGACCGCCGACCTCCGTGCTGTCTCAATTACTAATTTAGATCCATGCGCCTATTTCAATGTGGCTAT
>sample1_135 random sequence
GTAGATCCCATTGGAGCGAGTGGCAGCACAGTCTGCGTTGAAGGCAGGTGTCGGCGACCTGGTCGAAAGGCATGCCGCATCTGCGGCTAGAGGGGGCCAATTCGCTGTTTAGCCCGCAACGAGCTCAGAAATTAGAGTAAATGGTAAGAAGAAAGTTAGCCGAGGCGCATCTTAATGCGCATGGATACGGGCGTACGTGCTCGGTGTGTGGTTATTGGCTCGTATTCCTTTAGCGGGGGGCTAGTACTCA
>sample1_136 random sequence
TGTGTCGGAAACATTGTGCCTAAGTCCAAATGACCAACCCCCTCAAGCTCTTTAATTTCATTGACTGCACACGGCGGTCTCTCCGCCAAACCTTTCGAATCGCAATCAGCAAGTTCTGTAGCGCACGATGAGGTAGTATGATCGGACGGGTCCCTGAAAAAGTTAGTTGTCGGCCTTCGTAGTCTTTTCAAAACTTTCCATGGGTGGACTGCATGCGCGAGATCTGGACACCGTAACCCAAATCGCAGCG
>sample1_137 random sequence
GAAAGGTAGAGTCGGGAGTTGTGATGTCGTGTCGGCCATGAGTTGACATCCAGAGAAACGCTCTACTGCTCCTCGGGCATCGATGTGTTTCACTAGGAATCGATTTGGCGCACAGTATCAACAGCCCACCTGTAAAATTCGAATAGCCTGAAACCGCCCTACCTGTTGGGACGATTTGGCAGATACTCGAGAAAGGGAACTTCGGTTTTTTCCGGGGATTTGTATTCAGCACAGAGATAAACAGCTTGGC
>sample1_138 random sequence
TACACAGCGGGATCAACAAACCTTCCTGACTCCGTAGCAGCCCAAGGGTGCAGTGTGGTTGGTATAATTTATCCTCACTTCTTCGGCGCCCGGGGATTATCCCCGCCTGGACGACAGGGTTTCAGGGACGCGTAATAACCCTTAGAAACGAGAAGGAGAAGATTCACTACACAGATGGCTCTACCTAAAACAACTAGCACAACTGATCGCCTCTCTTCGTCGAAAGGCGAGAATGAACGAGACTCCCATT
>sample1_139 random sequence
GGCACCGCGGGCAAATGTCCCCTGCGTAGTCAATACAGCAAAACTCTCATGGCCCTTAGCGAATGAGTACCACATAGAACGCGCTAGTAGGACTTTATACCTCATAAGTAATAACTTCATTCTCACTCTACTCGGTGACCTATCCGGTGCGTTAACCGCTGAGGCGGCGCGTGCAAAGCGGTGGTTAAGGCGACAACGTGCCTGAGAGGTCACAATACAACTTCCAATACTGATGGGGCGTAGTTTAAAG
>sample1_140 random sequence
ACCGGTAGGTTATGGGTGTAAATCCTTTGCCAAGAAACCTCCGTAGTTTTTTATCGATCTGCCCTGATTACATTGAGTAGTTCCCTTTAGTCCCCAATTCTCCTAGCGATGGGTGCCTATTTGCCCAGAAGGGTCGTGTCTGTCTTCTAGTATGAGTTTAGACATTCGAACGCCAACTCCGGTGGCTGCTTGCCGCCTGCAGCAAGTATATCCAAAGCGCACCCGGAGCTGACTACCAACTTTTGCAGTT
>sample1_141 random sequence
TTCATGAATTTTTTCCGCCCATTCTAGATTAGCACTATACTGCCCAGGTTTACTCAATCAGATGTACGCAGGCATTAAAATTGCTAGACTCACCTATAGCGTTCTCCTGTAGATGAGTAAGGGGCTTCCGTATCGACTGCAGTAATTCCGTTTCACACATGAAGAAAAGGAGTTTACAGCGAAGGCCTGGTAAGGCGGGGAGGCGCTCATCGTAAAAGTGCGAATGAGCGGTAGGAGAACCAAGGTAAAC
>sample1_142 random sequence
AAGCTGATTAGATTTATTTGTCGAACCATGCCTCATCCCTTGACGCGCGATCCCCCGGCGATGCCTGTTGTCGGACTCACGTCTTAGTTATCTAAACTAGAAACCCCGCAGCTGGGACCTCAGGGTTCCCCTTTAGACACCCCAGGTGCTCACATGTCAGGAAACTTATTACTAGCCGACTCTTTAGCTGTGGTTTCTACACTGAGGGTGCAAAAAAGCAACCGACAGACTCGGGGGCGAGGATCATGTG
>sample1_143 random sequence
AGCAGAACGCGCATCCTGAATTTATAAAGCTAAAAACTGCGTTGGGCGGTTGCTCAAATGGTAGACTATGGGGCAGCCCGCTAAACCTAAGCGACGTCCTAAAAGTATTTGGCACATGGTGCGCAAGTTTTAGCTTAAACCTTAGAACCTTAAAGAGACGGTACCTGTCATGCTGCACACCTTACGCCCTAAATGCGGTGCATTCAAGGGCAAGCACTGGGTATTAGGGAGAGATATAACTCGTCAGAGC
>sample1_144 random sequence
CGTGTGAGCGCTACGGCTATACTCCCGTAACCACACTAGAGAGTGCAGGACCGCCGTATAATTGATTTGTCAGATAGGCTATGCTTGTCGTAGGCAGCTACCACAATAACGGTCGACTGGGGTGAGGAGGTATGTCTTTCAAGCTGCAGTCCTATCGCCGTACGAAAGTCATGATGATTGCAGGTTTAGACCCGACGGAACTCCGGTGCAACTCCAGGCGCCCCTTCTAAGGGCGCGCACTAGTCAAGCA
>sample1_145 random sequence
TCTGCTGGCTGCCCACTCGCTTGAGCTAAGATGCGAAAAGAATTTAAAGTTATAATTATAAGCAAGTTGGCTAACAGCCATGTCAGCAGGCGCGACGGGCGCCCACCTACTGACTCTTTACAACGTTATTGGTGTAGCTTCCTTTCAAAAAGCGGAGGCTCAACCAGGTGTATGACACCATCGTAACCACACGGACTGAGCCCGGTGGCGACAAACTGCACACGTGAACAACACACGTCATTTCGCAAGT
>sample1_146 random sequence
ACGCTCATAGAACTATTATGCTCAGACAATCAAACTTCTCTCAAAGTCGAGCGTAGTATCAGCGCTCTGGATGCTTTTGTAATCGTAGGGCGATGGAAATGAACCTGGATCCAAACCCCAACCTATATTTCGGCCCCACTGACAGGTACTTGTCTGTTGTATCCTTTCCCGCCACGGTGTCCTAAGAGTAGTGGAGTTTTGGTTCAAACATGGCCGATGCGGGGGGGAATCCGCCGTTCAATAATACTAC
>sample1_147 random sequence
GCCAGATGTCGTCAGCGAGTGAAAGGTTACACTTTGGAGGACCACTAAACACTCCTCGTTTGGTGTACGGAGGCGAGGCAACAAGGCTCTCGGCGCGGGACAATGAATGGGGGCACTAGCGTGCGAGGATATAGCTTTTTGTAGCAGGAAGGCTTTGAATGGTTTGGCAGGGCTACAAACGCTTAGTAGTACAGAGACTCTACTACGTACAACTCTTTCGGGGCGTACGAAGCGCTGGCTAAGTATAACT
>sample1_148 random sequence
TTCGGCTATACATTTAGCTATAGCTACATCGTCGGGTCGGACTTCGATATTGACGCATTGCCGCCCCTTAATGAACACTACTCGGTGCTGGGAGGGTACTCTAGGTACCCTTCGTGGTAATCCTTTGCTTTTACGCGGTAATTAGTACGCGAAGAGGAATCCTGTGCCCAACCGTGATCCGCAAGGGGCATGCCCGAACCACGTTCTACAAGGAAAGGATGACTCGAAGTGCAGCCTATTTAAACCCTCA
>sample1_149 random sequence
GGACTATTTGCAAGCCTTAGGACTGTCGCGTCTTAGGTGCTCGATTCGGCATTATTCTTGAAAGCTGCAGCCGAATTTGATGAGAAACAGTCCCCGAGAAACGAACCTCACGCAAGCAACCGGGCGACGTCACCTTGTTGACCCCAAGCCAACTTAGTATGAGGTCCTGAACTCTTTATACGCTACCTGCTCCCATACACTCCTCGGGAGGGGGCAGATGTAGACACGTGTTTTTGCTTCGCCCATATGA
>sample1_150 random sequence
TATAGCAAATATAGCTACGCGTTGACTAGTCGCCATATAGGGGCGGCGCCATCAATGTGTAGAGTTCGTAAAGCGCGGCCCAGAGTTTGTAGGACCCGGTCCCATGGCCAGCATCCAACGCCGGTAGTCACGCCCTTTCCCGCTTAGCGAGTACAAACGCAGGGAAGCATGCGGAAGCCTGCCATGCAGCACGTTCTGCCCTACTGTGAAGCAAGGTTTAGAAGCGACCAACTGGTCAGCGCACATGAGG
>sample1_151 random sequence
CCAGCTGTATGCTTGCGACGCATACAATGGACTTGCATGGAAGATTCCTATGACGCCTTATTGTGTCGAGGGGAAGGAGGTCGATTCAGGGCACACACTTTTAGGCATCACCCAGGTCACCTTTAATCGTCTTCTGGCAACAGCGGGCGGATTTTTTGGGCAGCTCAAGCTTGTGGATAAATAATGTGCCAACGTCGAGAAACGAGACTCTACGATTACGAGGGACTAGGAAGCGGATGAAAGATCGGTG
>sample1_152 random sequence
CGACGGCGAAATGGTCGACGAAAACAAGCGAAGCAGAGGGTTACCGCTTGGGCGACAGGAGCCTGGCAACTATCACAGCAATATTGTTGGTCATTTTAGCGGCAACTGCAACCACACTCCCGTATTAGGGTAAACCAATCCCTGACTTCTTAGCGTCATTTCGTGCAGCATCAGTAACGCGCGAACCGCGGGGACATGCAGCGTAGAAAATTATGGTAGGTACCTACTGCCACCGAGAGTTTACCCGCTG
>sample1_153 random sequence
TACTTACGTGTTAATAGATCTAAGTAGATGCCAGTAGCGGTGCTAATGCTTCTCCATAGACTCTCTAATCTCCTATGCGGGTATTCCAGACGTGATTTGGACCATCGCTGTTGACCGTCCCTCCTCGGCTGGGGACTTTCATCCGCATGCCGAAGGCGGAAGTGTCTTACCCGACCGTTAATGTGCTCCGACATGCACCTGTGTGGCCGACGACTCTTTGGGTAAGATCGGCGGTCTCGAGGATCTGACG
>sample1_154 random sequence
TGCAGTTCAAAAGAGAAATTTACCATCTTTTACATGGAAATCCCCTCTGTACTGGTAATCATACAGCGCAATCTCAGTTACCAATATGGTCTCGGTGCCATAGTATTAACATGTAAGCAAATCGGGCTGATCTCCCTTTAGAAACCTAGTCGGCATGGCTCTTATGAGTCGGGGGGAGTTCCGTGAATTTCTAACTACACAGGCTCGCGCCTGAAGAGCAGTCCCCCGACACAACCCGACCCTCATATAA
>sample1_155 random sequence
CAGGACTAAAAGCGTCAACTGTAAGGTACTCAGATACCAACAGCCACCTAGGCTATGGTTATCTTAGTCGGTTAGAACACACGCGTGGGTACTACGGTGGGGATATACGACATTCAATAGCCCGCCCTTACAGATGGATTCCCGCACTGCGACATCAACCTTACATGACCGGGCTGTACTACAGCCTCTTCATTCCTTATCCAAAAATGTGTGGCAGCATAATGCGGTATCGTCGGCCGCGGCACAGCGA
>sample1_156 random sequence
GTACGCATACAACTAAGGGCTCTCGACGACACACCTGTTTCACTCCGGGCCGACAAGTCCCCATAAACGCCTTATCCAGGTTCGTGACTCAATGCTGGCTGTAAAGCATCTGGCTTTTCTGCACTCAATCTGGGTAGTGTCAGGCATGTAGAATCAATATATTAGCATCCGCACCTCTGACCCTGAGTCAACCTCCGCTTGATCGTCGAAGAAAGATGACTTCTGTACAGTGACGAGGACTATTCAAACA
>sample1_157 random sequence
TGCTTTTTCTATCCGAGGAGCTTGCGATCCATGGGCATACAATCCATGTGAGTGCCCCAAGTCCTATTTCGTGGTACTCGGAGATGGGCTAACGCAGCGGGGTGTTGCGTCGTCCTAGGGGTTAGCATGACTAATTGCTTCAATGACGTTTGCCATCGTGATTCACACACCGCTGATCTTGGCAGCTTCAGGCCCAAATGAAAAACAGTCTCTCTATCCGCGCTCCGGACTAATATTGCAACTCCGATGC
>sample1_158 random sequence
ACCACACTGGCTGCCCGGAGGCTTAATCAAGATTCAACCGGCTGACGGCTTTACCAAGACAATCAGGGTGGCGCACGCTCCAGCCGGCCTATCAACTTATCTGAATATTGTCTTCTTATACACCCAAACTCTTTTGTCTTCCTGTTTACGTTGGATGGATTGGCGTTTCGTTGGAAATTCCAGCATCCGGGATGACAACGCACCTTCACATTCCTGATATTCACTGCAAACAAGCGTGCGATAGCGCCTT
>sample1_159 random sequence
GCTCGACAAGGAGTTCTATTTGAGATCCATGTACACCCAGGAATCACCCCTGTCGAGCAAATCGGTGTGCAGCACACTCTAGCCTCGGCTCGACAGTGACAGAAAGTCACCCGTCTTAATATCAATGTAGTCAAAGCCGCGGAGCGCCGCCAGCCGGCGGCCAAGGCCCCTTGGGTGACGCGGAACACCGAGGGGTTCATGCCGCCCCTTTCGATCTCATTTTTGATACGTACTGTAGGAATCAGAGGTA
>sample1_160 random sequence
AGAAGCAACACCCGACTAGGCGTTTCTTAACAGGGCGAACCGCGTCTTAAGATCAGAGTACCTTACCTTAACACCTCTGGGCCGGCAGACCCGCTTAGTCGAATCAGTATAGTTGAAAACATTTCTGTTCGGTCGCACCATGTAATGCTTGCGGCTGCCCATTAAGGAGACGTAAATTATGGAGATTCCCGTTAGAATTTTCTTCGATAGTCATATGCCCGGTTCACGTCAAAGCACCCGAGGACCTGTA
>sample1_161 random sequence
CCCGAGTCACCCTGTGGAGTGTGTAAGTTTAGCCTTTTTCTGCATCTTCCCCTCTCTTACGCGGTAATCTAGGGATTAAGTCAGTTCGTGAATAGTAGCCGCGCAAGTACCTTGGTCTTGATACAAGACCTGTATCACTGGAGTGCAAGCCGCCTAAGAAGAGTAAGGAAACTTCCCCAAAGTGGTACCCATAGAATCTCGCGACATACTACGACAAGCAAATGACCTACCTGCGGGTTATACGTTTACA
>sample1_162 random sequence
CGACCCCGGTTAGAACTTCAATCGTTCAACAGGGCTAATCATACCCTTCTAAACAACGGGCCGCAGACTTACTGTCCACCACCGTTGCAGTTTAGGGTAGGCACCGCCCCACGTCCCATCCCACACTTAAACCCATACGCATACTTTCCCTGCGATGGGACGGGTTACACCACCGAGCAGGGAGCTAATGCTTTCGATCTCCACGACCTCAAACCTGAAGGGGCCACTGACTGTGGCAGTGGGAGTATTC
>sample1_163 random sequence
GATCGTGAGCGGGAACTATCTCAGGTCAAGGCAAAAATCCTCATACTTCCGCATCTATTAAGTCGTCGTTCTGGATCACGGCGGCCTGCTCGCGCTCTACCCTCGCTCATGATTAACCAGAGATTACTTTATCTTAAGGAGTACCAGGTCTAGGTTGAATGAAAAGGTCGTATGCTTAAGAATAAATTCGCTAAGTGGGCGGCCACACTACCAGCTGGTGTCCACCCACTGGTTGGCGGACTCTGGATGA
>sample1_164 random sequence
GCGTCGCTACTCCGTATGAAAACGGAGTTATTGCGGTTTGTGACTAGACTCCTTTCAGACTAGCATGAAATTCGAGTCTAGTCGGTATCAATCATCACCCATTGATACGGAGCAGTCCACTATTTCACTCTATCTCAAGTATCCAACATGCCGAGTAATACCTAGGAGGTACATATGGCGCTGGGCATTTTCTGTTGGCAGGAAAGACACGGACTGCGACCTGTAACTGAATAAACGTGTAAACTGTAGA
>sample1_165 random sequence
CTGACGAAAAGATAATGTATCGCGTTAGGAGTTTTTCGTATGTGCGCGCTAGCTGCGAGTCACGGGTGCCTAGATGACGAAAGCATCTCAGAAATTGAGCGTGAGCGGTACAGCGTAATCGCGGTCAAACGTCTACCACCACATCACATTGACATATACGTGAAGAGAGTTATGTTGATAATGATAAGAAATCCCCGGCAGTATCTCGCAAGGTTACAAGAAGTGACTAATACCACTAGGATGAGCCCAG
>sample1_166 random sequence
AGTGCCTTTATGTAATGAGCCTGAATGTCAGGCTTTGCCCACTACAGGGCATTCACTCTACGCAACGACGCGTGCCTACAAGCAGGCCGAGATTAGCGTATCTTAGTTGACCCCGTCTTCCCTTTACGGTATCGATCGCTTCGGCTTGTTAACTACAGGCCGCTTGATCGAGCTGCAGTGAACGATACGACATCAGTGTGCCTTTTCTCCTGTGTGTTGATTTCGGTCCTGGCTAACGTTTAATCTACAT
>sample1_167 random sequence
CTATGTATTTACTAAATGCTAGGTATATGACAGAACGAGGAGGTCGAGGCTCGACTCGAGCTGCGATGGCTAGCCATTCTCCAAACCGCGATCGTGAGTACTAAAAGCAAGAAGAGATCGGTCTTAATCGCTCGTGAAACACAATGACTCACAGTCTGTTGCGCCAAGGAATCTTGGATGCAATAGGCTCCGATGCTGTTGAGAATAACAGGCGTTACTAACTACAGCCTATGTGTCCCGGCCCCACTTA
>sample1_168 random sequence
CATACCGGCCACTAGTGAGTCTCGCAACCCACAAACCGGCGCCGCACGTCTCGGGGCATGTAACGAGACTCTGTCGTTTACTCAACATTTTTGCCATCCGCAGGTCTGACCATGTCGGTGCGAAAAATAGCTCAGTTGAATAAAGACCAGATTTTTAGTTCTCTGCCCAACCTCCGGGACCGGGCTCGTTATGAGTGAGGATGGATGTGACCGTACTGGAACTTAATTCCTGCTGTTAGATACGTCACGT
>sample1_169 random sequence
GAGAAATTGACGGGTCAGTACGTTATGACAGTCGTAGCTCGAGACTTGACACGACATGTCTTTACATACCTAATAGGACTTTCGAACCACAGAATCAGTCCTAGTGTAAACCGGGTACCACTAATGGCTATGGCCGCCGTCCTAGCGACAAACCCAGTCGGCGGTAGCGCGTGGTTACTGAACCAGCGCTTTGATGGCTAGTCAACCCATGTCCCAGTATAATTCTTTTTGTCTGCCTGGGGGAATACAC
>sample1_170 random sequence
TAACCGGAGATAGGACTCACCGGTCCTAACTGACCAGGGGCATAACATCAGGTTTAACACACTCATAGGCTTACCCGTGGTGCGTCCGGCAGCAATAGATTCCAGCAAATCCTGTACTACCTGAAAAACTTCATCCGAGGGCCACTTGACACGCCCAACCCCTACGAACTCGCTGGTCTTAACCGCATGCCCACAAGATAATAACATTGATTGACTCATTGCGCTAAAGACCTCGTTTACGCTAAATTCG
>sample1_171 random sequence
TGGCAAGAGTTACTCACTTCCATTGAGTTCGCTTACTCAGCGGCTAAAATGGATTTCTGACTTCGCCGTTGGGTGCCGTTCAATCATCAATTCTAAGCTTGCAAGCGTCATTCCTTCTAAAGACGGTCTCCTGACCAGATTTATACGAGGGCGAATCTCGTCATCCCGATAACTATCGGGCCTTGAATGAAGATGAGAACGTACGGTCGTAATGCATAATTATCGATTCAAAGAGGACCTGCTGGGGGTG
>sample1_172 random sequence
CCGCTCCACGTTGATACGTCCAAATGGGGGCAATGAAGCGTTCCCTGAGTCCACAACACGAACTACAGCTAGAATAAACTTACACTTGATTCGCTTACGAAGCCAAATTCGACTTTACACGCTTCCGCTCGGAATGTTTGCGCATAAACGTCTATAGACAGCGTGTATTAGGTCACTTATAGTATGCCCGAAAATAACTACGTATCTTCCGTTCTTCAACTTTGGGGGACACTCCTGAGGACGTTATGAC
>sample1_173 random sequence
TACTTCCCTAGAGGTCTGAAACAGGGGCGCAGATTTGTTTCGGGCAACCTGCACCCGGGCCATAATCTTGGGGGACGGTGCATATTTGGAGTGGATGGAGGGGTCGGGGGCAATGCAGTTAAATTGCTGTATCGATAAGACACCCTAGCCTGCGGGTTGGGGTGGAAGGCTAGAAAACCGCTATCTGAAGGGGGATGGTTCTAGCTCTGGGGTCCCAGACTGCATCTGGCAGCGGATCTTTAGATTACGA
>sample1_174 random sequence
ATGTTAATACCTACATACCTTTTGCCCTCCAAGAACTTCCTCACCTCGAAAGGATCCGTTTCGTCTATAGGAGCAGTGTTCGGTCCGAAAAATGAATTTGTGAGTGACGAGCGAGCCACGACGGGCCACCCCCGCCCATAGATCAGGCAGCTATGTTCCAGAACTCGTTACGAACGCTTACCTCGAAGCTTTACGACCCTGAACCTACCGTAGTGGTTGCTAGATCGCACGGATTATTGACCATTTCCGA
>sample1_175 random sequence
GTAGTTGCCCAACCCTCATGGAAGCTGAGTCACTCCTCTGCATCTGTCGGCCCGTCTGCGATGATGGCCTATCTACGCTTGAGCCTCGTGGCAGGCACTATTGGAGGTGTGCTCGATTGCACGTAGGAGTTGACTTCATTTTGTCCATGTCGATATAGGAAATGCACCAGTTGATAAGGTGTGTTTCCATTCTACCCGACTATATAAAGTGCGCGTCTGAGGCCCGTAGATCTCACACCCACTTTTCAAT
>sample1_176 random sequence
CGTTGGTCCCCAAACTGTAATAGACTCCGCGCTGTGTGGCATCATCCGGCAGTAGGTTATTGAGAGCTTGAATGACTACAGTGTAGGAGATCTGTAGCCCATCGTGCTCGACTTTGCTCACTCCCCTATTTACCCGACATGGAAACGACCCCAGCTTTATGACAAGCGCGCTGGGCACGACACGGCCATTAAATACTCTCGGTGTAAGAAATTGTCTAATATCCGATGTGTGCGAGAACTAAATATAGTC
>sample1_177 random sequence
CCAAATTGTTTCAGCTGGCGATCATGACGAGGATACTAGAAGGCCGCTATGGCCATTTCCTAGCTCGGCATCCGGGTCCCACATGGAGTAGTCGATCACTCTAACTAACATTAGCCACCGTTTGCTCAGCCGTTCCCTGAGGACTCCGAGCAATTTCCGCTGAAGTTGGTGAGGAGCTTCATACTACGGCGAGCTTAATGATTTATGTTACTCCTGGCATCTCGTCCTTTGCTGACTTAGGTCCGTCCGT
>sample1_178 random sequence
CGAAACCCGCATGGTATTGGTATCTGGCCTATACAATATAGGCTCGATGTGAGGCTTTGGGTTACTGCACAGCCTCTGTTCCTCGTACGGCAACACAGGCTCTTCATACATCCATGTGAGTTGATTCATTCTCCTCCTGTGTTTAACTCATCACTAGAATCTACGGTCACCCATCCCGGCTCTTCGGTTTCTGCGTGTCCTCCATTCTAGCTGGTTGCAGCGTCGCGCAATGACCCAACCCAACTTCATC
>sample1_179 random sequence
CAACATCACGGCCCTATCCAGGACTAAGGATGCTAATGGCTAGTCCCTTCACGACAGAGGTAGTATGAGATTAGGAACTTGCTCAGCCAGGCGAGCGACCGGGATGGGCTCATCCTGCAGCCTGGAACCTTGACCCTCCTAGTTGATCACATAAATTGGGTCTAGAGTACGAGGCAACACCACTGGCTGGCCCGGCGCGGTTCACTCGGGCTCTATACCGGATATTGATATATTAGTGAGTTCCCGACTG
>sample1_180 random sequence
TCTGCTTGATTCACCCAGTGGGACAGATAGCCTAAGTGCCTCTGTTGACTTGCGGGGAAGTCCAATGGGTCAAGTTCGAAGTGGCCCGCATTACTACTCCTCTGTTTCTCCGCGCGGACCGTTAGATTCTCGGTTGAACTTGCCACGCGGTTCATCGAAACGCAACATATGATTCCCCGTCATGCGGAGAGGCCATTGCCAGGCCATGTCGTAACCCCTTTCAGTAGAACTTAAAGTTAACAACACGGGG
>sample1_181 random sequence
TGTATGAATTGGCACGACCGTGGTCAGAGTTGGTACGGCAAAAGAATGGGCGTAAGAAAGCGGGATCATAGTGTCCTTCAGGCTGCCAGGACTAAGCCACCCTCGGGTACCAAACACCAATAGTTAGAATGTCGCAGTTAGTTGGCGCTTAATCTAGCTGGCCCGACTTCGAATATCCTCGGGATTGGGCGAAAGCTTTCCGGGCTGGAGCTCTGCACCATACATGAATACGAGCCGTCACGTTATACTT
>sample1_182 random sequence
AGCCGTGGAGTATCTGCTACCTGATGAACTTTAATACTGACTTCGATGCCGCAGGTGGCTTAGGGGTACTTGACTATCTGCGATCTAAAGTCCGGAAGCACATGTTGATAGGAGCGCAGCGTCTTTACGGGATCGAAACCTCCTCCACCCCAGAGGATAGTTAGGGGCAGGTGGAATCTCCGTCTAGTCGTAGACCGCACTTAATTTCCTTGTGATCTGGCACACTTTAAATGTGTCTGGGTCATAACCG
>sample1_183 random sequence
ATAGCCTACTGAGCTCTTCCTGAGGATAGGACTGTGCTACGTTTATGATCTAGGCTCGCCGGAAATTAGTGTGGAAGTTCTAAGGGCCCCGGCAATGCGCCTACCTGGACCCGGCGTCGCACAAGGAGACTACCACCAACGCGTTGAGAGAGAGGCTTCGACACCTACGTTTAGGACGTCTCAGCCAGTCTATGATTTCTGGGGGCCTAGGCATCTGGTGTAGGGACTCTTGCACCTTAAACGTGACCTT
>sample1_184 random sequence
GCTTTCGCAATGAACTCCGTGAAAGCCTCGAAGGCGATTGTTTGTCTGTACTCTAGAGGCTGTACACCTTTCGTGTAGGATCGGACCTCCCACCCGGTGGCTGCCGTATTCGTCTCCCCGTGTATCGGATGCCTAGATGTAATTAATATCCTCTGAACCCAGTTGTCGCAAAAAATGCTAGTACGCTGTCCTCGGCAGGTCTACAGCTATCGAACGTTCTCATTGACTTCTGAATTTTGTTTAGCTCATA
>sample1_185 random sequence
CATTAGACTTTCCCCGTGTGTGGTTGTATGAGCTGCCGGGAATCGAGACATCCGACCCCTGCGGATTTGGCCATATTGTGCGAGGAATGTACACCATAGGGTCCAACACCCCACCGTAGCCTTAGTCCAGGTTACCAACGATCTCGACCTAGTAGCGATGATCCTCCACGCTTGTTAGCGGACAAGTTCTATGCTAATGAGAGTTGCGGTCCAATACACCTTCCCAGCCCGTCCGCTATATGCGGTCGAA
>sample1_186 random sequence
CCAGTGAAGGAAGCAATTCTTCTGGGATTAGTCTCGATCTGCTTGCCACAATTGTAGCCGCTCGCAGTCCAATGAACGTGCGTTTTAAAATCGCTTAGGTGGCAATGCAGGCGCGGGGCGATTGAACTACCGACACGTGTTATTCCCGGTATCTCCTGGCCTTGACATCATGTGGTCATGCTGGTTGGTACCAGAAGGCTGCATGGGGGCGGAGCAGAAGCCTTCTAAGGCGCCGACAGTAATATCCCAT
>sample1_187 random sequence
TCGCGTTTACTTTTGAACAGGGGGTTAGGGATCGTACTTGTTTTTTTAGTGCCCTCGCCTTGGCCCCTCGTCTTAGGCCCAAACAACTATGCAGCCTCCGCACCGTACGTAGATCATCTAGGTCCACGACGGAAGGGGGAGATACCTCGTAAATCATGGCAACTCACAAGTATACTGGGCGCTTCACAAGCGGGATGAAAGCAATGCGTACTTCGGTTCCGGAGCGCAGTTATCGGGTGCGTCCCGAAGT
>sample1_188 random sequence
GCCTAGCGTATTTTGACATTCTACCAATCGATCATACAGACGATGCGAGAAGTTAAATTCTAGTTACAAACAGAACTGGAAAGGGAGCGGTGACGAGTCAGTGGCACGTAGAAAGCTCGGTATAGGCGAAGACTTGGGAACGAACCCCTAGCGGGGAAAGGCACTCACAATAGCCCATCAAGAACGAGCCCGTGTCATTGAATAAAAATGTTCACACGAGTAGGAATAGGAATCCTGTTGAATGCAAACC
>sample1_189 random sequence
TCTGGTGACCACTGTAAATCTCTCGCTCGAGAGCCTTGAAAAGTAATAGATACAGTTCTTCGAACCTCCGGACCCCGAGGCGTGGAGCGTACTCACAATCCTCGATGGTCTGTTAGGACTATCGATCACGACATATCCCATTTACGCGATGGCTCCTCTATGTTGTCCTTTCGCCTCGCGTGAGCCCTATAGAATGATCTGGTGTCAAGGCGAGCACAGATCTCGCAGGCTCTACGACGAGCACACTGTT
>sample1_190 random sequence
GTTAAAATGTACTTAGGGTGCGGCAACTGGTTCCGGCGGTACGAAATAACGCACACAGCCGGCGTCGACATGCATAGCGCAGGCCCCCCAGCATGCGATTTCACGGAGCGTTTACATAAACGCCGTGTTACGCTGAATGTGGTGGCAGGTTAGCCTACATGGACAGAATCAACATTAATTCGGCATAGTGTGATTGGTACTGGCGCGATTACACATAGTTCCCTTCCTCAAAAGTCGAACGAGACGGCGG
>sample1_191 random sequence
ATCACGTGGACCGCACGGAGTGGAACAATTCGCTCCCGTTTTTTTCGAGCAACTAGACGAAGCATAGTCGAGTCGGACTGATCAGCCATCTTCAAAATCATGTAGAGTAAATTGCTCCCTATCGATTGGCCTCGAGCATTGAGAATCTCTCGAAGAGATTACAAATCACCCCACCGGATAAGGAAGCCAGTCGATTATTCAAGATCGTTGCAGGGACTTTCGATCTAGGTTAAGCATATAGCCTCCTTGC
>sample1_192 random sequence
ATGGTTCTAGCAGTTGTCATTCTTACGGTGTCGGATCATGCACGTACGCATGGCGCGGTCAGAAACAGGTGCGGACATTTCCCCGAGCGACAGGATAGAGTGCGAGGCGAACCTGCAAAAATTACAACGTGCACAGCCCTTGCGTACCGTCCATGTTTCGATCCGCTTTTTGTGCTTCAAGTCTGGAGCTAGTCCAGAACGGATCGCGTCCGCCAGACTGCTCACTGGAGTAGGGCTAATTACCTACCGC
>sample1_193 random sequence
GTCACCGCTAATGGTTGGTCGTGGCTGGAAAGCGACCCGCAAGCGCGGTAATTAATTAGCCAGCGGGCGTTAACTGGCTAAGTCTTGAAGGCCGCACTTGTCTGTGAGTGCATAGGCTGTGTATCGATCAAGATAAATTACACAGTCATTGCGGGTGAATACTCGGATAAGATCTCCGACTGGAGGATCGATCTTGGCCCGCGCGTGAGCCATTGCGTGCTAGCGTTCCTTTAAGCGTATCTCGTCCAAT
>sample1_194 random sequence
CAGAGGTTCCAGAGTTAGCCCTCCGAGTTCACTTAGGACAAGGTAAAGCTTGTTTCCTAAAGTATAGTGGAGCAACGCAGCCGGATGGGTCACGTGGCCTACGTCAAGACAGGTCCGACCTTTGCAGGTACTGTAGACTGCGAGGCCCGAAGTTGCTTTCCCGTTTCCGCTATGAGATTCTAGGACGAGAATAGAGAGTGAACCCGGGTCCTGGGTTGGGCTGTCAATCCACACATCTTTCACAACCTGT
>sample1_195 random sequence
TACGGCTTACTCCGAGGCTCCAACGTATGGCGTGGTTACATAGGCCATTGGTTGACGCCTACCAGGACAGCTCCGCGGTCCAGGAAGGGACCCTAGACACTGCCCAGGGGATCAATCCCGGATGTTATGCCGATATGCCCAGAATACCGCGGCGGGAAGAGGCAAGTGAAGTTTTGGCCGCTGGGGTAACGTCCAGGCTACGTTCAGCGGTCTCGGCTTTATACAGCAACTCTTCGTTCGCAACGTTCCA
>sample1_196 random sequence
ATGAGGGACCATCGCACGCCAGAGGTACTGTCCGTTGTGCGAATGGTAAGGAGTAACCTTCAATTCAACCCAATTACCCCAATGAAGACGCTAAAGCATTAAAGCAGTAACCTATGTCGAGTATGAAGTCTGTATATGCTCGTTCCCCCCTGACTCAGCAGGCGGTAAACACACGTGGATGCAATGCAGGGACATCTGCAAAACTTCATTGTCTTCCCGCGAGCCGTAATTACCGATTAAGGACTTAACA
>sample1_197 random sequence
TGCATAGTCTATCCCGAATACCACATTTAGTGGCCCGTGGTGTAAACCCGTTCTTGTTGTTTAGTCTTGCCAGGGGGCGATGTAGACGGGTCATCACCAGTATCCATCCTCCACAGGTCAGCCTCAGTTTGCCCGATGCAGACACTGTGAGCATCCGTGCGATAAACCCTTCGCTATGACGACAGTCTACGCCCAGGCTGCGACGACCCAAGTCAATGGGAATACGCTGTTGAATAGCTTATCGTTACAT
>sample1_198 random sequence
CGCCCTAATATGAGAGCGATCAGACCGAAACTACCTCCGTCAGAGCAGGCTAACAGCAGTTACCACAACCCTCAGCGAGTATCTTTTTGCCTGTTGGTAATGCCTGATATCTAGGTACTAGCAATACTGTTAAACTAAGGTCAAGGTAAGATTGCCGAGACTAATTAGTGCGACTGTTCGACGCGCAGCCTCTAGTATCGGAAGATAGCGTATCATGAGCGGGGGACTGGCATCTGCCCGTCCGGGCAGA
>sample1_199 random sequence
ACAAAGCTAGGGCGTTTGCAACTAAGTGACGGTGGCGACCTTCTTTATCAAATCTCATATTAAACAGCTACACATATACGCAGTGAGTCGAAGCCCGGACGCTTGGCGAGTGAAGAATGTAGCCACTCAGAAACCTTATTTGTGGTATCGCTAACTTTGCTTTCCCTTTGTTTTGGAAAGCATGTTACCGTGGGGCCGTCCGATTTTACATGAGACTCGGAGCATTCTAAGATCCATCCATGCCATCTCC
>sample1_200 repeat of previous random sequence
TACGGCTTACTCCGAGGCTCCAACGTATGGCGTGGTTACATAGGCCATTGGTTGACGCCTACCAGGACAGCTCCGCGGTCCAGGAAGGGACCCTAGACACTGCCCAGGGGATCAATCCCGGATGTTATGCCGATATGCCCAGAATACCGCGGCGGGAAGAGGCAAGTGAAGTTTTGGCCGCTGGGGTAACGTCCAGGCTACGTTCAGCGGTCTCGGCTTTATACAGCAACTCTTCGTTCGCAACGTTCCA
>sample1_201 repeat of previous random sequence
TACGGCTTACTCCGAGGCTCCAACGTATGGCGTGGTTACATAGGCCATTGGTTGACGCCTACCAGGACAGCTCCGCGGTCCAGGAAGGGACCCTAGACACTGCCCAGGGGATCAATCCCGGATGTTATGCCGATATGCCCAGAATACCGCGGCGGGAAGAGGCAAGTGAAGTTTTGGCCGCTGGGGTAACGTCCAGGCTACGTTCAGCGGTCTCGGCTTTATACAGCAACTCTTCGTTCGCAACGTTCCA
>sample1_202 repeat of previous random sequence
TACGGCTTACTCCGAGGCTCCAACGTATGGCGTGGTTACATAGGCCATTGGTTGACGCCTACCAGGACAGCTCCGCGGTCCAGGAAGGGACCCTAGACACTGCCCAGGGGATCAATCCCGGATGTTATGCCGATATGCCCAGAATACCGCGGCGGGAAGAGGCAAGTGAAGTTTTGGCCGCTGGGGTAACGTCCAGGCTACGTTCAGCGGTCTCGGCTTTATACAGCAACTCTTCGTTCGCAACGTTCCA
>sample1_203 repeat of previous random sequence
TACGGCTTACTCCGAGGCTCCAACGTATGGCGTGGTTACATAGGCCATTGGTTGACGCCTACCAGGACAGCTCCGCGGTCCAGGAAGGGACCCTAGACACTGCCCAGGGGATCAATCCCGGATGTTATGCCGATATGCCCAGAATACCGCGGCGGGAAGAGGCAAGTGAAGTTTTGGCCGCTGGGGTAACGTCCAGGCTACGTTCAGCGGTCTCGGCTTTATACAGCAACTCTTCGTTCGCAACGTTCCA
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


pick_ref_otus_refseqs2 = """>r0 PC.635_779
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTTGACTTGGTGGGCCGTTACCCCGCCAACTATCTAATGGAACGCATCCCCATCGATAACCGAAATTCTTTAATAGTGAAACCATGCGGAAATACTATACTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTTATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCATCAA
>r1 PC.636_263
CTGGGCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTCCGCCACTAGGGCG
>r10 PC.355_740
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATTGTCGCTTTGGTAGGCCGTTACCCCTGCCAACTGGCTAATCAGACGCGGGTCCATCCCATACCACCGGAGTTTTTCACACAGCAACATGCGTTGCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACT
>r100 PC.356_1196
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCTGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACAAGCTAATCAGACGCGGGTCCATCTCGCACCACCGGAGTTTTCAGGGCGGGGGCATGCGCCCCCCTCCCGTTATGCGGTGTTAGCACCTATTTCTGGGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCG
>r101 PC.634_99
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGAAAGGAACATGCGTCCCCTTCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTTCGCCACTCGGGCAC
>r102 PC.634_91
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACCTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAATGATTTAATTATTAAAAGATGCCTTCAAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>r103 PC.481_49
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACCGGGGCATGCGCCCCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATATGGCAGGTTACCCACGCGTTACTCACCCG
>r104 PC.607_151
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCCTGTCCGGCCGAAGCCTTTCCTGCCTCCGGGATGCCCCGGTGGCATGTACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>r105 PC.607_1176
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTTGCCCTCTCAGGCCGGCTACTGATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCAGACGCGGGCCCATCCCATACCGCGAAAGCTTTCCATGCAGGGGCATGCGCCCCCGCATGTGCATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCTCACGCGTTACTC
>r106 PC.636_850
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGTCCACGCGTTACTCACCCGTCCGCCGCTAAGCC
>r107 PC.355_301
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCCGCGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r108 PC.355_1283
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGGTCCATCTTACACCACTAATGTTTTTCACTCTGTCCCATGCGGGACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCGTCCG
>r109 PC.593_1236
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTTCGGCTACGCATCATTGCCTTGGTAAGCCTTTACCCCACCAACTAGCTAATGCGCCGCGGGCCCATCCAAAAGCGGTAGCATAGCCACCTTTTACATAGTTACCATGCGGTAACTATGGTTATGCGGTATTAGCACCTGTTTCCAAGTGTTATCCCCCTCTTTTGGGCAGGTTGCCCACGTGTTACT
>r11 PC.354_1171
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACCGGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTGTCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r110 PC.634_170
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAGCTTTCCCCACAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r111 PC.634_173
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTTTAACCACCGGAGTTTTTAACCCAAGAAGATGCCTTCTCGAATGTTATGGGGAATTAGTACCAATTTCTCAGTGTTATGCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCACCAAAG
>r112 PC.607_992
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACCGATCGTCGCTTTGGTGGGCTTCTACCCCGCCAACTGGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACTATACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAACTGTTTCCAGTTGGTATCCCCCTGCATACGGCAGGTTACC
>r113 PC.355_399
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACCGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCGGGCGCGGGCCCATCCCGTGCCGCCGGAGCTTTCCGCATACACCCATGCGGCTGTATGCGCTTATGCGGTATTAGCGCCTGTTTCCAGGCGGTATCCCCCGGCACGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>r114 PC.593_644
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGTCCAGCTAATGCGCCGCATGCCCATCCTTCCCCGATGAATCTTTCCTCCTCCGGAGATGCCTCCGGAGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGGGAAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>r115 PC.354_470
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGATCCATCCTGTACCACCGGAGTTTTTCACACCGGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCAATCATTTCTAATTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r116 PC.636_351
TTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTC
>r117 PC.481_1112
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCTCAGTTTTTCACACTGTCCCATGCGAGACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTGAGGCAGGTTATCCACG
>r118 PC.355_1046
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGCTCTCTCAAGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r119 PC.481_1111
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCTTTTACCTCACCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCGGAGTTTTTCGCACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTTTCCACGCGTTACT
>r12 PC.635_570
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGATATATCTGCCATGCGACAAACATATATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTCCT
>r120 PC.354_1008
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCTTTTATCTCACCAACTAGCTAATCAGACGCGGGTCGCTCTCCCAGCGGCAGCTGCTCTGCTCTGCCACCTTTCTTTCTTCCATCATGCGATGGTTGAACCTCATGCGGTATTAGCTGTGATTTCTCACGGTTATT
>r121 PC.355_546
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTGGACTTGGTGAGCCGTTACCTCACCAACTATCTAATGGAACGCATGCCCATCTATCAGCGATAATTCTTTAACAAATATTTCCATGTGGAACCCCTGTTTTATGCGGTATTAGTCCGACTTTCGCCGGGTTATTCCCTCTGATAGGTAGGTTGCATACGCGTTACTTCACCCGTCGCG
>r122 PC.355_545
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTTGACTTGGTAGGCCATTACCCCACCAACTATCTAATCAGACGCAAGCCCATCTATTAGCGGATTTCTCCTTTCCCACTAGAACCATGTGGTTCCTATGGCATATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTAATAGGCAGGTTGCTTACG
>r123 PC.355_544
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCATACCGATAAATCTTTTCCGTCCGTATCATGCGATACTAGCGGGTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTATGAGGTAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r124 PC.635_1034
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGCCGGCTGCTGATCGTCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTCCAGCGATAAATCTTTGACGTCGGAGTCATGCGGCTCCAACGCATCATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCT
>r125 PC.354_600
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACCGATCGTTGACTTGGTGGGCCGTTACCCCGCCAACTATCTAATCGGACGCGAGCCCATCTCCCAGCGGATTGCTCCTTTGATATATCTACCATGTGATAAATATATGTTATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGGAAGGCAGGTTGCTCACGTGTTACT
>r126 PC.635_1039
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCGCCCTCTCAGGCCGGCTGCTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACGAGCTAATCAGACGCGGGCCCATCCCGTACCACCGGAGTTTTTCACACCAGGGCATGCGCCCCTGTGCGCTTATGCGGTATTACCAGCCGTTTCCAGCTGCTATCCCCCTGTATGGGCCAGGTTATCCACGCGTTACTCACCCG
>r127 PC.354_692
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACTGTCCCATGCAGGACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTGTCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r128 PC.607_636
CTGGGCCGTGTCTCAGTCCCAGTGTGTCCGATCAGCCTCTCAGCTCGGATATGCATCGTTGCCTTGGTAGGCCATTACCCCACCAACTAGCTAATACAACGCAAGCCCATCTCAAAGTGAAGCAAATGCTCCTTTTAACGTATCTAGTTATCTAATTACGTATCATCTGGTATTAGCGTTCATTTCTAAACGTTATCCCAGTCTCTGAGGTAGGTTACCCACG
>r129 PC.354_696
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCTGTGTAAGGCAGGTTACCCACGCGTTACT
>r13 PC.635_574
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACCGATCGTCGCTTTGGTGGGCCTCTACCCCGCCAACTGGCTAATCGGACGCGGGTCCATCCCATACCACCGGAGTTTTTCACACTGTATCATGCGTTACTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGGTATCCCCCGGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>r130 PC.593_414
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACTCTCTCAAGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCTTACCAACTAGCTAATCAGACGCGGGTCCATCCCATAGCGATAAATCTTTGAAGAAAAAGACATGCGTCTCTTTCCTGCTATGCGGTATTAGCAGCCGTTTCCAGCTGTTGTCCCCCTCTATGGGGCAGGTTACCCACGCGTT
>r131 PC.354_987
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCGGCACCGGGATGCCCCGATGCCGCGTCACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTGGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCCGG
>r132 PC.354_984
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTAGGCCGTCACCCTTCCAACTAGCTAATCAGACGCGGGTCCATCCTTTTGCGATAGCATCTTCAGAGGCCATCTTTCCTTAAAGAATCATGCGATCCTCTATTATTATGCGGTATTAGCTCCGATTTCTCGAAGTTGTTCCCCACAAAAGGGCAGGTT
>r133 PC.356_1148
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACAAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGCTTCATGCAAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCG
>r134 PC.634_214
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGATGTCTTGGTGGGCCGTTACCCCGCCAACAAACTAATGGAACGCATCCCCATCGATTATCGAAATTCTTTAATAACAAGAAGATGCCTTCTCGTTATGCTATCCAGTATTAATCTTTCTTTCGAAAGGCTATCCCGGAATAATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG
>r135 PC.634_213
CTGGTCCGTGTCTCAGTACCAGTGTGGGGACCTTCCTCTCAGAACCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCGCCAACCAGCTAATGCGCCGCATGCCCATCCTTTACCGGATCGCTCCTTTGACATACCGATCATGCGGTCGGTATGTGTTATGAGGTATTAGACGGAATTTTCTCCCGATTATCCCTCTGTAAAGGGCAGGTCGCATACGTGTTACTCACCCG
>r136 PC.354_462
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACCGGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCGTCCGCCACT
>r137 PC.607_529
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCATTACCTCACCAACTAGCTAATCAGACGCGAGGCCATCTTTCGGCGATAAATCTTTGATGTCCAGTGCCATGCGACACAGACATATTATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCACTCCAAGGCAGG
>r138 PC.354_1097
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGTCCATCTCACACCACCGGAGTTTTTCACACCAGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTGAGGCAGGTTATCCACGCGTTACTCACCCG
>r139 PC.354_1091
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCGCACCACCGGAGTTTTTCACACCGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCAGCTGTTGTCCCCCAGTGCGGGGCAGGTTGCCCACGCGTTACTCACCCG
>r14 PC.607_541
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTTCACCCTCTCAGGTCGGCTATGTATCGTCGCCTAGGTGAGCCTTTACCTCACCTACTAGCTAATACAACGCAGGTCCATCTGGTAGTGATGCAGTTGCATCTTTCAAACTTCAATCATGCAATCAAAGATATTATGCGGTATTAGCTATCGTTTCCAATAGTTATCCCCCGCTACCAGGCAGG
>r140 PC.607_274
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACTATACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAACTGTTTCCAGTTGGTATCCCCCTGCATACGGCAGGTTATCCACGCGTTACTCACCCG
>r141 PC.481_1069
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACAAGACCATGCAGTCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTGAGGCAGGTTATCCACGCGTTACT
>r142 PC.607_278
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCAGAGTTTTTCACACCGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGTGTTACTCACCCG
>r143 PC.356_828
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGCGCCGCATGCCCAGCCGCCACCGGATTCCTCCTTTCGCCCGGTCCGGATGCCCGGTCCGGGCGGCATGGGGTATTAGGCCGGGTTTCCCCGGGTTATCCCCCTGTGGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAGGTG
>r144 PC.635_951
CTGGGCTGTATCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATGCGCCGCGAGCCCATCCTATACCGGATTTCTCCTTTTACCCCTCAGCCATGCGACTGCGTGGTCTTATACCGTATTAACTACAGTTTCCCGTAGGTATCCCGTTGTATAGGGCAGGTTGCTCACGTGTTACT
>r145 PC.635_320
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTTGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCGTATACCTCCGGAAATTTTCACACTCTGTCATGCGACAGTGTGCGCTTATGCGGTATTAGCAGTTGTTTCCAACTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCG
>r146 PC.607_50
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGACCTCTCAGTCCGGCTACCGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCGGACGCGAGCCCATCTTCGAGCGATAAAATCTTTGATACCAGCAGGATGTCCTCCCGGTATGTTATGCGGTATTAGCGACCGTTTCCAGCCGTTATTCCCCTCTCGAAGGCAGGTTGCTCACGTGTTACTCACCCGTCCG
>r147 PC.636_558
CTGGTCCGTGTCTCAGTACCAGCGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAACGG
>r148 PC.355_122
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>r149 PC.635_931
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTGCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGAACCATCGTATACCACCGGAGTTTTTCACACTTTATCATGCGATCTCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTT
>r15 PC.356_1131
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGGCCCATCCCATACCACCGGAGTTTTTCACACAATCCCATGCGGGATTGGGCGCTTATGCGGTTTTAGCACCTGTTTCCAGGTGTTATCCCCCGGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTAGGTT
>r150 PC.481_1018
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTACCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCCGCGGTCTTATGCGGTATTAGCAGTCATTTCTAGCTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r151 PC.636_1032
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGAAAACCTCTCAGTCCCCCTATGTATCGTAGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATACAACGCATGCCCATCTCTAACCATCGGAATTTTCAACCATCTAAGATGCCTTAAATGATATTATGGGGTATTAGACCAAATTTCTTCGGATTATCCCCCTGTTAAAGGTAGGTTGCATACGCGTTACGCACCCGTGCGCCGGTCGCCGG
>r152 PC.635_939
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCATACCACCGGAGTTTTTCACACTGTTCCATGCGAAACCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r153 PC.634_58
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCATCAACCACCGGAGTTTTCAACCCCATGCGATGCCGCATGAGATATTATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTTGATGGTAGGTTGCATACGCGTTACGCACCCG
>r154 PC.636_955
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCTCACCAACTATCTAATGGAACGCATGCCCATCTATCAGCGATAATTCTTTAACAAATATTTCCATGTGGAACCCCTGTTTTATGCGGTATTAGTCCGACTTTCGCCGGGTTAGTTCCCCTCGTACAGGCCAGGTTGGCCCACCGCGTTACTCACCCGTCCGCCACTCA
>r155 PC.607_1036
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGACGCGAGCCCATCTCATACCGGATTTCTCCTTTTACCCCTCTACCATGCGGTAGTGTGGTCTTATACCGTATTAACCGAAGTTTCCCTCGGGTATCCCGTTGTATGAGGCAGGTTGCTCACGCGTTACTCACCCGTTCG
>r156 PC.634_52
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACAAACTAATGGAACGCATCCCCATCGATTATCGAAATTCTTTAATAACAAGAAGATGCCTTCTCGTTATGCTATCCAGTATTAATCTTTCTTTCGAAAGGCTATCCCGGAATAATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCAT
>r157 PC.634_57
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGG
>r158 PC.354_838
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTATCATGCGATACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCG
>r159 PC.355_756
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCTCAGAACCATGCGGTTCTGGGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r16 PC.634_203
CTGGGCCGTATCTCAGTCCCAATGGTGGCCGATCGCCCTCCTCAGGCCGGCTACCCATCGCAGGCTAGGTGGGCCGTTGCCCCGCCTACTACCTAATGGGCCGCGACCCCATCCCGCACCGTCGGAGCTTTCCTCCGTGGCGCATGCGCGCCTCGGAGAGTATCCGGTATTAGCCGCCGTTTCCGGCGGTTGTCCCGGGGTGCGGGGCAGGTTGGTCACGTGTTACTCACCCGTTCGCC
>r160 PC.355_971
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGGCAGTCAGAGCCATGCGGCCCAACTGCATTATGCGGTATTAGCACATCTTTCGATGTGTTATTCCCCACTCCAAGGCAGG
>r161 PC.636_824
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTCGGTGGGCCGTTACCCCGCCGACTGCCTAATGCGCCGCATGCCCATCCTCCACCACCGGAGTTTTCCTCCCACGGAGATGCCTCCGCGGGATTTACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGGCAG
>r162 PC.355_13
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>r163 PC.356_787
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTAACTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTTCCTGCTGTGCCATGCGGCACCGCAGGCTTATGCGGTATTAGCAGCCGTTTCCAGCTGTTGTCCCCCTGTACAGGGCAGGTTGCCCACACGTTACTCACCCG
>r164 PC.636_417
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATGCGCCGCATGCCCATCCTTGTCCGGATAAATCCTTTGATCGAATTCTCATGCGAGAACCCGATGTCACGCGGTATTAGACCGGATTTCTCCGGCTTATCCCCCTGACAAGGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>r165 PC.356_918
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCCGCACCGGATCGCTCCTTTGACCGCTCCCCCATGCAGAGGAACGGTGTCATGCCGTATTAGTCCGGATTTCTCCGGGTTATCCGTCTGTGCGGGGCAGGTCGCATACGCGTTACT
>r166 PC.356_917
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTTGCCTTGGTGGGCTGTTATCTCACCAACTAGCTAATCAGATGCGGGCCCATCTTTTACCGAATTTCTCCTTTCCTTCTCAGAAGATGCCTCCTAAGAATATTATGCGGTATTAGTCACCGTTTCCAGTGATTATTCCCCAGTAAAAGGCAGGTTGCCCACACGTTACTTACCCGTCCG
>r167 PC.356_915
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTCTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCGTATACCACCGGAGTTTTTCACACCGGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r168 PC.635_508
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCTCTCTCAAGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGACCATCCTGTACCACCGGAGTTTTCCCTGCGCTGCCATGCGGCAGCGCAGGATTATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCACTCCAAGGCAGGTTACTCACG
>r169 PC.481_338
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTACCTAATCAGACGCGGGCCCATCTTACACCACCGGAGTTTTTACCACCAAACCATGCGGTTTTGTGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCG
>r17 PC.634_200
TTGGACCGTGTCTCAGTTCCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCATTGCCTTGGTGGGCCGTTACCCCCCAACTAGCTAATAGGACGCATGCCCATCTGATACCTCGAACGATTTAATTATTATAAGATGCCTTACAATAATATTATGGGGTGTTAATCCACGTTTCCATGGGCTATCCCCCTGTATCAGCCAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>r170 PC.593_241
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGGTCCATCTTATACCACCGGAGTTTTTTCACACTGTACCATGCAGTACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGGTTACCCACGCGTTACTCACCCG
>r171 PC.635_119
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGAAAGAAACATGCGTCCCCTTCTGTTATGCGGTATTAGCACCTGTTTCCAGGTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCCACTAACTCATACAT
>r172 PC.593_1268
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG
>r173 PC.593_1305
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGGCCCATCTCATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTGTTTCCAAGTGTTATCCCCCTCTTTTGGGCAGGTTGCCCACGTGTTACTCACCCGTTCG
>r174 PC.636_604
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATTGGACGCGAGCCCATCTCTCAGCGGATTTCTCCTTTGGTCCAGACAGGATGTCCCATCCGGACGTTATGCGGTATTAGCGATCGTTTCCAATCGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTTACTCACCCGTCCG
>r175 PC.634_133
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTATGGATCGTTGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCCAACGCGAGCCCATCCTTCGGCACCTCAGTTTTGGTATTTCTATCATGCGGTAAAAATACGTTATGCGGTATTACTGTCCGTTTCCAGACACTATCCCCCTCCGAAGGGCAGGTTGCTCACGCGTTACTCACCCGTCCGCC
>r176 PC.634_134
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCGTCGTCGGCTAGGTGAGCTGCTACCTCACCTACTACCTGATGGGCCGCGACCCCATCCCAGACCGCAAGAGCTTTCATCCCTCCGCCATGCGGTGGAAGGATAGTATCCGGTATTAGCTGCCGTTTCCGGCAGTTATCCCGATGTCTGGGGCAGGTTGGTCACGTGTT
>r177 PC.593_450
TTGGACCGTGTCTCAGTTCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATCAGACGCGGGTCCATCCTGTACCGCAAAAGCTTTGATACTTCTACCATGCGATAAAAGTATTTTATCTCGTATTAGCATACCTTTCGGTATGTTATCCGTGTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r178 PC.634_139
CTGGACCGTCTCTCAGTTCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACATGAGCCCATCTTTCAGCGGATTGCTCCTTTGATAACAGGATCATGCGATCCCGTTCATTTCATTGCGGTATTAGCACACCTTTCGATGTGTTATCCCCCTCTGAAAGG
>r179 PC.607_1213
CTGGTCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGGCTTGGTGAGCCGTTACCTCACCAACTACCTAATCAGACGCGAGGCCATCTTCCAGCGATAAATCTTTGGCGTCGGAGTCATGCGGCTCCAACGCATCATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCTCACGCG
>r18 PC.634_201
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCGTCGTCGGCTAGGTGGGCCACTGCCCCGCCTACTACCTGATGGGCCGCGACCCCATCCTGCACCGCCGGAGCTTTCATCCGCTCCCCATGCGGGGTGCGGATAGTATCCGGTATTAGCTGCCGTTTCCGGCAGTTATCCCGATGTGCAGGGCAGGTTGGTCACGTGTTACT
>r180 PC.356_965
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACAAGCTAATCAGACGCGGGTCCATCTCGCACCACCGGAGTTTTCAGGGCGGGGGCATGCGCCCCCTCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>r181 PC.593_362
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTCACCGGATTGCTCCTTTGACCGCTCCGGGATGCCCCGGAATGGTGTTACGCGGAATTAGTCGGAATTTCTTCCGGTTATTCCCCTGTGACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>r182 PC.481_876
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACCGATTGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCGGACGCGAGCCCATCTCTCAGCGGATTGCTCCTTTGATATATCTACCATGCGATACATATATTTCATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAGAGGCAGGTTGCTC
>r183 PC.481_872
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCGCCACCGGAATCAACCTTTGGCACCAACAGGATGTCCCGTCGATGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCCCCTGTACAGGGCAGGTTCCCCACGCGTT
>r184 PC.481_4
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>r185 PC.607_330
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATATGGCAGGTTACCCACGCGTTACTCACCCG
>r186 PC.355_1227
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCACCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTTCTAGCTATATCATGCGATACTACTAGCTTATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTGATAGGCAGG
>r187 PC.635_750
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTACCAACTAGCTAATCAGACGCGAGGCCATCTCTCAGCGATAAATCTTTGGTACCAGTACCATGCGATACCCGTACGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCACTGAAAGGCAGGTTCCTCACGCGTTACTCACCCG
>r188 PC.635_615
TTGGACCGTGTCTCAGTTCCAATGTGGCCGTTCACCCTCTCAGGTCGGCTACCGATCGTCGCCTTGGTGGGCCGTCACCCCGCCAACTAACTAATCGGACGCAGGCCCATCCCTCGCCGTTCCCTTTGGCCCCTCGAGCATGCGCTCCAAAGGCTTCATGGGGTCTTAATCCCGGTTTCCCGAGGCTGTCCCCCAGCGAGGGGCAGGTTGCCTACGCGTTACTCACCCGTCCGCCGCTTTCCGCCATCG
>r189 PC.356_863
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGCTGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGAACCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCAGTCTCTTGGGCAGGTTACCCACGTGTTACT
>r19 PC.481_775
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGAGCTCATCTTTGGACGATAAATCTTTGGCAGCAAAATCATGCGGTCTCGCTGCATCATGCGGTATTAGCAGTCGTTTCCGGCTGTTATCCCCCATCCAAAGGCAGATTG
>r190 PC.635_923
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACCGATCGTTGCTTTGGTGGGCCGTTGCCCCGCCAACTGGCTAATCGGACGCGGATCCATCGTATGCCGTTAAACCTTTTCACACTGTGTCATGCGACACCGTGCGCTTATGCGGTATTAGCAACTGTTTCCAGTTGGTATCCCCCTGCATACGGCAGGTTATCCACGCGTTACTCACCCG
>r191 PC.356_1113
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGCCCATCCATTACCACCGGAGTTTTCAACCCAAGAAGATGCCTCCCTGGATGTTATGGGGTATTAGTACCGATTTCTCAGTGTTATCCCCCTGTAATGGGTAGATTGCATACGCGTTACTCACCCGTGCGCCGGTCG
>r192 PC.593_84
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATGCGCCATAAGTCCATCCTCTACCAGTGCCTTGCAGCACTTTTAATACGGTCACCATGCAGTGTCCCTACCTATGCGGTCTTAGCTGCCGTTTCCAGCAGTTATCCCCCTGTAAAGGCCAGG
>r193 PC.356_138
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGATATCCAGAGCCATGCGACCCAGATATATTATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCATCCAAGGCAGGTT
>r194 PC.607_419
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCATCCCATACCGCGAAAGCTTTCCATGCAGGGGCATGCGCCCCCGCATGTGCATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTACCCACGCGTT
>r195 PC.356_431
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGCTGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGAACCATCATATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCGGTATACGGCAGGTTTTCCACGCGTTACTCACCCG
>r196 PC.356_435
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTTGACCATGCAGTCCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTT
>r197 PC.635_566
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTCCGGCTACCGATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTAGCTAATCGGCCGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGCTGCATCGGGATGCCCCATCGCAGCGTTATGCGGTATTAGCAGTCGTTTCCAACTGTTATCCCCCTCTGAAAGGCAGGTTGCTCACG
>r198 PC.481_1058
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCCACTGATCGTCGCCTTGGTGGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGGAAAGCCATGCGGCTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTGTCCCCCTGTACAGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTGGG
>r199 PC.607_81
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATCAGACGCGGGTCGCTCTATCAGCGATAGCCTCTCTCGAGTCCGCCACCTTTCCTTCTGCCATCATGCGATGACTGAACCTTATGCGGTATTAGCACTGATTTCTCATTGTTATTCCCCT
>r2 PC.636_264
CTGGACCGTGTCTCAGTTCCAGTGTGGCCGATCACCCTCTCAGGTCGGCTACGTATCGTCGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATACGGCGCGGGTCCATCTATAAGTGATAGCAAAACCATCTTTCACTTTAGAACCATGCGGTTCTAAATGTTATCCGGTATTAGCTCCGGTTTCCCGAAGTTATCCCAGTCTTATAGGTAGGTTACCCACGTGTTACTCACCCGTCCGCCGCTAAG
>r20 PC.607_530
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCCTCCTATACCACTATCGTTTTTCACACAGGGCCATGCGGCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCC
>r21 PC.355_361
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACTCTCTCAAGCCGGCTACTGATCGTTGCTTTGGTAGGCCATTACCCTGCCAACTGGCTAATCAGACGCGGGGCCATCGTATGCCGATAACTCTTTTCACACCATGCCATGCAGCATTGTGTGCTTATGCGGTATTAGCACCTATTTCTAACTGTTATCCCCCTGTGTAAGGCAGG
>r22 PC.607_267
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATGCGCCGCATGCCCATCCTCCGCCGGAATCGCTTCCTTTTAACACCCGTGGGATGCCCCACAGATGTACCACGCGGTATTAGTCCGGATTTCTCCGGGTTATCCCCCTGCGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCG
>r23 PC.635_818
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCATCCTCTCAGACCAGCTACCCATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATGGGACGCGGACTCATCGAAAAGCGGAAGCATATGCAGAGGCCCCCTTTCCTGCATAGACATATATGCAGCGTATGCGGTATTAGCAGCCGTTTCCAACTGTTATTCCCCACTCCCCGGTAGATTAT
>r24 PC.481_1072
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATAAATCTTTTCCGGCCGGATCATGCGATCCTTCCGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r25 PC.635_15
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTAGGCCGTTACCCCACCAACTATCTAATCAGACGCAAGCCCATCTATCAGCGGATTGCTCCTTTCCCATTTATATCATGTGATATTCATGGCATATGCGGTATTAGCAGTCATTTCTAACTGTTGTTCCCCTCTGATAGG
>r26 PC.593_1326
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCTTTACCCCGCCAACAAACTAATGCACCGCAGGTCCATCCGCGCCCCATCCCCTAAAGGATGTTTCACAGAAAGAAGATGCCTCCTTCCTGTACATCGGGATTTGTTCTCCGTTTCCAGAGCGTATTCCCGGTGCGCGGGCAGGTTCCCTACGTGTTACTCACCCGTTCGCC
>r27 PC.634_46
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCATATACCACCGGAGTTTTTCACACAGAAACATGCGTCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGATTACCCACGTGTTACTCACCCG
>r28 PC.634_44
CTGGGCCGTGTCTCAGTCCCAGTCTGGATGATCATCCTCTCAAACCATCTAACGATCGTCGACTTGGTGAGCCTTTACCTCACCAACTATCTAATCGTACGCAGGCCATTCCTAAAGCGCATAAATGCTTTAATCCGAAGATCATATGCGGTATTAGCCACTCTTTCGAGTAGTTATCCCTCACTTTAGGGTATGTTCCCACGCGTTACTCAGCCGTCCG
>r29 PC.634_41
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCAGACGCGGGGCCATCTCATACCACCGGAGTTTTTCACACCGTGTCATGCGACACTGTGCGCTTATGCGGCATTAGCAGTCATTTCTAACTGTTATTCCCCTGTATGAGGCAGGTTCCCCACGCGTTACT
>r3 PC.481_1280
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCGTATGCCACCGGAGTTTTTCACACAGGAGCATGCGCTCCCGTGCGCTTATGCGGTGTTAGCACCTATTTCTAAGTGTTATCCCCCTGCATACGGCAGGTTACCCACGCGTTACT
>r30 PC.635_131
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCTTCCACCACAAAAAGCTTTCAACCCAGAGAGATGCCTCTCCGAATTATATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTGAAAGGTAGGTTGCATACGCGTTACGCACCCGTCCGCCGGTCG
>r31 PC.481_71
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCATCGCATACCACCGGAGTTTTTCACACCAAGCCATGCGGCTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCGGTATGCGGCAGGTTGCCCACGCG
>r32 PC.607_1004
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTATGCATCATCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCTTCACCGAAATTCTTTCCCCCTCGGAAGATGCCTCCCAAGGAGTATATGCGGTATTAGGCGAAATTTCTTCCGGTTATCCCGCTGTAAAGGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGG
>r33 PC.356_1009
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGGATCATGCAATCCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTGTAAGGCAGGTTCTCCACGCGTT
>r34 PC.635_645
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCGTCGTCGGCTAGGTAGGCCACTGCCCCACCTACTACCTGATGGGCCGCGACCCCATCTCACGCCGCGAAAGCTTTCATCCGCTCCTCATGCGAGGTGCGGATAATATCCGGTATTAGCTGCCGTTTCCGGCAGTTATCCCGATGCGTGAGGCAGGTTGGTC
>r35 PC.607_1002
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTCAACTCGGCTACGTATCATCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATACGCCGCAGGTCCATCCATGTCCACACCCGAAGGTCCTTTAATATAGTCAAGATGCCTTGTCTATACCTATCCGGTTTTAGCAGACGTTTCCATCTGTTATCTCGGGTACATGGGCAGGTTCCCTACGTGTTACTCACCCG
>r36 PC.607_1003
CTGGGCCGTGTCTCAGTCCCAATCTGGCCGGTCGGTCTCTCAACCCGGCTACCCATAGAAGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAACAGGCCGCGGGCCCATCCCTCTCCGCCAGAGCTTTCCCGACGCTTCCATGCGAAAGCGTCGGAGTATCCGGTATTATCCCCGGTTTCCCGAGGCTATCCCGATGAGAGGGGCAGGTTGCCCACGTGTTACTCAGCCGTTCGCCACTATCC
>r37 PC.356_668
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCTCCAACTAGCTAATCAGACGCGGGTCCATCTCATACCGTCTCGGCTTTTCACCCCGAACCATGCGGTTCTGTGTGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r38 PC.636_396
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACCAGCTAATCAGACGCGGGCTCATCTTATACTACCGGAGTTTTTCACACAGAAACATGCGTCCCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAAGGCAGATTACCCACGTGTTACTCACCCGTCCGCCG
>r39 PC.607_302
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGTCTTGGTGGTCCGTTACACCGCCAACTAACTAATGCGACGCATGCCCATCCTTCACCGAAATTCTTTCCCCCTCGGAAGATGCCTCCCAAGGAGTATATGCGGTATTAGGCGAAATTTCTTCCGGTTATCCCGCTGTAAAGGGTAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAG
>r4 PC.354_610
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTCTATCCGTCTCGGCTTTTCACACCGTGCCATGCGGCACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACT
>r40 PC.607_752
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCACACCGCCCGAAGGCTTTTCACACCGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTGTGGGGCAGGTTGCCCACGCGTTACT
>r41 PC.607_631
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATCAGACGCGAGCCCATCTCAGAGCGATAAATCTTTGGCAGTCAGGGAGATGTCTCCCAACTGCTTCATGCGGTATTAGCAGATGTTTCCATCTGTTATTCCCCACTCCAAGGCAGGTT
>r42 PC.607_1189
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTATACCACCGGAGTTTTTCACACTGGGCCATGCAGCCCCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>r43 PC.354_518
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTAGGTGGGCCGTTACCCCGCCTACCAGCTAATCAGACGCGAGGCCATCTTCCAGCGATAAATCTTTGGTATCAGGGTCATGCGGCCCCAATACTTTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCTCACGCGTTACTCACCAG
>r44 PC.607_1186
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCTGCCAACTAGCTAATCAGACGCGAGGCCATCTTCCAGCGATAAATCTTTGACATTCAAACCATGCGATTCAAATGTGTCATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGGAAGGCAGGTTCCTCACGCG
>r45 PC.593_1053
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCAACCTCTCAGTCCGGCTACCGATCGTCGGCTTGGTGGGCCGTTACCCCACCAACTACCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATTATCCTACCATGCGATAAGATAATGTCATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGCGTTACT
>r46 PC.607_1133
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGGCAATCAGTACCATGCGATACAATTGCATTATGCGGTATTAGCACAGCTTTCGTTGTGTTATTCCCCACTCCAAGGCAGGTT
>r47 PC.354_334
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGTTCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTTACTCACCCGTCCG
>r48 PC.354_729
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGCCGGCTGCTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAAAGGATTACCATGCGGTATTCCCTATTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCGGTACAGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTAAG
>r49 PC.636_402
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGGCCACCCTCTCAGGCCGGCTACCCATCGTCGCCTTGGTAGGCCATTACCCTACCAACTAGCTAATGGGACGCGAGTCCATCTTTCAGCGTCAGGTCTCCCCAACTTTTCCTATATAAACATGCACTTATATAACTCATCCGGCATTAGCTCACCTTTCGGTAAGTTGTTCCAGTCTAAAAGGCAGGTCACTCACGTGTTACT
>r5 PC.593_1329
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTCAAGTAAAAGACATGCGTCTCCTACTGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r50 PC.481_481
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCACCCTCTCAGGTCGGCTACTGATCGCTGCCTTGGTGGGCCGTTACCCCGCCAACTAGCTAATCAGACGCGGATCCATCCCATACCACCGGAGTTTTTCACACTGCGCCATGCGGCGCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCGGTATACGGCAGGTTCTCCACGCGTTACT
>r51 PC.481_994
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACCGTATCATGCGATACTGTGCGCTTATGCGGTTTTAGCAGCCGTTTCCGGCTGTTATCCCCCGGTACGGGGCAGGTTCCCCACGCGTTACTCACCCGTCCGCC
>r52 PC.593_1290
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCCCACCAACAAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGAGCCATGCAGCTCCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACG
>r53 PC.354_1210
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACCGGACCATGCGGTCCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGAGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC
>r54 PC.593_1292
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCAGGTCCATCTTACACCGATAAAATCTTTTCCGTCCGGGCCATGCGGCCCTAGCGGGTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTACCCACGCGTTACTCACCCG
>r55 PC.593_400
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGTGCTTTTACCACTCCAACTAGCTAATCAGACGCGGGTCCATCTCATACCAGATTTCTCCTTTTCACACGAATACATGCGTCTTCGTGCGCTTATGCGGTTTTACCACCAGTTTCCCGATGCTATCCCCCTGTATGAGGCAGGTTG
>r56 PC.607_197
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTATACCACCGGAGTTTTTCACACCGGAGCATGCGCTCCTGTGCGCTTATGCGGTATTAACAGTCGTTTCCAACTGTTATCCCCCTGTATAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT
>r57 PC.636_389
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTTGCCTTGGTGGGCCGTTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTCAAGTAAAAGACATGCGTCTCCTACTGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTAAAGGCCAGGTTACTTATGTATTACTCACCCGTTCGCCACTCGGGC
>r58 PC.354_1037
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTAGGCCGTTACCCTACCAACTAGCTAATGCGCCATAAGTCCATCCTTCAGCTATCCCGTGGGATATTTAATATGTCGATCATGCGATCTTCATACCTATGCGGTCTTAGCTATCATTTCTAATAGTTATCCCCCTCTTAAGGCCAGGTTACTTA
>r59 PC.636_222
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG
>r6 PC.481_1214
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCCGTACCACCGGAGTTTTCAAGGGGTCCCCATGCAGGGTCCCCTGTTATGCGGTATTAGCACCTGTTTCCAAGTGTTATCCCCCGGTACGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACTAAAACAGTCCGGGG
>r60 PC.607_949
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCGGGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAG
>r61 PC.607_927
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACTAGCTAATCAGACGCGGGCCCATCCCACACCGCCGGAGCTTTCCGCACCGCCCCATGCGGGGCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTGTCCCCCTGTGTGGGGTAGGTTGCCCACGCGTTACTCACCCGTCCG
>r62 PC.593_1279
CTGGGCCGTATCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTAGGCCTTTACCCCACCAACTAGCTAATCAGACGCAAGCCCATCTGTAAGTGGATTGCTCCTTTCCTCACTCTAACATGTGTCAGTGTGAGCATATGCGGTATTAGCAATGATTTCTCACTGTTATTCCCCTCTTACAGG
>r63 PC.636_148
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGTTAGGTGGGCCGTTACCCCGCCTACTGCCTAATGCGCCGCATGCCCATCCTCCACCGGTAATCCTTTCCTCCCCCAAGGATGCCCCCAAGGGATATACGCGGGATTAGCCTCCCTTTCGGAAGGTTGTCCCCCTGTGGAGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCAGTCGCCGG
>r64 PC.607_395
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCAGGACATGCATCCCTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCAGCTGTTGTCCCCCAGTGTGGGGCAGGTTGCCCACGCGTTACTCACCCG
>r65 PC.636_635
CTGGTCCGTGTCTCAGTACCAGTGTGGCCGGCCACCCTCTCAGGTCGGCTACCGATCGTCGCCTTGGTGGGCCTCTCCCCCGCCAACCAGCTAATCAGACGCGGATCCATCCCATACCACCGGAGTTTTTCACACAGGAGCATGCGCTCCCGTGCGCTTATGCGGTATTAGCACCCGTTTCCAGGTGGTATCCCCCGGTATGGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r66 PC.356_1256
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTCTTACGCGGTATTAGACGGAATTTCTCCCGATTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTTACTC
>r67 PC.355_293
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGCCACCGGTAATCCCTTTGGCACCAACAGGATGTCCCATAGGTGCATTATGCCGTATTAGACGGAATTTCTCCCGATTATCCGGCTGTGGCAGGCAGGTTGCATACGTGTTACT
>r68 PC.607_891
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCGTATACCACCTCGGTTTTTCACACTGTCCCATGCGGAACCGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG
>r69 PC.607_1068
TTGGACCGTGTCCCAGTCCCAGTGTGGCCGATCACCCTCTCAGGTCGGCTACCCATCGTTGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGGGTCGCGAGCCCATCCTATACCGATAAATCTTTTACCCCTGTATCATGCGATACTGTGGTCTTATACCGTATTAGCACAAATTTCTCTGTGTTATCCCGTTGTATAGGGCAGGTTGCTCACGTGTTACTCACCCGTTCG
>r7 PC.593_1320
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAACCTCTCAGTTCGGCTACGCATCATTGCCTTGGTAAGCCCTTTACCCCACCAACTAGCTAATGCGCCGCGGGCCCATCCAAAAGCGGTAGCATAGCCACCTTTTACATAGTTACCATGCGGTAACTATGGTTATGCGGTATTAGCACCTGTTTCCAAGTGTTACTCCCCCTTCTTTTGGGCAGGGTTCGCCCCACGTTGTTACT
>r70 PC.593_526
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCG
>r71 PC.481_860
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGGCTTGGTGAGCCGTTACCCCACCAACTACCTAATCAGACGCGGGTCCATCTTACACCACCTCAGTTTTTACCACTGTACCATGCAGTACTGTGGTCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACGCGTTACTCACCCG
>r72 PC.607_573
CTGGTCCGTGTCTCAGTCCCAATGTGGCCGGTTGCCCTCTCAGGCCGGCTACTGATCGTCGGTTTGGTGGGCCGTTACCCCGCCAACTGCCTAATCGGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGATATCAGTGCCATGCGACTCCGATATTTTATGCGGTATTAGCGTTCGTTTCCAAACGTTATCCCCCTCTGAAAGGCAGGTTGCTCACGTGTTACTCACCCG
>r73 PC.607_673
TTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCAACGG
>r74 PC.635_72
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCACGGTGGGCCGTTACCCCGCCGTCAAGCTAATGTCACGCGAGCCTATCCTCATCCGACGAATCTTTAGATGGATTCAGATGCCTGATTCCATCACCATGGGGCGATTAGACGCCGTTTCCTAGCGATTATTCCCCTCGATGAGGGCAAGTTGCTCACG
>r75 PC.635_647
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCCCACCAACCAGCTAATCAGACGCGAGTCCATCTCAGAGCGATAAATCTTTGGCAGTCAGAGCCATGCGGCCCAACTGCATCATGCGGTATTAGCAGCTGTTTCCAGCTGTTATTCCCCACTCCAAGG
>r76 PC.355_1239
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r77 PC.355_1238
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCG
>r78 PC.635_742
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGGCTTGGTGGGCCGTTACCTCACCAACTACCTAATCAGACGCGGGTCCCTCCTATACCACTATCGTTTTTCACTCAGCGCCATGCGGCTCCGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTCCCCACGCGTTACTCACCCGTCCG
>r79 PC.354_1044
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCTTGCACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT
>r8 PC.354_665
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAAGCCGGCTACTGATCGTTGCTTTGGTAGGCCATTACCCTGCCAACTGGCTAATCAGACGCGGGGCCATCGTATGCCGATAACTCTTTTCACACCATGCCATGCAGCATTGTGTGCTTATGCGGTATTAGCAGTCATTTCTGACTGTTGTCCCCCTGCATACGG
>r80 PC.635_621
TTGGGCCGTGTCTCAGTCCCAATGTGGGCCGTCTACCCTCTCAGGCCGGCTACGTATCATCGCTTTGGTGGGCCGTTACCCGCCAACTGGCTAATACGCCGCAGGTCCATCCATGTCCACACCCGAAGGTCCTTTAATATAGTCAAGATGCCTTGTCTATACCTATCCGGTTTTAGCAGACGTTTCCATCTGTTATCCCGGGTACATGGGCAGGTTCCCTACGTGTTACTCACCCGTTCG
>r81 PC.356_1168
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATGTCACGCATGCCCATCCGTTGCCGGAATCGCTTCCTTTGGCCCCAGGGCCATGCAGCCCCGCGGCATTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCG
>r82 PC.356_1161
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACCGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCGGACGCGGACCCATCGCATGCCGGGATCGCTCCCTTTCCGCACTGCGCCATGCGGCGCCGTGCGCATATGCGGTATTAGCGGTTGTTTCCAGCCGGTATCCCCCTGCATGCGGCAGGTTG
>r83 PC.356_1166
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTATACCACCTCAGTTTTTCACACCGAACCATGCGGTTTTAGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTATAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCC
>r84 PC.636_6
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG
>r85 PC.356_249
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGCCGGCTGCTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACT
>r86 PC.481_1022
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACTCTCTCAAGCCGGCTACTGATCGTCGCTTTGGTAGGCCATTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCTATAAATATTTGATAGCAATATCATGCGGTATCGCTATGTTATGCGGTGTTAGCAGCCGTTTCCAGCTGTTATCCCCCTGTACAGGGCAGGTTACCCACG
>r87 PC.481_1025
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCATCGCATACCACCGGAGTTTTTCACACCAAGCCATGCGGCTCTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCCTGTGTAAGGCAGGTTGCCCACG
>r88 PC.481_738
CTGGGCCGTATCTCAGTCCCAATGTGGCCGATCGACCTCTCAGTCCGGCTACCCATCGCAGGCTCGGTGGGCCGTTACCCCGCCGACTACCTAATGGGCCGCGACCCCATCTCACACCGCGAGAGCTTTCATCCGGCTAACATGTGTTGGCCGGATAGTATCCGGTATTAGCTACCGTTTCCAGTAGTTATCCCGGTGTGTGAGGCAGGTTGGTCACGTGTTACTCACCCGTTCG
>r89 PC.635_968
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGGCTTGGTGGGCCGTTGCCCCGCCAACTACCTAATGTCGCGCGTGCCCGTCCCGTACCACCGGAATTTTAAATCGAGAGCCATGCGGCTCTCGAGTATCATGGGATGTTAGTCCACGTTTCCGCGGGTTATCTCCCGGTACGGGGTTGGTTGCACACGTGTTACTCACCCGTGCGCCGGTCGCCGGCG
>r9 PC.481_1193
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCTTACCAACCAGCTAATCAGACGCGGGCCCATCCTGTACCACCGTGGTTTTCCCTGCTGTTCCATGCGGCACAGCAGGCTTATGCGGTATTAGCAGCCATTTCTGGCTGTTGTCCCCCGGTACAGGGCAGGTTGCCCACGCGTTACT
>r90 PC.635_966
CTGGACCGTGTCTCAGTTCCAGTGTGGGGGACCTTCCTCTCAGAACCCCTAGACATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGCGTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGGCGG
>r91 PC.635_364
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTGGGCTTGGTGGGCCGTTACCCCGCCAACTACCTAATCAGACGCGGACCCATCGTGTACCGTACTAGATAAGATCTAGGCTTTCCGCCCTGTGCCATGCGGCACTGTGCGCATATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTACACGGCAGGTTG
>r92 PC.356_233
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGCCCATCCTGTACCACCGGAGTTTTCAGGGAAAAGCCATGCGGCTTCCCCCGTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCGTCCGCCACT
>r93 PC.635_367
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCAACCTCTCAGTCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTGCCCCGCCAACTAGCTAATCAGACGCGAGCTCATCTCAGAGCGATAAATCTTTGGCGTCCAGAGAGATGCCTCCCAGACGCATCATGCGGTATTAGCGGCTGTTTCCAACCGTTATTCCCCACTCCAAGG
>r94 PC.593_932
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACGCATCGTCGCCGTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCG
>r95 PC.356_238
TTGGACCGTATCTCAGTTCCAATGTGGCCGATCAGCCTCTCAGCTCGGCTATGCATCGTTGCCTTGGTAGGCCATTGCCCCACCAACTAGCTAATACACCGCAAGCTCATCCTAAGGTGAAGCAAACGCTTCTTTTAACATATCTAGTTATCTAGTTATGTATCATCCGGTATTAGCGTTCGTTTCCAAACGTTATCCCAGTCCCTAGGGTAGATTACCC
>r96 PC.481_283
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGAACCATGCGGTTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r97 PC.354_773
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCAGACGCGGGTCCATCTTATACCACCGGAGTTTTTCACACACTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>r98 PC.636_592
CTGGTCCGTGTCTCAGTACCAGCGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGTGGACTTGGTGAGCCGTTACCTCACCAACTATCTAATGGAACGCATGCCCATCTATCAGCGATAATTCTTTAACAAATATTTCCATGTGGAACCCCTGTTTTATGCGGTATTAGTCCGACTTTCGCCGGGTTATTCCCTCTGATAGGTAGGTTGCATACGCGTTACTCACCCGTGCGCCGG
>r99 PC.356_336
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGCCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTCTACCCCGCCAACAAGCTAATCAGACGCGGGTCCATCGTATACCACCGGAGTTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGCCAGGTTGCCCACGCGTTACTCACCCG
>r1000 PC.607_1000
GCCGTCAGTTCGTGCCGTAAGGTGTTCTGTTAAGTCAGATAACGAACGAGACCCATGCCATTAGTTGCTATCTGTTCCTTCGGGGACAGAGCACTCTAATGGGACCGCTGCTGCTAAAGCAGAGGAAGGTGTGGGCAACGGTAGGTCAGTATGCCCCGAATCTCCCGGGCTACACGCGGACTACAATGGTTGAAACAATGGGCTGCTACGCCGAGAGGCGACGCTAATCTCTAAATTCAATCCAAGTTCGGA
"""

if __name__ == "__main__":
    main()