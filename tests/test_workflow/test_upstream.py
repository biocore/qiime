#!/usr/bin/env python
# File created on 20 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from shutil import rmtree
from glob import glob
from os.path import join, exists, getsize, split, splitext
from cogent import LoadTree, LoadSeqs
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.util import get_tmp_filename
from qiime.util import (load_qiime_config,
                        get_qiime_temp_dir,
                        create_dir)
from qiime.parse import (parse_qiime_parameters)
from biom.parse import parse_biom_table
from qiime.test import (initiate_timeout,
                        disable_timeout,
                        get_test_data_fps)
from qiime.workflow.util import (call_commands_serially,
                                 no_status_updates)
from qiime.workflow.upstream import (run_pick_closed_reference_otus,
                                     run_pick_de_novo_otus)

class UpstreamWorkflowTests(TestCase):
    
    def setUp(self):
        """ """
        self.test_data = get_test_data_fps()
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
        
        self.qiime_config = load_qiime_config()
        self.params = parse_qiime_parameters([])
        
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

    def test_run_pick_closed_reference_otus(self):
        """run_pick_closed_reference_otus generates expected results"""
        run_pick_closed_reference_otus(
         self.test_data['seqs'][0],
         self.test_data['refseqs'][0],
         self.test_out,
         self.test_data['refseqs_tax'][0],
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)

        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_ref_picked_otus',
         '%s_otus.txt' % input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2']
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)

        # Number of OTUs matches manually confirmed result
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,3)

        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # otu ids are as expected
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        
        # expected number of sequences in OTU table
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,117)
        
        # One tax assignment per otu
        self.assertEqual(len(otu_table.ObservationMetadata),3)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
        
    def test_run_pick_closed_reference_otus_parallel(self):
        """run_pick_closed_reference_otus generates expected results in parallel
        """
        run_pick_closed_reference_otus(
         self.test_data['seqs'][0],
         self.test_data['refseqs'][0],
         self.test_out,
         self.test_data['refseqs_tax'][0],
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=True,
         status_update_callback=no_status_updates)

        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_ref_picked_otus',
         '%s_otus.txt' % input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2']
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)

        # Number of OTUs matches manually confirmed result
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,3)

        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # otu ids are as expected
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        
        # expected number of sequences in OTU table
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,117)
        
        # One tax assignment per otu
        self.assertEqual(len(otu_table.ObservationMetadata),3)

        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_pick_de_novo_otus_rdp_tax_assign(self):
        """run_pick_de_novo_otus generates expected results with rdp tax assignment
        """
        self.params['assign_taxonomy'] = \
         {'id_to_taxonomy_fp':self.test_data['refseqs_tax'][0],
          'reference_seqs_fp':self.test_data['refseqs'][0],
          'assignment_method':'rdp'}
        self.params['align_seqs'] = \
         {'template_fp':self.test_data['refseqs_aligned'][0]}
        self.params['filter_alignment'] = \
         {'lane_mask_fp':self.test_data['refseqs_aligned_lanemask'][0]}
        actual_tree_fp, actual_otu_table_fp = run_pick_de_novo_otus(
         self.test_data['seqs'][0], 
         self.test_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_picked_otus',
         '%s_otus.txt' % input_file_basename)
        alignment_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_aligned.fasta' % 
          input_file_basename)
        failures_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_failures.fasta' % 
          input_file_basename)
        taxonomy_assignments_fp = join(self.test_out,
         'rdp_assigned_taxonomy','%s_rep_set_tax_assignments.txt' %
         input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        tree_fp = join(self.test_out,'rep_set.tre')
        
        self.assertEqual(actual_tree_fp,tree_fp)
        self.assertEqual(actual_otu_table_fp,otu_table_fp)
        
        input_seqs = LoadSeqs(self.test_data['seqs'][0],
                              format='fasta',
                              aligned=False)
        
        # Number of OTUs falls within a range that was manually 
        # confirmed
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,14)
        
        # all otus get taxonomy assignments
        taxonomy_assignment_lines = list(open(taxonomy_assignments_fp))
        self.assertEqual(len(taxonomy_assignment_lines),num_otus)
        
        # number of seqs which aligned + num of seqs which failed to
        # align sum to the number of OTUs
        aln = LoadSeqs(alignment_fp)
        failures = LoadSeqs(failures_fp,aligned=False)
        self.assertTrue(aln.getNumSeqs() + failures.getNumSeqs(),num_otus)
         
        # number of tips in the tree equals the number of sequences that
        # aligned
        tree = LoadTree(tree_fp)
        self.assertEqual(len(tree.tips()),aln.getNumSeqs())
        
        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2','not16S.1']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # otu ids are as expected
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        # number of sequences in the full otu table equals the number of
        # input sequences
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,input_seqs.getNumSeqs())
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)

    def test_run_pick_de_novo_otus_uclust_tax_assign(self):
        """run_pick_de_novo_otus generates expected results with uclust tax assignment
        """
        self.params['assign_taxonomy'] = \
         {'id_to_taxonomy_fp':self.test_data['refseqs_tax'][0],
          'reference_seqs_fp':self.test_data['refseqs'][0],
          'assignment_method':'uclust'}
        self.params['align_seqs'] = \
         {'template_fp':self.test_data['refseqs_aligned'][0]}
        self.params['filter_alignment'] = \
         {'lane_mask_fp':self.test_data['refseqs_aligned_lanemask'][0]}
        actual_tree_fp, actual_otu_table_fp = run_pick_de_novo_otus(
         self.test_data['seqs'][0], 
         self.test_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_picked_otus',
         '%s_otus.txt' % input_file_basename)
        alignment_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_aligned.fasta' % 
          input_file_basename)
        failures_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_failures.fasta' % 
          input_file_basename)
        taxonomy_assignments_fp = join(self.test_out,
         'uclust_assigned_taxonomy','%s_rep_set_tax_assignments.txt' %
         input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        tree_fp = join(self.test_out,'rep_set.tre')
        
        self.assertEqual(actual_tree_fp,tree_fp)
        self.assertEqual(actual_otu_table_fp,otu_table_fp)
        
        input_seqs = LoadSeqs(self.test_data['seqs'][0],
                              format='fasta',
                              aligned=False)
        
        # Number of OTUs falls within a range that was manually 
        # confirmed
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,14)
        
        # all otus get taxonomy assignments
        taxonomy_assignment_lines = list(open(taxonomy_assignments_fp))
        self.assertEqual(len(taxonomy_assignment_lines),num_otus)
        
        # number of seqs which aligned + num of seqs which failed to
        # align sum to the number of OTUs
        aln = LoadSeqs(alignment_fp)
        failures = LoadSeqs(failures_fp,aligned=False)
        self.assertTrue(aln.getNumSeqs() + failures.getNumSeqs(),num_otus)
         
        # number of tips in the tree equals the number of sequences that
        # aligned
        tree = LoadTree(tree_fp)
        self.assertEqual(len(tree.tips()),aln.getNumSeqs())
        
        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2','not16S.1']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # otu ids are as expected
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        # number of sequences in the full otu table equals the number of
        # input sequences
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,input_seqs.getNumSeqs())
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)


    def test_run_pick_de_novo_otus_parallel(self):
        """run_pick_de_novo_otus generates expected results in parallel
        """
        self.params['assign_taxonomy'] = \
         {'id_to_taxonomy_fp':self.test_data['refseqs_tax'][0],
          'reference_seqs_fp':self.test_data['refseqs'][0]}
        self.params['align_seqs'] = \
         {'template_fp':self.test_data['refseqs_aligned'][0]}
        self.params['filter_alignment'] = \
         {'lane_mask_fp':self.test_data['refseqs_aligned_lanemask'][0]}
        actual_tree_fp, actual_otu_table_fp = run_pick_de_novo_otus(
         self.test_data['seqs'][0], 
         self.test_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=True,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_picked_otus',
         '%s_otus.txt' % input_file_basename)
        alignment_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_aligned.fasta' % 
          input_file_basename)
        failures_fp = join(self.test_out,
         'pynast_aligned_seqs','%s_rep_set_failures.fasta' % 
          input_file_basename)
        taxonomy_assignments_fp = join(self.test_out,
         'uclust_assigned_taxonomy','%s_rep_set_tax_assignments.txt' %
         input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        tree_fp = join(self.test_out,'rep_set.tre')
        
        self.assertEqual(actual_tree_fp,tree_fp)
        self.assertEqual(actual_otu_table_fp,otu_table_fp)
        
        input_seqs = LoadSeqs(self.test_data['seqs'][0],
                              format='fasta',
                              aligned=False)
        
        # Number of OTUs falls within a range that was manually 
        # confirmed
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,14)
        
        # all otus get taxonomy assignments
        taxonomy_assignment_lines = list(open(taxonomy_assignments_fp))
        self.assertEqual(len(taxonomy_assignment_lines),num_otus)
        
        # number of seqs which aligned + num of seqs which failed to
        # align sum to the number of OTUs
        aln = LoadSeqs(alignment_fp)
        failures = LoadSeqs(failures_fp,aligned=False)
        self.assertTrue(aln.getNumSeqs() + failures.getNumSeqs(),num_otus)
         
        # number of tips in the tree equals the number of sequences that
        # aligned
        tree = LoadTree(tree_fp)
        self.assertEqual(len(tree.tips()),aln.getNumSeqs())
        
        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2','not16S.1']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # otu ids are as expected
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        # number of sequences in the full otu table equals the number of
        # input sequences
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,input_seqs.getNumSeqs())
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
    def test_run_pick_de_novo_otus_muscle(self):
        """run_pick_de_novo_otus w muscle generates expected results
        """
        self.params['assign_taxonomy'] = \
         {'id_to_taxonomy_fp':self.test_data['refseqs_tax'][0],
          'reference_seqs_fp':self.test_data['refseqs'][0]}
        self.params['align_seqs'] = {'alignment_method':'muscle'}
        self.params['filter_alignment'] = \
         {'suppress_lane_mask_filter':None,
          'entropy_threshold':'0.10'}
        
        run_pick_de_novo_otus(
         self.test_data['seqs'][0], 
         self.test_out, 
         call_commands_serially,
         self.params, 
         self.qiime_config, 
         parallel=False,
         status_update_callback=no_status_updates)
         
        input_file_basename = splitext(split(self.test_data['seqs'][0])[1])[0]
        otu_map_fp = join(self.test_out,'uclust_picked_otus',
         '%s_otus.txt' % input_file_basename)
        alignment_fp = join(self.test_out,
         'muscle_aligned_seqs','%s_rep_set_aligned.fasta' % 
          input_file_basename)
        taxonomy_assignments_fp = join(self.test_out,
         'uclust_assigned_taxonomy','%s_rep_set_tax_assignments.txt' %
         input_file_basename)
        otu_table_fp = join(self.test_out,'otu_table.biom')
        tree_fp = join(self.test_out,'rep_set.tre')
        
        input_seqs = LoadSeqs(self.test_data['seqs'][0],
                              format='fasta',
                              aligned=False)
        
        # Number of OTUs falls within a range that was manually 
        # confirmed
        otu_map_lines = list(open(otu_map_fp))
        num_otus = len(otu_map_lines)
        otu_map_otu_ids = [o.split()[0] for o in otu_map_lines]
        self.assertEqual(num_otus,14)
        
        # all otus get taxonomy assignments
        taxonomy_assignment_lines = list(open(taxonomy_assignments_fp))
        self.assertEqual(len(taxonomy_assignment_lines),num_otus)
        
        # all OTUs align
        aln = LoadSeqs(alignment_fp)
        self.assertTrue(aln.getNumSeqs(),num_otus)
        
        # all OTUs in tree
        tree = LoadTree(tree_fp)
        self.assertEqual(len(tree.tips()),num_otus)
         
        # check that the two final output files have non-zero size
        self.assertTrue(getsize(tree_fp) > 0)
        self.assertTrue(getsize(otu_table_fp) > 0)
        
        # Check that the log file is created and has size > 0
        log_fp = glob(join(self.test_out,'log*.txt'))[0]
        self.assertTrue(getsize(log_fp) > 0)
        
        # parse the otu table
        otu_table = parse_biom_table(open(otu_table_fp,'U'))
        expected_sample_ids = ['f1','f2','f3','f4','p1','p2','t1','t2','not16S.1']
        # sample IDs are as expected
        self.assertEqualItems(otu_table.SampleIds,expected_sample_ids)
        # expected OTUs
        self.assertEqualItems(otu_table.ObservationIds,otu_map_otu_ids)
        # number of sequences in the full otu table equals the number of
        # input sequences
        number_seqs_in_otu_table = sum([v.sum() for v in otu_table.iterSampleData()])
        self.assertEqual(number_seqs_in_otu_table,input_seqs.getNumSeqs())

if __name__ == "__main__":
    main()
