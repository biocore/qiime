#!/usr/bin/env python
# File created on 02 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from glob import glob
from os import chdir, getcwd
from os.path import exists
from shutil import rmtree
from tempfile import mkdtemp
from unittest import TestCase, main

from skbio.tree import TreeNode
from skbio.util import remove_files
from biom import load_table

from qiime.util import load_qiime_config, count_seqs, get_qiime_temp_dir
from qiime.workflow.util import (call_commands_serially, no_status_updates,
                                 WorkflowError)
from qiime.parse import parse_qiime_parameters, fields_to_dict
from qiime.test import initiate_timeout, disable_timeout, get_test_data_fps
from qiime.workflow.pick_open_reference_otus import (
    pick_subsampled_open_reference_otus,
    iterative_pick_subsampled_open_reference_otus,
    final_repset_from_iteration_repsets,
    convergent_pick_subsampled_open_reference_otus)
from bfillings.sortmerna_v2 import build_database_sortmerna


allowed_seconds_per_test = 120


class PickSubsampledReferenceOtusThroughOtuTableTests(TestCase):

    """ """
    def setUp(self):
        """
        """
        self.start_dir = getcwd()
        self.qiime_config = load_qiime_config()
        self.dirs_to_remove = []
        self.files_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()
        self.test_data = get_test_data_fps()

        # override default reference data for reduced runtime and better
        # defined tests
        self.params = parse_qiime_parameters(
            ["assign_taxonomy:id_to_taxonomy_fp	%s"
             % self.test_data['refseqs_tax'][0],
             "assign_taxonomy:reference_seqs_fp	%s"
             % self.test_data['refseqs'][0],
             "align_seqs:template_fp	%s"
             % self.test_data['refseqs_aligned'][0],
             "filter_alignment:lane_mask_fp	%s"
             % self.test_data['refseqs_aligned_lanemask'][0]])

        self.wf_out = mkdtemp(dir=self.tmp_dir,
                              prefix='qiime_wf_out', suffix='')
        self.dirs_to_remove.append(self.wf_out)

        initiate_timeout(allowed_seconds_per_test)

    def tearDown(self):
        """ """
        # turn off the alarm
        disable_timeout()
        # change back to the start dir - some workflows change directory
        chdir(self.start_dir)
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_pick_subsampled_open_reference_otus_step_1_and_4(self):
        """pick_subsampled_open_reference_otus functions as expected
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=100000)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU0) is not in "
                        "the final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_step_1_through_4(self):
        """pick_subsampled_open_reference_otus functions as expected
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in "
                        "the final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_rdp_tax_assign(self):
        """pick_subsampled_open_reference_otus fns when assigning tax with rdp
        """
        self.params.update(
            parse_qiime_parameters(
                ['assign_taxonomy:assignment_method rdp',
                 'assign_taxonomy:reference_seqs_fp %s'
                 % self.test_data['refseqs'][0],
                 'assign_taxonomy:id_to_taxonomy_fp %s'
                 % self.test_data['refseqs_tax'][0]]))
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        rdp_tax_fp = ('%s/rdp_assigned_taxonomy/rep_set_tax_assignments.txt'
                      % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")
        self.assertTrue(
            exists(rdp_tax_fp),
            "rdp taxonomy assignment result doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in "
                        "the final OTU map.")

        # confirm that number of tips in the tree is the same as the number
        # of sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_usearch(self):
        """pick_subsampled_open_reference_otus functions as expected with
        usearch
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            denovo_otu_picking_method='usearch61',
            reference_otu_picking_method='usearch61_ref',
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in "
                        "the final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        num_tree_tips = len(list(TreeNode.from_newick(open(tree_fp)).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_sortmerna_sumaclust(self):
        """pick_subsampled_open_reference_otus functions as expected with
        sortmerna and sumaclust (input reference database not indexed),
        the minimum_failure_threshold=100000 so steps 2 and 3 are skipped
        and all failures.fasta from step 1 are directly passed for de novo
        clustering in step 4 """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=None,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            denovo_otu_picking_method='sumaclust',
            reference_otu_picking_method='sortmerna',
            status_update_callback=no_status_updates,
            minimum_failure_threshold=100000)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        reads_for_denovo_fp = glob(
            '%s/step1_otus/*_failures.txt' %
            self.wf_out)[0]
        reads_for_otumap_fp = glob(
            '%s/step1_otus/*_otus.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(reads_for_denovo_fp),
            "Denovo reads file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")

        # number of OTUs in final OTU table equals the number of sequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])

        # the correct sequences failed the prefilter (E-value)
        # reads_for_denovo_ids all passed E-value and are considered
        # as true rRNA (having <97% id and/or <97% query coverage)
        reads_for_denovo_ids = [s.strip()
                                for s in open(reads_for_denovo_fp, 'U')]
        self.assertEqual(len(reads_for_denovo_ids), 59)

        self.assertTrue('t1_1' not in reads_for_denovo_ids)
        self.assertTrue('p1_2' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_130' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_151' not in reads_for_denovo_ids)

        # reads_for_otumap_ids all passed E-value and are considered
        # as true rRNA (having >=97% id and >=97% query coverage)
        reads_for_otumap_ids = [seq_id for line in
                                open(reads_for_otumap_fp, 'U')
                                for seq_id in line.strip().split('\t')]

        self.assertEqual(len(reads_for_otumap_ids), 120)

        self.assertTrue('t1_1' not in reads_for_otumap_ids)
        self.assertTrue('p1_2' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_130' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_151' not in reads_for_otumap_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU0) is not in "
                        "the final OTU map.")
        self.assertTrue('wf.test.otu.CleanUp.ReferenceOTU5' in otu_map,
                        "Failure OTU (wf.test.otu.CleanUp.ReferenceOTU5) is not in "
                        "the final OTU map.")

        # confirm size of full OTU map
        self.assertTrue(len(otu_map), 9)

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        num_tree_tips = len(list(TreeNode.from_newick(open(tree_fp)).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_sortmerna_sumaclust_step2_3(self):
        """pick_subsampled_open_reference_otus functions as expected with
        sortmerna and sumaclust (input reference database not indexed),
        minimum_failure_threshold=0 and percent_subsample=0.7, so a non-empty
        subsampled_failures.fasta file is created in step 1 and passed to steps 2/3/4
        for analysis. """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.7,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=None,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            denovo_otu_picking_method='sumaclust',
            reference_otu_picking_method='sortmerna',
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        reads_for_denovo_fp = glob(
            '%s/step1_otus/*_failures.txt' %
            self.wf_out)[0]
        reads_for_otumap_fp = glob(
            '%s/step1_otus/*_otus.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(reads_for_denovo_fp),
            "Denovo reads file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")

        # number of OTUs in final OTU table equals the number of sequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])

        # the correct sequences failed the prefilter (E-value)
        # reads_for_denovo_ids all passed E-value and are considered
        # as true rRNA (having <97% id and/or <97% query coverage)
        reads_for_denovo_ids = [s.strip()
                                for s in open(reads_for_denovo_fp, 'U')]
        self.assertEqual(len(reads_for_denovo_ids), 59)

        self.assertTrue('t1_1' not in reads_for_denovo_ids)
        self.assertTrue('p1_2' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_130' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_151' not in reads_for_denovo_ids)

        # reads_for_otumap_ids all passed E-value and are considered
        # as true rRNA (having >=97% id and >=97% query coverage)
        reads_for_otumap_ids = [seq_id for line in
                                open(reads_for_otumap_fp, 'U')
                                for seq_id in line.strip().split('\t')]

        self.assertEqual(len(reads_for_otumap_ids), 120)

        self.assertTrue('t1_1' not in reads_for_otumap_ids)
        self.assertTrue('p1_2' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_130' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_151' not in reads_for_otumap_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is in "
                        "the final OTU map.")

        # confirm size of full OTU map
        self.assertTrue(len(otu_map), 11)

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        num_tree_tips = len(list(TreeNode.from_newick(open(tree_fp)).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_sortmerna_sumaclust_db(self):
        """pick_subsampled_open_reference_otus functions as expected with
        sortmerna and sumaclust (input reference database indexed) """
        
        # Create sortmerna indexed database
        sortmerna_db, files_to_remove = \
            build_database_sortmerna(self.test_data['refseqs'][0],
                                     max_pos=10000,
                                     output_dir=self.tmp_dir)

        self.files_to_remove.extend(files_to_remove)

        self.params['pick_otus']['sortmerna_db'] = sortmerna_db

        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=None,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            denovo_otu_picking_method='sumaclust',
            reference_otu_picking_method='sortmerna',
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        reads_for_denovo_fp = glob(
            '%s/step1_otus/*_failures.txt' %
            self.wf_out)[0]
        reads_for_otumap_fp = glob(
            '%s/step1_otus/*_otus.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(reads_for_denovo_fp),
            "Denovo reads file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")

        # number of OTUs in final OTU table equals the number of sequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])

        # the correct sequences failed the prefilter (E-value)
        # reads_for_denovo_ids all passed E-value and are considered
        # as true rRNA (having <97% id and/or <97% query coverage)
        reads_for_denovo_ids = [s.strip()
                                for s in open(reads_for_denovo_fp, 'U')]
        self.assertEqual(len(reads_for_denovo_ids), 59)

        self.assertTrue('t1_1' not in reads_for_denovo_ids)
        self.assertTrue('p1_2' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_130' not in reads_for_denovo_ids)
        self.assertTrue('not16S.1_151' not in reads_for_denovo_ids)

        # reads_for_otumap_ids all passed E-value and are considered
        # as true rRNA (having >=97% id and >=97% query coverage)
        reads_for_otumap_ids = [seq_id for line in
                                open(reads_for_otumap_fp, 'U')
                                for seq_id in line.strip().split('\t')]

        self.assertEqual(len(reads_for_otumap_ids), 120)

        self.assertTrue('t1_1' not in reads_for_otumap_ids)
        self.assertTrue('p1_2' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_130' not in reads_for_otumap_ids)
        self.assertTrue('not16S.1_151' not in reads_for_otumap_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in "
                        "the final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        num_tree_tips = len(list(TreeNode.from_newick(open(tree_fp)).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_suppress_assign_tax(self):
        """pick_subsampled_open_reference_otus functions without assign tax
        step
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
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

        otu_table_fp = '%s/otu_table_mc2_no_pynast_failures.biom' % self.wf_out
        otu_table_w_tax_fp = '%s/otu_table_mc2_w_tax.biom' % self.wf_out

        self.assertFalse(exists(otu_table_w_tax_fp),
                         "OTU table w tax exists (it shouldn't).")
        self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")

    def test_pick_subsampled_open_reference_otus_suppress_align_and_tree(self):
        """pick_subsampled_open_reference_otus functions without align/tree
        step
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
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

        otu_table_fp = '%s/otu_table_mc2_w_tax.biom' % self.wf_out
        tree_fp = '%s/rep_set.tre' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
        self.assertFalse(exists(tree_fp),
                         "Tree exists (it shouldn't).")
        self.assertFalse(exists(pynast_failures_fp),
                         "PyNAST failures file exists (it shouldn't).")

    def test_pick_subsampled_open_reference_otus_no_prefilter(self):
        """pick_subsampled_open_reference_otus functions as expected without
        prefilter
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=None,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)

        prefilter_output_directory = '%s/prefilter_otus/'
        self.assertFalse(exists(prefilter_output_directory))

        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6 (note: this is the same as without the prefilter
        # because these reads are getting filtered from the final otu table because
        # they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # sequences that are ordinarily prefiltered are in the OTU map
        otu_map_seq_ids = []
        for l in open(otu_map_w_singletons_fp, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertTrue('t1_1' in otu_map_seq_ids)
        self.assertTrue('p1_2' in otu_map_seq_ids)
        self.assertTrue('not16S.1_130' in otu_map_seq_ids)
        self.assertTrue('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_parallel(self):
        """pick_subsampled_open_reference_otus functions as expected in
        parallel
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=True,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)
        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertFalse(exists(final_failure_fp),
                         "Final failures file shouldn't exist, but it does")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 6
        self.assertEqual(len(otu_table.ids(axis='observation')), 6)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 6)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment (can't confirm the )
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_suppress_step4(self):
        """pick_subsampled_open_reference_otus functions as expected wo step 4
        """
        pick_subsampled_open_reference_otus(
            input_fp=self.test_data['seqs'][0],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            qiime_config=self.qiime_config,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=True,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)

        step4_output_dir = '%s/step4_otus/' % self.wf_out
        self.assertFalse(exists(step4_output_dir),
                         "Step 4 output directory exists - it shouldn't.")

        otu_map_w_singletons_fp = '%s/final_otu_map.txt' % self.wf_out
        final_failure_fp = '%s/final_failures.txt' % self.wf_out
        final_repset_fp = '%s/rep_set.fna' % self.wf_out
        new_refseqs_fp = '%s/new_refseqs.fna' % self.wf_out
        prefilter_failures_fp = glob(
            '%s/prefilter_otus/*_failures.txt' %
            self.wf_out)[0]
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)

        self.assertTrue(
            exists(otu_map_w_singletons_fp),
            "OTU map doesn't exist")
        self.assertTrue(
            exists(final_failure_fp),
            "Final failures file doesn't exist")
        self.assertTrue(
            exists(final_repset_fp),
            "Final representative set doesn't exist")
        self.assertTrue(
            exists(new_refseqs_fp),
            "New refseqs file doesn't exist")
        self.assertTrue(
            exists(prefilter_failures_fp),
            "Prefilter failures file doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(otu_table_fp), "Final BIOM table doesn't exist")
        self.assertTrue(
            exists(pynast_failures_fp),
            "PyNAST failures file doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of sequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is either 4, 5 or 6 (it can vary due to the random
        # subsampling)
        num_obs_ids = len(otu_table.ids(axis='observation'))
        self.assertTrue(4 <= num_obs_ids <= 6,
                        "Obtained %d observations, but expected 4, 5 or 6."
                        % num_obs_ids)

        # the correct sequences failed the prefilter
        prefilter_failure_ids = [s.strip()
                                 for s in open(prefilter_failures_fp, 'U')]
        self.assertEqual(len(prefilter_failure_ids), 24)
        self.assertTrue('t1_1' in prefilter_failure_ids)
        self.assertTrue('p1_2' in prefilter_failure_ids)
        self.assertTrue('not16S.1_130' in prefilter_failure_ids)
        self.assertTrue('not16S.1_151' in prefilter_failure_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        #
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(otu_map_w_singletons_fp))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        self.assertTrue('wf.test.otu.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment (and we already checked that we're happy
        # with that number above when it was compared to the number of OTUs)
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        otu_table = load_table(otu_table_fp)
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_pick_subsampled_open_reference_otus_invalid_input(self):
        """pick_subsampled_open_reference_otus raises error on refseqs in
        params file
        """
        self.assertRaises(WorkflowError,
                          pick_subsampled_open_reference_otus,
                          input_fp=self.test_data['seqs'][0],
                          refseqs_fp=self.test_data['refseqs'][0],
                          output_dir=self.wf_out,
                          percent_subsample=0.5,
                          new_ref_set_id='wf.test.otu',
                          command_handler=call_commands_serially,
                          params=parse_qiime_parameters(
                              ['pick_otus:refseqs_fp %s'
                               % self.test_data['refseqs'][0]]),
                          prefilter_refseqs_fp=None,
                          qiime_config=self.qiime_config,
                          step1_otu_map_fp=None,
                          step1_failures_fasta_fp=None,
                          parallel=False,
                          suppress_step4=False,
                          logger=None,
                          status_update_callback=no_status_updates,
                          minimum_failure_threshold=0)

    def test_iterative_pick_subsampled_open_reference_otus_no_prefilter(self):
        """pick_subsampled_open_reference_otus functions as expected without
        prefilter
        """
        iterative_pick_subsampled_open_reference_otus(
            input_fps=[self.test_data['seqs'][0],
                       self.test_data['extra_seqs'][0]],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            prefilter_refseqs_fp=None,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)

        for i in (0, 1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out, i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out, i)
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, i)
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out, i)
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out, i)

            self.assertTrue(
                exists(final_otu_map_fp),
                "Final OTU map doesn't exist")
            self.assertFalse(
                exists(final_failure_fp),
                "Final failures file shouldn't exist, but it does")
            self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
            self.assertTrue(exists(repset_fp), "Rep set doesn't exist.")
            self.assertTrue(
                exists(new_refseqs_fp),
                "New refseqs file doesn't exist.")

        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        iter0_otu_map_w_singletons = '%s/0/final_otu_map.txt' % self.wf_out
        iter1_otu_map_w_singletons = '%s/1/final_otu_map.txt' % self.wf_out

        self.assertTrue(exists(otu_table_fp), "Final OTU table doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(pynast_failures_fp),
                        "PyNAST failurs file doesn't exist")
        self.assertTrue(exists(iter0_otu_map_w_singletons),
                        "Iteration 0 OTU map with singletons doesn't exist")
        self.assertTrue(exists(iter1_otu_map_w_singletons),
                        "Iteration 1 OTU map with singletons doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 7 (note: this is the same as without the prefilter
        # because these reads are getting filtered from the final otu table because
        # they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 7)

        # sequences that are ordinarily prefiltered are in the OTU map
        otu_map_seq_ids = []
        for l in open(iter0_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertTrue('t1_1' in otu_map_seq_ids)
        self.assertTrue('p1_2' in otu_map_seq_ids)
        self.assertTrue('not16S.1_130' in otu_map_seq_ids)
        self.assertTrue('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(iter0_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        # OTU from first iteration is in final map
        self.assertTrue('wf.test.otu.0.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")
        # OTU from second iteration is in final map
        otu_map = fields_to_dict(open(iter1_otu_map_w_singletons))
        self.assertTrue('wf.test.otu.1.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")

        # OTUs from each iteration are in the final OTU table
        self.assertTrue('295053' in otu_table.ids(axis='observation'))
        self.assertTrue(
            'wf.test.otu.1.ReferenceOTU0' in otu_table.ids(axis='observation'))

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 7)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_iterative_pick_subsampled_open_reference_otus(self):
        """pick_subsampled_open_reference_otus functions as expected with
        prefilter
        """
        iterative_pick_subsampled_open_reference_otus(
            input_fps=[self.test_data['seqs'][0],
                       self.test_data['extra_seqs'][0]],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)

        for i in (0, 1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out, i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out, i)
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, i)
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out, i)
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out, i)

            self.assertTrue(
                exists(final_otu_map_fp),
                "Final OTU map doesn't exist")
            self.assertFalse(
                exists(final_failure_fp),
                "Final failures file shouldn't exist, but it does")
            self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
            self.assertTrue(exists(repset_fp), "Rep set doesn't exist.")
            self.assertTrue(
                exists(new_refseqs_fp),
                "New refseqs file doesn't exist.")

        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        iter0_otu_map_w_singletons = '%s/0/final_otu_map.txt' % self.wf_out
        iter1_otu_map_w_singletons = '%s/1/final_otu_map.txt' % self.wf_out

        self.assertTrue(exists(otu_table_fp), "Final OTU table doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(pynast_failures_fp),
                        "PyNAST failurs file doesn't exist")
        self.assertTrue(exists(iter0_otu_map_w_singletons),
                        "Iteration 0 OTU map with singletons doesn't exist")
        self.assertTrue(exists(iter1_otu_map_w_singletons),
                        "Iteration 1 OTU map with singletons doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 7 (note: this is the same as without the prefilter
        # because these reads are getting filtered from the final otu table because
        # they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 7)

        # non-16S sequences are prefiltered, so not in the OTU map
        otu_map_seq_ids = []
        for l in open(iter0_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertFalse('t1_1' in otu_map_seq_ids)
        self.assertFalse('p1_2' in otu_map_seq_ids)
        self.assertFalse('not16S.1_130' in otu_map_seq_ids)
        self.assertFalse('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(iter0_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        # OTU from first iteration is in final map
        self.assertTrue('wf.test.otu.0.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map (first iteration).")
        # OTU from second iteration is in final map
        otu_map = fields_to_dict(open(iter1_otu_map_w_singletons))
        self.assertTrue('wf.test.otu.1.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map (second iteration).")

        # OTUs from each iteration are in the final OTU table
        self.assertTrue('295053' in otu_table.ids(axis='observation'))
        self.assertTrue(
            'wf.test.otu.1.ReferenceOTU0' in otu_table.ids(axis='observation'))

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 7)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_iterative_pick_subsampled_open_reference_otus_parallel(self):
        """pick_subsampled_open_reference_otus functions as expected in
        parallel
        """
        iterative_pick_subsampled_open_reference_otus(
            input_fps=[self.test_data['seqs'][0],
                       self.test_data['extra_seqs'][0]],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=True,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0)

        for i in (0, 1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out, i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out, i)
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, i)
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out, i)
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out, i)

            self.assertTrue(
                exists(final_otu_map_fp),
                "Final OTU map doesn't exist")
            self.assertFalse(
                exists(final_failure_fp),
                "Final failures file shouldn't exist, but it does")
            self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
            self.assertTrue(exists(repset_fp), "Rep set doesn't exist.")
            self.assertTrue(
                exists(new_refseqs_fp),
                "New refseqs file doesn't exist.")

        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        iter0_otu_map_w_singletons = '%s/0/final_otu_map.txt' % self.wf_out
        iter1_otu_map_w_singletons = '%s/1/final_otu_map.txt' % self.wf_out

        self.assertTrue(exists(otu_table_fp), "Final OTU table doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(pynast_failures_fp),
                        "PyNAST failurs file doesn't exist")
        self.assertTrue(exists(iter0_otu_map_w_singletons),
                        "Iteration 0 OTU map with singletons doesn't exist")
        self.assertTrue(exists(iter1_otu_map_w_singletons),
                        "Iteration 1 OTU map with singletons doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 7 (note: this is the same as without the prefilter
        # because these reads are getting filtered from the final otu table because
        # they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 7)

        # non-16S sequences are prefiltered, so not in the OTU map
        otu_map_seq_ids = []
        for l in open(iter0_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertFalse('t1_1' in otu_map_seq_ids)
        self.assertFalse('p1_2' in otu_map_seq_ids)
        self.assertFalse('not16S.1_130' in otu_map_seq_ids)
        self.assertFalse('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(iter0_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        # OTU from first iteration is in final map
        self.assertTrue('wf.test.otu.0.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")
        # OTU from second iteration is in final map
        otu_map = fields_to_dict(open(iter1_otu_map_w_singletons))
        self.assertTrue('wf.test.otu.1.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.ReferenceOTU0) is not in the"
                        " final OTU map.")

        # OTUs from each iteration are in the final OTU table
        self.assertTrue('295053' in otu_table.ids(axis='observation'))
        self.assertTrue(
            'wf.test.otu.1.ReferenceOTU0' in otu_table.ids(axis='observation'))

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 7)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053',axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

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

        exp = [("o1", "ACCGT"), ("o2", "AGG"), ("o3", "ACCGTT"),
               ("o4", "TACCGT"), ("o5", "TAGG"), ("o6", "TACCGTT")]
        actual = list(final_repset_from_iteration_repsets([repset1, repset2]))
        self.assertEqual(actual, exp)

        exp = [("o1", "ACCGT"), ("o2", "AGG"), ("o3", "ACCGTT"),
               ("o4", "TACCGT"), ("o5", "TAGG"), ("o6", "TACCGTT"),
               ('o7', 'AAAA')]
        actual = list(
            final_repset_from_iteration_repsets([repset1, repset2, repset3]))
        self.assertEqual(actual, exp)

    def test_convergent_pick_subsampled_open_reference_otus_no_prefilter(self):
        """pick_subsampled_open_reference_otus functions as expected without
        prefilter in convergent mode
        """
        convergent_pick_subsampled_open_reference_otus(
            input_fps=[self.test_data['seqs'][0],
                       self.test_data['extra_seqs'][0]],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            prefilter_refseqs_fp=None,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0,
            num_seqs=90)

        for i in (0, 1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out, i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out, i)
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, i)
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out, i)
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out, i)

            self.assertTrue(
                exists(final_otu_map_fp),
                "Final OTU map doesn't exist")
            self.assertFalse(
                exists(final_failure_fp),
                "Final failures file shouldn't exist, but it does")
            self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
            self.assertTrue(exists(repset_fp), "Rep set doesn't exist.")
            self.assertTrue(
                exists(new_refseqs_fp),
                "New refseqs file doesn't exist.")

        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        iter0_otu_map_w_singletons = '%s/0/final_otu_map.txt' % self.wf_out
        iter1_otu_map_w_singletons = '%s/1/final_otu_map.txt' % self.wf_out

        self.assertTrue(exists(otu_table_fp), "Final OTU table doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(pynast_failures_fp),
                        "PyNAST failures file doesn't exist")
        self.assertTrue(exists(iter0_otu_map_w_singletons),
                        "Iteration 0 OTU map with singletons doesn't exist")
        self.assertTrue(exists(iter1_otu_map_w_singletons),
                        "Iteration 1 OTU map with singletons doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 7 (note: this is the same as without the
        # prefilter because these reads are getting filtered from the final
        # otu table because they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 7)

        # sequences that are ordinarily prefiltered are in the OTU map
        otu_map_seq_ids = []
        for l in open(iter0_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertTrue('p2_31582' in otu_map_seq_ids)
        self.assertTrue('f1_2278' in otu_map_seq_ids)

        otu_map_seq_ids = []
        for l in open(iter1_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertTrue('not16S.1_130' in otu_map_seq_ids)
        self.assertTrue('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(iter0_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        # OTU from first iteration is in final map
        self.assertTrue('wf.test.otu.Iter.0.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.Iter.0.ReferenceOTU0) is not"
                        " in the final OTU map.")
        # OTU from second iteration is in final map
        otu_map = fields_to_dict(open(iter1_otu_map_w_singletons))
        self.assertTrue('wf.test.otu.Iter.1.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.Iter.1.ReferenceOTU0) is not"
                        " in the final OTU map.")

        # in the test dataset, the OTUs from the second iteration are filtered
        # out due to pynast failures. Thus, we will test that OTUs from each
        # iteration are in the unfiltered table
        unf_table = load_table('%s/otu_table_mc2.biom' % self.wf_out)
        self.assertTrue('295053' in unf_table.ids(axis='observation'))
        self.assertTrue(
            'wf.test.otu.Iter.1.ReferenceOTU0' in
            unf_table.ids(axis='observation'))

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 7)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

    def test_convergent_pick_subsampled_open_reference_otus(self):
        """pick_subsampled_open_reference_otus functions as expected with
        prefilter in convergent mode
        """
        convergent_pick_subsampled_open_reference_otus(
            input_fps=[self.test_data['seqs'][0],
                       self.test_data['extra_seqs'][0]],
            refseqs_fp=self.test_data['refseqs'][0],
            output_dir=self.wf_out,
            percent_subsample=0.5,
            new_ref_set_id='wf.test.otu',
            command_handler=call_commands_serially,
            params=self.params,
            qiime_config=self.qiime_config,
            prefilter_refseqs_fp=None,
            prefilter_percent_id=0.60,
            step1_otu_map_fp=None,
            step1_failures_fasta_fp=None,
            parallel=False,
            suppress_step4=False,
            logger=None,
            status_update_callback=no_status_updates,
            minimum_failure_threshold=0,
            num_seqs=150)

        for i in (0, 1):
            final_otu_map_fp = '%s/%d/final_otu_map.txt' % (self.wf_out, i)
            final_failure_fp = '%s/%d/final_failures.txt' % (self.wf_out, i)
            otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, i)
            repset_fp = '%s/%d/rep_set.fna' % (self.wf_out, i)
            new_refseqs_fp = '%s/%d/new_refseqs.fna' % (self.wf_out, i)

            self.assertTrue(
                exists(final_otu_map_fp),
                "Final OTU map doesn't exist")
            self.assertFalse(
                exists(final_failure_fp),
                "Final failures file shouldn't exist, but it does")
            self.assertTrue(exists(otu_table_fp), "OTU table doesn't exist.")
            self.assertTrue(exists(repset_fp), "Rep set doesn't exist.")
            self.assertTrue(
                exists(new_refseqs_fp),
                "New refseqs file doesn't exist.")

        otu_table_fp = ('%s/otu_table_mc2_w_tax_no_pynast_failures.biom'
                        % self.wf_out)
        tree_fp = '%s/rep_set.tre' % self.wf_out
        aln_fp = '%s/pynast_aligned_seqs/rep_set_aligned.fasta' % self.wf_out
        pynast_failures_fp = ('%s/pynast_aligned_seqs/rep_set_failures.fasta'
                              % self.wf_out)
        iter0_otu_map_w_singletons = '%s/0/final_otu_map.txt' % self.wf_out
        iter1_otu_map_w_singletons = '%s/1/final_otu_map.txt' % self.wf_out

        self.assertTrue(exists(otu_table_fp), "Final OTU table doesn't exist")
        self.assertTrue(exists(tree_fp), "Final tree doesn't exist")
        self.assertTrue(exists(aln_fp), "Final alignment doesn't exist")
        self.assertTrue(exists(pynast_failures_fp),
                        "PyNAST failurs file doesn't exist")
        self.assertTrue(exists(iter0_otu_map_w_singletons),
                        "Iteration 0 OTU map with singletons doesn't exist")
        self.assertTrue(exists(iter1_otu_map_w_singletons),
                        "Iteration 1 OTU map with singletons doesn't exist")

        # all OTUs in final OTU table occur more than once
        otu_table = load_table(otu_table_fp)
        for row in otu_table.iter_data(axis='observation'):
            self.assertTrue(sum(row) >= 2,
                            "Singleton OTU detected in OTU table.")
        # number of OTUs in final OTU table equals the number of seequences in
        # the alignment...
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         count_seqs(aln_fp)[0])
        # ... and that number is 7 (note: this is the same as without the
        # prefilter because these reads are getting filtered from the final
        # otu table because they fail to align with PyNAST)
        self.assertEqual(len(otu_table.ids(axis='observation')), 7)

        # non-16S sequences are prefiltered, so not in the OTU map
        otu_map_seq_ids = []
        for l in open(iter0_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertFalse('t1_1' in otu_map_seq_ids)
        self.assertFalse('p1_2' in otu_map_seq_ids)
        self.assertFalse('not16S.1_130' in otu_map_seq_ids)
        self.assertFalse('not16S.1_151' in otu_map_seq_ids)

        otu_map_seq_ids = []
        for l in open(iter1_otu_map_w_singletons, 'U'):
            otu_map_seq_ids.extend(l.strip().split()[1:])
        self.assertFalse('t1_1' in otu_map_seq_ids)
        self.assertFalse('p1_2' in otu_map_seq_ids)
        self.assertFalse('not16S.1_130' in otu_map_seq_ids)
        self.assertFalse('not16S.1_151' in otu_map_seq_ids)

        # confirm that the new reference sequences is the same length as the
        # input reference sequences plus the number of new non-singleton otus
        self.assertEqual(count_seqs(new_refseqs_fp)[0],
                         count_seqs(self.test_data['refseqs'][0])[0] +
                         len([o for o in otu_table.ids(axis='observation')
                              if o.startswith('wf.test.otu')]) +
                         count_seqs(pynast_failures_fp)[0])

        # spot check a few of the otus to confirm that we're getting reference
        # and new otus in the final otu map. This is done on the OTU map
        # singletons get filtered before building the otu table
        otu_map = fields_to_dict(open(iter0_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map.")
        self.assertTrue('42684' in otu_map,
                        "Failure OTU (42684) is not in the final OTU map.")
        # OTU from first iteration is in final map
        self.assertTrue('wf.test.otu.Iter.0.ReferenceOTU0' in otu_map,
                        "Failure OTU (wf.test.otu.Iter.0.ReferenceOTU0) is not"
                        " in the final OTU map (first iteration).")
        # OTU from second iteration is in final map
        # Due to prefiltering, only the reference OTU map is kept
        otu_map = fields_to_dict(open(iter1_otu_map_w_singletons))
        self.assertTrue('295053' in otu_map,
                        "Reference OTU (295053) is not in the final OTU map "
                        "(second iteration).")

        # In the test dataset, the only OTU present in the second OTU table is
        # a reference OTU that is also present in the first iteration. To check
        # that the results form the second iteration are also present in the
        # final OTU table we are going to get the counts of that OTU
        iter0_otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, 0)
        iter1_otu_table_fp = '%s/%d/otu_table_mc2.biom' % (self.wf_out, 1)

        i0_counts = load_table(iter0_otu_table_fp).data(
            '295053', axis='observation').sum()
        i1_counts = load_table(iter1_otu_table_fp).data(
            '295053', axis='observation').sum()
        final_counts = load_table('%s/otu_table_mc2.biom' % self.wf_out).data(
            '295053', axis='observation').sum()
        self.assertEqual(final_counts, i0_counts + i1_counts)

        # confirm that number of tips in the tree is the same as the number of
        # sequences in the alignment
        with open(tree_fp) as f:
            num_tree_tips = len(list(TreeNode.from_newick(f).tips()))
        num_align_seqs = count_seqs(aln_fp)[0]
        self.assertEqual(num_tree_tips, num_align_seqs)
        self.assertEqual(num_tree_tips, 7)

        # OTU table without singletons or pynast failures has same number of
        # otus as there are aligned sequences
        self.assertEqual(len(otu_table.ids(axis='observation')),
                         num_align_seqs)

        # Reference OTUs have correct taxonomy assignment
        obs_idx = otu_table.index('295053', axis='observation')
        self.assertEqual(
            otu_table.metadata(axis='observation')[obs_idx]['taxonomy'],
            ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__", "s__"])
        # All observations have 'taxonomy' metadata, and are at least assigned
        # to 'bacteria'
        for o in otu_table.iter(axis='observation'):
            self.assertTrue(
                o[2]['taxonomy'][0] in ['k__Bacteria', 'Unassigned'])

if __name__ == "__main__":
    main()
