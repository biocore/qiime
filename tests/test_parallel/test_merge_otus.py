#!/usr/bin/env python

from unittest import TestCase, main
from qiime.parallel.merge_otus import mergetree, mergeorder, \
    initial_nodes_to_merge, initial_has_dependencies, job_complete, \
    torque_job, local_job, start_job, JobError, reset_internal_count
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


class MergeTests(TestCase):

    def setUp(self):
        reset_internal_count()

    def test_mergetree(self):
        """construct a merge subtreetree with various properties set"""
        exp = "(A,B)0;"
        obs = mergetree(['A.biom'], ['B.biom'], 'foo')
        self.assertEqual(obs.getNewick(escape_name=False), exp)

        self.assertEqual(obs.Children[0].Name, 'A')
        self.assertEqual(obs.Children[0].FilePath, 'A.biom')
        self.assertEqual(obs.Children[0].Processed, False)
        self.assertEqual(obs.Children[0].PollPath, None)
        self.assertEqual(obs.Children[0].FullCommand, None)

        self.assertEqual(obs.Children[1].Name, 'B')
        self.assertEqual(obs.Children[1].FilePath, 'B.biom')
        self.assertEqual(obs.Children[1].Processed, False)
        self.assertEqual(obs.Children[1].PollPath, None)
        self.assertEqual(obs.Children[1].FullCommand, None)

        self.assertEqual(obs.Name, '0')
        self.assertEqual(obs.FilePath, 'foo/0.biom')
        self.assertEqual(obs.Processed, False)
        self.assertEqual(obs.PollPath, 'foo/0.biom.poll')
        self.assertEqual(obs.FullCommand, None)

    def test_mergeorder(self):
        """recursively build and join all the subtrees"""
        exp = "((A,B)0,(C,(D,E)1)2)3;"
        obs = mergeorder(['A', 'B', 'C', 'D', 'E'], 'foo')
        self.assertEqual(obs.getNewick(escape_name=False), exp)

    def test_initial_nodes_to_merge(self):
        """determine the first nodes to merge"""
        t = mergeorder(['A', 'B', 'C', 'D', 'E'], 'foo')
        exp = set([t.Children[0], t.Children[1].Children[1]])
        obs = initial_nodes_to_merge(t)
        self.assertEqual(obs, exp)

    def test_initial_has_dependencies(self):
        """determine initial has_dependencies"""
        t = mergeorder(['A', 'B', 'C', 'D', 'E'], 'foo')
        exp = [t, t.Children[1]]
        obs = initial_has_dependencies(t, initial_nodes_to_merge(t))
        self.assertEqual(obs, exp)

    def test_job_complete(self):
        """check if a job is complete"""
        t = mergeorder(['A', 'B', 'C', 'D', 'E'], 'foo')
        self.assertFalse(job_complete(t))
        self.assertFalse(job_complete(t.Children[0]))
        self.assertFalse(job_complete(t.Children[1].Children[1]))

        self.assertRaises(JobError, job_complete, t.Children[0].Children[0])

        f = 'test_parallel_merge_otus_JOB_COMPLETE_TEST.poll'
        self.assertFalse(os.path.exists(f))

        testf = open(f, 'w')
        testf.write('0\n')
        testf.close()
        t.PollPath = f
        t.StartTime = 10

        self.assertTrue(job_complete(t))
        self.assertNotEqual(t.EndTime, None)
        self.assertNotEqual(t.TotalTime, None)

        testf = open(f, 'w')
        testf.write('1\n')
        testf.close()

        self.assertRaises(JobError, job_complete, t)
        t.Processed = False
        self.assertRaises(JobError, job_complete, t)

        os.remove(f)

    def test_torque_job(self):
        """wrap a torque job"""
        exp = 'echo "abc; echo $? > xyz" | qsub -k oe -N MOTU -q queue'
        obs = torque_job('abc', 'xyz', '123', 'queue')
        self.assertEqual(obs, exp)

    def test_start_job(self):
        """start a job"""
        exp = 'echo "y -i A.biom,B.biom -o foo/0.biom; echo $? > foo/0.biom.poll" | qsub -k oe -N MOTU -q ignored'
        t = mergeorder(['A.biom', 'B.biom', 'C', 'D', 'E'], 'foo')
        start_job(t.Children[0], 'y', 'ignored', torque_job, False)
        self.assertEqual(t.Children[0].FullCommand, exp)

    def test_local_job(self):
        """fire off a local job"""
        exp = "abc; echo $? > xyz"
        obs = local_job('abc', 'xyz', 'notused', 'notused')
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
