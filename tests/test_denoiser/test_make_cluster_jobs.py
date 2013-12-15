#!/usr/bin/env python

"""Tests for submitting jobs via qsub."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, environ
from os.path import exists
from time import sleep

from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename


from qiime.denoiser.make_cluster_jobs import QSUB_TEXT, make_jobs, \
    submit_jobs

class Test_make_cluster_jobs(TestCase):

    def setUp(self):
        
        self.home = environ['HOME']
        self.queue = "friendlyq"

        self.tmp_result_file = get_tmp_filename(tmp_dir= self.home,
                                                prefix = "/test_hello_",
                                                suffix =".txt")
        self.command = "echo hello > %s\n" % self.tmp_result_file
        self.tmp_name = get_tmp_filename(tmp_dir="/tmp",
                                         prefix="make_cluster_jobs_test_",
                                         suffix = ".txt")
        fh = open(self.tmp_name,"w")
        fh.write(self.command)
        fh.close()

    def tearDown(self):
        remove(self.tmp_name)
        if exists(self.tmp_result_file):
            remove(self.tmp_result_file)
            
    def test_make_jobs(self):
        """writing the job files works"""
        #no commands should make no jobs files
        self.assertEqual(make_jobs([], "test", self.queue), [])
     
        #one job file should be created
        filenames = make_jobs([self.command], "test_qsub", self.queue)   
        self.assertTrue(len(filenames)==1)
        observed_text= list(open(filenames[0]))

        self.assertEqual("".join(observed_text),
                         QSUB_TEXT % ("72:00:00", 1, 1, self.queue,
                                      "test_qsub", "oe",
                                      self.command))


if __name__ == "__main__":
    main()
