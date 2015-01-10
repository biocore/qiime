#!/usr/bin/env python

"""Tests for submitting jobs via qsub."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, environ, close
from os.path import exists
from time import sleep
from tempfile import mkstemp

from unittest import TestCase, main


from qiime.denoiser.make_cluster_jobs import QSUB_TEXT, make_jobs, \
    submit_jobs


class Test_make_cluster_jobs(TestCase):

    def setUp(self):

        self.home = environ['HOME']
        self.queue = "friendlyq"

        fd, self.tmp_result_file = mkstemp(dir=self.home,
                                           prefix="test_hello_",
                                           suffix=".txt")
        close(fd)
        self.command = "echo hello > %s\n" % self.tmp_result_file
        fd, self.tmp_name = mkstemp(dir="/tmp",
                                    prefix="make_cluster_jobs_test_",
                                    suffix=".txt")
        close(fd)
        fh = open(self.tmp_name, "w")
        fh.write(self.command)
        fh.close()

    def tearDown(self):
        remove(self.tmp_name)
        if exists(self.tmp_result_file):
            remove(self.tmp_result_file)

    def test_make_jobs(self):
        """writing the job files works"""
        # no commands should make no jobs files
        self.assertEqual(make_jobs([], "test", self.queue), [])

        # one job file should be created
        filenames = make_jobs([self.command], "test_qsub", self.queue)
        self.assertTrue(len(filenames) == 1)
        observed_text = list(open(filenames[0]))

        self.assertEqual("".join(observed_text),
                         QSUB_TEXT % ("72:00:00", 1, 1, self.queue,
                                      "test_qsub", "oe",
                                      self.command))


if __name__ == "__main__":
    main()
