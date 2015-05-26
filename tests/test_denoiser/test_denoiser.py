#!/usr/bin/env python

"""tests for denoiser functions
"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

import signal
import os
from os import remove, rmdir

from shutil import rmtree
from subprocess import Popen, PIPE, STDOUT

from unittest import TestCase, main
from skbio.parse.sequences import parse_fasta

from qiime.util import get_qiime_project_dir
from qiime.denoiser.utils import check_flowgram_ali_exe
from qiime.denoiser.preprocess import make_tmp_name

PROJECT_HOME = get_qiime_project_dir()

# timeout handling taken from test_workflow.py


class TimeExceededError(Exception):
    pass

allowed_seconds_per_test = 240


def timeout(signum, frame):
    raise TimeExceededError("Test failed to run in allowed time (%d seconds)."
                            % allowed_seconds_per_test)


class DenoiserTests(TestCase):

    def setUp(self):
        # abort all tests without the alignment binary
        check_flowgram_ali_exe()

        signal.signal(signal.SIGALRM, timeout)
        # set the 'alarm' to go off in allowed_seconds seconds
        signal.alarm(allowed_seconds_per_test)

        self.test_dir = "denoiser_main_test" + make_tmp_name() + "/"
        self.expected = ">FS8APND01D3TW3 | cluster size: 94 \nCTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACCCGTCAAAGCCTTGGTAAGCCACTACCCCACCAACAAGCTGATAAGCCGCGAGTCCATCCCCAACCGCCGAAACTTTCCAACCCCCACCATGCAGCAGGAGCTCCTATCCGGTATTAGCCCCAGTTTCCTGAAGTTATCCCAAAGTCAAGGGCAGGTTACTCACGTGTTACTCACCCGTTCGCC\n"
        self.expected_map_string = """FS8APND01EWRS4:
FS8APND01BSTVP:
FS8APND01DXG45:
FS8APND01D3TW3:\tFS8APND01CSXFN\tFS8APND01DQ8MX\tFS8APND01DY7QW\tFS8APND01B5QNI\tFS8APND01CQ6OG\tFS8APND01C7IGN\tFS8APND01DHSGH\tFS8APND01DJ17E\tFS8APND01CUXOA\tFS8APND01EUTYG\tFS8APND01EKK7T\tFS8APND01D582W\tFS8APND01B5GWU\tFS8APND01D7N2A\tFS8APND01BJGHZ\tFS8APND01D6DYZ\tFS8APND01C6ZIM\tFS8APND01D2X6Y\tFS8APND01BUYCE\tFS8APND01BNUEY\tFS8APND01DKLOE\tFS8APND01C24PP\tFS8APND01EBWQX\tFS8APND01ELDYW\tFS8APND01B0GCS\tFS8APND01D4QXI\tFS8APND01EMYD9\tFS8APND01EA2SK\tFS8APND01DZOSO\tFS8APND01DHYAZ\tFS8APND01C7UD9\tFS8APND01BTZFV\tFS8APND01CR78R\tFS8APND01B39IE\tFS8APND01ECVC0\tFS8APND01DM3PL\tFS8APND01DELWS\tFS8APND01CIEK8\tFS8APND01D7ZOZ\tFS8APND01CZSAI\tFS8APND01DYOVR\tFS8APND01BX9XY\tFS8APND01DEWJA\tFS8APND01BEKIW\tFS8APND01DCKB9\tFS8APND01EEYIS\tFS8APND01DDKEA\tFS8APND01DSZLO\tFS8APND01C6EBC\tFS8APND01EE15M\tFS8APND01ELO9B\tFS8APND01C58QY\tFS8APND01DONCG\tFS8APND01DVXX2\tFS8APND01BL5YT\tFS8APND01BIL2V\tFS8APND01EBSYQ\tFS8APND01CCX8R\tFS8APND01B2YCJ\tFS8APND01B1JG4\tFS8APND01DJ024\tFS8APND01BIJY0\tFS8APND01CIA4G\tFS8APND01DV74M\tFS8APND01ECAX5\tFS8APND01DC3TZ\tFS8APND01EJVO6\tFS8APND01D4VFG\tFS8APND01DYYYO\tFS8APND01D1EDD\tFS8APND01DQUOT\tFS8APND01A2NSJ\tFS8APND01DDC8I\tFS8APND01BP1T2\tFS8APND01DPY6U\tFS8APND01CIQGV\tFS8APND01BPUT8\tFS8APND01BDNH4\tFS8APND01DOZDN\tFS8APND01DS866\tFS8APND01DGS2J\tFS8APND01EDK32\tFS8APND01EPA0T\tFS8APND01CK3JM\tFS8APND01BKLWW\tFS8APND01DV0BO\tFS8APND01DPNXE\tFS8APND01B7LUA\tFS8APND01BTTE2\tFS8APND01CKO4X\tFS8APND01DGGBY\tFS8APND01C4NHX\tFS8APND01DYPQN
FS8APND01EFK0W:
FS8APND01DCIOO:
FS8APND01CKOMZ:
"""
        self.expected_titanium_map_string = """FS8APND01EWRS4:
FS8APND01BSTVP:
FS8APND01DXG45:
FS8APND01D3TW3:\tFS8APND01CSXFN\tFS8APND01DQ8MX\tFS8APND01DY7QW\tFS8APND01B5QNI\tFS8APND01CQ6OG\tFS8APND01C7IGN\tFS8APND01DHSGH\tFS8APND01DJ17E\tFS8APND01CUXOA\tFS8APND01EUTYG\tFS8APND01EKK7T\tFS8APND01D582W\tFS8APND01B5GWU\tFS8APND01D7N2A\tFS8APND01BJGHZ\tFS8APND01D6DYZ\tFS8APND01C6ZIM\tFS8APND01D2X6Y\tFS8APND01BUYCE\tFS8APND01BNUEY\tFS8APND01DKLOE\tFS8APND01C24PP\tFS8APND01EBWQX\tFS8APND01ELDYW\tFS8APND01B0GCS\tFS8APND01D4QXI\tFS8APND01EMYD9\tFS8APND01EA2SK\tFS8APND01DZOSO\tFS8APND01DHYAZ\tFS8APND01C7UD9\tFS8APND01BTZFV\tFS8APND01CR78R\tFS8APND01B39IE\tFS8APND01ECVC0\tFS8APND01DM3PL\tFS8APND01DELWS\tFS8APND01CIEK8\tFS8APND01D7ZOZ\tFS8APND01CZSAI\tFS8APND01DYOVR\tFS8APND01BX9XY\tFS8APND01DEWJA\tFS8APND01BEKIW\tFS8APND01DCKB9\tFS8APND01EEYIS\tFS8APND01DDKEA\tFS8APND01DSZLO\tFS8APND01C6EBC\tFS8APND01EE15M\tFS8APND01ELO9B\tFS8APND01C58QY\tFS8APND01DONCG\tFS8APND01DVXX2\tFS8APND01BL5YT\tFS8APND01BIL2V\tFS8APND01EBSYQ\tFS8APND01CCX8R\tFS8APND01B2YCJ\tFS8APND01B1JG4\tFS8APND01DJ024\tFS8APND01BIJY0\tFS8APND01CIA4G\tFS8APND01DV74M\tFS8APND01ECAX5\tFS8APND01DC3TZ\tFS8APND01EJVO6\tFS8APND01D4VFG\tFS8APND01DYYYO\tFS8APND01D1EDD\tFS8APND01DQUOT\tFS8APND01A2NSJ\tFS8APND01DDC8I\tFS8APND01BP1T2\tFS8APND01DPY6U\tFS8APND01CIQGV\tFS8APND01BPUT8\tFS8APND01BDNH4\tFS8APND01DOZDN\tFS8APND01DS866\tFS8APND01DGS2J\tFS8APND01EDK32\tFS8APND01EPA0T\tFS8APND01CK3JM\tFS8APND01BKLWW\tFS8APND01DV0BO\tFS8APND01DPNXE\tFS8APND01B7LUA\tFS8APND01BTTE2\tFS8APND01CKO4X\tFS8APND01C4NHX\tFS8APND01DYPQN\tFS8APND01DGGBY
FS8APND01EFK0W:
FS8APND01DCIOO:
FS8APND01CKOMZ:
"""

    def tearDown(self):
        """remove tmp files and stop alarm"""

        # turn off the alarm
        signal.alarm(0)

        if hasattr(self, "result_dir") and self.result_dir:
            try:
                rmtree(self.result_dir)
            except OSError:
            # directory probably not empty (e.g. stale nfs files)
                pass
        try:
            rmdir(self.test_dir)
        except OSError:
            # directory probably not empty, better not remove
            pass

    def test_main(self):
        """Denoiser should always give same result on test data"""

        expected = ">FS8APND01D3TW3 | cluster size: 94 \nCTCCCGTAGGAGTCTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACCCGTCAAAGCCTTGGTAAGCCACTACCCCACCAACAAGCTGATAAGCCGCGAGTCCATCCCCAACCGCCGAAACTTTCCAACCCCCACCATGCAGCAGGAGCTCCTATCCGGTATTAGCCCCAGTTTCCTGAAGTTATCCCAAAGTCAAGGGCAGGTTACTCACGTGTTACTCACCCGTTCGCC\n"

        expected_map = """FS8APND01EWRS4:
FS8APND01DXG45:
FS8APND01D3TW3:\tFS8APND01CSXFN\tFS8APND01DQ8MX\tFS8APND01DY7QW\tFS8APND01B5QNI\tFS8APND01CQ6OG\tFS8APND01C7IGN\tFS8APND01DHSGH\tFS8APND01DJ17E\tFS8APND01CUXOA\tFS8APND01EUTYG\tFS8APND01EKK7T\tFS8APND01D582W\tFS8APND01B5GWU\tFS8APND01D7N2A\tFS8APND01BJGHZ\tFS8APND01D6DYZ\tFS8APND01C6ZIM\tFS8APND01D2X6Y\tFS8APND01BUYCE\tFS8APND01BNUEY\tFS8APND01DKLOE\tFS8APND01C24PP\tFS8APND01EBWQX\tFS8APND01ELDYW\tFS8APND01B0GCS\tFS8APND01D4QXI\tFS8APND01EMYD9\tFS8APND01EA2SK\tFS8APND01DZOSO\tFS8APND01DHYAZ\tFS8APND01C7UD9\tFS8APND01BTZFV\tFS8APND01CR78R\tFS8APND01B39IE\tFS8APND01ECVC0\tFS8APND01DM3PL\tFS8APND01DELWS\tFS8APND01CIEK8\tFS8APND01D7ZOZ\tFS8APND01CZSAI\tFS8APND01DYOVR\tFS8APND01BX9XY\tFS8APND01DEWJA\tFS8APND01BEKIW\tFS8APND01DCKB9\tFS8APND01EEYIS\tFS8APND01DDKEA\tFS8APND01DSZLO\tFS8APND01C6EBC\tFS8APND01EE15M\tFS8APND01ELO9B\tFS8APND01C58QY\tFS8APND01DONCG\tFS8APND01DVXX2\tFS8APND01BL5YT\tFS8APND01BIL2V\tFS8APND01EBSYQ\tFS8APND01CCX8R\tFS8APND01B2YCJ\tFS8APND01B1JG4\tFS8APND01DJ024\tFS8APND01BIJY0\tFS8APND01CIA4G\tFS8APND01DV74M\tFS8APND01ECAX5\tFS8APND01DC3TZ\tFS8APND01EJVO6\tFS8APND01D4VFG\tFS8APND01DYYYO\tFS8APND01D1EDD\tFS8APND01DQUOT\tFS8APND01A2NSJ\tFS8APND01DDC8I\tFS8APND01BP1T2\tFS8APND01DPY6U\tFS8APND01CIQGV\tFS8APND01BPUT8\tFS8APND01BDNH4\tFS8APND01DOZDN\tFS8APND01DS866\tFS8APND01DGS2J\tFS8APND01EDK32\tFS8APND01EPA0T\tFS8APND01CK3JM\tFS8APND01BKLWW\tFS8APND01DV0BO\tFS8APND01DPNXE\tFS8APND01B7LUA\tFS8APND01BTTE2\tFS8APND01CKO4X\tFS8APND01DGGBY\tFS8APND01C4NHX\tFS8APND01DYPQN
FS8APND01BSTVP:
FS8APND01EFK0W:
FS8APND01DCIOO:
FS8APND01CKOMZ:
"""

        command = " ".join(["denoiser.py",
                            "--force", "-o", self.test_dir, "-i",
                            "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        observed = "".join(list(open(self.result_dir + "centroids.fasta")))
        self.assertEqual(observed, expected)

        self.assertEqual(
            len(list(parse_fasta(open(self.result_dir + "singletons.fasta")))),
            6)

        observed = "".join(
            list(open(self.result_dir + "denoiser_mapping.txt")))
        self.assertEqual(observed, expected_map)

    def test_main_with_fasta(self):
        """Denoiser with fasta file should always give same result on test data"""

        command = " ".join(["denoiser.py",
                            "--force", "-o", self.test_dir,
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        observed = "".join(list(open(self.result_dir + "centroids.fasta")))
        self.assertEqual(observed, self.expected)

        observed = "".join(
            list(open(self.result_dir + "denoiser_mapping.txt")))
        self.assertEqual(observed, self.expected_map_string)

    def test_main_with_titanium_error(self):
        """Denoiser with titanium error should always give same result on test data"""

        command = " ".join(["denoiser.py",
                            "--force", "-o", self.test_dir,
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME,
                            "-e", "%s/qiime/support_files/denoiser/Data/Titanium_error_profile.dat" % PROJECT_HOME])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        observed = "".join(list(open(self.result_dir + "centroids.fasta")))
        self.assertEqual(observed, self.expected)

        observed = "".join(
            list(open(self.result_dir + "denoiser_mapping.txt")))
        self.assertEqual(observed, self.expected_titanium_map_string)

    def test_main_on_cluster(self):
        """Denoiser works in a cluster environment"""

        command = " ".join(["denoiser.py",
                            "--force", "-o", self.test_dir, "-c", "-n", "2",
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        observed = "".join(list(open(self.result_dir + "centroids.fasta")))
        self.assertEqual(observed, self.expected)

    def test_main_low_mem(self):
        """Denoiser works using low_memory"""

        command = " ".join(["denoiser.py",
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME,
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-o", self.test_dir, "--low_memory"])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        observed = "".join(list(open(self.result_dir + "centroids.fasta")))
        self.assertEqual(observed, self.expected)

    def test_main_split(self):
        """Denoiser in split mode should always give same result on test data"""

        command = " ".join(["denoiser.py",
                            "-S", "--force",
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME,
                            "-o", self.test_dir])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        for subdir in ["0/", "1/"]:
            observed = "".join(
                list(open(self.result_dir + subdir + "centroids.fasta")))
            self.assertEqual(observed, expected_centroids[subdir])

            observed = "".join(
                list(open(self.result_dir + subdir + "denoiser_mapping.txt")))
            self.assertEqual(observed, expected_map_string[subdir])

    def test_main_split_cluster(self):
        """Denoiser on cluster in split mode should always give same result on test data"""

        command = " ".join(["denoiser.py",
                            "-S", "--force", '-c', '-n 2',
                            "-i", "%s/qiime/support_files/denoiser/TestData/denoiser_test_set.sff.txt" % PROJECT_HOME,
                            "-f", "%s/qiime/support_files/denoiser/TestData/test_set_seqs.fna" % PROJECT_HOME,
                            "-o", self.test_dir])

        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        self.result_dir = self.test_dir

        for subdir in ["0/", "1/"]:
            observed = "".join(
                list(open(self.result_dir + subdir + "centroids.fasta")))
            self.assertEqual(observed, expected_centroids[subdir])

            observed = "".join(
                list(open(self.result_dir + subdir + "denoiser_mapping.txt")))
            self.assertEqual(observed, expected_map_string_on_cluster[subdir])

expected_centroids = {
    '0/':
    ">FS8APND01D3TW3 | cluster size: 49 \nCTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACCCGTCAAAGCCTTGGTAAGCCACTACCCCACCAACAAGCTGATAAGCCGCGAGTCCATCCCCAACCGCCGAAACTTTCCAACCCCCACCATGCAGCAGGAGCTCCTATCCGGTATTAGCCCCAGTTTCCTGAAGTTATCCCAAAGTCAAGGGCAGGTTACTCACGTGTTACTCACCCGTTCGCC\n",
    '1/': ">FS8APND01DQ8MX | cluster size: 45 \nCTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGCCGGCTACCCGTCAAAGCCTTGGTAAGCCACTACCCCACCAACAAGCTGATAAGCCGCGAGTCCATCCCCAACCGCCGAAACTTTCCAACCCCCACCATGCAGCAGGAGCTCCTATCCGGTATTAGCCCCAGTTTCCTGAAGTTATCCCAAAGTCAAGGGCAGGTTACTCACGTGTTACTCACCCGTTCGCC\n"}

expected_map_string = {'0/': """FS8APND01EFK0W:
FS8APND01D3TW3:\tFS8APND01CSXFN\tFS8APND01DY7QW\tFS8APND01C7IGN\tFS8APND01EUTYG\tFS8APND01EKK7T\tFS8APND01B5GWU\tFS8APND01D7N2A\tFS8APND01C6ZIM\tFS8APND01BUYCE\tFS8APND01BNUEY\tFS8APND01EBWQX\tFS8APND01D4QXI\tFS8APND01EMYD9\tFS8APND01EA2SK\tFS8APND01DHYAZ\tFS8APND01BTZFV\tFS8APND01B39IE\tFS8APND01ECVC0\tFS8APND01DELWS\tFS8APND01CIEK8\tFS8APND01DEWJA\tFS8APND01BEKIW\tFS8APND01DCKB9\tFS8APND01C6EBC\tFS8APND01EE15M\tFS8APND01C58QY\tFS8APND01DVXX2\tFS8APND01EBSYQ\tFS8APND01CCX8R\tFS8APND01B1JG4\tFS8APND01BIJY0\tFS8APND01DV74M\tFS8APND01DC3TZ\tFS8APND01D1EDD\tFS8APND01DQUOT\tFS8APND01A2NSJ\tFS8APND01DDC8I\tFS8APND01BP1T2\tFS8APND01DPY6U\tFS8APND01CIQGV\tFS8APND01BDNH4\tFS8APND01DOZDN\tFS8APND01DS866\tFS8APND01EDK32\tFS8APND01BKLWW\tFS8APND01DPNXE\tFS8APND01DGGBY\tFS8APND01C4NHX
""",
                       '1/': """FS8APND01EWRS4:
FS8APND01DQ8MX:\tFS8APND01B5QNI\tFS8APND01CQ6OG\tFS8APND01DHSGH\tFS8APND01DJ17E\tFS8APND01CUXOA\tFS8APND01D582W\tFS8APND01BJGHZ\tFS8APND01D6DYZ\tFS8APND01D2X6Y\tFS8APND01DKLOE\tFS8APND01C24PP\tFS8APND01ELDYW\tFS8APND01B0GCS\tFS8APND01DZOSO\tFS8APND01C7UD9\tFS8APND01CR78R\tFS8APND01DM3PL\tFS8APND01D7ZOZ\tFS8APND01CZSAI\tFS8APND01DYOVR\tFS8APND01BX9XY\tFS8APND01EEYIS\tFS8APND01DDKEA\tFS8APND01DSZLO\tFS8APND01ELO9B\tFS8APND01DONCG\tFS8APND01BL5YT\tFS8APND01BIL2V\tFS8APND01B2YCJ\tFS8APND01DJ024\tFS8APND01CIA4G\tFS8APND01ECAX5\tFS8APND01EJVO6\tFS8APND01D4VFG\tFS8APND01DYYYO\tFS8APND01BPUT8\tFS8APND01DGS2J\tFS8APND01EPA0T\tFS8APND01CK3JM\tFS8APND01DV0BO\tFS8APND01B7LUA\tFS8APND01BTTE2\tFS8APND01CKO4X\tFS8APND01DYPQN
FS8APND01DXG45:
FS8APND01DCIOO:
FS8APND01BSTVP:
FS8APND01CKOMZ:
"""}

# Except for the ordering this is the same as expected_map_string
# If we'd have a more clever way of comparing otu maps, we could use it here...
expected_map_string_on_cluster = {'0/': """FS8APND01EFK0W:
FS8APND01D3TW3:\tFS8APND01CSXFN\tFS8APND01DY7QW\tFS8APND01C7IGN\tFS8APND01EUTYG\tFS8APND01EKK7T\tFS8APND01B5GWU\tFS8APND01D7N2A\tFS8APND01C6ZIM\tFS8APND01BUYCE\tFS8APND01BNUEY\tFS8APND01EBWQX\tFS8APND01D4QXI\tFS8APND01EMYD9\tFS8APND01EA2SK\tFS8APND01DHYAZ\tFS8APND01BTZFV\tFS8APND01B39IE\tFS8APND01ECVC0\tFS8APND01DELWS\tFS8APND01CIEK8\tFS8APND01DEWJA\tFS8APND01BEKIW\tFS8APND01DCKB9\tFS8APND01C6EBC\tFS8APND01EE15M\tFS8APND01C58QY\tFS8APND01DVXX2\tFS8APND01EBSYQ\tFS8APND01CCX8R\tFS8APND01B1JG4\tFS8APND01BIJY0\tFS8APND01DV74M\tFS8APND01DC3TZ\tFS8APND01D1EDD\tFS8APND01DQUOT\tFS8APND01A2NSJ\tFS8APND01DDC8I\tFS8APND01BP1T2\tFS8APND01DPY6U\tFS8APND01CIQGV\tFS8APND01BDNH4\tFS8APND01DOZDN\tFS8APND01DS866\tFS8APND01EDK32\tFS8APND01BKLWW\tFS8APND01DPNXE\tFS8APND01DGGBY\tFS8APND01C4NHX
""",
                                  '1/': """FS8APND01EWRS4:
FS8APND01DXG45:
FS8APND01DCIOO:
FS8APND01BSTVP:
FS8APND01DQ8MX:	FS8APND01B5QNI	FS8APND01CQ6OG	FS8APND01DHSGH	FS8APND01DJ17E	FS8APND01CUXOA	FS8APND01D582W	FS8APND01BJGHZ	FS8APND01D6DYZ	FS8APND01D2X6Y	FS8APND01DKLOE	FS8APND01C24PP	FS8APND01ELDYW	FS8APND01B0GCS	FS8APND01DZOSO	FS8APND01C7UD9	FS8APND01CR78R	FS8APND01DM3PL	FS8APND01D7ZOZ	FS8APND01CZSAI	FS8APND01DYOVR	FS8APND01BX9XY	FS8APND01EEYIS	FS8APND01DDKEA	FS8APND01DSZLO	FS8APND01ELO9B	FS8APND01DONCG	FS8APND01BL5YT	FS8APND01BIL2V	FS8APND01B2YCJ	FS8APND01DJ024	FS8APND01CIA4G	FS8APND01ECAX5	FS8APND01EJVO6	FS8APND01D4VFG	FS8APND01DYYYO	FS8APND01BPUT8	FS8APND01DGS2J	FS8APND01EPA0T	FS8APND01CK3JM	FS8APND01DV0BO	FS8APND01B7LUA	FS8APND01BTTE2	FS8APND01CKO4X	FS8APND01DYPQN
FS8APND01CKOMZ:
"""}

if __name__ == "__main__":
    main()
