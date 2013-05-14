#!/usr/bin/env python

"""test the denoiser settings
"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Release"

from os import access, X_OK, R_OK
from os.path import exists
from subprocess import Popen, PIPE, STDOUT
from cogent.util.misc import app_path
from cogent.util.unit_test import TestCase, main
from qiime.util import get_qiime_scripts_dir, load_qiime_config
from qiime.denoiser.utils import get_flowgram_ali_exe

class DenoiserTests(TestCase):

    def test_cluster_jobs_script(self):
        """cluster_jobs_fp is set to a good value"""

        qiime_config = load_qiime_config()
        submit_script = qiime_config['cluster_jobs_fp']
        
        if (submit_script):
            full_path = app_path(submit_script)
            if full_path:
                submit_script = full_path
            self.assertTrue(exists(submit_script),
                            "cluster_jobs_fp is not set to a valid path in qiime config: %s" % submit_script)
            #check if executable
            self.assertTrue(access(submit_script, X_OK),
                            "cluster_jobs_fp is not executable: %s" % submit_script)
        else:
            #Can't run in parallel, but not a critical error
            pass

    def test_denoiser_min_per_core(self):
        """denoiser_min_per_core is set to a good value"""

        qiime_config = load_qiime_config()
        min_per_core = qiime_config['denoiser_min_per_core']
        if (min_per_core):
            self.assertTrue( int(min_per_core) > 0)
        else:
            self.fail('denoiser_min_per_core not defined in qiime_config.')
 
    def test_denoise_worker(self):
        """denoiser_worker.py is where it belongs and is callable."""

        qiime_config = load_qiime_config()
        PYTHON_BIN = qiime_config['python_exe_fp']

        DENOISE_WORKER = get_qiime_scripts_dir()+"/denoiser_worker.py"

        self.assertTrue(exists(DENOISE_WORKER),
                        "DENOISER_WORKER is not where it's supposed to be: %s"
                        % DENOISE_WORKER)
                       
        #test if its callable and actually works
        command = "%s %s -h"% (PYTHON_BIN, DENOISE_WORKER)
        proc = Popen(command,shell=True,universal_newlines=True,\
                       stdout=PIPE,stderr=STDOUT)

        if (proc.wait() != 0):
            self.fail("Calling %s failed. Check permissions and that it is in fact an executable." \
                          % DENOISE_WORKER)

        result = proc.stdout.read()     
        #check that the help string looks correct
        self.assertTrue(result.startswith("Usage"))

    def test_flowgramAli_bin(self):
        """Check if we have a working FlowgramAligner"""

        ali_fp = get_flowgram_ali_exe()

        self.assertTrue(exists(ali_fp),
                        "The alignment program is not where it's supposed to be: %s"
                        %  ali_fp)
                       
        #test if its callable and actually works
        command = "%s -h" % ali_fp
        proc = Popen(command, shell=True, universal_newlines=True,\
                       stdout=PIPE,stderr=STDOUT)

        if (proc.wait() != 0):
            self.fail("Calling %s failed. Check permissions and that it is in fact an executable."
                      % ali_fp)

        result = proc.stdout.read()     
        #check that the help string looks correct
        self.assertTrue(result.startswith("Usage"))

if __name__ == "__main__":
    main()
