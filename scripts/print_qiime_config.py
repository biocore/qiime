#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder","Dan Knights"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from optparse import make_option
from os import access, X_OK, R_OK, W_OK, getenv
from os.path import isdir, exists, split
from sys import platform, version as python_version, executable
from subprocess import Popen, PIPE, STDOUT

from cogent.util.unit_test import TestCase, main as test_main
from cogent.util.misc import app_path
from cogent.app.util import ApplicationNotFoundError

from qiime.parse import parse_qiime_config_file
from qiime.util import load_qiime_config, get_qiime_project_dir, parse_command_line_parameters
from qiime import __version__ as qiime_lib_version
from cogent import __version__ as pycogent_lib_version
from numpy import __version__ as numpy_lib_version
try:
    from matplotlib import __version__ as matplotlib_lib_version
except ImportError:
    matplotlib_lib_version = "Not installed."
try:
    from pynast import __version__ as pynast_lib_version
except ImportError:
    pynast_lib_version = "Not installed."

try:
    from Denoiser import __version__ as denoiser_version
except ImportError:
    denoiser_version = "Not installed."

script_info = {}
script_info['brief_description']= """Print out the qiime config settings."""
script_info['script_description'] = """A simple scripts that prints out the qiime config settings and does some sanity checks."""
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("Example 1","""Print qiime config settings:""","""print_qiime_config.py"""))
script_info['script_usage'].append(
    ("Example 2","""Print and check qiime config settings for sanity:""",
     """print_qiime_config.py -t"""))

script_info['output_description'] = """This prints the qiime_config to stdout."""
script_info['version'] = __version__
script_info['help_on_no_arguments'] = False

script_info['optional_options'] = [\
    make_option('-t','--test', action='store_true',
                dest='test', default = False,
                help='Test the qiime config for sanity '
                +'[default: %default]')]

class Qiime_config(TestCase):
    
    def setUp(self):
        self.config = load_qiime_config()
    
    def test_python_exe_fp(self):
        """python_exe_fp is set to a working python env"""
        
        python = self.config["python_exe_fp"]
        command = "%s --version" % python
        proc = Popen(command,shell=True,universal_newlines=True,\
                         stdout=PIPE,stderr=STDOUT)
        #Check if callable
        if proc.wait() !=0:
            self.fail("Calling python failed. Check you python_exe_fp:%s" %python)
        
        # Does it give its version string?
        out_string = proc.stdout.read()
        if not out_string:
            self.fail("Something is wrong with your python\n." \
                          +" Check you python_exe_fp:%s" %python)

        #get the version number out of "python 2.6.2\n"
        version = out_string.rstrip().split(" ")[1].split(".")
        version = tuple(map(int, version))
        self.assertTrue(version >= (2, 6, 0),
                        "Qiime requires at least python version 2.6.0. "+\
                        "You are running " + out_string)
    
    def test_cluster_jobs_fp(self):
        """cluster_jobs_fp is set to a valid path and is executable"""       
        
        test_qiime_config_variable("cluster_jobs_fp", self.config, self, X_OK,)
   
    def test_blastmat_dir(self):
        """blastmat_dir is set to a valid path."""
        
        test_qiime_config_variable("blastmat_dir", self.config, self)
        
    def test_blastall_fp(self):
        """blastall_fp is set to a valid path"""
        
        blastall = self.config["blastall_fp"]
        if not self.config["blastall_fp"].startswith("/"):
            #path is relative, figure out absolute path
            blast_all = app_path(blastall)
            if not blast_all:
                raise ApplicationNotFoundError("blastall_fp set to %s, but is not in your PATH. Either use an absolute path to or put it in your PATH." % blastall)
            self.config["blastall_fp"] = blast_all

        test_qiime_config_variable("blastall_fp", self.config, self, X_OK)

    def test_rdp_classifier_fp(self):
        """rdp_classifier is set to a valid path"""
            
        test_qiime_config_variable("rdp_classifier_fp", self.config, self)
        
    def test_pynast_template_alignment_fp(self):
        """pynast_template_alignment is set to a valid path"""
            
        test_qiime_config_variable("pynast_template_alignment_fp",
                                   self.config, self)
            
    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb is set to a valid path"""
            
        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)
    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb is set to a valid path"""
        
        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)
        
    def test_template_alignment_lanemask_fp(self):
        """template_alignment_lanemask is set to a valid path"""
            
        test_qiime_config_variable("template_alignment_lanemask_fp",
                                   self.config, self)
    
    def test_qiime_scripts_dir(self):
        """qiime_scripts_dir is set to a valid path"""

        scripts_dir = self.config["qiime_scripts_dir"]
        
        if scripts_dir:
            self.assertTrue(exists(scripts_dir),
                            "qiime_scripts_dir does not exist: %s" % scripts_dir)
            self.assertTrue(isdir(scripts_dir),
                            "qiime_scripts_dir is not a directory: %s" % scripts_dir)
        else:
            pass
            #self.fail("scripts_dir is not set.")
                        
    def test_working_dir(self):
        """working_dir is set to a valid path"""

        working_dir = self.config["working_dir"]
        
        if working_dir:
            self.assertTrue(exists(working_dir), 
                            "working dir does not exist: %s" % working_dir)
            self.assertTrue(isdir(working_dir),
                            "working_dir is not a directory: %s" % working_dir)        
            self.assertTrue(access(working_dir, W_OK),
                            "working_dir not writable: %s" % working_dir)
        else:
            pass
            #self.fail("working_dir is not set.")
        
    #we are not testing these values from the qiime_config:
    # jobs_to_start   1
    # seconds_to_sleep        60


    def test_for_obsolete_values(self):
        """local qiime_config has no extra params"""
        
        qiime_project_dir = get_qiime_project_dir()
        orig_config = parse_qiime_config_file(open(qiime_project_dir +
                                             '/qiime/support_files/qiime_config'))
        
        #check the env qiime_config
        qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
        if qiime_config_env_filepath:
            qiime_config_via_env= parse_qiime_config_file(open(qiime_config_env_filepath))
            extra_vals = []
            for key in qiime_config_via_env:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The qiime_config file set via QIIME_CONFIG_FP"+
                          "enviroment variable contains obsolete parameters:\n"+
                          ", ".join(extra_vals))
        # check the qiime_config in $HOME/.qiime_config
        home_dir = getenv('HOME')        
        if (exists(home_dir+"/.qiime_config")):
            qiime_config_home = parse_qiime_config_file(open(home_dir+"/.qiime_config"))
            extra_vals = []
            for key in qiime_config_home:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The .qiime_config in your HOME contains obsolete "+
                          "parameters:\n" + ", ".join(extra_vals))

    def test_chimeraSlayer_install(self):
        """no obvious problems with ChimeraSlayer install """

        #The ChimerSalyer app requires that all its components are installed
        # relative to the main program ChimeraSlayer.pl.
        # We therefore check that at least one the files is there.
        # However, if the directory structure of ChimeraSlayer changes, this test will most
        # likely fail as well and need to be updated.
        # Tested with the version of microbiomeutil_2010-04-29

        chim_slay = app_path("ChimeraSlayer.pl")
        self.assertTrue(chim_slay,"ChimeraSlayer was not found in your $PATH")
        dir, app_name = split(chim_slay)
        self.assertTrue(exists(dir+"/ChimeraParentSelector/chimeraParentSelector.pl"),
         "ChimeraSlayer depends on external files in directoryies relative to its "
         "install directory. Thesedo not appear to be present.")

def test_qiime_config_variable(variable, qiime_config, test,
                               access_var=R_OK, fail_on_missing=False):
    """test if a variable is set and set to a readable path."""
    
    fp = qiime_config[variable]
    
    if not fp:
        if fail_on_missing:
            test.fail("%s not set."%variable)
        else:
            # non-essential file, so do not fail
            return
    #test if file exists    
    test.assertTrue(exists(fp), "%s set to an invalid file path: %s" %\
                        (variable,fp))

    modes = {R_OK:"readable",
             W_OK:"writable",
             X_OK:"executable"}
    #test if file readable    
    test.assertTrue(access(fp, access_var),
                    "%s is not %s: %s" % (variable, modes[access_var], fp))

if __name__ == "__main__":
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    qiime_config = load_qiime_config()

    system_info = [
     ("Platform", platform),
     ("Python version",python_version.replace('\n', ' ')),
     ("Python executable",executable)]
    max_len =  max([len(e[0]) for e in system_info])
    print "\nSystem information"
    print  "==================" 
    for v in system_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])

    version_info = [
     ("PyCogent version", pycogent_lib_version),
     ("NumPy version", numpy_lib_version),
     ("matplotlib version", matplotlib_lib_version),
     ("QIIME library version", qiime_lib_version),
     ("QIIME script version", __version__),
     ("PyNAST version (if installed)", pynast_lib_version),
     ("Denoiser version (if installed)", denoiser_version)]
    max_len =  max([len(e[0]) for e in version_info])
    print "\nDependency versions"
    print  "===================" 
    for v in version_info:
        print "%*s:\t%s" % (max_len,v[0],v[1])
    
    print "\nQIIME config values"
    print  "==================="    
    max_len =  max([len(key) for key in qiime_config])
    for key,value in  qiime_config.items():
        print "%*s:\t%s"%(max_len,key,value)

    #run the Testcase.main function to do the tests
    # need to mess with the arg string, otherwise TestCase complains
    if (opts.test):
        print "\n\nrunning checks:\n"
        test_main(argv=["","-v"])
