#!/usr/bin/env python
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Dan Knights", "Antonio Gonzalez Pena",
               "Justin Kuczynski", "Jai Ram Rideout", "Greg Caporaso",
               "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import re
from os import access, X_OK, R_OK, W_OK, getenv, environ, remove, devnull
from os.path import isdir, exists, split, join
from sys import platform, version as python_version, executable, stdout
from unittest import TestLoader, TextTestRunner, TestCase
from shutil import rmtree
from subprocess import Popen, PIPE, STDOUT
from optparse import SUPPRESS_HELP

core_dependency_missing_msg = "See the QIIME Installation Guide: http://qiime.org/install/install.html"

try:
    from numpy import __version__ as numpy_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from scipy import __version__ as scipy_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from tempfile import mkdtemp, remove_files
    from skbio.app.util import ApplicationNotFoundError, ApplicationError
    from cogent import __version__ as pycogent_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from qiime.parse import parse_qiime_config_file
    from qiime.util import (load_qiime_config,
                            get_qiime_project_dir,
                            parse_command_line_parameters,
                            get_qiime_library_version,
                            get_qiime_scripts_dir,
                            get_rdp_jarpath,
                            get_java_version,
                            get_pynast_version,
                            make_option,
                            qiime_system_call,
                            get_qiime_temp_dir)
    from qiime.denoiser.utils import check_flowgram_ali_exe
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from biom import __version__ as biom_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from qcli import __version__ as qcli_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from pyqi import __version__ as pyqi_lib_version
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from skbio import __version__ as skbio_lib_version
    from skbio.app.util import which
except ImportError as e:
    raise ImportError("%s\n%s" % (e, core_dependency_missing_msg))

try:
    from matplotlib import __version__ as matplotlib_lib_version
except ImportError:
    matplotlib_lib_version = "Not installed."

try:
    from emperor import __version__ as emperor_lib_version
except ImportError:
    emperor_lib_version = "Not installed."


pynast_lib_version = get_pynast_version()
if pynast_lib_version is None:
    pynast_lib_version = "Not installed."

script_info = {}
script_info['brief_description'] = """Print out the qiime config settings."""
script_info[
    'script_description'] = """A simple scripts that prints out the qiime config settings and does some sanity checks."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example 1", """Print qiime config settings:""", """print_qiime_config.py"""))
script_info['script_usage'].append(
    ("Example 2", """Print and check qiime config settings for sanity:""",
     """print_qiime_config.py -t"""))

script_info[
    'output_description'] = """This prints the qiime_config to stdout."""
script_info['version'] = __version__
script_info['help_on_no_arguments'] = False
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-t', '--test',
                action='store_true',
                default=False,
                help='Test the QIIME install and configuration [default: %default]'),
    make_option('-b',
                '--qiime_base_install',
                action='store_true',
                default=False,
                help='If passed, report only on dependencies required for the QIIME base install [default: %default]'),
    make_option('--haiku',
                action='store_true',
                default=False,
                help=SUPPRESS_HELP)
]


class QIIMEConfig(TestCase):

    def setUp(self):
        self.config = load_qiime_config()

    def test_cluster_jobs_fp(self):
        """cluster_jobs_fp is set to a valid path and is executable"""
        fp = self.config["cluster_jobs_fp"]

        if fp:
            full_path = which(fp)
            if full_path:
                fp = full_path

            # test if file exists or is in $PATH
            self.assertTrue(exists(fp),
                            "cluster_jobs_fp set to an invalid file path or is not in $PATH: %s" % fp)

            modes = {R_OK: "readable",
                     W_OK: "writable",
                     X_OK: "executable"}
            # test if file readable
            self.assertTrue(access(fp, X_OK),
                            "cluster_jobs_fp is not %s: %s" % (modes[X_OK], fp))

    def test_blastmat_dir(self):
        """blastmat_dir is set to a valid path."""

        test_qiime_config_variable("blastmat_dir", self.config, self)

    def test_pynast_template_alignment_fp(self):
        """pynast_template_alignment, if set, is set to a valid path"""

        test_qiime_config_variable("pynast_template_alignment_fp",
                                   self.config, self)

    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb, if set, is set to a valid path"""

        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)

    def test_pynast_template_alignment_blastdb_fp(self):
        """pynast_template_alignment_blastdb, if set, is set to a valid path"""

        test_qiime_config_variable("pynast_template_alignment_blastdb_fp",
                                   self.config, self)

    def test_template_alignment_lanemask_fp(self):
        """template_alignment_lanemask, if set, is set to a valid path"""

        test_qiime_config_variable("template_alignment_lanemask_fp",
                                   self.config, self)

    def test_get_qiime_scripts_dir(self):
        """Test that we can find the directory containing QIIME scripts."""
        # get_qiime_scripts_dir will raise an error if it can't find a scripts
        # directory.
        scripts_dir = get_qiime_scripts_dir()
        self.assertTrue(isdir(scripts_dir), "The QIIME scripts directory does "
                        "not exist: %s" % scripts_dir)

    def test_temp_dir(self):
        """temp_dir is set to a valid path"""
        temp_dir = get_qiime_temp_dir()

        self.assertTrue(exists(temp_dir),
                        "temp_dir does not exist: %s" % temp_dir)
        self.assertTrue(isdir(temp_dir),
                        "temp_dir is not a directory: %s" % temp_dir)
        self.assertTrue(access(temp_dir, W_OK),
                        "temp_dir is not writable: %s" % temp_dir)

    # we are not testing these values from the qiime_config:
    # jobs_to_start   1
    # seconds_to_sleep        60

    def test_for_unrecognized_values(self):
        """qiime_config has no extra values"""
        error_msg_fragment = (" contains unrecognized values:\n%s\nYou can "
                              "safely remove these values from your QIIME "
                              "config file as they will be ignored by QIIME.")

        qiime_project_dir = get_qiime_project_dir()
        orig_config = parse_qiime_config_file(open(qiime_project_dir +
                                                   '/qiime/support_files/qiime_config'))

        # check the env qiime_config
        qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
        if qiime_config_env_filepath:
            qiime_config_via_env = parse_qiime_config_file(
                open(qiime_config_env_filepath))
            extra_vals = []
            for key in qiime_config_via_env:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The QIIME config file set via the QIIME_CONFIG_FP "
                          "environment variable" +
                          error_msg_fragment % ", ".join(extra_vals))
        # check the qiime_config in $HOME/.qiime_config
        home_dir = getenv('HOME')
        if (exists(home_dir + "/.qiime_config")):
            qiime_config_home = parse_qiime_config_file(
                open(home_dir + "/.qiime_config"))
            extra_vals = []
            for key in qiime_config_home:
                if key not in orig_config:
                    extra_vals.append(key)
            if extra_vals:
                self.fail("The .qiime_config in your HOME" +
                          error_msg_fragment % ", ".join(extra_vals))


class QIIMEDependencyBase(QIIMEConfig):

    def test_uclust_supported_version(self):
        """uclust is in path and version is supported """
        acceptable_version = (1, 2, 22)
        self.assertTrue(which('uclust'),
                        "uclust not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = 'uclust --version'
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split('v')[-1].strip('q')
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported uclust version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_python_supported_version(self):
        """python is in path and version is supported """
        min_acceptable_version = (2, 7, 0)
        min_unacceptable_version = (3, 0, 0)

        command = 'python --version'
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()

        version_str_matches = re.findall('Python\s+(\S+)\s*', stdout.strip())
        self.assertEqual(len(version_str_matches), 1,
                         "Could not determine the Python version in '%s'." %
                         stdout)
        version_string = version_str_matches[0]

        try:
            if version_string[-1] == '+':
                version_string = version_string[:-1]
            version = tuple(map(int, version_string.split('.')))
            if len(version) == 2:
                version = (version[0], version[1], 0)
            pass_test = (version >= min_acceptable_version and
                         version < min_unacceptable_version)
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported Python version. Must be >= %s and < %s, "
                        "but running %s."
                        % ('.'.join(map(str, min_acceptable_version)),
                           '.'.join(map(str, min_unacceptable_version)),
                           version_string))

    def test_numpy_supported_version(self):
        """numpy version is supported """
        min_acceptable_version = (1, 5, 1)
        max_acceptable_version = (1, 7, 1)
        try:
            from numpy import __version__ as numpy_lib_version
            version = tuple(map(int, numpy_lib_version.split('.')))
            pass_test = (version >= min_acceptable_version and
                         version <= max_acceptable_version)
            version_string = str(numpy_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(pass_test,
                        "Unsupported numpy version. Must be >= %s and <= %s , but running %s."
                        % ('.'.join(map(str, min_acceptable_version)),
                           '.'.join(map(str, max_acceptable_version)),
                            version_string))

    def test_scipy_supported_version(self):
        """scipy version is supported """
        min_acceptable_version = (0, 13, 0)
        try:
            from scipy import __version__ as scipy_lib_version
            version = tuple(map(int, scipy_lib_version.split('.')))
            pass_test = version >= min_acceptable_version
            version_string = str(scipy_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(
            pass_test,
            "Unsupported scipy version. Must be >= %s, but running %s." %
            ('.'.join(map(str, min_acceptable_version)), version_string))

    def test_matplotlib_suported_version(self):
        """matplotlib version is supported """
        min_acceptable_version = (1, 1, 0)
        max_acceptable_version = (1, 3, 1)
        try:
            from matplotlib import __version__ as matplotlib_lib_version
            version = matplotlib_lib_version.split('.')
            if version[-1].endswith('rc'):
                version[-1] = version[-1][:-2]
            version = tuple(map(int, version))
            pass_test = (version >= min_acceptable_version and
                         version <= max_acceptable_version)
            version_string = str(matplotlib_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"
        self.assertTrue(pass_test,
                        "Unsupported matplotlib version. Must be >= %s and <= %s , but running %s."
                        % ('.'.join(map(str, min_acceptable_version)),
                           '.'.join(map(str, max_acceptable_version)), version_string))

    def test_pynast_suported_version(self):
        """pynast version is supported """
        min_acceptable_version = (1, 2)
        max_acceptable_version = (1, 2, 2)
        try:
            from pynast import __version__ as pynast_lib_version
            version = pynast_lib_version.split('.')
            if version[-1][-4:] == '-dev':
                version[-1] = version[-1][:-4]
            version = tuple(map(int, version))
            pass_test = (version >= min_acceptable_version and
                         version <= max_acceptable_version)
            version_string = str(pynast_lib_version)
        except ImportError:
            pass_test = False
            version_string = "Not installed"

        min_version_str = '.'.join(map(str, min_acceptable_version))
        max_version_str = '.'.join(map(str, max_acceptable_version))
        error_msg = ("Unsupported pynast version. Must be >= %s and <= %s, "
                     "but running %s." % (min_version_str, max_version_str,
                                          version_string))
        self.assertTrue(pass_test, error_msg)

    def test_FastTree_supported_version(self):
        """FastTree is in path and version is supported """
        acceptable_version = (2, 1, 3)
        self.assertTrue(which('FastTree'),
                        "FastTree not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")

        # If FastTree is run interactively, it outputs the following line:
        #     Usage for FastTree version 2.1.3 SSE3:
        #
        # If run non-interactively:
        #     FastTree Version 2.1.3 SSE3
        command = "FastTree 2>&1 > %s | grep -i version" % devnull
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read().strip()

        version_str_matches = re.findall('ersion\s+(\S+)\s+', stdout)
        self.assertEqual(len(version_str_matches), 1,
                         "Could not find FastTree version info in usage text "
                         "'%s'." % stdout)

        version_str = version_str_matches[0]

        try:
            version = tuple(map(int, version_str.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False

        acceptable_version_str = '.'.join(map(str, acceptable_version))
        self.assertTrue(pass_test,
                        "Unsupported FastTree version. %s is required, but "
                        "running %s." % (acceptable_version_str, version_str))


class QIIMEDependencyFull(QIIMEDependencyBase):

    def test_ampliconnoise_install(self):
        """ AmpliconNoise install looks sane."""
        url = "http://qiime.org/install/install.html#ampliconnoise-install-notes"

        pyro_lookup_file = getenv('PYRO_LOOKUP_FILE')
        self.assertTrue(pyro_lookup_file is not None,
                        "$PYRO_LOOKUP_FILE variable is not set. See %s for help." % url)
        self.assertTrue(exists(pyro_lookup_file),
                        "$PYRO_LOOKUP_FILE variable is not set to an existing filepath.")

        seq_lookup_file = getenv('SEQ_LOOKUP_FILE')
        self.assertTrue(seq_lookup_file is not None,
                        "$SEQ_LOOKUP_FILE variable is not set. See %s for help." % url)
        self.assertTrue(exists(seq_lookup_file),
                        "$SEQ_LOOKUP_FILE variable is not set to an existing filepath.")

        self.assertTrue(which("SplitKeys.pl"),
                        "Couldn't find SplitKeys.pl. " +
                        "Perhaps AmpliconNoise Scripts directory isn't in $PATH?" +
                        " See %s for help." % url)

        self.assertTrue(which("FCluster"),
                        "Couldn't find FCluster. " +
                        "Perhaps the AmpliconNoise bin directory isn't in $PATH?" +
                        " See %s for help." % url)

        self.assertTrue(which("Perseus"),
                        "Couldn't find Perseus. " +
                        "Perhaps the AmpliconNoise bin directory isn't in $PATH?" +
                        " See %s for help." % url)

    def test_sourcetracker_installed(self):
        """sourcetracker is installed"""

        sourcetracker_path = getenv('SOURCETRACKER_PATH')
        self.assertNotEqual(sourcetracker_path, None,
                            ("SOURCETRACKER_PATH is not set. This is "
                             "only important if you plan to use SourceTracker."))
        self.assertTrue(exists(sourcetracker_path),
                        "SOURCETRACKER_PATH is not set to a valid path: %s" %
                        sourcetracker_path)

    def test_chimeraSlayer_install(self):
        """no obvious problems with ChimeraSlayer install """

        # The ChimerSalyer app requires that all its components are installed
        # relative to the main program ChimeraSlayer.pl.
        # We therefore check that at least one the files is there.
        # However, if the directory structure of ChimeraSlayer changes, this test will most
        # likely fail as well and need to be updated.
        # Tested with the version of microbiomeutil_2010-04-29

        chim_slay = which("ChimeraSlayer.pl")
        self.assertTrue(chim_slay, "ChimeraSlayer was not found in your $PATH")
        dir, app_name = split(chim_slay)
        self.assertTrue(
            exists(dir + "/ChimeraParentSelector/chimeraParentSelector.pl"),
            "ChimeraSlayer depends on external files in directoryies relative to its "
            "install directory. These do not appear to be present.")

    def test_blastall_fp(self):
        """blastall_fp is set to a valid path"""

        blastall = self.config["blastall_fp"]
        if not self.config["blastall_fp"].startswith("/"):
            # path is relative, figure out absolute path
            blast_all = which(blastall)
            if not blast_all:
                raise ApplicationNotFoundError(
                    "blastall_fp set to %s, but is not in your PATH. Either use an absolute path to or put it in your PATH." %
                    blastall)
            self.config["blastall_fp"] = blast_all

        test_qiime_config_variable("blastall_fp", self.config, self, X_OK)

    def test_blast_supported_version(self):
        """blast is in path and version is supported """
        acceptable_version = (2, 2, 22)
        self.assertTrue(which('blastall'),
                        "blast not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = 'blastall | grep blastall'
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[1].strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported blast version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_cdbtools_supported_version(self):
        """cdbtools is in path and version is supported """
        acceptable_version = (0, 99)
        self.assertTrue(which('cdbfasta'),
                        "cdbtools not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "cdbfasta -v"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported cdbtools version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_INFERNAL_supported_version(self):
        """INFERNAL is in path and version is supported """
        acceptable_version = (1, 0, 2)
        self.assertTrue(which('cmbuild'),
                        "Infernal not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "cmbuild -h | grep INF"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported INFERNAL version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_muscle_supported_version(self):
        """muscle is in path and version is supported """
        acceptable_version = (3, 8, 31)
        self.assertTrue(which('muscle'),
                        "muscle not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "muscle -version"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[1].strip('v')
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported muscle version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_mothur_supported_version(self):
        """mothur is in path and version is supported """
        acceptable_version = (1, 25, 0)
        self.assertTrue(which('mothur'),
                        "mothur not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        # mothur creates a log file in cwd, so create a tmp and cd there first
        log_file = join(get_qiime_temp_dir(), 'mothur.log')
        command = "mothur \"#set.logfile(name=%s)\" | grep '^mothur v'" % log_file
        stdout, stderr, exit_Status = qiime_system_call(command)

        # remove log file
        remove_files([log_file], error_on_missing=False)

        version_string = stdout.strip().split(' ')[1].strip('v.')
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported mothur version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_denoiser_supported_version(self):
        """denoiser aligner is ready to use """

        pass_test = True
        try:
            check_flowgram_ali_exe()
        except (ApplicationNotFoundError, ApplicationError):
            pass_test = False

        self.assertTrue(pass_test,
                        "Denoiser flowgram aligner not found or not "
                        "executable. This may or may not be a problem "
                        "depending on which components of QIIME you plan to "
                        "use.")

    def test_raxmlHPC_supported_version(self):
        """raxmlHPC is in path and version is supported """
        acceptable_version = [(7, 3, 0), (7, 3, 0)]
        self.assertTrue(which('raxmlHPC'),
                        "raxmlHPC not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "raxmlHPC -v | grep version"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[4].strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version in acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported raxmlHPC version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_clearcut_supported_version(self):
        """clearcut is in path and version is supported """
        acceptable_version = (1, 0, 9)
        self.assertTrue(which('clearcut'),
                        "clearcut not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "clearcut -V"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip().split(' ')[2].strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version == acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported clearcut version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_cdhit_supported_version(self):
        """cd-hit is in path and version is supported """
        self.assertTrue(which('cd-hit'),
                        "cd-hit not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        # cd-hit does not have a version print in their program

    def test_rtax_supported_version(self):
        """rtax is in path and version is supported """
        acceptable_version = [(0, 984)]
        self.assertTrue(which('rtax'),
                        "rtax not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "rtax 2>&1 > %s | grep Version | awk '{print $2}'" % devnull
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version in acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported rtax version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_pplacer_supported_version(self):
        """pplacer is in path and version is supported """
        acceptable_version = [(1, 1), (1, 1)]
        self.assertTrue(which('pplacer'),
                        "pplacer not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "pplacer --version"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip()[1:4]
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version in acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported pplacer version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_ParsInsert_supported_version(self):
        """ParsInsert is in path and version is supported """
        acceptable_version = ["1.04"]
        self.assertTrue(which('ParsInsert'),
                        "ParsInsert not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "ParsInsert -v | grep App | awk '{print $3}'"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()

        # remove log file generated
        remove_files(['ParsInsert.log'], error_on_missing=False)

        version_string = stdout.strip()
        try:
            pass_test = version_string in acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported ParsInsert version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_usearch_supported_version(self):
        """usearch is in path and version is supported """
        acceptable_version = [(5, 2, 236), (5, 2, 236)]
        self.assertTrue(which('usearch'),
                        "usearch not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "usearch --version"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.split('v')[1]
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = version in acceptable_version
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported usearch version. %s is required, but running %s."
                        % ('.'.join(map(str, acceptable_version)), version_string))

    def test_R_supported_version(self):
        """R is in path and version is supported """
        minimum_version = (2, 12, 0)
        self.assertTrue(which('R'),
                        "R not found. This may or may not be a problem depending on " +
                        "which components of QIIME you plan to use.")
        command = "R --version | grep 'R version' | awk '{print $3}'"
        proc = Popen(command, shell=True, universal_newlines=True,
                     stdout=PIPE, stderr=STDOUT)
        stdout = proc.stdout.read()
        version_string = stdout.strip()
        try:
            version = tuple(map(int, version_string.split('.')))
            pass_test = False
            if version[0] == minimum_version[0]:
                if version[1] == minimum_version[1]:
                    if version[2] >= minimum_version[2]:
                        pass_test = True
                elif version[1] > minimum_version[1]:
                    pass_test = True
            elif version[0] > minimum_version[0]:
                pass_test = True
        except ValueError:
            pass_test = False
            version_string = stdout
        self.assertTrue(pass_test,
                        "Unsupported R version. %s or greater is required, but running %s."
                        % ('.'.join(map(str, minimum_version)), version_string))

    def test_gdata_install(self):
        """gdata is installed"""
        # We currently can't programmatically find the version of gdata. An
        # issue has been created alerting the gdata devs.
        pass_test = True
        try:
            import gdata
        except ImportError:
            pass_test = False
        self.assertTrue(pass_test, "gdata is not installed.")


def test_qiime_config_variable(variable, qiime_config, test,
                               access_var=R_OK, fail_on_missing=False):
    """test if a variable is set and set to a readable path."""

    fp = qiime_config[variable]

    if not fp:
        if fail_on_missing:
            test.fail("%s not set." % variable)
        else:
            # non-essential file, so do not fail
            return
    # test if file exists
    test.assertTrue(exists(fp), "%s set to an invalid file path: %s" %
                    (variable, fp))

    modes = {R_OK: "readable",
             W_OK: "writable",
             X_OK: "executable"}
    # test if file readable
    test.assertTrue(access(fp, access_var),
                    "%s is not %s: %s" % (variable, modes[access_var], fp))


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.haiku:
        print "QIIME provides insight\nmicrobial in nature\nto ecology"
        exit(0)

    qiime_config = load_qiime_config()
    test = opts.test
    qiime_base_install = opts.qiime_base_install

    rdp_jarpath = get_rdp_jarpath()
    if rdp_jarpath is None:
        rdp_version = "Not installed."
    else:
        rdp_version = split(rdp_jarpath)[1]

    java_version = get_java_version()
    if java_version is None:
        java_version = "Not installed."

    system_info = [
        ("Platform", platform),
        ("Python version", python_version.replace('\n', ' ')),
        ("Python executable", executable)]
    max_len = max([len(e[0]) for e in system_info])
    print "\nSystem information"
    print "=================="
    for v in system_info:
        print "%*s:\t%s" % (max_len, v[0], v[1])

    version_info = [
        ("PyCogent version", pycogent_lib_version),
        ("NumPy version", numpy_lib_version),
        ("SciPy version", scipy_lib_version),
        ("matplotlib version", matplotlib_lib_version),
        ("biom-format version", biom_lib_version),
        ("qcli version", qcli_lib_version),
        ("pyqi version", pyqi_lib_version),
        ("scikit-bio version", skbio_lib_version),
        ("QIIME library version", get_qiime_library_version()),
        ("QIIME script version", __version__),
        ("PyNAST version (if installed)", pynast_lib_version),
        ("Emperor version", emperor_lib_version)]
    if not qiime_base_install:
        version_info += [
            ("RDP Classifier version (if installed)", rdp_version),
            ("Java version (if installed)", java_version)]

    max_len = max([len(e[0]) for e in version_info])
    print "\nDependency versions"
    print "==================="
    for v in version_info:
        print "%*s:\t%s" % (max_len, v[0], v[1])

    print "\nQIIME config values"
    print "==================="
    max_len = max([len(key) for key in qiime_config])
    for key, value in qiime_config.items():
        print "%*s:\t%s" % (max_len, key, value)

    if test:
        if qiime_base_install:
            suite = TestLoader().loadTestsFromTestCase(QIIMEDependencyBase)
        else:
            suite = TestLoader().loadTestsFromTestCase(QIIMEDependencyFull)
        if opts.verbose:
            verbosity = 2
        else:
            verbosity = 1
        TextTestRunner(stream=stdout, verbosity=verbosity).run(suite)

if __name__ == "__main__":
    main()
