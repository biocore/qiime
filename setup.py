#!/usr/bin/env python
# File created on 17 Feb 2010
from __future__ import division
from setuptools import setup
from stat import S_IEXEC
from os import (chdir, getcwd, listdir, chmod, walk, rename, remove, chmod,
                stat, devnull, environ)
from os.path import join, abspath
from sys import platform, argv, exc_info
from subprocess import call, Popen, PIPE
from glob import glob
from urllib import FancyURLopener
from tempfile import mkdtemp
from shutil import rmtree, copy
import re

__author__ = "QIIME development team"
__copyright__ = "Copyright (c) 2011--, %s" % __author__
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout",
               "Yoshiki Vazquez Baeza", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Greg Caporaso"
__email__ = "qiime.help@gmail.com"


long_description = """QIIME: Quantitative Insights Into Microbial Ecology
http://www.qiime.org

QIIME Allows Integration and Analysis of High-Throughput Community Sequencing Data
J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K. Goodrich, Jeffrey I. Gordon, Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight.
Nature Methods, 2010.
"""

# don't build any of the non-Python dependencies if the following modes are
# invoked
build_stack = all([e not in argv for e in 'egg_info', 'sdist', 'register'])

def build_denoiser():
    """ Build the denoiser code binary """
    if not app_available('ghc'):
        print ("GHC not installed, so cannot build the denoiser binary "
               "'FlowgramAli_4frame'.")
        return

    cwd = getcwd()
    denoiser_dir = join(cwd, 'qiime/support_files/denoiser/FlowgramAlignment')

    try:
        chdir(denoiser_dir)

        stdout, stderr, return_value = system_call('make clean')

        if return_value != 0:
            print ("Unable to clean denoiser build directory.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        stdout, stderr, return_value = system_call('make')

        if return_value != 0:
            print ("Unable to build denoiser.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        print "Denoiser built."
    finally:
        chdir(cwd)

# heavily based on lib.util.download_file from github.com/qiime/qiime-deploy


class URLOpener(FancyURLopener):
    def http_error_default(self, url, fp, errcode, errmsg, headers):
        raise IOError(
            'Could not download %s\nPlease ensure the URL is valid and that '
            'you have an active Internet connection.' % url)

# heavily based on lib.util.download_file from github.com/qiime/qiime-deploy


def download_file(URL, dest_dir, local_file, num_retries=4):
    """General file downloader

    Inputs:
    URL: string to download the file from
    dest_dir: directory where you want to download the file
    local_file: output filename of the download
    num_retries: number of times the function will try to download the file

    Output:
    return_code: exit status for the download 0 = success, 1 = fail
    """
    url_opener = URLOpener()
    localFP = join(dest_dir, local_file)
    tmpDownloadFP = '%s.part' % localFP

    return_code = 1
    while num_retries > 0:
        try:
            tmpLocalFP, headers = url_opener.retrieve(URL, tmpDownloadFP)
            rename(tmpDownloadFP, localFP)
            return_code = 0
        except IOError as msg:
            if num_retries == 1:
                print 'Download of %s failed.' % URL
            else:
                print 'Download failed. Trying again... %d tries remain.' % (
                    num_retries - 1)
            num_retries -= 1
        else:
            num_retries = 0
            print '%s downloaded successfully.' % local_file
    return return_code


def build_FastTree():
    """Download and build FastTree then copy it to the scripts directory"""
    if not app_available('gcc'):
        print "GCC not installed, so cannot build FastTree."
        return

    if download_file('http://www.microbesonline.org/fasttree/FastTree-2.1.3.c',
                     'scripts/', 'FastTree.c'):
        print 'Could not download FastTree, not installing it.'
        return

    cwd = getcwd()
    scripts = join(cwd, 'scripts')

    try:
        chdir(scripts)

        # as suggested by the docs in FastTree.c
        stdout, stderr, return_value = system_call(
            'gcc -Wall -O3 -finline-functions -funroll-loops -o FastTree '
            'FastTree.c -lm')

        if return_value != 0:
            print ("Unable to build FastTree.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        print "FastTree built."
    finally:
        # remove the source
        remove('FastTree.c')
        chdir(cwd)


def build_SortMeRNA():
    """Download and build SortMeRNA then copy it to the scripts directory"""
    # SortMeRNA's configure script doesn't correctly guess the C/C++ compilers
    # to use on OS X. Try to figure that out here.
    if platform.lower() in ['darwin', 'macos']:
        cxx = 'clang++'
        cc = 'clang'
    elif platform.lower() in ['linux', 'linux2']:
        cxx = 'g++'
        cc = 'gcc'
    else:
        print "Unknown or unsupported platform %r, so cannot build SortMeRNA." % platform
        return

    for compiler in cxx, cc:
        if not app_available(compiler):
            print "%r not installed, so cannot build SortMeRNA." % compiler
            return

    cwd = getcwd()
    scripts = join(cwd, 'scripts')
    cxx_old = environ.get('CXX', '')
    cc_old = environ.get('CC', '')

    try:
        tempdir = mkdtemp()
        if download_file('ftp://ftp.microbio.me/pub/sortmerna-2.0-no-db.tar.gz',
                         tempdir, 'sortmerna-2.0-no-db.tar.gz'):
            print "Could not download SortMeRNA, so cannot install it."
            return

        chdir(tempdir)

        stdout, stderr, return_value = system_call(
            'tar xzf sortmerna-2.0-no-db.tar.gz')

        if return_value != 0:
            print ("Unable to extract SortMeRNA archive.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        chdir('sortmerna-2.0')

        environ['CXX'] = cxx
        environ['CC'] = cc

        stdout, stderr, return_value = system_call('bash build.sh')

        if return_value != 0:
            print ("Unable to build SortMeRNA.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        copy('sortmerna', scripts)
        copy('indexdb_rna', scripts)
        print "SortMeRNA built."
    finally:
        environ['CXX'] = cxx_old
        environ['CC'] = cc_old

        # remove the source
        rmtree(tempdir)
        chdir(cwd)


def build_SUMACLUST():
    """Download and build SUMACLUST then copy it to the scripts directory"""
    cwd = getcwd()
    scripts = join(cwd, 'scripts')

    try:
        tempdir = mkdtemp()
        if download_file('ftp://ftp.microbio.me/pub/QIIME-v1.9.0-dependencies/suma_package_V_1.0.00.tar.gz',
                         tempdir, 'suma_package_V_1.0.00.tar.gz'):
            print "Could not download SUMACLUST, so cannot install it."
            return

        chdir(tempdir)

        stdout, stderr, return_value = system_call(
            'tar xzf suma_package_V_1.0.00.tar.gz')

        if return_value != 0:
            print ("Unable to extract SUMACLUST archive.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        chdir('suma_package_V_1.0.00/sumaclust')

        stdout, stderr, return_value = system_call('make')

        if return_value != 0:
            print ("Unable to build SUMACLUST.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        copy('sumaclust', scripts)
        print "SUMACLUST built."
    finally:
        # remove the source
        rmtree(tempdir)
        chdir(cwd)


def build_swarm():
    """Download and build swarm then copy it to the scripts directory"""
    cwd = getcwd()
    scripts = join(cwd, 'scripts')

    try:
        tempdir = mkdtemp()
        if download_file('https://github.com/torognes/swarm/archive/1.2.19.tar.gz',
                         tempdir, '1.2.19.tar.gz'):
            print "Could not download swarm, so cannot install it."
            return

        chdir(tempdir)

        stdout, stderr, return_value = system_call('tar xzf 1.2.19.tar.gz')

        if return_value != 0:
            print ("Unable to extract swarm archive.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        chdir('swarm-1.2.19')

        stdout, stderr, return_value = system_call('make')

        if return_value != 0:
            print ("Unable to build swarm.\nstdout:\n%s\n\nstderr:\n%s\n" %
                   (stdout, stderr))
            return

        copy('swarm', scripts)
        copy('scripts/amplicon_contingency_table.py', scripts)
        copy('scripts/swarm_breaker.py', scripts)
        print "swarm built."
    finally:
        # remove the source
        rmtree(tempdir)
        chdir(cwd)


def download_UCLUST():
    """Download the UCLUST executable and set it to the scripts directory"""
    if platform == 'darwin':
        URL = 'http://www.drive5.com/uclust/uclustq1.2.22_i86darwin64'
    elif platform == 'linux2':
        URL = 'http://www.drive5.com/uclust/uclustq1.2.22_i86linux64'
    else:
        print "Platform %r not supported by UCLUST" % platform
        return

    return_value = download_file(URL, 'scripts/', 'uclust')

    # make the file an executable file
    if not return_value:
        chmod('scripts/uclust', stat('scripts/uclust').st_mode | S_IEXEC)
        print "UCLUST installed."
    else:
        print "UCLUST could not be installed."


def app_available(app_name):
    """Check if a binary is available and on the user Path

    Inputs:
    app_name: Name of the binary, i. e. 'ls', 'gcc' etc.

    Output:
    False if the binary is not found, True if the binary is found
    """
    # redirect all output to /dev/null so nothing is seen on screen
    devnull_fd = open(devnull, 'w')
    output = True

    try:
        call([app_name], stdout=devnull_fd, stderr=devnull_fd)
    except OSError:
        output = False
    finally:
        devnull_fd.close()

    return output


def system_call(cmd, shell=True):
    """Call cmd and return (stdout, stderr, return_value).

    cmd can be either a string containing the command to be run, or a sequence
    of strings that are the tokens of the command.

    Please see Python's subprocess. Popen for a description of the shell
    parameter and how cmd is interpreted differently based on its value.

    This function is ported from qcli (previously qcli_system_call).
    """
    proc = Popen(cmd,
                 shell=shell,
                 universal_newlines=True,
                 stdout=PIPE,
                 stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return stdout, stderr, return_value


def catch_install_errors(install_function, name):
    try:
        install_function()
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        print ("Skipping installation of %s due to failure while downloading, "
               "building, or installing:\n%s\n" % (name, exc_info()[0]))


# do not compile and build any of these if running under pip's egg_info
if build_stack:
    catch_install_errors(build_denoiser, 'denoiser')
    catch_install_errors(download_UCLUST, 'UCLUST')
    catch_install_errors(build_FastTree, 'FastTree')
    catch_install_errors(build_SortMeRNA, 'SortMeRNA')
    catch_install_errors(build_SUMACLUST, 'SUMACLUST')
    catch_install_errors(build_swarm, 'swarm')

# taken from PyNAST
classes = """
    Development Status :: 5 - Production/Stable
    Environment :: Console
    License :: OSI Approved :: GPL License
    Topic :: Software Development :: Bioinformatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: UNIX
    Operating System :: MacOS X
    Operating System :: POSIX :: BSD
    Operating System :: POSIX :: Linux
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

# compile the list of all qiime_test_data files that need to be installed.
# these must be relative file paths, beginning after the qiime_test_data
# directory
qiime_test_data_files = []
for root, dnames, fnames in walk('qiime_test_data'):
    try:
        # strip 'qiime_test_data/' from the beginning of root
        root = root.split('/', 1)[1]
    except IndexError:
        # if there is no '/', then we're in qiime_test_data itself
        # so there is nothing to do
        continue
    else:
        # if there is a slash, we're in a script test data directory,
        # so compile all relative filepaths
        for fname in fnames:
            qiime_test_data_files.append(join(root, fname))

setup(name='qiime',
      version=__version__,
      description='Quantitative Insights Into Microbial Ecology',
      long_description=long_description,
      author=__author__,
      classifiers=classifiers,
      author_email=__email__,
      maintainer=__author__,
      maintainer_email=__email__,
      url='http://www.qiime.org',
      packages=['qiime', 'qiime/parallel', 'qiime/pycogent_backports',
                'qiime/denoiser', 'qiime/workflow', 'qiime_test_data'],
      scripts=glob('scripts/*py') + glob('scripts/FlowgramAli_4frame') +
      glob('scripts/FastTree') + glob('scripts/uclust') +
      glob('scripts/indexdb_rna') + glob('scripts/sortmerna') +
      glob('scripts/sumaclust') + glob('scripts/swarm'),
      package_data={'qiime':
                    ['support_files/qiime_config',
                     'support_files/css/*css',
                     'support_files/html_templates/*html',
                     'support_files/images/*png',
                     'support_files/js/*js',
                     'support_files/R/*r',
                     'support_files/denoiser/Data/*',
                     'support_files/denoiser/TestData/*',
                     'support_files/denoiser/FlowgramAlignment/*.lhs',
                     'support_files/denoiser/FlowgramAlignment/Makefile'],
                    'qiime_test_data': qiime_test_data_files},
      license=__license__,
      keywords=['bioinformatics', 'microbiome', 'microbiology', 'qiime'],
      platforms=['MacOS', 'Linux'],
      install_requires=['numpy >= 1.9.0',
                        'scipy >= 0.14.0',
                        'matplotlib >= 1.1.0',
                        'pynast == 1.2.2', 'qcli >= 0.1.1, < 0.2.0', 'gdata',
                        'biom-format >= 2.1.2, < 2.2.0',
                        'emperor >= 0.9.5, < 1.0.0',
                        'scikit-bio >= 0.2.2, < 0.3.0',
                        'burrito-fillings >= 0.1.0, < 0.2.0',
                        'pandas >= 0.13.1', 'burrito  < 1.0.0',
                        'qiime-default-reference >= 0.1.1, < 0.2.0'],
      extras_require={'all': ['ipython[all]', 'sphinx >= 0.3']}
      )
