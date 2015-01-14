#!/usr/bin/env python
# File created on 17 Feb 2010
from __future__ import division
from setuptools import setup
from stat import S_IEXEC
from os import (chdir, getcwd, listdir, chmod, walk, rename, remove, chmod,
                stat, devnull, environ)
from os.path import join, abspath
from sys import platform, argv
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

# if egg_info is passed as an argument do not build any of the dependencies
build_stack = 'egg_info' not in argv


def build_denoiser():
    """ Build the denoiser code binary """
    if not app_available('ghc'):
        print ("GHC not installed, so cannot build the denoiser binary "
               "'FlowgramAli_4frame'.")
        return

    cwd = getcwd()
    denoiser_dir = join(cwd, 'qiime/support_files/denoiser/FlowgramAlignment')
    chdir(denoiser_dir)
    call(["make", "clean"])
    call(["make"])
    chdir(cwd)
    print "Denoiser built."

# heavily based on lib.util.download_file from github.com/qiime/qiime-deploy


class URLOpener(FancyURLopener):

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        msg = 'ERROR: Could not download %s\nIs the URL valid?' % url
        raise IOError(msg)

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
    denoiser_dir = join(cwd, 'scripts')
    chdir(denoiser_dir)

    # as suggested by the docs in FastTree.c
    call(['gcc', '-Wall', '-O3', '-finline-functions', '-funroll-loops', '-o',
          'FastTree', 'FastTree.c', '-lm'])

    # remove the source
    remove('FastTree.c')
    chdir(cwd)
    print "FastTree built."


def build_SortMeRNA():
    """Download and build SortMeRNA then copy it to the scripts directory"""
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

    tempdir = mkdtemp()
    if download_file('ftp://ftp.microbio.me/pub/sortmerna-2.0-no-db.tar.gz',
                     tempdir, 'sortmerna-2.0-no-db.tar.gz'):
        print "Could not download SortMeRNA, so cannot install it."
        rmtree(tempdir)
        return

    cwd = getcwd()
    scripts = join(cwd, 'scripts')
    chdir(tempdir)

    call(['tar', 'xzf', 'sortmerna-2.0-no-db.tar.gz'])
    chdir('sortmerna-2.0')

    cxx_old = environ.get('CXX', '')
    cc_old = environ.get('CC', '')

    environ['CXX'] = cxx
    environ['CC'] = cc

    proc = Popen('bash build.sh',
                 universal_newlines=True,
                 shell=True,
                 stdout=PIPE,
                 stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    if return_value != 0:
        raise ValueError("Unable to build SortMeRNA")

    copy('sortmerna', scripts)
    copy('indexdb_rna', scripts)

    environ['CXX'] = cxx_old
    environ['CC'] = cc_old

    # remove the source
    rmtree(tempdir)
    chdir(cwd)
    print "SortMeRNA built."


def download_UCLUST():
    """Download the UCLUST executable and set it to the scripts directory"""

    if platform == 'darwin':
        URL = 'http://www.drive5.com/uclust/uclustq1.2.22_i86darwin64'
    elif platform == 'linux2':
        URL = 'http://www.drive5.com/uclust/uclustq1.2.22_i86linux64'
    else:
        raise SystemError(("Platform not supported by UCLUST"))

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

# do not compile and build any of these if running under pip's egg_info
if build_stack:
    build_denoiser()
    build_FastTree()
    build_SortMeRNA()
    download_UCLUST()

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
      glob('scripts/indexdb_rna') + glob('scripts/sortmerna'),
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
