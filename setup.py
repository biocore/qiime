#!/usr/bin/env python
# File created on 17 Feb 2010
from __future__ import division
from setuptools import setup
from stat import S_IEXEC
from os import (chdir, getcwd, listdir, chmod, walk, rename, remove, chmod,
                stat, devnull)
from os.path import join, abspath
from sys import platform, argv
from subprocess import call
from glob import glob
from urllib import FancyURLopener
import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout",
               "Yoshiki Vazquez Baeza", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


long_description = """QIIME: Quantitative Insights Into Microbial Ecology
http://www.qiime.org

QIIME Allows Integration and Analysis of High-Throughput Community Sequencing Data
J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K. Goodrich, Jeffrey I. Gordon, Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight.
Nature Methods, 2010.
"""

doc_imports_failed = False
try:
    import sphinx
except ImportError:
    doc_imports_failed = True

# if egg_info is passed as an argument do not build any of the dependencies
build_stack = 'egg_info' not in argv


def build_html():
    """ Build the sphinx documentation

    The code for building sphinx documentation is based on
    PyCogent's setup.py.

    """
    cwd = getcwd()
    doc_dir = join(cwd, 'doc')
    chdir(doc_dir)
    call(["make", "html"])
    chdir(cwd)
    index_html_path = join(abspath(doc_dir), '_build', 'html', 'index.html')
    print "Local documentation built with Sphinx. " +\
          "Open to following path with a web browser:\n%s" %\
        index_html_path


def build_denoiser():
    """ Build the denoiser code binary """
    cwd = getcwd()
    denoiser_dir = join(cwd, 'qiime/support_files/denoiser/FlowgramAlignment')
    chdir(denoiser_dir)
    # make sure we compile the executable
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
    return return_value


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
    if app_available('ghc'):
        build_denoiser()
    else:
        print "GHC not installed, so cannot build the Denoiser binary."

    if app_available('gcc'):
        build_FastTree()
    else:
        print "GCC not installed, so cannot build FastTree"

    if download_UCLUST():
        print "UCLUST could not be installed."

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
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://www.qiime.org',
      packages=['qiime', 'qiime/parallel', 'qiime/pycogent_backports',
                'qiime/denoiser', 'qiime/workflow', 'qiime_test_data'],
      scripts=glob('scripts/*py') + glob('scripts/ec2*') +
      glob('scripts/FlowgramAli_4frame') + glob('scripts/FastTree') +
      glob('scripts/uclust'),
      package_data={'qiime':
                    ['support_files/qiime_config',
                     'support_files/css/*css',
                     'support_files/html_templates/*html',
                     'support_files/images/*png',
                     'support_files/jar/*jar',
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
      install_requires=['numpy >= 1.7.1',
                        'click >= 1.0',
                        'ipython[all]',
                        'networkx',
                        'matplotlib >= 1.1.0, <= 1.3.1',
                        'pynast == 1.2.2', 'qcli', 'gdata',
                        'biom-format == 2.0.1-dev', 'emperor == 0.9.3-dev',
                        'scikit-bio == 0.1.4', 'brokit == 0.0.0-dev',
                        'pandas >= 0.13.1', 'future', 'h5py>=2.2.0'],
      dependency_links=[
          'https://github.com/biocore/brokit/archive/master.zip#egg=brokit-0.0.0-dev',
          'https://github.com/biocore/biom-format/archive/master.zip#egg=biom-format-2.0.1-dev',
          'https://github.com/biocore/emperor/archive/master.zip#egg=emperor-0.9.3-dev',
          'https://github.com/biocore/scikit-bio/archive/master.zip#egg=scikit-bio-0.1.1-dev'
      ],
      extras_require={'all': ['tornado', 'sphinx >= 0.3',
                              # the following are optional for pycogent, should
                              # remove when pycogent is no longer a dependency
                              'MySQL-python', 'SQLAlchemy', 'mpi4py']}
      )

if build_stack:
    if doc_imports_failed:
        print "Sphinx not installed, so cannot build local html documentation."
    else:
        build_html()
