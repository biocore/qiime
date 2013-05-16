#!/usr/bin/env python
# File created on 17 Feb 2010
from __future__ import division
from distutils.core import setup
from distutils.sysconfig import get_python_lib
from os import chdir, getcwd, listdir, chmod, walk
from os.path import join, abspath
from subprocess import call
from glob import glob
import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 
 
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

try:
    import cogent
except ImportError:
    print "PyCogent not installed but required. (Is it installed? Is it in the current user's $PYTHONPATH or site-packages?) See http://www.pycogent.org."
    exit(1)
from cogent.util.misc import app_path

def build_html():
    """ Build the sphinx documentation 
    
    The code for building sphinx documentation is based on 
    PyCogent's setup.py.
    
    """
    cwd = getcwd()
    doc_dir = join(cwd,'doc')
    chdir(doc_dir)
    call(["make", "html"])
    chdir(cwd)
    index_html_path = join(abspath(doc_dir),'_build','html','index.html')
    print "Local documentation built with Sphinx. "+\
          "Open to following path with a web browser:\n%s" %\
            index_html_path
            
def build_denoiser():
    """ Build the denoiser code binary """
    cwd = getcwd()
    denoiser_dir = join(cwd,'qiime/support_files/denoiser/FlowgramAlignment')
    chdir(denoiser_dir)
    call(["make"])
    chdir(cwd)
    print "Denoiser built."

pycogent_version = tuple([int(v) \
        for v in re.split("[^\d]", cogent.__version__) if v.isdigit()])
        
if pycogent_version < (1,5,3):
    print "PyCogent >= 1.5.3 required, but %s is installed." % cogent.__version__
    exit(1)

if app_path("ghc"):
    build_denoiser()
else:
    print "GHC not installed, so cannot build the Denoiser binary."

# compile the list of all qiime_test_data files that need to be installed. 
# these must be relative file paths, beginning after the qiime_test_data
# directory
qiime_test_data_files = []
for root, dnames, fnames in walk('qiime_test_data'):
    try:
        # strip 'qiime_test_data/' from the beginning of root
        root = root.split('/',1)[1]
    except IndexError:
        # if there is no '/', then we're in qiime_test_data itself
        # so there is nothing to do
        continue
    else:
        # if there is a slash, we're in a script test data directory,
        # so compile all relative filepaths
        for fname in fnames:
            qiime_test_data_files.append(join(root,fname))

setup(name='QIIME',
      version=__version__,
      description='Quantitative Insights Into Microbial Ecology',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://www.qiime.org',
      packages=['qiime','qiime/parallel','qiime/pycogent_backports',
                'qiime/denoiser','qiime/workflow','qiime_test_data'],
      scripts=glob('scripts/*py')+glob('scripts/ec2*')+
              glob('scripts/FlowgramAli_4frame'),
      package_data={'qiime':
                   ['support_files/qiime_config',
                    'support_files/css/*css',
                    'support_files/html_templates/*html',
                    'support_files/images/*png',
                    'support_files/jar/*jar',
                    'support_files/js/*js',
                    'support_files/R/*r',
                    'support_files/denoiser/Data/*',
                    'support_files/denoiser/TestData/*'],
                    'qiime_test_data':qiime_test_data_files},
      long_description=long_description
)

if doc_imports_failed:
    print "Sphinx not installed, so cannot build local html documentation."
else:
    build_html()
