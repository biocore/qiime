#!/usr/bin/env python
# File created on 17 Feb 2010
from __future__ import division
from distutils.core import setup
from os import chdir, getcwd, listdir
from os.path import join, abspath
from subprocess import call
from glob import glob
import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 
 
long_description = """QIIME: Quantitative Insights Into Microbial Ecology
http://qiime.sourceforge.net

QIIME Allows Integration and Analysis of High-Throughput Community Sequencing Data
J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K. Goodrich, Jeffrey I. Gordon, Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight. 
Manuscript under review. 2010.
"""

doc_imports_failed = False
try:
    import sphinx
except ImportError:
    doc_imports_failed = True

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

try:
    import cogent
except ImportError:
    print "PyCogent not installed but required. (Is it installed? Is it in the current users $PYTHONPATH or site-packages?) See http://pycogent.sourceforge.net."
    exit(1)

pycogent_version = tuple([int(v) \
        for v in re.split("[^\d]", cogent.__version__) if v.isdigit()])
        
if pycogent_version < (1,4):
    print "PyCogent >= 1.4.0 required, but %s is installed." % cogent.__version__
    exit(1)
    
setup(name='QIIME',
      version=__version__,
      description='Quantitative Insights Into Microbial Ecology',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://qiime.sourceforge.net',
      packages=['qiime','qiime/parallel'],
      scripts=glob('scripts/*py'),
      package_data={'qiime':\
                   ['support_files/qiime_config',\
                    'support_files/css/*css',\
                    'support_files/html_templates/*html',\
                    'support_files/images/*png',\
                    'support_files/jar/*jar',\
                    'support_files/js/*js',\
                    'support_files/sra_xml_templates/*xml']},
      long_description=long_description,
)

if doc_imports_failed:
    print "Sphinx not installed, so cannot build local html documentation."
else:
    build_html()