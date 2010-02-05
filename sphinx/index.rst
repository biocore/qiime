.. _introduction:
.. QIIME documentation master file, created by
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Welcome to QIIME!
#################
QIIME (canonically pronounced 'Chime') stands for "Quantitative Insights Into Microbial Ecology". 
It allows a range of community analyses suitable for microbiome data using traditional and 
high-throughput sequencing methods. 

Overview
========
QIIME takes as input FASTA-format sequence data, optional quality scores, and a metadata file that 
assigns samples to specific categories or numerical values. It supports a wide range of microbial 
community analyses that have been useful in recent high-profile publications. It is a series of 
command line scripts, built in Python and the open-source PyCogent_ toolkit It makes extensive use 
of unit testing to ensure the accuracy of the results, and is highly modular (i.e. components for 
specific stages such as choosing OTUs (Operational Taxonomic Units), sequence alignment, and 
phylogeny, including 3rd party applications, can be easily integrated) 

Helpful Links
=============
.. toctree::
   :maxdepth: 1

`Tutorial <http://qiime.sourceforge.net/tutorial/Qiime_Tutorial.pdf>`_

Getting QIIME: Stable pre-release version
=========================================
Currently the most stable version of QIIME is our 0.9 pre-release, which you can download 
`here <https://sourceforge.net/projects/qiime/files/releases/Qiime-0.9.tar.gz/download>`_.
You can also download the `install documents <http://qiime.sourceforge.net/QIIME_Install.pdf>`_,
the software `documentation <http://qiime.sourceforge.net/documentation/Qiime_Documentation.pdf>`_, 
a `data analysis tutorial <http://qiime.sourceforge.net/tutorial/Qiime_Tutorial.pdf>`_ 
and an `SRA submission tutorial <http://qiime.sourceforge.net/SRA_Submission_Tutorial.pdf>`_.

Getting QIIME: Latest development version
=========================================
To get the latest development version of QIIME, you should check it out of our Sourceforge 
repository. While this code is subject to minor changes in interface, it will provide access to 
the latest and greatest features. The official web documentation is likely to be out-of-date with 
respect to the development software. You should instead refer to the svn documentation in 
Qiime/doc. Check out the latest version of QIIME using svn with the command::

	svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

Support
=======
Add a `feature request <https://sourceforge.net/tracker/?atid=1157167&group_id=272178&func=browse>`_

Report a `bug <https://sourceforge.net/tracker/?atid=1157164&group_id=272178&func=browse>`_

Citing QIIME
============
If you use QIIME for any published research, please include the following citation:

**QIIME Allows Integration and Analysis of High-Throughput Community Sequencing Data**
J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, 
Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Peña, Julia K. Goodrich, Jeffrey I. Gordon, 
Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, 
Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, 
William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob 
Knight. Manuscript under review. 2009.

Contributors
============
J. Gregory Caporaso [1], Justin Kuczynski [2], Jesse Stombaugh [1], 
Kyle Bittinger [3], Frederic D. Bushman [3], Elizabeth K. Costello [1], Noah Fierer [4], 
Antonio Gonzalez Peña [5], Julia K. Goodrich [5], Jeff I. Gordon6, Gavin Huttley [7], 
Scott T. Kelley [8], Dan Knights [5], Jeremy E. Koenig [9], Ruth E. Ley [9], 
Cathy A. Lozupone [1], Daniel McDonald [1], Brian D. Muegge [6], Megan Pirrung [1], 
Jens Reeder [1], Joel R. Sevinsky [10], Peter J. Turnbaugh [6], Will Van Treuren [1], 
William A. Walters [2], Jeremy Widmann [1], Tanya Yatsunenko [6], Jesse Zaneveld [2] and 
Rob Knight [1],[11] 

* [1] Department of Chemistry and Biochemistry, UCB 215, University of Colorado, Boulder, CO 80309 
* [2] Department of Molecular, Cellular and Developmental Biology, UCB 347, University of Colorado, Boulder, CO 80309 
* [3] Department of Microbiology, Johnson Pavilion 425, University of Pennsylvania, Philadelphia, PA 19104 
* [4] Cooperative Institute for Research in Environmental Sciences, University of Colorado, Boulder, CO 80309, USA.; Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO 80309, USA. 
* [5] Department of Computer Science, University of Colorado, Boulder, Colorado, USA. 
* [6] Center for Genome Sciences, Washington University School of Medicine, St. Louis, MO 63108 
* [7] Computational Genomics Laboratory, John Curtin School of Medical Research, The Australian National University, Canberra, Australian Capital Territory, Australia.
* [8] Department of Biology, San Diego State University, San Diego CA 92182
* [9] Department of Microbiology, Cornell University, Ithaca NY 14853
* [10] Luca Technologies, 500 Corporate Circle, Suite C, Golden, Colorado 80401 
* [11] Howard Hughes Medical Institute

.. _PyCogent: http://pycogent.sourceforge.net