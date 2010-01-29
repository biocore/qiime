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

   tutorial.rst

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
J. Gregory Caporaso [#f1]_, Justin Kuczynski [#f2]_, Jesse Stombaugh [#f1]_, 
Kyle Bittinger [#f3]_, Frederic D. Bushman [#f3]_, Elizabeth K. Costello [#f1]_, Noah Fierer [#f4]_, 
Antonio Gonzalez Peña [#f5]_, Julia K. Goodrich [#f5]_, Jeff I. Gordon6, Gavin Huttley [#f7]_, 
Scott T. Kelley [#f8]_, Dan Knights [#f5]_, Jeremy E. Koenig [#f9]_, Ruth E. Ley [#f9]_, 
Cathy A. Lozupone [#f1]_, Daniel McDonald [#f1]_, Brian D. Muegge [#f6]_, Megan Pirrung [#f1]_, 
Jens Reeder [#f1]_, Joel R. Sevinsky [#f10]_, Peter J. Turnbaugh [#f6]_, Will Van Treuren [#f1]_, 
William A. Walters [#f2]_, Jeremy Widmann [#f1]_, Tanya Yatsunenko [#f6]_, Jesse Zaneveld [#f2]_ and 
Rob Knight [#f1]_, [#f11]_ 

.. [#f1] Department of Chemistry and Biochemistry, UCB 215, University of Colorado, Boulder, CO 80309 
.. [#f2] Department of Molecular, Cellular and Developmental Biology, UCB 347, University of Colorado, Boulder, CO 80309 
.. [#f3] Department of Microbiology, Johnson Pavilion 425, University of Pennsylvania, Philadelphia, PA 19104 
.. [#f4] Cooperative Institute for Research in Environmental Sciences, University of Colorado, Boulder, CO 80309, USA.; Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO 80309, USA. 
.. [#f5] Department of Computer Science, University of Colorado, Boulder, Colorado, USA. 
.. [#f6] Center for Genome Sciences, Washington University School of Medicine, St. Louis, MO 63108 
.. [#f7] Computational Genomics Laboratory, John Curtin School of Medical Research, The Australian National University, Canberra, Australian Capital Territory, Australia.
.. [#f8] Department of Biology, San Diego State University, San Diego CA 92182
.. [#f9] Department of Microbiology, Cornell University, Ithaca NY 14853
.. [#f10] Luca Technologies, 500 Corporate Circle, Suite C, Golden, Colorado 80401 
.. [#f11] Howard Hughes Medical Institute

.. _PyCogent: http://pycogent.sourceforge.net