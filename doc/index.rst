
.. QIIME documentation master file, created by
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Welcome to QIIME!
#################
QIIME (canonically pronounced 'Chime') stands for "Quantitative Insights Into Microbial Ecology". It allows a range of community analyses suitable for microbiome data using traditional and high-throughput sequencing methods.

About
========
QIIME takes as input FASTA-format sequence data, optional quality scores, and a metadata file that assigns samples to specific categories or numerical values. It supports a wide range of microbial community analyses that have been useful in recent high-profile publications. It is a series of command line scripts, built in Python and the open-source PyCogent_ toolkit It makes extensive use of unit testing to ensure the accuracy of the results, and is highly modular (i.e. components for specific stages such as choosing OTUs (Operational Taxonomic Units), sequence alignment, and phylogeny, including 3rd party applications, can be easily integrated)


Download
========

Virtual Box 
^^^^^^^^^^^^^^^^^^

The Virtual Box is pre-loaded with QIIME and its dependencies on a Ubuntu operating system.  To get the Virtual Box, please go here:

.. toctree::
   :maxdepth: 1

   virtual_box.rst


Stable Pre-Release
^^^^^^^^^^^^^^^^^^

Currently the most stable version of QIIME is our 0.92 pre-release, which you can download `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-0.92.tar.gz/download>`_. You can also view the `install documents <./install.html>`_, the software `documentation <./documentation.html>`_, a `data analysis tutorial <./tutorial.html>`_ and an `SRA submission tutorial <./doc_sra_submission.html>`_.

Development Version
^^^^^^^^^^^^^^^^^^^

To get the latest development version of QIIME, you should check it out of our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in Qiime/doc. Check out the latest version of QIIME using svn with the command::

	svn co http://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime



User Manuals
============

Installation
^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   install.rst

Tutorial
^^^^^^^^
.. toctree::
   :maxdepth: 1

   tutorial.rst

Documentation
^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   documentation.rst

Contact Us
===========

.. toctree::
   :maxdepth: 1

   developer/reporting_and_fixing_bugs.rst


- Add a `feature request <http://sourceforge.net/tracker/?atid=1157167&group_id=272178&func=browse>`_
* Contact `QIIME Support <qiime.help@colorado.edu>`_





References
============
If you use QIIME for any published research, please include the following citation:

\J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Peña, Julia K. Goodrich, Jeffrey I. Gordon, Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight. (Submitted) **QIIME allows integration and analysis of high-throughput community sequencing data.**

Contributors
============

.. note::

 \J. Gregory Caporaso :superscript:`1`, Justin Kuczynski :superscript:`2`, Jesse Stombaugh :superscript:`1`, Kyle Bittinger :superscript:`3`, Frederic D. Bushman :superscript:`3`, Elizabeth K. Costello :superscript:`1`, Noah Fierer :superscript:`4`, Antonio Gonzalez Peña :superscript:`5`, Julia K. Goodrich :superscript:`5`, Jeff I. Gordon :superscript:`6`, Gavin Huttley :superscript:`7`, Scott T. Kelley :superscript:`8`, Dan Knights :superscript:`5`, Jeremy E. Koenig :superscript:`9`, Ruth E. Ley :superscript:`9`, Cathy A. Lozupone :superscript:`1`, Daniel McDonald :superscript:`1`, Brian D. Muegge :superscript:`6`, Megan Pirrung :superscript:`1`, Jens Reeder :superscript:`1`, Joel R. Sevinsky :superscript:`10`, Peter J. Turnbaugh :superscript:`6`, Will Van Treuren :superscript:`1`, William A. Walters :superscript:`2`, Jeremy Widmann :superscript:`1`, Tanya Yatsunenko :superscript:`6`, Jesse Zaneveld :superscript:`2` and Rob Knight :superscript:`1,11`

 * :superscript:`1` Department of Chemistry and Biochemistry, UCB 215, University of Colorado, Boulder, CO 80309 
 * :superscript:`2` Department of Molecular, Cellular and Developmental Biology, UCB 347, University of Colorado, Boulder, CO 80309 
 * :superscript:`3` Department of Microbiology, Johnson Pavilion 425, University of Pennsylvania, Philadelphia, PA 19104 
 * :superscript:`4` Cooperative Institute for Research in Environmental Sciences, University of Colorado, Boulder, CO 80309, USA.; Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO 80309, USA. 
 * :superscript:`5` Department of Computer Science, University of Colorado, Boulder, Colorado, USA. 
 * :superscript:`6` Center for Genome Sciences, Washington University School of Medicine, St. Louis, MO 63108 
 * :superscript:`7` Computational Genomics Laboratory, John Curtin School of Medical Research, The Australian National University, Canberra, Australian Capital Territory, Australia.
 * :superscript:`8` Department of Biology, San Diego State University, San Diego CA 92182
 * :superscript:`9` Department of Microbiology, Cornell University, Ithaca NY 14853
 * :superscript:`10` Luca Technologies, 500 Corporate Circle, Suite C, Golden, Colorado 80401 
 * :superscript:`11` Howard Hughes Medical Institute

.. _PyCogent: http://pycogent.sourceforge.net