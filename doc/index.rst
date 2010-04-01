
.. QIIME documentation master file, created by
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

######################################################
QIIME: Quantitative Insights Into Microbial Ecology
######################################################
QIIME (canonically pronounced 'Chime') is a pipeline for performing microbial community analysis that integrates many third party tools which have become standard in the field. To stay up-to-date on what's new with QIIME, you should subscribe to the `blog <http://qiime.wordpress.com>`_.

QIIME can run on a laptop, a supercomputer, and everything in between. For ease of installation beginner users or users interested in testing QIIME should start with `QIIME virtual box <./virtual_box.html>`_. Linux and Mac OS X are supported natively, and Windows, Linux, and Mac OS X are supported via the virtual box.

Rather than reimplementing commonly used algorithms, QIIME wraps popular implementations of those algorithms. This allows us to make use of the many excellent tools available in this area, and allows faster integration of new tools. If you use tools that you think would be useful additions to QIIME, consider submitting a `feature request <http://sourceforge.net/tracker/?atid=1157167&group_id=272178&func=browse>`_.

A standard QIIME analysis begins with sequence data from one or more sequencing platforms, including Sanger, Roche/454, and Illumina GAIIx. QIIME can perform library de-multiplexing and quality filtering; denoising with PyroNoise; OTU and representative set picking with uclust, cdhit, mothur, BLAST, or other tools; taxonomy assignment with BLAST or the RDP classifier; sequence alignment with PyNAST, muscle, infernal, or other tools; phylogeny reconstruction with FastTree, raxml, clearcut, or other tools; alpha diversity and rarefaction, including visualization of results, using over 20 metrics including Phylogenetic Diversity, chao1, and observed species; beta diversity and rarefaction, including visualization of results, using over 25 metrics including weighted and unweighted UniFrac, Euclidean distance, and Bray-Curtis; summarization and visualization of taxonomic composition of samples using pie charts and histograms; and many other features.

QIIME includes parallelization capabilities for many of the computationally intensive steps. By default, these are configured to utilize a mutli-core environment, and are easily configured to run in a cluster environment. QIIME is built in Python using the open-source PyCogent_ toolkit. It makes extensive use of unit tests, and is highly modular to facilitate custom analyses.

Blog and Mailing List
======================
We recommend that all QIIME users keep an eye on the QIIME blog for important announcements. You can `subscribe to the RSS feed <http://qiime.wordpress.com/feed/>`_ or `sign up for e-mail notifications on the front page of the blog <http://qiime.wordpress.com>`_. This is a very low traffic list (typically around one message per month), and we will not share subscriber information with anyone.


Download
========

Virtual Box 
^^^^^^^^^^^^^^^^^^

The Virtual Box is pre-loaded with QIIME and its dependencies on a Ubuntu operating system.  To get the Virtual Box, please go here:

.. toctree::
   :maxdepth: 1

   virtual_box.rst


Stable Release
^^^^^^^^^^^^^^

Currently the most stable version of QIIME is our 0.92 pre-release, which you can download `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-0.92.tar.gz/download>`_. You can also view the `install documents <./install.html>`_, the software `documentation <./documentation.html>`_, a `data analysis tutorial <./tutorial.html>`_ and an `SRA submission tutorial <./doc_sra_submission.html>`_.

Development Version
^^^^^^^^^^^^^^^^^^^

QIIME is under very active development. To get the latest development version of QIIME, you access our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in Qiime/doc. Check out the latest version of QIIME using svn with the command::

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

\J. Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D. Bushman, Elizabeth K. Costello, Noah Fierer, Antonio Gonzalez Peña, Julia K. Goodrich, Jeffrey I. Gordon, Gavin A. Huttley, Scott T. Kelley, Dan Knights, Jeremy E. Koenig, Ruth E. Ley, Cathy A. Lozupone, Daniel McDonald, Brian D. Muegge, Meg Pirrung, Jens Reeder, Joel R. Sevinsky, Peter J. Turnbaugh, William van Treuren, William A. Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight. (Submitted) **QIIME allows analysis of high-throughput community sequencing data.**

The data presented in this paper can be found `here <http://tajmahal.colorado.edu/qiime/qiime_paper_data.zip>`_.

Contributors
============

.. note::

 \J. Gregory Caporaso :superscript:`1`, Justin Kuczynski :superscript:`2`, Jesse Stombaugh :superscript:`1`, Kyle Bittinger :superscript:`3`, Frederic D. Bushman :superscript:`3`, Elizabeth K. Costello :superscript:`1`, Noah Fierer :superscript:`4`, Antonio Gonzalez Peña :superscript:`5`, Julia K. Goodrich :superscript:`5`, Jeff I. Gordon :superscript:`6`, Gavin Huttley :superscript:`7`, Scott T. Kelley :superscript:`8`, Dan Knights :superscript:`5`, Jeremy E. Koenig :superscript:`9`, Ruth E. Ley :superscript:`9`, Cathy A. Lozupone :superscript:`1`, Daniel McDonald :superscript:`1`, Brian D. Muegge :superscript:`6`, Megan Pirrung :superscript:`1`, Jens Reeder :superscript:`1`, Joel R. Sevinsky :superscript:`10`, Peter J. Turnbaugh :superscript:`6`, William A. Walters :superscript:`2`, Jeremy Widmann :superscript:`1`, Tanya Yatsunenko :superscript:`6`, Jesse Zaneveld :superscript:`2` and Rob Knight :superscript:`1,11`

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