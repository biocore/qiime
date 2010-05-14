
.. QIIME documentation master file, created by
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

######################################################
QIIME: Quantitative Insights Into Microbial Ecology
######################################################
QIIME (canonically pronounced 'Chime') is a pipeline for performing microbial community analysis that integrates many third party tools which have become standard in the field. QIIME can run on a laptop, a supercomputer, and systems in between such as multicore desktops.  Linux and Mac OS X are supported natively, and Windows, Linux, and Mac OS X are supported via the `QIIME virtual box <./install/virtual_box.html>`_.

**The quickest way to start using QIIME is with the** `QIIME virtual box <./install/virtual_box.html>`_. To stay up-to-date on what's new with QIIME, you should subscribe to the `blog <http://qiime.wordpress.com>`_.

Rather than reimplementing commonly used algorithms, QIIME wraps popular implementations of those algorithms. This allows us to make use of the many excellent tools available in this area, and allows faster integration of new tools. If you use tools that you think would be useful additions to QIIME, consider submitting a `feature request <http://sourceforge.net/tracker/?atid=1157167&group_id=272178&func=browse>`_.

A standard QIIME analysis begins with sequence data from one or more sequencing platforms, including Sanger, Roche/454, and Illumina GAIIx. QIIME can perform library de-multiplexing and quality filtering; denoising with PyroNoise; OTU and representative set picking with uclust, cdhit, mothur, BLAST, or other tools; taxonomy assignment with BLAST or the RDP classifier; sequence alignment with PyNAST, muscle, infernal, or other tools; phylogeny reconstruction with FastTree, raxml, clearcut, or other tools; alpha diversity and rarefaction, including visualization of results, using over 20 metrics including Phylogenetic Diversity, chao1, and observed species; beta diversity and rarefaction, including visualization of results, using over 25 metrics including weighted and unweighted UniFrac, Euclidean distance, and Bray-Curtis; summarization and visualization of taxonomic composition of samples using pie charts and histograms; and many other features.

QIIME includes parallelization capabilities for many of the computationally intensive steps. By default, these are configured to utilize a mutli-core environment, and are easily configured to run in a cluster environment. QIIME is built in Python using the open-source PyCogent_ toolkit. It makes extensive use of unit tests, and is highly modular to facilitate custom analyses.

Blog and Mailing List
======================
We recommend that all QIIME users keep an eye on the QIIME blog for important announcements. You can `subscribe to the RSS feed <http://qiime.wordpress.com/feed/>`_ or `sign up for e-mail notifications on the front page of the blog <http://qiime.wordpress.com>`_. This is a very low traffic list (typically around one message per month), and we will not share subscriber information with anyone.


Download QIIME
===============

 * Virtual Box: The QIIME Virtual Box is an Ubuntu Linux virtual machine, pre-loaded with QIIME and its dependencies. This is the quickest way to start using QIIME. To get the Virtual Box, please `go here <./install/virtual_box.html>`_.

 * Stable Release: Currently the most stable version of QIIME is our 1.1.0 release, which you can `download here <http://sourceforge.net/projects/qiime/files/releases/Qiime-1.1.0.tar.gz/download>`_.

 * Development Version: QIIME is under very active development. To get the latest development version of QIIME, you access our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in Qiime/doc. Check out the latest version of QIIME using svn with the command::

	svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

Installing and using QIIME
==========================
New users should begin with the `QIIME installation guide <./install/install.html>`_ or the `QIIME virtual box <./install/virtual_box.html>`_. After installing QIIME, you should move on to the `QIIME overview tutorial <./tutorials/tutorial.html>`_ to analyze an example data set. As you begin using QIIME on your own data, you'll want to refer to the `QIIME documentation <./documentation/index.html>`_ and the general `QIIME tutorials <./tutorials/index.html>`_.

Contact Us
===========
For technical support, contact `QIIME Support <qiime.help@colorado.edu>`_. This is likely to be the fastest way to get help from the QIIME developers. As multiple developers monitor this address, response is likely to be a lot faster than contacting developers individually.

Users can also submit `bug reports <http://sourceforge.net/tracker/?group_id=272178&atid=1157164>`_ and `feature requests <http://sourceforge.net/tracker/?group_id=272178&atid=1157167>`_ using via Sourceforge.


QIIME Development
====================

QIIME is an open-source project, primarily developed in the Knight Lab at the University of Colorado at Boulder. If you are interested in getting involved, check out the `developer notes <./developer/index.html>`_.

Citing QIIME
============
If you use QIIME for any published research, please include the following citation:

	**QIIME allows analysis of high-throughput community sequencing data**
	
	J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010; doi:10.1038/nmeth.f.303


You can find the `QIIME paper here <http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.f.303.html>`_, and the data presented in this paper can be found `here <http://bmf.colorado.edu/QIIME/QIIME_NM_2010.tgz>`_.

.. _PyCogent: http://pycogent.sourceforge.net