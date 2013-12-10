.. _doc_install:
.. QIIME documentation master file, created by Jesse Stombaugh
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. index:: Installing QIIME

=========================
QIIME installation guide
=========================

QIIME consists of native python code and additionally wraps many external applications. This gives the user flexibility to easily build their own analysis pipelines, making use of popular microbial community analysis tools. QIIME handles the processing of input and output of these applications, so the user can spend time analyzing their data rather than parsing, writing, and converting file formats.

As a consequence of this 'pipeline' architecture, **QIIME has a lot of dependencies and can be very challenging to install**.


How to not install QIIME
========================

Because QIIME is hard to install, we have attempted to shift this burden to the QIIME development group rather than our users by providing virtual machines with QIIME and all of its dependencies pre-installed. We, and third-party developers, have also created several automated installation procedures. These alternatives (`summarized here <../index.html#downloading-and-installing-qiime>`_) allow you to bypass the complex installation procedure and have access to a full, working QIIME installation. 

**We highly recommend going with one of these solutions if you're new to QIIME, or just want to test it out to see if it will do what you want.**

How to install QIIME
====================

If you want to customize QIIME, work with QIIME in a multi-user environment (e.g., a Linux cluster), are interested in getting involved in QIIME development, or want to use the development version of QIIME, you may need to install QIIME manually (but you also may not - see `MacQIIME <http://www.wernerlab.org/software/macqiime>`_ and `qiime-deploy <https://github.com/qiime/qiime-deploy>`_).

Depending on the features of QIIME that you or your users are interested in, a *QIIME base install* may be sufficient. This is much easier than a *QIIME full install*. Both of these procedures are covered in the following sections.

To install most of following dependencies, either for the base or full install, you need to have a build environment on your machine. On OS X, this involves installing the `developer tools <http://developer.apple.com/technologies/xcode.html>`_. On Debian-based Linux (e.g., Ubuntu), this involves installing the ``build-essential`` package::

	sudo apt-get install build-essential

Shortcuts in this document
--------------------------
For simplicity throughout this document, we assume that you have downloaded QIIME in ``$HOME/``. You should consider all occurrences of ``$HOME/`` in the remainder of this document as references to the directory which contains the QIIME directory which you'll have after downloading and unpacking QIIME.

QIIME base install (standard QIIME workflow with default parameters)
--------------------------------------------------------------------

When getting started with QIIME, most users will want to begin with the QIIME base install. This allows users to run through the core QIIME workflows (including, but not limited to `validate_mapping_file.py <../scripts/validate_mapping_file.html>`_, `split_libraries.py <../scripts/split_libraries.html>`_, `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_, `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_, `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, and `core_diversity_analyses.py <../scripts/core_diversity_analyses.html>`_) while limiting the time and effort spent on installation.

You can install the QIIME base dependencies either via pip or by manually downloading and installing each package.

Installing QIIME via pip
++++++++++++++++++++++++

The easiest way to install the latest QIIME release and its base dependencies is by using pip::

	pip install numpy==1.7.1
	pip install qiime

**Note:** You may need to prefix the above commands with ``sudo`` if you do not have permission to place files in the default locations. You will also need an active Internet connection. If you do not have pip, the easiest way to install it is by running::

	easy_install pip

That's it!

Manually installing QIIME
+++++++++++++++++++++++++

If you decided not to install QIIME using pip, you can install it (and its dependencies) manually.

The following QIIME base install dependencies are grouped by installation method.

The first are the core scientific python dependencies. The easiest way to install these is by installing `Canopy Express <https://www.enthought.com/canopy-express/>`_ (formerly EPD Free), which contains core modules for python scientific computing, including those required by QIIME, but also packages such as IPython and Pandas, which QIIME users may also find useful.

* Python 2.7.3 (`src_python <http://www.python.org/ftp/python/2.7.3/Python-2.7.3.tgz>`_) (license: PSF)
* Numpy 1.5.1 - 1.7.1 (`src_numpy <http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/numpy-1.7.1.tar.gz/download>`_) (license: BSD)
* MatPlotLib 1.1.0 - 1.3.1 (`src_matplotlib <http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.1.0/matplotlib-1.1.0.tar.gz>`_) (license: PFS)

The next are python packages not included in Canopy Express. Each of these can be installed either via pip (e.g., ``pip install <package-name>``) or by downloading the package, unzipping with ``tar -xzf``, changing to the resulting directory, and running ``python setup.py install`` (see :ref:`Installing with setup.py <python-setup>` for some notes that may be useful).

* QIIME (see :ref:`Getting QIIME <getting-qiime>`)
* PyCogent 1.5.3 (`src_pycogent <https://pypi.python.org/packages/source/c/cogent/cogent-1.5.3.tgz>`_) (license: GPL)
* biom-format 1.3.1 (`src_biom <https://pypi.python.org/packages/source/b/biom-format/biom-format-1.3.1.tar.gz>`_) (license: BSD)
* qcli 0.1.0 (`src_qcli <https://pypi.python.org/packages/source/q/qcli/qcli-0.1.0.tar.gz>`_) (license: GPL)
* PyNAST 1.2.2 (`src_pynast <https://pypi.python.org/packages/source/p/pynast/pynast-1.2.2.tar.gz>`_) (license: BSD)
* Emperor 0.9.3 (`src_emperor <https://pypi.python.org/packages/source/e/emperor/emperor-0.9.3.tar.gz>`_) (license: BSD)

Next, there are two non-python dependencies required for the QIIME base package. These should be installed by following their respective install instructions. 

* uclust 1.2.22q (`src_uclust <http://www.drive5.com/uclust/downloads1_2_22q.html>`_) See :ref:`uclust install notes <uclust-install>`. (licensed specially for Qiime and PyNAST users)
* fasttree 2.1.3 (`src_fasttree <http://www.microbesonline.org/fasttree/FastTree-2.1.3.c>`_) See `FastTree install instructions <http://www.microbesonline.org/fasttree/#Install>`_ (license: GPL)

Data files necessary for the QIIME base installation
++++++++++++++++++++++++++++++++++++++++++++++++++++

After you've installed the base QIIME dependencies, there are several data files that are likely to be useful in your analyses. These can all be obtained using ``wget`` or ``curl``, and unzipping where necessary.

* greengenes core set data file (`fasta <http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta.imputed>`_)
* greengenes alignment lanemask file (`txt <http://greengenes.lbl.gov/Download/Sequence_Data/lanemask_in_1s_and_0s>`_)
* Marker gene reference OTUs, taxonomies, and trees  (follow the *Resources* link from the `QIIME homepage <http://www.qiime.org>`_)

You should next `write your QIIME config file <./qiime_config.html>`_.

Testing the QIIME base installation
-----------------------------------

After installing the QIIME base packages, you can test this for sanity by running::

	print_qiime_config.py -tb

You should see output that looks like the following::

	System information
	==================
	         Platform:	darwin
	   Python version:	2.7.1 (r271:86832, Aug 30 2012, 10:07:33)  [GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)]
	Python executable:	$HOME/.virtualenvs/qiime/bin/python

	Dependency versions
	===================
	             PyCogent version:	1.5.3
	                NumPy version:	1.5.1
	           matplotlib version:	1.1.0
	          biom-format version:	1.2.0
	                 qcli version:	0.1.0
	        QIIME library version:	1.7.0
	         QIIME script version:	1.7.0
	PyNAST version (if installed):	1.2.1
	              Emperor version:	0.9.2

	QIIME config values
	===================
	                     blastmat_dir:	None
	                         sc_queue:	all.q
	      topiaryexplorer_project_dir:	$HOME/code/TopiaryExplorer-0.9.1/
	     pynast_template_alignment_fp:	$HOME/data/greengenes_core_sets/core_set_aligned_imputed.fasta_11_8_07.no_dots
	                  cluster_jobs_fp:	start_parallel_jobs.py
	pynast_template_alignment_blastdb:	None
	assign_taxonomy_reference_seqs_fp:	$HOME/data/gg_13_5_otus/rep_set/97_otus.fasta
	                     torque_queue:	friendlyq
	   template_alignment_lanemask_fp:	$HOME/data/greengenes_core_sets/lanemask_in_1s_and_0s.txt
	                    jobs_to_start:	2
	                cloud_environment:	False
	                qiime_scripts_dir:	$HOME/code/Qiime/scripts
	            denoiser_min_per_core:	50
	                      working_dir:	None
	                    python_exe_fp:	python
	                         temp_dir:	$HOME/temp
	                      blastall_fp:	blastall
	                 seconds_to_sleep:	1
	assign_taxonomy_id_to_taxonomy_fp:	$HOME/data/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt
	................
	----------------------------------------------------------------------
	Ran 16 tests in 0.440s
	
	OK

This indicates that you have a complete QIIME base install. 

You should next :ref:`run QIIME's unit tests <run-test-suite>`. You will experience some test failures as a result of not having a full QIIME install. If you have questions about these failures, you should post to the `QIIME Forum <http://forum.qiime.org>`_.

QIIME full install (for access to advanced features in QIIME, and non-default processing pipelines)
---------------------------------------------------------------------------------------------------

The dependencies described below will support a full QIIME install. These are grouped by the features that each dependency will provide access to. Installation instructions should be followed for each individual package (e.g., from the project's website or README/INSTALL file). 

Alignment, tree-building, taxonomy assignment, OTU picking, and other data generation steps (required for non-default processing pipelines):

* jre1.6.0_05 (`src_jre <http://java.sun.com/javase/downloads/index.jsp>`_) (license: GPL2)
* rdp_classifier-2.2 (`src_rdp <http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.2.zip/download>`_) See :ref:`RDP install notes <rdp-install>`. (license: GPL)
* tax2tree 1.0.0 (`src_tax2tree <https://downloads.sourceforge.net/project/tax2tree/tax2tree-v1.0.tar.gz>`_)
* blast-2.2.22 (legacy BLAST from NCBI, *NOT* BLAST+) (`OS X <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-universal-macosx.tar.gz>`_ or `linux 32-bit <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-ia32-linux.tar.gz>`_) (license: GNU)
* cd-hit 3.1.1 (`src_cdhit <http://www.bioinformatics.org/download/cd-hit/cd-hit-2007-0131.tar.gz>`_) (license: Free access)
* ChimeraSlayer (via microbiomeutil_2010-04-29) (`src_chimeraslayer <http://sourceforge.net/projects/microbiomeutil/files/>`_) See :ref:`ChimeraSlayer install notes <chimeraslayer-install>`.
* mothur 1.25.0 (`src_mothur <http://www.mothur.org/w/images/6/6d/Mothur.1.25.0.zip>`_) (license: GPL)
* clearcut v1.0.9 (`src_clearcut <http://www.mothur.org/w/images/9/91/Clearcut.source.zip>`_)
* raxml 7.3.0 (`src_raxml <ftp://thebeast.colorado.edu/pub/QIIME-v1.5.0-dependencies/stamatak-standard-RAxML-5_7_2012.tgz>`_)
* infernal 1.0.2 (`src_infernal <ftp://selab.janelia.org/pub/software/infernal/infernal.tar.gz>`_) (license: GPL)
* cdbtools (`src_cdbtools <ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/cdbfasta/cdbfasta.tar.gz>`_)
* muscle 3.8.31 (`src_muscle <http://www.drive5.com/muscle/downloads.htm>`_) (Public domain)
* rtax 0.984 (`src_rtax <http://static.davidsoergel.com/rtax-0.984.tgz>`_) (license: BSD)
* pplacer 1.1 (`src_pplacer <http://matsen.fhcrc.org/pplacer/builds/pplacer-v1.1-Linux.tar.gz>`_) (license: GPL)
* ParsInsert 1.04 (`src_parsinsert <http://downloads.sourceforge.net/project/parsinsert/ParsInsert.1.04.tgz>`_) (license: GPL)
* usearch v5.2.236 and/or usearch v6.1 (`src_usearch <http://www.drive5.com/usearch/>`_) (license: see http://www.drive5.com/usearch/nonprofit_form.html) **At this stage two different versions of usearch are supported.** usearch v5.2.236 is referred to as ``usearch`` in QIIME, and usearch v6.1 is referred to as ``usearch61``.

Processing sff files:

* sfffile and sffinfo (optional, QIIME 1.2.0 and later contain built-in tools for processing sff files although they are about 10x slower than the tools from Roche) (license: proprietary - must be obtained from Roche/454)

Denoising 454 data:

* GNU Science Library (required by AmpliconNoise) (`src_gsl <ftp://ftp.gnu.org/gnu/gsl/gsl-1.9.tar.gz>`_)
* AmpliconNoise 1.27 (`src_ampliconnoise <http://ampliconnoise.googlecode.com/files/AmpliconNoiseV1.27.tar.gz>`_) See :ref:`AmpliconNoise install notes <ampliconnoise-install>`.
* ghc 6.8 (required by the QIIME denoiser) (`src_ghc <http://haskell.org/ghc>`_)

Visualization and plotting steps:

* cytoscape v2.7.0 (`src_cytoscape <http://www.cytoscape.org/>`_) (license: LGPL)

Supervised learning (``supervised_learning.py``) and ``compare_categories.py``:

* R 3.0.2 (`src_r <http://www.r-project.org/>`_) See :ref:`R install notes <R-install>`. (license: GPL2)

If you plan to build the QIIME documentation locally:

* Sphinx 1.0.4 (`src <http://pypi.python.org/pypi/Sphinx>`_) See :ref:`Building the QIIME documentation <build-qiime-docs>` (license: BSD)

If you plan to use remote mapping files (stored as Google Spreadsheets) with QIIME (see the tutorial `here <../tutorials/remote_mapping_files.html>`_):

* gdata 2.0.17 (`src <http://gdata-python-client.googlecode.com/files/gdata-2.0.17.tar.gz>`_) (license: Apache 2.0)

If you plan to use SourceTracker with QIIME:

* SourceTracker 0.9.5 (`src <http://downloads.sourceforge.net/project/sourcetracker/sourcetracker-0.9.5.tar.gz>`_) (license: GPL)

Testing the QIIME full installation
-----------------------------------

After installing the QIIME base packages, you can test this for sanity by running::

	print_qiime_config.py -t

You should see output that looks like the following::

	System information
	==================
	         Platform:	darwin
	   Python version:	2.7.1 (r271:86832, Aug 30 2012, 10:07:33)  [GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)]
	Python executable:	$HOME/.virtualenvs/qiime/bin/python

	Dependency versions
	===================
	                     PyCogent version:	1.5.3
	                        NumPy version:	1.5.1
	                   matplotlib version:	1.1.0
	                  biom-format version:	1.2.0-dev
	                         qcli version:	0.1.0
	                QIIME library version:	1.7.0-dev
	                 QIIME script version:	1.7.0-dev
	        PyNAST version (if installed):	1.2.1
	                      Emperor version:	0.9.2-dev
	RDP Classifier version (if installed):	rdp_classifier-2.2.jar
	          Java version (if installed):	1.6.0_43

	QIIME config values
	===================
	                     blastmat_dir:	/Applications/blast-2.2.22/data/
	                         sc_queue:	all.q
	      topiaryexplorer_project_dir:	$HOME/code/TopiaryExplorer-0.9.1/
	     pynast_template_alignment_fp:	$HOME/data/greengenes_core_sets/core_set_aligned_imputed.fasta_11_8_07.no_dots
	                  cluster_jobs_fp:	start_parallel_jobs.py
	pynast_template_alignment_blastdb:	None
	assign_taxonomy_reference_seqs_fp:	$HOME/data/gg_13_5_otus/rep_set/97_otus.fasta
	                     torque_queue:	friendlyq
	   template_alignment_lanemask_fp:	$HOME/data/greengenes_core_sets/lanemask_in_1s_and_0s.txt
	                    jobs_to_start:	2
	                cloud_environment:	False
	                qiime_scripts_dir:	$HOME/code/Qiime/scripts
	            denoiser_min_per_core:	50
	                      working_dir:	None
	                    python_exe_fp:	python
	                         temp_dir:	$HOME/temp
	                      blastall_fp:	blastall
	                 seconds_to_sleep:	1
	assign_taxonomy_id_to_taxonomy_fp:	$HOME/data/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt
	...................................
	----------------------------------------------------------------------
	Ran 35 tests in 0.641s

	OK

You should next :ref:`run QIIME's unit tests <run-test-suite>`. All tests should pass if you have a working full QIIME installation. If you have questions about these failures, you should post to the `QIIME Forum <http://forum.qiime.org>`_.

==========================================
QIIME installation guide: Additional notes
==========================================

The following sections are referenced from the installation guide above.

.. _getting-qiime:

Getting QIIME
=============

First, change to the directory where you would like to download QIIME::

	cd $HOME

Stable Release
--------------
Currently the most stable version of QIIME is our |release| release, which you can download from `here <https://pypi.python.org/pypi/qiime>`_.

Latest Development Version
--------------------------
To get the latest development version of QIIME, you should check it out of our git repository, which is hosted on GitHub. While this code is subject to changes in interface and hasn't been as extensively tested as the release version, it will provide access to the latest and greatest QIIME features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the documentation in ``Qiime/doc``. Check out the latest version of QIIME using git with the command::

	git clone git://github.com/qiime/qiime.git Qiime

If you are using the latest development version of QIIME, you should periodically update your checkout by running the following command (from within your checkout)::

	git pull

Unpacking QIIME (release only)
------------------------------
After downloading the QIIME release tar file you'll need to unpack the code. For simplicity in this document, we will assume that you have downloaded QIIME to the directory ``$HOME/``.

Unpack the release .tar.gz file with the commands::

	cd $HOME
	tar -xvzf qiime-1.7.0.tar.gz
	ln -s $HOME/qiime-1.7.0 $HOME/Qiime

If you have downloaded the development version from GitHub, QIIME is already unpacked.

Installing QIIME
----------------
QIIME consists of library code (in ``Qiime/qiime``), test code (in ``Qiime/tests``), example script input and output (in ``Qiime/qiime_test_data``), documentation (in ``Qiime/doc``), and scripts (in ``Qiime/scripts``). Installing QIIME consists installing the library code in a place where python knows where to find it, and installing the scripts in a place where the shell looks for executable files, and running the tests (optional, but highly recommended).

.. _python-setup:

Installing with setup.py
------------------------

Using ``Qiime/setup.py`` (and thereby python's ``distutils`` package) is the recommended way of installing the Qiime library code and scripts. You can optionally specify where the library code and scripts should be installed -- depending on your setup, you may want to do this. By default, the QIIME library code will be placed under python's ``site-packages``, and the QIIME scripts will be place in ``/usr/local/bin/``. You may need to run ``setup.py`` using ``sudo`` if you do not have permission to place files in the default locations.

First, ensure that you are in the top-level QIIME directory::

	cd $HOME/Qiime

By default the QIIME scripts will be installed in ``/usr/local/bin``. As there are a lot of QIIME scripts, we highly recommend customizing the script directory to keep your system organized. This can be customized with the ``--install_scripts`` option. You also can specify an alternate directory for the library files with ``--install-purelib``. An example command is::

	python setup.py install --install-scripts=$HOME/bin/ --install-purelib=$HOME/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/release/2.7.1/install/index.html#alternate-installation>`_.

If you used default values for ``--install-scripts`` and ``--install-purelib`` (by not specifying them), your installation should be complete. If you specified an alternate value for ``--install-scripts``, you'll need to ensure that the shell knows where to look for the scripts. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::

	echo "export PATH=$HOME/bin/:$PATH" >> $HOME/.bashrc

If you specified an alternate value for ``--install-purelib``, you'll need to be sure that python knows where to look for Qiime. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::

	echo "export PYTHONPATH=$HOME/lib/:$PYTHONPATH" >> $HOME/.bashrc

The source your ``.bashrc``::

	source $HOME/.bashrc

.. _set-script-dir:

Finally, you'll need to create and edit a custom ``qiime_config`` file to tell QIIME where to look for the QIIME scripts. Create a custom ``qiime_config`` file by copying the default ``qiime_config`` packaged with QIIME::

	cp $HOME/Qiime/qiime/support_files/qiime_config $HOME/.qiime_config

Open the new file, ``$HOME/.qiime_config``, in a text editor such as TextEdit (on Mac), gedit (on Linux), vim, or emacs (but not Microsoft Word, which is a `word processor <http://en.wikipedia.org/wiki/Word_processor>`_, not a `text editor <http://en.wikipedia.org/wiki/Text_editor>`_!). Find the line beginning ``qiime_scripts_dir`` and add a tab, followed by the QIIME scripts directory. If you've used the default value (i.e., you didn't specify ``--install-scripts``) the value you add will be ``/usr/local/bin/``. Otherwise, specify the value that you provided for ``--install-scripts``. In the example above, this would look like::

	qiime_scripts_dir	$HOME/bin/

Note that the delimiter between the key and the value here is a tab, not a space! For additional information on the qiime_config file, `see this document <./qiime_config.html>`_.

.. _run-test-suite:

Running the test suite
----------------------
Next you should run the test suite. Execute the following commands::

	cd $HOME/Qiime/tests/
	python all_tests.py

You will see test output on the terminal indicating test successes and failures. Some failures are OK. The ``all_tests.py`` command will complete with a summary of test failures. Some tests may fail due to missing external applications -- these will be noted separately from other test failures. If these are related to features of QIIME that you are not using, this is acceptable. Otherwise, you'll need to ensure that you have the external applications installed correctly (and the correct versions), and re-run the tests.

License information for external dependencies
=============================================
We have attempted to provide accurate licensing information for the above dependencies for the convenience of our users. This information is by no means definitive and may contain errors. Any questions about licenses or the legality of specific uses of these software packages should be directed to the authors of the software. Do not rely solely on the license information presented above!

External application install notes
==================================

PATH Environment Variable
-------------------------

External applications used by QIIME need to be visible to the shell by existing in executable search path (i.e., listed in the ``$PATH`` environment variable). For example, if you plan to use cd-hit, and have the cd-hit executables installed in ``$HOME/bin`` you can add this directory to your system path with the commands::

	echo "export PATH=$HOME/bin/:$PATH" >> $HOME/.bashrc
	source $HOME/.bashrc

PYTHONPATH Environment Variable
-------------------------------

Qiime, PyCogent, and NumPy must be visible to python for all features of QIIME. matplotlib must be visible to python if you plan to use graphics features of QIIME; PyNAST must be visible to python if you plan to use PyNAST for multiple sequence alignment; and Denoiser must be visible to python if you plan to denoise 454 data. With the exception of Denoiser, all of these packages come with setup.py scripts. If you have used these, you should not need to modify your PYTHONPATH to make the library code visible. If you haven't used the respective setup.py scripts, or if you specified an alternate value for ``--install-purelib``, you may need to add the locations of these libraries to your PYTHONPATH environment variable.

For example, if you've installed PyNAST in ``$HOME/PyNAST`` you can add this to your PYTHONPATH with the commands::

	echo "export PYTHONPATH=$HOME/PyNAST/:$PYTHONPATH" >> $HOME/.bashrc
	source $HOME/.bashrc

.. _rdp-install:

RDP_JAR_PATH Environment Variable
---------------------------------

If you plan to use the RDP classifier for taxonomy assignment you must define an ``RDP_JAR_PATH`` environment variable. If you downloaded and unzipped the RDP classifier folder in ``$HOME/app/``, you can do this with the following commands::

	echo "export RDP_JAR_PATH=$HOME/app/rdp_classifier_2.2/rdp_classifier-2.2.jar" >> $HOME/.bashrc
	source $HOME/.bashrc

Note that you will need the contents inside ``rdp_classifier_2.2`` for the program to function properly.

.. _uclust-install:

uclust Install Notes
--------------------

The uclust binary must be called ``uclust``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``uclust1.2.21q_i86linux64`` as ``$HOME/bin/uclust1.2.21q_i86linux64``, we recommend creating a symbolic link to this file::

	ln -s $HOME/bin/uclust1.2.21q_i86linux64 $HOME/bin/uclust

.. _usearch-install:

usearch Install Notes
---------------------

The usearch binary must be called ``usearch``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``usearch5.2.236_i86linux32`` as ``$HOME/bin/usearch5.2.236_i86linux32``, we recommend creating a symbolic link to this file::

	ln -s $HOME/bin/usearch5.2.236_i86linux32 $HOME/bin/usearch

.. _chimeraslayer-install:

ChimeraSlayer Install Notes
---------------------------

ChimeraSlayer can only be run from the directory where it was unpacked and built as it depends on several of its dependencies being in specific places relative to the executable (``ChimeraSlayer/ChimeraSlayer.pl``). Carefully follow the ChimeraSlayer install instructions. Then add the directory containing ``ChimeraSlayer.pl`` to your ``$PATH`` environment variable. If your ``ChimeraSlayer`` folder is in ``$HOME/app/`` you can set the ``$PATH`` environment variable as follows::

	echo "export PATH=$HOME/app/ChimeraSlayer:$PATH" >> $HOME/.bashrc
	source $HOME/.bashrc

If you're having trouble getting ChimeraSlayer to work via QIIME, you should first check to see if you can run it directly from a directory other than its install directory. For example, try running ``ChimeraSlayer.pl`` from your home directory.

Once you have configured Qiime, you can test your ChimeraSlayer install by running::

	print_qiime_config.py -t

This includes a check for obvious problems with your ChimeraSlayer install, and should help you determine if you have it installed correctly.

.. _R-install:

R Install Notes
---------------

To install R visit http://www.r-project.org/ and follow the install instructions. Once R is installed, run R and excecute the following commands::

	install.packages('randomForest')
	install.packages('optparse')
	install.packages('vegan')
	install.packages('ape')
	install.packages('MASS')
	install.packages('gtools')
	install.packages('klaR')
	install.packages('RColorBrewer')
	q()

.. _ampliconnoise-install:

AmpliconNoise Install Notes
---------------------------

AmpliconNoise requires that several environment variables are set. After you've installed AmpliconNoise, you can set these with the following commands (assuming your AmpliconNoise install directory is ``$HOME/AmpliconNoiseV1.27/``)::

	echo "export PATH=$HOME/AmpliconNoiseV1.27/Scripts:$HOME/AmpliconNoiseV1.27/bin:$PATH" >> $HOME/.bashrc

	echo "export PYRO_LOOKUP_FILE=$HOME/AmpliconNoiseV1.27/Data/LookUp_E123.dat" >> $HOME/.bashrc
	echo "export SEQ_LOOKUP_FILE=$HOME/AmpliconNoiseV1.27/Data/Tran.dat" >> $HOME/.bashrc

QIIME Denoiser Install Notes
----------------------------

If you do not install QIIME using ``setup.py`` and you plan to use the QIIME Denoiser, you'll need to compile the FlowgramAlignment program. To do this you'll need to have ``ghc`` installed. Then from the ``Qiime/qiime/support_files/denoiser/FlowgramAlignment/`` directory, run the following command::

	make ; make install

.. _build-qiime-docs:

Building The QIIME Documentation
================================

If you are using the development version of QIIME, you may want to build the documentation locally for access to the latest version. You can change to the ``Qiime/doc`` directory and run::

	make html

We try to update the documentation as we update the code, but development version users may notice some discrepancies. After building the documentation, you can view it in a web browser by opening the file ``Qiime/doc/_build/html/index.html``. You may want to bookmark that page for easy access.
