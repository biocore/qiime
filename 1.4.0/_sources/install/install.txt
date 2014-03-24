.. _doc_install:
.. QIIME documentation master file, created by Jesse Stombaugh
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. index:: Installing QIIME

===========================
Installing QIIME natively
===========================
QIIME consists of native code and additionally wraps many external applications. This gives the user flexibility to easily build their own analysis pipelines, making use of popular microbial community analysis tools. QIIME handles the processing of input and output of these applications, so the user can spend time analyzing their data rather than parsing, writing, and converting file formats. 

As a consequence of this 'pipeline' architecture, depending on the features of QIIME that you plan to use, you may or may not need all of QIIME dependencies. Getting all of these applications working correctly can be difficult, which is why we distribute the QIIME virtual box. If you are a beginner user, or just testing QIIME to see if it will meet your needs, you may want to begin with the `virtual box <./virtual_box.html>`_ which will greatly ease installation.

The following programs and datasets were used to generate the QIIME tutorial. QIIME has been tested with the same versions listed below. Other versions are not guaranteed to work. Download links are provided for convenience, but these are subject to change. We'll do our best to keep these up-to-date.

You should follow instructions provided by the package developers to install the dependencies.

Dependencies required for all features of QIIME
-----------------------------------------------

To install most of following dependencies you need to have a build environment on your machine. On OS X, this involves installing the `developer tools <http://developer.apple.com/technologies/xcode.html>`_. On Debian-based Linux (e.g., Ubuntu), this involves installing the ``build-essential`` package::

	sudo apt-get install build-essential

The following are required by QIIME:

* Python 2.7.1 (`src <http://www.python.org/ftp/python/2.7.1/Python-2.7.1.tgz>`_) (license: PSF)
* PyCogent 1.5.1 (`src <http://sourceforge.net/projects/pycogent/files/PyCogent/1.5.1/PyCogent-1.5.1.tgz/download>`_) (license: GPL)
* Numpy 1.5.1 (`src <http://sourceforge.net/projects/numpy/files/NumPy/1.5.1/numpy-1.5.1.tar.gz/download>`_) (license: BSD)

Dependencies required for a subset of QIIME's features
------------------------------------------------------


PyNAST alignment, tree-building, taxonomy assignment, OTU picking, and other data generation steps (required in default pipeline):

* uclust 1.2.22q (`binaries <http://www.drive5.com/uclust/downloads1_2_22q.html>`_) See :ref:`uclust install notes <uclust-install>`. (licensed specially for Qiime and PyNAST users)
* PyNAST 1.1 (`src  <http://sourceforge.net/projects/pynast/files/PyNAST%20releases/PyNAST-1.1.tgz/download>`_) (license: GPL)
* greengenes core set data file (`fasta <http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta.imputed>`_)
* greengenes alignment lanemask file (`txt <http://greengenes.lbl.gov/Download/Sequence_Data/lanemask_in_1s_and_0s>`_)
* fasttree 2.1.3 (`src <http://www.microbesonline.org/fasttree/FastTree-2.1.3.c>`_) (license: GPL)
* jre1.6.0_05 (`link <http://java.sun.com/javase/downloads/index.jsp>`_) (license: GPL2)
* rdp_classifier-2.2 (`zip <http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.2.zip/download>`_) See :ref:`RDP install notes <rdp-install>`. (license: GPL)

Alignment, tree-building, taxonomy assignment, OTU picking, and other data generation steps (required for alternative pipelines):

* blast-2.2.22 (legacy BLAST from NCBI, *NOT* BLAST+) (`OS X <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-universal-macosx.tar.gz>`_ or `linux 32-bit <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-ia32-linux.tar.gz>`_) (license: GNU)
* cd-hit 3.1.1 (`src <http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131.tar.gz>`_) (license: Free access)
* ChimeraSlayer (via microbiomeutil_2010-04-29) (`src <http://sourceforge.net/projects/microbiomeutil/files/>`_) See :ref:`ChimeraSlayer install notes <chimeraslayer-install>`.
* mothur v.1.6.0 (`web <http://www.mothur.org/w/images/e/e8/Mothur.1.6.0.zip>`_) (license: GPL)
* clearcut v1.0.9 (`src <http://www.mothur.org/w/images/9/91/Clearcut.source.zip>`_)
* raxml v7.0.3 (`src <http://wwwkramer.in.tum.de/exelixis/r703-source.php>`_)
* infernal 1.0.2 (`src <ftp://selab.janelia.org/pub/software/infernal/infernal.tar.gz>`_) (license: GPL)
* cdbtools (`src <ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/cdbfasta/cdbfasta.tar.gz>`_)
* muscle 3.8.31 (`link <http://www.drive5.com/muscle/downloads.htm>`_) (Public domain)

Processing sff files:

* sfffile and sffinfo (optional, QIIME 1.2.0 and later contain built-in tools for processing sff files although they are about 10x slower than the tools from Roche) (license: proprietary - must be obtained from Roche/454)

Denoising 454 data:

* GNU Science Library (required by AmpliconNoise) (`src <ftp://ftp.gnu.org/gnu/gsl/gsl-1.9.tar.gz>`_)
* AmpliconNoise 1.25 (`src <http://ampliconnoise.googlecode.com/files/AmpliconNoiseV1.25.tar.gz>`_) See :ref:`AmpliconNoise install notes <ampliconnoise-install>`.
* ghc 6.8 (required by the QIIME denoiser) (`src <http://haskell.org/ghc>`_)


Visualization and plotting steps:

* MatPlotLib 1.1.0 (`src <http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.1.0/matplotlib-1.1.0.tar.gz>`_) (license: PFS)
* cytoscape v2.7.0 (`web <http://www.cytoscape.org/>`_) (license: LGPL)

Supervised learning (``supervised_learning.py``):

* R 2.12.0 (`source <http://www.r-project.org/>`_) See :ref:`R install notes <R-install>`. (license: GPL2)

Assigning taxonomy using BLAST or picking OTUs against Greengenes filtered at 97% identity:

* Greengenes 97% OTUs, taxonomies, and tree (`zip <http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010.zip>`_)

If you plan to build the QIIME documentation locally:

* Sphinx 1.0.4 (`src <http://pypi.python.org/pypi/Sphinx>`_) See :ref:`Building the QIIME documentation <build-qiime-docs>` (license: BSD)

License information for external dependencies
---------------------------------------------
We have attempted to provide accurate licensing information for the above dependencies for the convenience of our users. This information is by no means definitive and may contain errors. Any questions about licenses or the legality of specific uses of these software packages should be directed to the authors of the software. Do not rely solely on the license information presented above!

Shortcuts in this document
--------------------------
For simplicity throughout this document, we assume that you have downloaded QIIME in ``/home/qiime/``. You should consider all occurrences of ``/home/qiime/`` in the remainder of this document as references to the directory which contains the QIIME directory which you'll have after downloading and unpacking QIIME.

Getting QIIME
----------------
First, change to the directory where you would like to download QIIME::

	cd /home/qiime_user

Stable Pre-Release
^^^^^^^^^^^^^^^^^^
Currently the most stable version of QIIME is our |release| release, which you can download from `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-1.4.0.tar.gz/download>`_.

Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^
To get the latest development version of QIIME, you should check it out of our Sourceforge repository. While this code is subject to changes in interface and hasn't been as extensively tested as the release version, it will provide access to the latest and greatest QIIME features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in ``Qiime/doc``. Check out the latest version of QIIME using svn with the commands::

	svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

svn users should periodically update QIIME by using the following command::

	svn update /home/qiime/Qiime/


Unpacking QIIME (release only)
---------------------------------------
After downloading the QIIME release tar file you'll need to unpack the code. For simplicity in this document, we will assume that you have downloaded QIIME to the directory ``/home/qiime/``. 

Unpack the release Qiime tar file with the commands::

	cd /home/qiime_user
	tar -xvzf Qiime-1.4.0.tar.gz
	ln -s /home/qiime/Qiime-1.4.0 /home/qiime/Qiime
	
If you have downloaded from svn, QIIME is already unpacked.
	
Installing QIIME
----------------
QIIME consists of library code (in ``Qiime/qiime``), test code (in ``Qiime/tests``), documentation (in ``Qiime/doc``), and scripts (in ``Qiime/scripts``). Installing QIIME consists of running the tests (optional, but highly recommend), installing the library code in a place where python knows where to find it, and installing the scripts in a place where the shell looks for executable files.



Installing the library code and scripts with setup.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using ``Qiime/setup.py`` (and thereby python's ``distutils`` package) is the recommended way of installing the Qiime library code and scripts. You can optionally specify where the library code and scripts should be installed -- depending on your setup, you may want to do this. By default, the QIIME library code will be placed under python's ``site-packages``, and the QIIME scripts will be place in ``/usr/local/bin/``. You may need to run ``setup.py`` using ``sudo`` if you do not have permission to place files in the default locations. 

First, ensure that you are in the top-level QIIME directory::
	
	cd /home/qiime/Qiime

By default the QIIME scripts will be installed in ``/usr/local/bin``. As there are a lot of QIIME scripts, we highly recommend customizing the script directory to keep your system organized. This can be customized with the ``--install_scripts`` option. You also can specify and alternate directory for the library files with ``--install-purelib``, but if you do so you must also specify ``--install-data`` as the same directory. Failure to do this will result in a broken QIIME install. An example command is::
	
	python setup.py install --install-scripts=/home/qiime/bin/ --install-purelib=/home/qiime/lib/ --install-data=/home/qiime/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/release/2.6.6/install/index.html#alternate-installation-the-home-scheme>`_.

If you used default values for ``--install-scripts`` and ``--install-purelib`` (by not specifying them), your installation should be complete. If you specified an alternate value for ``--install-scripts``, you'll need to ensure that the shell knows where to look for the scripts. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PATH=/home/qiime/bin/:$PATH" >> /home/qiime/.bashrc

If you specified an alternate value for ``--install-purelib``, you'll need to be sure that python knows where to look for Qiime. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PYTHONPATH=/home/qiime/lib/:$PYTHONPATH" >> /home/qiime/.bashrc
	
The source your ``.bashrc``::

	source /home/qiime/.bashrc

.. _set-script-dir:

Finally, you'll need to create and edit a custom ``qiime_config`` file to tell QIIME where to look for the QIIME scripts. Create a custom ``qiime_config`` file by copying the default ``qiime_config`` packaged with QIIME::

	cp /home/qiime/Qiime/qiime/support_files/qiime_config /home/qiime/.qiime_config
	
Open the new file, ``/home/qiime/.qiime_config``, in a text editor such as TextEdit (on Mac), gedit (on Linux), vim, or emacs (but not Microsoft Word, which is a `word processor <http://en.wikipedia.org/wiki/Word_processor>`_, not a `text editor <http://en.wikipedia.org/wiki/Text_editor>`_!). Find the line beginning ``qiime_scripts_dir`` and add a tab, followed by the QIIME scripts directory. If you've used the default value (i.e., you didn't specify ``--install-scripts``) the value you add will be ``/usr/local/bin/``. Otherwise, specify the value that you provided for ``--install-scripts``. In the example above, this would look like::

	qiime_scripts_dir	/home/qiime/bin/
	
Note that the delimiter between the key and the value here is a tab, not a space! For additional information on the qiime_config file, `see this document <./qiime_config.html>`_.

Running the test suite
----------------------
Next you should run the test suite. Execute the following commands::
	
	cd /home/qiime/Qiime/tests/
	python all_tests.py

You will see test output on the terminal indicating test successes and failures. Some failures are OK. The ``all_tests.py`` command will complete with a summary of test failures. Some tests may fail due to missing external applications -- these will be noted separately from other test failures. If these are related to features of QIIME that you are not using, this is acceptable. Otherwise, you'll need to ensure that you have the external applications installed correctly (and the correct versions), and re-run the tests. 

Testing your QIIME installation
-------------------------------
If QIIME is installed correctly, you should be able to run the QIIME scripts. Try the following::
	
	cd
	align_seqs.py -h
	
This should give you help text describing the interface to the align_seqs.py script. (Note that if you do not have a /home/qiime/.bashrc you may get an error at the ``source`` step. If you did not specify alternate values for ``--install-purelib`` or ``--install-scripts`` this shouldn't be a problem.)

External application install notes
----------------------------------

PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^

External applications used by QIIME need to be visible to the shell by existing in executable search path (i.e., listed in the ``$PATH`` environment variable). For example, if you plan to use cd-hit, and have the cd-hit executables installed in ``/home/qiime/bin`` you can add this directory to your system path with the commands::
	
	echo "export PATH=/home/qiime/bin/:$PATH" >> /home/qiime/.bashrc
	source /home/qiime/.bashrc

PYTHONPATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Qiime, PyCogent, and NumPy must be visible to python for all features of QIIME. matplotlib must be visible to python if you plan to use graphics features of QIIME; PyNAST must be visible to python if you plan to use PyNAST for multiple sequence alignment; and Denoiser must be visible to python if you plan to denoise 454 data. With the exception of Denoiser, all of these packages come with setup.py scripts. If you have used these, you should not need to modify your PYTHONPATH to make the library code visible. If you haven't used the respective setup.py scripts, or if you specified an alternate value for ``--install-purelib``, you may need to add the locations of these libraries to your PYTHONPATH environment variable. 

For example, if you've installed PyNAST in ``/home/qiime/PyNAST`` you can add this to your PYTHONPATH with the commands::
	
	echo "export PYTHONPATH=/home/qiime/PyNAST/:$PYTHONPATH" >> /home/qiime/.bashrc
	source /home/qiime/.bashrc


RDP_JAR_PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _rdp-install:

If you plan to use the RDP classifier for taxonomy assignment you must also define an RDP_JAR_PATH variable. If you have the RDP classifier jar file (``rdp_classifier-2.0.1.jar``) in ``/home/qiime/app`` you can do this with the following command::

	echo "export RDP_JAR_PATH=/home/qiime/app/rdp_classifier-2.0.1.jar" >> /home/qiime/.bashrc
	source /home/qiime/.bashrc
	
uclust Install Notes
^^^^^^^^^^^^^^^^^^^^^^^

.. _uclust-install:

The uclust binary must be called ``uclust``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``uclust1.2.21q_i86linux64`` as ``/home/qiime/bin/uclust1.2.21q_i86linux64``, we recommend creating a symbolic link to this file::
	
	ln -s /home/qiime/bin/uclust1.2.21q_i86linux64 /home/qiime/bin/uclust
	
ChimeraSlayer Install Notes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _chimeraslayer-install:

ChimeraSlayer can only be run from the directory where it was unpacked and built as it depends on several of its dependencies being in specific places relative to the executable (``ChimeraSlayer/ChimeraSlayer.pl``). Carefully follow the ChimeraSlayer install instructions. Then add the directory containing ``ChimeraSlayer.pl`` to your ``$PATH`` environment variable. If your ``ChimeraSlayer`` folder is in ``/home/qiime/app/`` you can set the ``$PATH`` environment variable as follows::

	echo "export PATH=/home/qiime/app/ChimeraSlayer:$PATH" >> /home/qiime/.bashrc
	source /home/qiime/.bashrc

If you're having trouble getting ChimeraSlayer to work via QIIME, you should first check to see if you can run it directly from a directory other than its install directory. For example, try running ``ChimeraSlayer.pl`` from your home directory.

Once you have configured Qiime, you can test your ChimeraSlayer install by running::

	print_qiime_config.py -t
	
This includes a check for obvious problems with your ChimeraSlayer install, and should help you determine if you have it installed correctly.

R Install Notes
^^^^^^^^^^^^^^^

.. _R-install:

To install R visit http://www.r-project.org/ and follow the install instructions. Once R is installed, run R and excecute the command::

	install.packages('randomForest')
	q()

AmpliconNoise Install Notes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ampliconnoise-install:

AmpliconNoise requires that several environment variables are set. After you've installed AmpliconNoise, you can set these with the following commands (assuming your AmpliconNoise install directory is ``/home/qiime/AmpliconNoiseV1.24/``)::

	echo "export PATH=/home/qiime/AmpliconNoiseV1.24/Scripts:/home/qiime/AmpliconNoiseV1.24/bin:$PATH" >> /home/qiime/.bashrc
	
	echo "export PYRO_LOOKUP_FILE=/home/qiime/AmpliconNoiseV1.24/Data/LookUp_E123.dat" >> /home/qiime/.bashrc
	echo "export SEQ_LOOKUP_FILE=/home/qiime/AmpliconNoiseV1.24/Data/Tran.dat" >> /home/qiime/.bashrc

QIIME Denoiser Install Notes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you do not install QIIME using ``setup.py`` and you plan to use the QIIME Denoiser, you'll need to compile the FlowgramAlignment program. To do this you'll need to have ``ghc`` installed. Then from the ``Qiime/qiime/support_files/denoiser/FlowgramAlignment/`` directory, run the following command::

	make ; make install


Building The QIIME Documentation
---------------------------------

.. _build-qiime-docs:

If you are using the svn version of QIIME, you may want to build the documentation locally for access to the latest version. You can change to the ``Qiime/doc`` directory and run::

	make html
	
We try to update the documentation as we update the code, but svn users may notice some discrepancies. After building the documentation, you can view it in a web browser by opening the file ``Qiime/doc/_build/html/index.html``. You may want to bookmark that page for easy access. 
