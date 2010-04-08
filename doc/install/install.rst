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

* Python 2.6 (`src <http://www.python.org/ftp/python/2.6.4/Python-2.6.4.tgz>`_)
* PyCogent 1.4.1 (`src <http://sourceforge.net/projects/pycogent/files/PyCogent/1.4.1/PyCogent-1.4.1.tgz/download>`_)
* Numpy 1.3.0 (`src <http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download>`_)

Dependencies required for a subset of QIIME's features
------------------------------------------------------


PyNAST alignment, tree-building, taxonomy assignment, OTU picking, and other data generation steps:

* uclust 1.1.579 (`src and binaries <http://www.drive5.com/uclust/downloads1_1_579.html>`_) See :ref:`uclust install notes <uclust-install>`.
* PyNAST 1.1 (`src  <http://sourceforge.net/projects/pynast/files/PyNAST%20releases/PyNAST-1.1.tgz/download>`_)
* greengenes core set data file (`fasta <http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta.imputed>`_)
* greengenes alignment lanemask file (`txt <http://greengenes.lbl.gov/Download/Sequence_Data/lanemask_in_1s_and_0s>`_)
* blast-2.2.22 (legacy BLAST from NCBI, *NOT* BLAST+) (`OS X <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-universal-macosx.tar.gz>`_ or `linux 32-bit <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-ia32-linux.tar.gz>`_)
* fasttree 2.1.0 (`src <http://www.microbesonline.org/fasttree/FastTree-2.1.0.c>`_)
* cd-hit 3.1 (`src <http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131.tar.gz>`_)
* jre1.6.0_05 (`link <http://java.sun.com/javase/downloads/index.jsp>`_)
* rdp_classifier-2.0.1 (`src <http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.0.1/rdp_classifier_2.0.1.tar.gz>`_) See :ref:`RDP install notes <rdp-install>`.

Denoising 454 data:

* Denoiser 0.83 (`src <http://www.microbio.me/denoiser/>`_ -- remember to add the top-level ``Denoiser_0.83`` directory to your ``$PYTHONPATH``)

Graphics and other data analysis steps:

* MatPlotLib 0.98.5.2 (`src  <http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.2.tar.gz>`_)

If you plan to build the QIIME documentation locally:

* Sphinx 0.6.3 (`src <http://pypi.python.org/pypi/Sphinx>`_) See :ref:`Building the QIIME documentation <build-qiime-docs>`


Shortcuts in this document
--------------------------
For simplicity throughout this document, we assume that you have downloaded QIIME in ``/home/qiime_user/``. You should consider all occurrences of ``/home/qiime_user/`` in the remainder of this document as references to the directory which contains the QIIME directory which you'll have after downloading and unpacking QIIME.

Getting QIIME
----------------
First, change to the directory where you would like to download QIIME::

	cd /home/qiime_user

Stable Pre-Release
^^^^^^^^^^^^^^^^^^
Currently the most stable version of QIIME is our 1.0.0 release, which you can download from `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-1.0.0.tar.gz/download>`_.

Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^
To get the latest development version of QIIME, you should check it out of our Sourceforge repository. While this code is subject to changes in interface and hasn't been as extensively tested as the release version, it will provide access to the latest and greatest QIIME features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in ``Qiime/doc``. Check out the latest version of QIIME using svn with the commands::

	svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

svn users should periodically update QIIME by using the following command::

	svn update /home/qiime_user/Qiime/


Unpacking QIIME (release only)
---------------------------------------
After downloading the QIIME release tar file you'll need to unpack the code. For simplicity in this document, we will assume that you have downloaded QIIME to the directory ``/home/qiime_user/``. 

Unpack the release Qiime tar file with the commands::

	cd /home/qiime_user
	tar -xvzf Qiime-1.0.0.tar.gz
	ln -s /home/qiime_user/Qiime-1.0.0 /home/qiime_user/Qiime
	
If you have downloaded from svn, QIIME is already unpacked.
	
Installing QIIME
----------------
QIIME consists of library code (in ``Qiime/qiime``), test code (in ``Qiime/tests``), documentation (in ``Qiime/doc``), and scripts (in ``Qiime/scripts``). Installing QIIME consists of running the tests (optional, but highly recommend), installing the library code in a place where python knows where to find it, and installing the scripts in a place where the shell looks for executable files.



Installing the library code and scripts with setup.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using ``Qiime/setup.py`` (and thereby python's ``distutils`` package) is the recommended way of installing the Qiime library code and scripts. You can optionally specify where the library code and scripts should be installed -- depending on your setup, you may want to do this. By default, the QIIME library code will be placed under python's ``site-packages``, and the QIIME scripts will be place in ``/usr/local/bin/``. You may need to run ``setup.py`` using ``sudo`` if you do not have permission to place files in the default locations. 

First, ensure that you are in the top-level QIIME directory::
	
	cd /home/qiime_user/Qiime

By default the QIIME scripts will be installed in ``/usr/local/bin``. As there are a lot of QIIME scripts, we recommend customizing the script directory to keep your system organized. This can be customized with the ``--install_scripts`` option::
	
	python setup.py install --install-scripts=/home/qiime_user/bin/
	
You can similarly install the library code in an alternate location using the ``--install-purelib`` option::
	
	python setup.py install --install-purelib=/home/qiime_user/lib/


Combine these options as follows::
	
	python setup.py install --install-scripts=/home/qiime_user/bin/ --install-purelib=/home/qiime_user/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/install/index.html#alternate-installation-the-home-scheme>`_.

If you used default values for ``--install-scripts`` and ``--install-purelib`` (by not specifying them), your installation should be complete. If you specified an alternate value for ``--install-scripts``, you'll need to ensure that the shell knows where to look for the scripts. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PATH=/home/qiime_user/bin/:$PATH" >> /home/qiime_user/.bashrc

If you specified an alternate value for ``--install-purelib``, you'll need to be sure that python knows where to look for Qiime. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PYTHONPATH=/home/qiime_user/lib/:$PYTHONPATH" >> /home/qiime_user/.bashrc
	
The source your ``.bashrc``::

	source /home/qiime_user/.bashrc

.. _set-script-dir:

Finally, you'll need to create and edit a custom ``qiime_config`` file to tell QIIME where to look for the QIIME scripts. Create a custom ``qiime_config`` file by copying the default ``qiime_config`` packaged with QIIME::

	cp /home/qiime_user/Qiime/qiime/support_files/qiime_config /home/qiime_user/.qiime_config
	
Open the new file, ``/home/qiime_user/.qiime_config``, in a text editor such as TextEdit (on Mac), gedit (on Linux), vim, or emacs (but not Microsoft Word, which is a `word processor <http://en.wikipedia.org/wiki/Word_processor>`_, not a `text editor <http://en.wikipedia.org/wiki/Text_editor>`_!). Find the line beginning ``qiime_scripts_dir`` and add a tab, followed by the QIIME scripts directory. If you've used the default value (i.e., you didn't specify ``--install-scripts``) the value you add will be ``/usr/local/bin/``. Otherwise, specify the value that you provided for ``--install-scripts``. In the example above, this would look like::

	qiime_scripts_dir	/home/qiime_user/bin/
	
Note that the delimiter between the key and the value here is a tab, not a space! For additional information on the qiime_config file, `see this document <./qiime_config.html>`_.

Running the test suite
----------------------
Next you should run the test suite. Execute the following commands::
	
	cd /home/qiime_user/Qiime/tests/
	python all_tests.py

You will see test output on the terminal indicating test successes and failures. Some failures are OK. The ``all_tests.py`` command will complete with a summary of test failures. Some tests may fail due to missing external applications -- these will be noted separately from other test failures. If these are related to features of QIIME that you are not using, this is acceptable. Otherwise, you'll need to ensure that you have the external applications installed correctly (and the correct versions), and re-run the tests. 

Testing your QIIME installation
-------------------------------
If QIIME is installed correctly, you should be able to run the QIIME scripts. Try the following::
	
	cd
	align_seqs.py -h
	
This should give you help text describing the interface to the align_seqs.py script. (Note that if you do not have a /home/qiime_user/.bashrc you may get an error at the ``source`` step. If you did not specify alternate values for ``--install-purelib`` or ``--install-scripts`` this shouldn't be a problem.)

External application install notes
----------------------------------

PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^

External applications used by QIIME need to be visible to the shell by existing in executable search path (i.e., listed in the ``$PATH`` environment variable). For example, if you plan to use cd-hit, and have the cd-hit executables installed in ``/home/qiime_user/bin`` you can add this directory to your system path with the commands::
	
	echo "export PATH=/home/qiime_user/bin/:$PATH" >> /home/qiime_user/.bashrc
	source /home/qiime_user/.bashrc

PYTHONPATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Qiime, PyCogent, and NumPy must be visible to python for all features of QIIME. matplotlib must be visible to python if you plan to use graphics features of QIIME; PyNAST must be visible to python if you plan to use PyNAST for multiple sequence alignment; and Denoiser must be visible to python if you plan to denoise 454 data. With the exception of Denoiser, all of these packages come with setup.py scripts. If you have used these, you should not need to modify your PYTHONPATH to make the library code visible. If you haven't used the respective setup.py scripts, or if you specified an alternate value for ``--install-purelib``, you may need to add the locations of these libraries to your PYTHONPATH environment variable. 

For example, if you've installed PyNAST in ``/home/qiime_user/PyNAST`` you can add this to your PYTHONPATH with the commands::
	
	echo "export PYTHONPATH=/home/qiime_user/PyNAST/:$PYTHONPATH" >> /home/qiime_user/.bashrc
	source /home/qiime_user/.bashrc


RDP_JAR_PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _rdp-install:

If you plan to use the RDP classifier for taxonomy assignment you must also define an RDP_JAR_PATH variable. If you have the RDP classifier jar file (``rdp_classifier-2.0.1.jar``) in ``/home/qiime_user/app`` you can do this with the following command::

	echo "export RDP_JAR_PATH=/home/qiime_user/app/rdp_classifier-2.0.1.jar" >> /home/qiime_user/.bashrc
	
uclust Install Notes
^^^^^^^^^^^^^^^^^^^^^^^

.. _uclust-install:

The uclust binary must be called ``uclust``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``uclust1.1.579_i86linux32`` as ``/home/qiime_user/bin/uclust1.1.579_i86linux32``, we recommend creating a symbolic link to this file::
	
	ln -s /home/qiime_user/bin/uclust1.1.579_i86linux32 /home/qiime_user/bin/uclust
	
Building The QIIME Documentation
---------------------------------

.. _build-qiime-docs:

If you are using the svn version of QIIME, you may want to build the documentation locally for access to the latest version. You can change to the ``Qiime/doc`` directory and run::

	make html
	
We try to update the documentation as we update the code, but svn users may notice some discrepancies. After building the documentation, you can view it in a web browser by opening the file ``Qiime/doc/_build/html/index.html``. You may want to bookmark that page for easy access. 
