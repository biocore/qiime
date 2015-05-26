.. _doc_install:
.. QIIME documentation master file, created by Jesse Stombaugh
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. index:: Installing QIIME

=================
Installing QIIME
=================

.. toctree::
   :maxdepth: 2


To use QIIME, the user should install the following programs and set the appropriate environment variables.

Software Dependencies
---------------------
The following programs and datasets were used to generate this tutorial and it is recommended to use the same versions as shown below.

* Python 2.6 - Download `Python-2.6.4.tgz <http://www.python.org/ftp/python/2.6.4/Python-2.6.4.tgz>`_
* PyCogent 1.4.0 - Download `PyCogent-1.4.tgz  <http://sourceforge.net/projects/pycogent/files/PyCogent/1.4/PyCogent-1.4.tgz/download>`_
* PyNAST 1.0 - Download `PyNAST-1.0-tar.gz  <http://sourceforge.net/projects/pynast/files/PyNAST%20releases/PyNAST-1.0.tar.gz/download>`_
* Numpy 1.3.0 - Download `numpy-1.3.0-tar.gz  <http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download>`_
* MatPlotLib 0.98.5.2 - Download `matplotlib-0.98.5.2.tar.gz  <http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.2.tar.gz>`_
* jre1.6.0_05 - Download `JavaForMacOSX10.5Update1.dmg <http://wsidecar.apple.com/cgi-bin/nph-reg3rdpty2.pl/product=18844&cat=59&platform=osx&method=sa/JavaForMacOSX10.5Update1.dmg>`_
* rdp_classifier-2.0.1 - Download `rdp_classifier_2.0.1.tar.gz <http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.0.1/rdp_classifier_2.0.1.tar.gz>`_
* greengenes core set data file - Download `core_set_aligned.fasta.imputed <http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta.imputed>`_
* greengenes alignment lanemask file - Download `lanemask_in_1s_and_0s <http://greengenes.lbl.gov/Download/Sequence_Data/lanemask_in_1s_and_0s>`_
* blast-2.2.22 (please note that this refers to legacy BLAST from NCBI, not BLAST+) - Download `blast-2.2.22-universal-macosx.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/blast-2.2.22-universal-macosx.tar.gz>`_
* fasttree 2.1.0 - Download `FastTree-2.1.c <http://www.microbesonline.org/fasttree/FastTree-2.1.0.c>`_
* cd-hit 3.1 - Download `cd-hit-2007-0131.tar.gz <http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131.tar.gz>`_

Installing QIIME
----------------

Stable Pre-Release
^^^^^^^^^^^^^^^^^^
Currently the most stable version of QIIME is our 0.92 pre-release, which you can download from `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-0.92.tar.gz/download>`_.

Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^
To get the latest development version of QIIME, you should check it out of our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in "Qiime/doc". Check out the latest version of QIIME using svn with the command::

	svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

The user can update QIIME, by using the following command::

	svn update /path/to/QIIME/

Environment Variables
---------------------
Make sure the following environment variable are set:

PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^

Your $PATH environment variable should contain the following:

* /path/to/python/bin
* /path/to/cd-hit-est
* /path/to/FastTree
* /path/to/blast-2.2.21/bin/
* /path/to/PyNAST/scripts/
* /path/to/jre1.6_0_16/

For all path description throughout this tutorial, the "/path/to/" refers to the physical location of each program on your local computer. For instance, the "/path/to/python/bin/" refers to "/Library/Frameworks/Python.framework/Versions/2.5/bin" on Mac OS X version 10.5.

In the bash shell, you can use the following command (example only shows the path to python, so other softwares can be added to the PATH, using a similar approach): ::

	export PATH=/path/to/python/bin/:$PATH

PYTHONPATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your PYTHONPATH should contain the following:

* /path/to/PyCogent
* /path/to/QIIME
* /path/to/PyNAST

In the bash shell, you can use the following command (example only shows the path to QIIME, so other softwares can be added to the PYTHONPATH, using a similar approach): ::

	export PYTHONPATH=/path/to/QIIME/:$PYTHONPATH

RDP_JAR_PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user should also define an RDP_JAR_PATH variable, since this tutorial uses the RDP Classifier:

* /path/to/rdp_classifier-2.0.jar

In the bash shell, you can use the following command: ::

	export RDP_JAR_PATH=/path/to/rdp_classifier-2.0.jar

Testing QIIME Install
---------------------
Once the source code is downloaded, the user should test QIIME to be sure all essential software is properly installed and the correct environment variables are set.

In a terminal window the user, should cd to their qiime/test directory using the following command: ::

	cd /path/to/QIIME/tests/

Then run the following test command: ::

	python all_tests.py -v

If all tests run properly, then QIIME was properly installed. Some test scripts may fail, due to optional third party applications not being installed. Test failures which contain the text ApplicationNotFoundError can be safely ignored. The following scripts may generate errors:

1. test_align_seqs.py - will fail if PyNAST, blast, muscle or infernal are not installed. PyNAST requires BLAST, and is the default sequence aligner. If you are not using PyNAST, you can safely ignore fails regarding blast and/or PyNAST. If you are not using muscle or infernal, you can safely disregard those error messages.
2. test_pick_otus.py – will fail if cd-hit and/or blast are not installed. If you are not planning to use these methods for OTU picking, you do not need to worry about these failures.
3. test_process_sff.py - will fail if you do not have sffinfo installed.
4. test_pycogent_backports/test_mothur.py - will fail is mothur is not installed.
5. test_pycogent_backports/test_uclust.py - will fail is uclust is not installed.
6. test_pyronoise.py - will fail in PyroNoise is not installed.

