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

* Python 2.6
* PyCogent 1.4.0 or later with ReportLab
* PyNAST 1.0
* Numpy 1.3
* MatPlotLib 0.98.5.2
* jre1.6.0_16
* rdp_classifier-2.0
* greengenes core set data file (http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/core_set_aligned.fasta.imputed)
* greengenes alignment lanemask file (http://greengenes.lbl.gov/Download/Sequence_Data/lanemask_in_1s_and_0s)
* blast-2.2.22 (please note that this refers to legacy BLAST from NCBI, not BLAST+)
* fasttree 2.0
* cd-hit 3.1

Installing QIIME
----------------

Stable Pre-Release
^^^^^^^^^^^^^^^^^^
Currently the most stable version of QIIME is our 0.91 pre-release, which you can download from `here <http://sourceforge.net/projects/qiime/files/releases/Qiime-0.91.tar.gz/download>`_.

Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^
To get the latest development version of QIIME, you should check it out of our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in "Qiime/doc". Check out the latest version of QIIME using svn with the command::

	$ svn co https://qiime.svn.sourceforge.net/svnroot/qiime/trunk Qiime

The user can update QIIME, by using the following command::

	$ svn update /path/to/QIIME/

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

	$ export PATH=/path/to/python/bin/:$PATH

PYTHONPATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Your PYTHONPATH should contain the following:

* /path/to/PyCogent
* /path/to/QIIME
* /path/to/PyNAST

In the bash shell, you can use the following command (example only shows the path to QIIME, so other softwares can be added to the PYTHONPATH, using a similar approach): ::

	$ export PYTHONPATH=/path/to/QIIME/:$PYTHONPATH

RDP_JAR_PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user should also define an RDP_JAR_PATH variable, since this tutorial uses the RDP Classifier:

* /path/to/rdp_classifier-2.0.jar

In the bash shell, you can use the following command: ::

	$ export RDP_JAR_PATH=/path/to/rdp_classifier-2.0.jar

Testing QIIME Install
---------------------
Once the source code is downloaded, the user should test QIIME to be sure all essential software is properly installed and the correct environment variables are set.

In a terminal window the user, should cd to their qiime/test directory using the following command: ::

	$ cd /path/to/QIIME/tests/

Then run the following test command: ::

	$ python all_tests.py -v

If all tests run properly, then QIIME was properly installed. Some test scripts may fail, due to optional third party applications not being installed. Test failures which contain the text ApplicationNotFoundError can be safely ignored. The following scripts may generate errors:

1. test_align_seqs.py - will fail if PyNAST, BLAST, muscle or infernal are not installed. PyNAST requires BLAST, and is the default sequence aligner. If you are not using PyNAST, you can safely ignore fails regarding BLAST and/or PyNAST. If you are not using muscle or infernal, you can safely disregard those error messages.
2. test_pick_otus.py – will fail if cd-hit and/or blast are not installed. If you are not planning to use these methods for OTU picking, you do not need to worry about these failures.
3. test_pyronoise.py - Will fail in PyroNoise is not installed.

