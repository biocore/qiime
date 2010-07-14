.. _doc_sra_submission:

.. index:: SRA Submission

========================= 
SRA Submission 
=========================

QIIME's SRA submission code is changing rapidly in the ``1.1.0-dev`` versions of QIIME. We're making every effort to keep the documentation up-to-date with the code, but if you are interested in using this code you'll most likely want to get on the ``qiime-sra-users`` mailing list to stay current on the state of the code. Head `here <http://groups.google.com/group/qiime-sra-users?hl=en>`_ for information on joining that public group.

Introduction 
------------

This document discusses using QIIME to submit your own barcoded 454 16S community sequencing data to the SRA. When you're ready to start this process you should have your sff files and your metadata mapping file. The following steps are covered: 

	1. Getting an SRA study accession number from SRA.
	1. Generation of template files using a web service provided by the QIIME team. 
	1. Fill in template files to describe your study and data.
	1. Applying the ``process_sra_submission.py`` workflow in QIIME to generate per-sample sff files; screen and remove human contaminants from your sequence set (if desired); and generate the SRA XML files.
	1. Submit your data to the SRA.

In addition to QIIME's standard dependencies, you will need several of QIIME's optional dependencies to complete this process:

	* ``uclust``
	* 454 off-instrument tools (``sffinfo``, ``sfffile``). These must be obtained from Roche 454.

You should refer to the `QIIME install pages <../install/index.html>`_ for all information related to getting up and running with QIIME. 

If you are interested in runner the latter steps with a test data set, you should review the `QIIME SRA submission tutorial <../tutorials/doc_sra_submission.html>`_.


Step 1. Getting an SRA study accession number from SRA.
-------------------------------------------------------
Greg will fill in.

Step 2. Generation of template files using a web service.
---------------------------------------------------------
Doug will fill in.

Step 3. Fill in template files to describe your study and data.
---------------------------------------------------------------
Kyle will fill in.

Step 4. Applying the ``process_sra_submission.py`` workflow in QIIME.
---------------------------------------------------------------------
Greg will fill in.

Step 5. Submit your data to the SRA.
------------------------------------
Need someone to fill in exactly what files are submitted, URLs, how to get an upload account, etc.


Troubleshooting
---------------


Standard sra_parameters.txt file for barcoded 16S community sequencing on 454
-----------------------------------------------------------------------------

Currently our standard parameters files looks like the following. You can copy and paste this to a text file, and pass it with ``-p`` to ``process_sra_submission.py``. The ``pick_otus:similarity`` value has been carefully chosen to exclude human sequences but include bacterial/archaeal 16S sequences, so it's not a good idea to change that without exploring the affect it will have.

::
	
	# split_libraries parameters
	split_libraries:min-qual-score	5
	split_libraries:min-seq-length	30
	split_libraries:max-seq-length	1000
	split_libraries:barcode-type	12
	split_libraries:max-homopolymer	1000
	split_libraries:max-primer-mismatch	100
	split_libraries:max-ambig	1000

	# pick_otus parameters
	pick_otus:similarity	0.70
	pick_otus:enable_rev_strand_match	True