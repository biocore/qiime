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
	2. Generation of template files using a web service provided by the QIIME team. 
	3. Fill in template files to describe your study and data.
	4. Applying the ``process_sra_submission.py`` workflow in QIIME to generate per-sample sff files; screen and remove human contaminants from your sequence set (if desired); and generate the SRA XML files.
	5. Submit your data to the SRA.

In addition to QIIME's standard dependencies, you will need several of QIIME's optional dependencies to complete this process:

	* ``uclust``
	* 454 off-instrument tools (``sffinfo``, ``sfffile``). These must be obtained from Roche 454.

You should refer to the `QIIME install pages <../install/index.html>`_ for all information related to getting up and running with QIIME. 

If you are interested in running the latter steps with a test data set, you should review the `QIIME SRA submission tutorial <../tutorials/doc_sra_submission.html>`_.


Step 1. Getting an SRA study accession number from SRA.
-------------------------------------------------------
Go to the `SRA homepage <http://www.ncbi.nlm.nih.gov/Traces/sra>`_ and click the Submit/Submissions tag. Register (only the first time, of course) and log in via the NCBI PDA route. Now, you can create a submission. The Submission ID can be choosen freely (SRA says it should make sense to the submitter). This ID should/must/can (??) be used later as SAMPLE_ALIAS in the submission process. You will receive a SRA accession ID and a new submission without any data loaded.

If you registered as a new user and belong to a center, you might want to consider to be added to the center's user list by contacting sra@ncbi.nlm.nih.gov. 



Step 2. Generation of template files using a web service.
---------------------------------------------------------
To generate your template files, please visit the following website:

`http://microbio.me/qiime <http://microbio.me/qiime>`_

	* If you already have an account you can use that to log into the website. If not, you can use the links at the bottom to create a new account for yourself.
	* Once logged in, click the "Create a New Study" link and fill out the form.
	* After creating a new study, new links will be available on the main page. Click the "Generate a MIENS-compliant metadata template" link to begin the process of generating your template files.
	* On the new template page you will be presented with all required fields for your MIENS study type(s). You will also have the option of selecting any additional fields that apply to your study. If you do not see a field you'd like to add, scroll to the bottom of the page and use the search utility. If you still do not find the field you wish to add, you may then manually add a column to the appropriate template file after you've generated your template files (make sure these fields are less than 25 characters in length, start with a letter, and contain only letters, numbers, and the "_" character).
	* Once you've made your field selections, click the "Generate Templates" button at the bottom of the page. You may download each file individually or download the entire bundle as a .zip archive.
	
Working with Ontologies:

Within the generated template files is a field reference which describes each field and what the expected values are. For fields that require an ontology value it may be helpful to use the Ontology Lookup utility located on the main page (after creating or selecting a study). You may also wish to browse the ontology online. Links to the referenced ontologies are included in the field reference documentation generated when you created your templates. Look up the field in question and read the description to find the URL of the ontology of interest.

Step 3. Fill in template files to describe your study and data.
---------------------------------------------------------------
Kyle will fill in.

Step 4. Applying the ``process_sra_submission.py`` workflow in QIIME.
---------------------------------------------------------------------
Greg will fill in.

(That should be mainly covered in the other tutorial, right?)

Step 5. Submit your data to the SRA.
------------------------------------

SRA distinguishes between two types of submitters: individuals or centers. E.g. the Knight lab has registered a center name with the NCBI called "Center for Comparative Microbial Ecology" (CCME). All projects overseen by this center must be submitted via the project's account (Greg has the credentials for CCME). If you submit as an individual write an email to trace@ncbi.nlm.nih.gov to request the curent ftp address of the anonymous ftp server.

The actual submission consists of several XML files and the demultiplexed sff files, usually produced with make_sra_submission.py or process_sra_submission.py.

Metadata files:

- Study: XML file specifying sequencing study
- Sample: XML file specifying the target of sequencing
- Experiment: XML file specifying experimental organization and parameters 
- Run: One of more XML descriptors linking run data to their experiments
- Submission XML file specifying submission session


To check your files against the XML schema::

   xmllint --schema  http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA/SRA.run.xsd?view=co run.xml

Replace "run" in the URL with "sample", "study", or "experiment" to validate the other files.

The sff files are already tared and zipped by process_sra_submission.py, but the xml files should be collected in one directory and then be tared and zipped. These two files will be uploaded to the SRA ftp site, usually to the short_read subdirectory. At this point, it's a good idea to send an email to your contact at the SRA (can we give a general email sdress here?) to inform them of your upload.


Troubleshooting
---------------


