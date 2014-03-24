.. _doc_sra_submission:

.. index:: SRA Submission

========================= 
SRA Submission 
=========================

QIIME's SRA submission code is changing rapidly in the ``1.2.0-dev`` versions of QIIME. We're making every effort to keep the documentation up-to-date with the code, but if you are interested in using this code you'll most likely want to get on the ``qiime-sra-users`` mailing list to stay current on the state of the code. Head `here <http://groups.google.com/group/qiime-sra-users?hl=en>`_ for information on joining that public group.

Introduction 
--------------

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
The SRA submission workflow involves four input files:

1. Study input file
2. Sample input file
3. Experiment input file
4. Submission input file

The SRA input files are formatted as tab-delimited text.  The first
non-blank line is the header line, and specifies the field labels for
the input file.  Field labels may appear in any order, but all
required field labels must be present for each file type.  Subsequent
lines specify a set of values for the fields.  SRA input files may
contain comment lines, which must start with '#'.

A spreadsheet program may be used to edit the SRA input files.  After
editing, they must be exported to TSV format for use with QIIME.  On
Google Docs, the menu item *File -> Download as -> Text (current sheet)*
will generate a file in the correct format.

The `SRA Field List <sra_field_list.html>`_ describes all available
fields for each type of input file.

Step 4. Applying the ``process_sra_submission.py`` workflow in QIIME.
---------------------------------------------------------------------
After filling in your templates, you will run the ``process_sra_submission.py`` workflow script. This code generates per-library sff files, as required by the SRA, and performs human screen if it is necessary to remove human contaminant sequences from your sff files. (If you're not sure whether you need to human screen, perform the human screen.) 

There are two types of human screens that can be applied: a positive screen and a negative screen. A positive screen is where all sequences are searched against a database, and those found to match the database are retained. This is used, for example, when you have 16S reads and you want to search them against greengenes to discard all sequences that don't look like they are 16S and therefore represent possible human contaminants. This is a conservative screen, but it is typically a better idea to err on the side of discarding too many reads than to accidentally submit human sequences. A negative screen, on the other hand, is where all sequences are searched against a database, and those found not to match the database are retained. This is less conservative, but necessary for example when you have metagenomic reads, and searching against a single sequence collection therefore won't work (because you'd need full genome sequences of all microbes). In this case, you can screen against the full human genome, and retain only sequences that do not match. When using a negative screen, be sure that your reference set contains the full genome, not just (e.g.) the coding regions. To bypass the human screening step, do not provide a reference set via ``-r``.

You can run ``process_sra_submission.py`` as follows::

	process_sra_submission.py -s sff_files/ -e experiment.txt -r 16S_reference_set.fasta -u submission.txt -p sra_parameters.txt -o sra_out/


The options are:

	* -s : directory containing sff files
	* -e : the experiment.txt file generated in Step 2 and filled in in Step 3
	* -r : reference set fasta file for human screening (optional: to bypass the human screen do not pass -r, but be certain that this is what you want to do)
	* -u : the submission.txt file generated in Step 2 and filled in in Step 3
	* -p : the parameter file -- the standard parameter file used here can be copied from `here <../tutorials/doc_sra_submission.html#standard-sra-parameters-txt-file-for-barcoded-16s-community-sequencing-on-454>`_
	* -o : the directory where the output should be written

A tutorial is provided `here <../tutorials/doc_sra_submission.html>`_ that illustrates how to run ``process_sra_submission.py`` using an example data set.


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


Resubmitting and changing your data in SRA
--------------------------------------------

As long as your data in not succesfully loaded and was assigned an SRA ID, you can simply re-upload new XML files. The new ones will replace the old ones. Once at least part of your data is accepted you have to use the MODIFY flag in the submission.xml like this::

   <ACTION>
	<MODIFY source="experiment_modified.xml" target="SRX025053 " notes="Modifying SRX025053"/>
   ...


