.. _doc_sra_submission:

.. index:: SRA Submission

========================= 
SRA Submission 
=========================

Introduction 
------------

This tutorial illustrates how to use QIIME for SRA submission of barcoded 16S community sequencing data generated on the Roche 454. To run this tutorial, you'll need to have QIIME installed, as well as the 454 off-instrument tools (``sffinfo``, ``sfffile``).

Future versions of QIIME will support similar workflows and tutorials for SRA submission of meta-genomics data, and data generated on the Illumina platforms. 

This tutorial uses as an example the `Fierer et al. 2008 hand dataset <http://www.pnas.org/content/105/46/17994.long>`_. This represents the first barcoded SRA submission, and the finished product can be found via SRA accession #: `SRS001216 <http://www.ncbi.nlm.nih.gov/sites/entrez?db=sra&term=SRS001216>`_. 

To get started, download the data from: `SRA Tutorial Data <http://bmf.colorado.edu/QIIME/knight_handstudy_demo-v1.1.0.zip>`_

Philosophy of QIIME's SRA submission code
-----------------------------------------

This document and code does not seek to capture the full complexity of what is possible in the SRA. Instead, the goal is to provide a simple pathway for submission to the SRA of the most common types of data. For example, more normalization could be provided through requiring additional tables, but this additional complexity is not justified when the appropriate rows can be copied/pasted/updated in Excel in a few seconds. More flexibility and better normalization is certainly possible and may be supported in future: the focus here is on getting something that works now.

Originally, the plan was to do this as a series of standalone scripts, but this ended up with an unacceptable level of reimplementation of things that are already in QIIME. We have therefore contributed the code into QIIME.

In addition to QIIME's standard dependencies, you will need several of QIIME's optional dependencies to complete this tutorial:

	* `BLAST <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_ 
	* `cdhit <http://www.bioinformatics.org/cd-hit/>`_ 
	* 454 off-instrument tools (``sffinfo``, ``sfffile``). These must be obtained from Roche 454.

You should refer to the `QIIME install pages <../install/index.html>`_ for all information related to getting up and running with QIIME. 

Overview of the Submission Process 
----------------------------------

An SRA submission can consist of metadata only, or of metadata together with sequence data. The SRA has recently added support for submission of barcoded pyrosequencing runs. This document describes how to prepare such submissions in a two-stage process:

	1. Submission of the study and sample metadata, including clinical metadata. 
	2. Submission of the experiment and run data, and associated technical metadata.

Updates are handled by regenerating the submission xml files (and optionally the .sff files), at which point they will be reloaded into the SRA. The current pipeline does not provide a mechanism for doing partial updates (e.g. adding an experiment or a run to an existing submission) -- all the metadata and/or data for the submission must be replaced with clean files. It may be possible in the future to allow incremental submissions but this is not presently supported. 


Preparing a submission using the QIIME workflow scripts
-------------------------------------------------------

A submission consists of a ``submission.xml``, metadata file, which references other xml metadata files and optionally tar files of sequence data files.

Set up your data for the tutorial run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download and unpack the tutorial data, and change to the resulting directory::

	wget http://bmf.colorado.edu/QIIME/knight_handstudy_demo-v1.1.0.zip
	unzip knight_handstudy_demo-v1.1.0.zip
	cd knight_handstudy_demo-v1.1.0

Generate the study and sample metadata submissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	Input: 
		a. ``study.txt`` - tabular metadata about the study (this is used to accession the study). 
		b. ``sample.txt`` - tabular metadata about each sample (this is used to accession samples). 
		c. ``submission.txt`` - a two-column file of tabular metadata about the submission.
	
	Output: 
		a. ``study.xml`` - xml-format metadata about the study. 
		b. ``sample.xml`` - xml-format metadata about each sample 
		c. ``submission.xml`` - xml-format metadata about the study and sample submission

Run the following command::

	make_sra_submission.py -a sample.txt -t study.txt -u submission.txt

This produces ``sample.xml``, ``study.xml``, and ``submission.xml`` from the corresponding tab-delimited text files. 

Generate the experiment and run metadata submissions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	**Input:** 
		a. ``experiment.txt`` - tabular metadata about the contents of each combination of library and sff file. 
		b. ``sff_files`` - a directory of multiple sff files containing the actual sequence data. 
		c. ``submission_second_stage.txt`` - a two-column file of tabular metadata about the submission.
		d. ``sra_parameters.txt`` - a qiime parameters file, defining what parameters should be passed to the individual component scripts
		e. ``greengenes_unaligned.fasta-OTUs_at_0.05.fasta`` - a FASTA file of known 16S sequences to screen for non-16S contaminants (here, we use a subset of the unaligned greengenes database, filtered at 95% sequence identity)  

	**Output:** 
		a. ``experiment.xml`` - xml-format metadata about the set of experiments described in ``library.txt`` 
		b. ``run.xml`` - xml-format metadata about each run, i.e. the association between a specific member of a pool and a specific xml file. 
		c. ``data.tgz`` - a gzipped tar archive containing individual sff files for each SRA RUN (see the Questions above if you are unclear on the distinction between the SRA RUN concept and the concept of an instrument run). 

The workflow script `process_sra_submission.py <../scripts/process_sra_submission.html>`_ may be used to create a submission of experiment and run metadata in one step.  

Run the following command::

	process_sra_submission.py -s sff_files -e experiment.txt -r greengenes_unaligned.fasta-OTUs_at_0.05.fasta -u submission_second_stage.txt -p sra_parameters.txt -o sra_out

This produces a tar archive of per-sample SFF files, :file:`experiment.xml`, :file:`run.xml`, and :file:`submission_second_stage.xml` from the input files. The list of commands that were actually run is available in the log file that in the top-level ``sra_out/`` directory.


Notes regarding individual steps of the SRA submission process
--------------------------------------------------------------

The `process_sra_submission.py <../scripts/process_sra_submission.html>`_ workflow script combines many separate QIIME commands. This section provides a discussion of the key components of the SRA submission workflow.

Print the commands to be run by the workflow without actually running them
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is sometimes useful to get the individual commands that will be run, but not actually run them. This is useful, for example, if you want to tweak one or more of the commands and then run them all via a bash script. To get the commands, but not run them, you can append ``-w`` to the ``process_sra_submission.py`` call::

	process_sra_submission.py -s sff_files -e experiment.txt -r greengenes_unaligned.fasta-OTUs_at_0.05.fasta -u submission_second_stage.txt -p sra_parameters.txt -o sra_out -w

Get fasta and qual from sff files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step converts the sff files into text formats that are more usable. 

**Output:** Makes .fna and .qual files for each sff file.

Produce valid mapping file for library demultiplexing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step converts the input experiment file into separate mapping files for each combination of STUDY and RUN_PREFIX (separating by run prefix is necessary when the same barcodes are used in different runs). This allows demultiplexing of the separate studies, which will then be sent in as separate submissions, and of the different barcoded plates, which will be demultiplexed separately.

Note: the LINKER field is no longer required in the spreadsheet.

**Output:** Produces valid mapping files per 454 plate: :file:`fierer_hand_study_E86FECS.map` and :file:`fierer_hand_study_FA6P1OK.map`

Demultiplex libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step assigns each sequence to a library, dropping low-quality sequences and producing a log explaining why specific sequences were dropped.

NOTE: The SRA requests that you deposit ALL your sequence data, including bad reads, unless there is an IRB reason not to do so (i.e. human contamination). Therefore the quality and length filtering should be turned off. We do this by setting high values for the quality and length filtering, that in practice are not exceeded. For details on the quality and length filtering options run::

	split_libraries.py -h

**Output:** Produces two files: :file:`seqs.fna` with valid sequences assigned to samples via barcodes, and :file:`split_libraries_log.txt` with info about which sequences failed QC. The parameters used are essentially turning off the default quality filters. You can turn off the quality filtering steps if you want to make sure that all the sequences appear in the output. You should do this by editing the appropriate values in your ``sra_parameters.txt`` file.

Briefly, the relevant settings in ``sra_parameters.txt`` require an average qual score of at least 5; a minimum sequence length of 30 (basically just the primer_barcode); a maximum sequence length of 1000; max homopolymer run of 1000; up to 100 errors in the primer; etc. In this run, we specify that we are using 12-base barcodes, (turning off the Golay error-correction which would be specified with ``split_libraries:barcode-type golay_12``). These parameters are specified in ``sra_parameters.txt`` as::

	split_libraries:min-qual-score	5
	split_libraries:min-seq-length	30
	split_libraries:max-seq-length	1000
	split_libraries:barcode-type	12
	split_libraries:max-homopolymer	1000
	split_libraries:max-primer-mismatch	100
	split_libraries:max-ambig	1000

Reduce sequence complexity by picking OTUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step reduces the number of sequences to do the human screen by picking OTUs with cdhit at 95%. We make the simplifying assumption that sequences that are identical over the first 100 bases will fall into the same OTU. These parameters are specified in ``sra_parameters.txt`` as::

	# pick_otus parameters
	pick_otus:otu_picking_method	cdhit
	pick_otus:max_cdhit_memory	4000
	pick_otus:prefix_prefilter_length	100
	pick_otus:similarity	0.95

**Output:** Produces two files: :file:`E86FECS_demultiplex/seqs_otus.txt` and :file:`E86FECS_demultiplex/seqs_otus.log` (which have the OTUs and the log file describing the analysis respectively).

The same procedure is applied to each library.

We then pick a representative sequence for each OTU by choosing the most abundant sequence from each OTU.

Blast the representative set sequences against 95% OTUs in greengenes to eliminate sequences that aren't really 16S rRNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step performs a human/contaminant screen the "safe" way by identifying and excluding sequences that aren't 16S rRNA. We are using blastn with a word size of 10, requiring 25% coverage of the sequence, and an E-value of 1e-20. Our tests suggest that this is sufficient to screen out human genomic reads (the human 18S sequence hits bacterial 16S with E-value between 1e-18 and 1e-10 depending on lineage). These parameters are specified in ``sra_parameters.txt`` as::

	# exclude_seqs_by_blast parameters
	exclude_seqs_by_blast:word_size	10
	exclude_seqs_by_blast:percent_aligned	0.25
	exclude_seqs_by_blast:e_value	1e-20

**Output:** This produces a bunch of log files and output; the file of screened seqs (i.e. that failed to hit a known 16S rRNA with even relaxed criteria).

The same procedure is applied to each library.

Make per-library files of "good" ids to pass to sfffile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step maps the ids of the representative set back onto the ids of the OTUs they came from so that we can get all the members of the OTUs that had a representative that matched a known 16S rRNA.

**Output:** This makes a new directory called :file:`E86FECS_demultiplex/per_lib_idlists`, which contains a separate file with an id list for each library.

The same procedure is applied to each library.

Use sfffile to make per-library sff files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step takes the good lists of ids from step 7 and extracts a separate sff file for each of those lists.

Use sfffile to quality-trim the barcodes, primers and linkers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SRA requires that the user reset the left trim in the sff file to eliminate the technical reads (barcode, primer, linker if present). This means figuring out the length of the technical parts of the read, the length of the current read, writing out a text file with the per-id info, and running sfffile to reset the read lengths.

Optional post-processing: modifying the second-stage submission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `make_sra_submission.py <../scripts/make_sra_submission.html>`_ script has the ability to include per-experiment attributes or links.  The attributes and links should be specified in separate, tab-delimited files. For example, a file named :file:`attributes.txt` can be created with the following contents:

::

  #EXPERIMENT_ALIAS	Attribute	Value
  fierer_hand_study_FA6P1OK	library strategy	targeted-locus
  fierer_hand_study_FA6P1OK	gene	16S rRNA V1-V2 region
  fierer_hand_study_E86FECS	library strategy	targeted-locus
  fierer_hand_study_E86FECS	gene	16S rRNA V1-V2 region

The following command will then add "gene" and "library strategy" attributes to both experiments in the resulting XML. (The experiment alias is specified in :file:`experiment.txt`, under the field 'EXPERIMENT_ALIAS'.) ::

  make_sra_submission.py -u submission_second_stage.txt -e experiment.txt -s per_run_sff --experiment_attribute_fp=attributes.txt

Links may be added to the experiments in a similar manner. After the `make_sra_submission.py <../scripts/make_sra_submission.html>`_ script has been run, the resulting XML files are ready to submit to the SRA.

Note: SRA prefers you give the individual files more meaningful names than the defaults, so suggest not just using generic names like experiment etc.

Questions about the Submission Process
--------------------------------------

	Q1. Can I submit multiplexed pyrosequnecing runs now?

	A1. Yes.

	Q2. Can I combine mock and clinical samples on the same 454 plate? (or, more generally, can I combine samples from different studies on the same 454 plate?)

	A2. Yes, but you must specify in the :file:`library.txt` input file which samples go with which study.

	Q3. Can I associate the same sample (and thus reads) with more than one study?

	A3. No.

	Q4. Can I combine samples that use different primers on the same run?

	A4. Yes, but you must specify in the :file:`library.txt` input file which primer was used for each "member" of the pooled library.

	Q5. Who will submit what?

	A5. At this stage, we expect the DACC to submit both the sample/study metadata and the experiment/library metadata and sequence data for the pilot. Later, the Centers will have the capacity to submit their own data. Centers will be credited with their data appropriately regardless of the mechanism by which the submission is actually performed. The submission will be a two-stage process: (1) the creation of study and sample records by the DACC, (2) the submission of sequence data and associated metadata by the DACC and/ or the Centers.

	Q6. Can I associate the same sample with more than one barcode and/or primer?

	A6. Yes, but you must specify a unique identifier for each "member" of the pool that associates the sample, primer and barcode.

	Q7. What is the distinction between a STUDY, an EXPERIMENT, and a RUN?

	A7. As SRA uses the terms, a STUDY is a collection of EXPERIMENTS. An EXPERIMENT is a LIBRARY (potentially a library of many samples that form a POOL, if multiplexing was used -- each MEMBER of a pool is associated with a sample, a primer, and a barcode) that was sequenced using one or more instrument runs. A RUN is the sequencing of a particular MEMBER of a pooled library on a particular instrument at a particular time. Thus, a single instrument run gives rise to many RUN entries in SRA.

	Q8. Is there an intermediate level between STUDY and EXPERIMENT?

	A8. Not for practical purposes. SRA will eventually allow a hierarchy of STUDY entries but this is not yet implemented.

	Q9. Do I really have to make a separate sff file for every MEMBER of every POOL for every instrument run?

	A9. Yes, and you also have to reset the quality trimming to correspond to the primer that was used for that particular member. The SRA will, in future, provide the demultiplexing service, but for now requires that the submissions be demultiplexed in advance. Fortunately, the accompanying scripts assist with this process.

	Q10. Is it OK for primers to be different lengths on the same 454 run?

	A10. Yes, but not within the same MEMBER of a library (i.e. if you have primers of different lengths, the different lengths are considered different MEMBER entries and should be marked as such in :file:`library.txt`).

	Q11. How should degenerate primers be handled?

	A11. All possible sequences that match the degenerate primer should be allowed using the EXPECTED_BASECALL_TABLE mechanism in :file:`experiment.xml` (see example).

Standard sra_parameters.txt file for barcoded 16S community sequencing on 454
-----------------------------------------------------------------------------

Currently our standard parameters files looks like the following. You can copy and paste this to a text file, and pass it with ``-p`` to ``process_sra_submission.py``.

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
	pick_otus:otu_picking_method	cdhit
	pick_otus:max_cdhit_memory	4000
	pick_otus:prefix_prefilter_length	100
	pick_otus:similarity	0.95

	# exclude_seqs_by_blast parameters
	exclude_seqs_by_blast:word_size	10
	exclude_seqs_by_blast:percent_aligned	0.25
	exclude_seqs_by_blast:e_value	1e-20