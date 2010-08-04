.. _doc_sra_field_list:

.. index:: SRA Field List

========================= 
SRA Field List
=========================

The make_sra_submission.py QIIME script uses tab-delimited input files
to generate valid XML files for an SRA submission.  This document
lists the fields in each type of input file, and gives details on how
to fill them in.

If a field is labeled as *optional*, it may be omitted from the
header.  For fields labeled *derived automatically*, a value will be
generated for each entry in the input file if the field is missing or
empty.

Study Input File
----------------

STUDY_ALIAS

  One study per publication (i.e. the STUDY is supposed to be about
  the same amount of info as in a paper). 

  The STUDY_ALIAS is used as a unique identifier for the study (it
  should be unique among all studies submitted from your institution).
  It should contain only alphanumeric characters and underscores.

  *Example*: hand_study

STUDY_TITLE 

  Expected (or actual) title of the paper that will be published about
  the study. Free text.

STUDY_TYPE

  Should be "Metagenomics" for 16S surveys (regrettably). Other
  choices relate to whole-genome studies. Controlled vocabulary.

STUDY_ABSTRACT

  Abstract, e.g. of the publication. Free text.

STUDY_DESCRIPTION

  Use "Targeted Gene Survey" for 16S or other target gene studies

CENTER_NAME

  NCBI-approved name of the center coordinating the overall study,
  e.g. WUGSC. If you don't have a center name, you need to get NCBI to
  define one and then use that - you can't do this as free text. This
  will often be the same as other center names.

CENTER_PROJECT_NAME

  Name of project as used by the center responsible for the study,
  NULL if none.

PMID

  PubMed ID of paper describing project, if supplied will write out
  STUDY_LINK block, can be multiple (comma-delimited), can be absent
  if no linked publication yet.


Submission Input File
---------------------

ACCESSION 

  Accession number for the submission: leave blank if not assigned
  yet, e.g. if new submission rather than replacing XML for an
  existing submission.

SUBMISSION_ID

  Internally unique id for the submission: this is used as an ALIAS
  elsewhere so needs to be short, alphanumeric+underscores, no special
  characters.

CENTER_NAME

  Name of the center (e.g. sequencing center) preparing the
  submission: assigned by NCBI so you need to get a name for your
  institution rather than making something up here.

SUBMISSION_COMMENT

  Free-text comments regarding submission

LAB_NAME

  Name of lab preparing submission, can differ from center (usually
  refers to the PI's info, not the sequencing center's)
  
SUBMISSION_DATE

  Timestamp of submission.
  
CONTACT

  Use semicolon to separate email address from name, can be multiple
  contacts. (Example: Rob Knight;Rob.Knight@Colorado.edu) 

FILE

  Leave blank if not submitting sequence data, put in filename
  otherwise.  (Example: fierer_hand_study.seqs.tgz)

Sample Input File
-----------------

SAMPLE_ALIAS

  Sample name that you make up for each sample (note: for dbGAP submissions, you need to get the name of each specimen from dbGAP). Used as identifier so must be alphanumeric plus underscores, no special characters.

TITLE

  (Example: "human hand microbiome") Arbitrary title for your sample that you make up.

TAXON_ID

  Taxon Id is what is getting sequenced: i.e txid 539655 = human skin
  metagenome, species, metagenomes

COMMON_NAME

  Common name of what is being sequenced, should match taxon id's
  name, e.g. human skin metagenome. (Example: "human skin metagenome")

ANONYMIZED_NAME

  Anonymized name of the subject, if applicable (e.g. deidentified
  subject IDs from dbGAP, deidentified subject ids from your
  study). Only applies to human studies, leave blank if not
  applicable.  (Example: "subject 1")

DESCRIPTION

  Free-text description of this specific sample.  (Example: "female
  right palm")

HOST_TAXID

  Required if there is a host (skip otherwise): taxon id that refers
  to the host. e.g. 9606 = Human.

Experiment Input File
---------------------

EXPERIMENT_TITLE

  Title of the experiment. Must be the same for every member of a
  given pool. Free text.

  *Example*: ``Sampling and pyrosequencing methods for quantifying
  bacterial communities in the human gut``

  *Output*: This field is used as the text of the <TITLE> element in
  the SRA Experiment XML file.

STUDY_REF

  Official alias of the study registered with SRA.  Must be the same
  for every member of a given pool but can be different for different
  pools. If you put items from multiple STUDY records (e.g. clinical
  and mock) on the same run, create separate pools but have them
  reference the same RUN_PREFIX so they can pull sffs from the same
  files.

  This field is used as the *refname* attribute of the
  <STUDY_REF> element in the SRA Experiment XML file.  It is also used
  to derive several optional fields.

STUDY_CENTER

  Name of the center associated with the overall STUDY, i.e. whoever
  is designated as having overall responsibility for the STUDY (this
  is a controlled vocabulary, assigned by NCBI). Needs to be the same
  for every member of a pool.

  This field is used as the *refcenter* attribute of the
  <EXPERIMENT_REF> element in the SRA Run XML file.

SAMPLE_ALIAS

  Unique (within the STUDY referenced) ID for each sample. You can use
  the same sample in multiple pools referenced in the same
  EXPERIMENT. If you mixed samples from more than one STUDY in the
  same EXPRIMENT, the components from each STUDY need to be registered
  as a separate EXPERIMENT.

  This field is used as the *refname* attribute of the
  <MEMBER> element in the SRA Experiment XML file.  It is also used to
  derive several optional fields.

POOL_PROPORTION

  Floating-point number representing the fraction of the pool that was
  intended to come from that library member.

  This field is used as the *proportion* attribute of the
  <MEMBER> element in the SRA Experiment XML file.

BARCODE

  Barcode sequence used for each pool member.  Each combination of
  barcode, primer and plate region must be unique.

  This field is used as the text of the <BASECALL> element in
  the SRA Experiment XML file.  It is also used to derive several
  optional fields.

RUN_PREFIX

  The 454 instrument usually produces more than one sff file. This
  should be the prefix of the sff file name that was produced by a
  given run (usually these will have 01, 02, etc. sufixes). This
  allows you to designate a pool as per-library rather than per sff
  file (otherwise you would need to duplicate all the info per run for
  each sff file).

  This field is used as the *name* attribute of the
  <DATA_BLOCK> element in the SRA Run XML file.  It is also used to
  derive several optional fields.

EXPERIMENT_DESIGN_DESCRIPTION

  Free text description of the overall motivation for the experiment
  (i.e. pool) - why those samples were mixed together, what it was
  for, etc.  Needs to be the same for every member of a pool.

LIBRARY_CONSTRUCTION_PROTOCOL

  Free-text description of how the library was put together (e.g. from
  the methods section of a paper).  Needs to be the same for
  everything in a given pool.

SAMPLE_CENTER *

  Name of the center that provided the sample, can be separate for
  each sample.  If sample information is stored in dbGAP, the
  SAMPLE_CENTER should be set to "NCBI".

PLATFORM *

  This is the sequencing platform, e.g. FLX or Titanium.  If the
  platform value is not found in a table of supported platforms, a
  KeyError is raised.

KEY_SEQ *

  This is a technical aspect of the 454 platform, is usually TCAG, can
  be obtained from the sff file using the sfftools.

REGION *

  Region of the plate that was sequenced (in cases where there was a
  split run and the same primer/barcode means different things in
  different parts of the plate).

RUN_CENTER *

  Name of the institution that performed the run, assigned by
  NCBI. You can use the center name for your lab for this even if you
  had the sequencing done elsewhere according to SRA.

EXPERIMENT_CENTER *

  Official abbreviation for the sequencing center associated with the
  experiment, i.e. who made the pool. Needs to be the same for every
  member of a given pool. This is your center name as assigned by NCBI
  and is often the same as the STUDY center.

EXPERIMENT_ALIAS (optional, derived automatically)

  Unique id (within the submission) for the experiment.  Needs to be
  the same for everything in a given pool.  If absent, the value will
  be derived as <STUDY_REF>_<RUN_PREFIX>.

RUN_ALIAS (optional, derived automatically)

  Alias for the run.  Presently, this should be different for every
  pool member, since each pool member gets a unique RUN element in the
  run XML.  In the future, we plan to change this behavior, and create
  only a single RUN element of multiple pool members share the same
  RUN_ALIAS. Needs to be a short identifier, alphanumeric and
  underscores only (no special characters).  If absent, this field is
  automatically derived as <STUDY_REF>_<SAMPLE_ALIAS>_<RUN_PREFIX>.

RUN_ACCESSION (optional)

  Optional accession number for the run. Leave blank if not already
  assigned.

STUDY_ACCESSION (optional)

  Optional accession number for study. You should already have created
  the study in SRA in the first stage submission and may reuse that id
  here.

EXPERIMENT_ACCESSION (optional)

  Optional accession number for the experiment. If you already created
  the Experiment accession in SRA, use it -- otherwise, leave blank.

DEFAULT_SAMPLE_CENTER (optional)

  Optional default sample center.  If absent, the value from the
  SAMPLE_CENTER field is used instead.

DEFAULT_SAMPLE_ACCESSION (optional)

  Optional default sample accession number, if available (leave blank
  if you don't have e.g. an accession assigned by dbGAP).

DEFAULT_SAMPLE_NAME (optional, derived automatically)

  Optional reference name for the default sample.  If this field is
  not present, and no DEFAULT_SAMPLE_ACCESSION is provided, the name
  is automatically derived as <STUDY_REF>_default.  Otherwise, the
  default sample is specified by the accession number alone, and this
  attribute is not inserted into the XML output.

POOL_MEMBER_ACCESSION (optional)

  Optional accession number for pool member. This field should be
  blank or not present if an SRA accession number is not already
  assigned.

POOL_MEMBER_NAME (optional, derived automatically)

  Unique (within the pool) id for each pool member. In the hand
  example, we only used V2 primers, so I am calling the pool members
  S1_V2 etc. If you mixed primers, a reasonable thing to do would be
  to use sample_primer codes; if you did replicates doing different
  barcodes you might want to use sample_primer_barcode or
  sample_primer_replicate, if you used different numbers of PCR cycles
  you might want to use sample_numcycle, etc.

  If absent, the derived value of this field depends on the primer.
  If the PRIMER field is not blank, it is derived as
  <RUN_PREFIX>_<SAMPLE_ALIAS>_<PRIMER_READ_GROUP_TAG>.  Otherwise, a
  value of <RUN_PREFIX>_<SAMPLE_ALIAS> is used.

POOL_MEMBER_FILENAME (optional, derived automatically)

  Filename for SFF file containing sequences from this pool member.
  The SFF files are searched for in a subdirectory of the sff_dir
  named after the RUN_PREFIX.  If the field is blank or absent, a
  default value of <POOL_MEMBER_NAME>.sff is used.

BARCODE_READ_GROUP_TAG (optional, derived automatically)

  Pool that a sample will be assigned to based on the barcode.  If
  this field is absent, a value of <RUN_PREFIX>_<BARCODE> is derived
  automatically.

LINKER (optional)

  Linker sequence between the primer and the barcode (to reduce
  differences in hybridization based on the barcode).  This field may
  be empty.

PRIMER (optional)

  Primer sequence that was used for this particular library member. If
  you used more than one primer for a given pool member (which is
  allowed) you need to duplicate the whole row with the additional
  primer information. This needs to be the actual sequence of the
  primer, not the name of the primer (i.e. not V2).

PRIMER_READ_GROUP_TAG (optional, derived automatically)

  Read group that samples will be assigned to based on the primer,
  e.g. V2 for the V2 primers. By default, multidimensional
  demultiplexing on the barcode and primer is performed.  If it is not
  present, this field will be derived using a table of standard primer
  read group tags.  If the primer is not found in the table, a
  KeyError is raised.

LIBRARY_STRATEGY (optional, default value provided automatically)

  Sequencing technique intended for this library (optional
  field). This will usually be AMPLICON (default) or METAGENOMIC.

LIBRARY_SOURCE (optional, default value provided automatically)

  Type of source material that is being sequenced (optional
  field). This will usually be GENOMIC (default) or METAGENOMIC.

LIBRARY_SELECTION (optional, default value provided automatically)

  Whether any method was used to select and/or enrich the material
  being sequenced (optional field). This is used in cases where
  e.g. the cells were sorted, if PCR was used to make a specific
  amplicon, if fractionation for viruses was done, etc.  The default
  value is PCR.

RUN_DATE (optional, currently unused)

  Date the run was performed: this can be obtained from the sff file.

INSTRUMENT_NAME (optional, currently unused)

  This field is used if the specific machine used has a name or label
  (i.e. a label on that specific piece of equipment, not the type of
  instrument). Some sequencing centers assign names to specific
  instruments."

SAMPLE_ACCESSION (**DEPRECATED**)

  Please use DEFAULT_SAMPLE_ACCESSION instead.  If the new field is
  blank or absent, this valie is used.  This field will continue to
  work, but will produce a warning.

