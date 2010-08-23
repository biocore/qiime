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

One study per publication (i.e. the STUDY is supposed to be about the
same amount of info as in a paper).



STUDY_ALIAS

  The STUDY_ALIAS is used as a unique identifier for the study (it
  should be unique among all studies submitted from your institution).
  It should contain only alphanumeric characters and underscores.

  *Example*: ``hand_study``

STUDY_TITLE 

  Expected (or actual) title of the paper that will be published about
  the study. Free text.

  *Example*: ``The influence of sex, handedness, and washing on the diversity of hand surface bacteria``

STUDY_TYPE

  Should be "Metagenomics" for 16S surveys (regrettably). Other
  choices relate to whole-genome studies. Controlled vocabulary.

  *Example*: ``Metagenomics``

STUDY_ABSTRACT

  Abstract, e.g. of the publication. Free text.

  *Example*: ``This project aims to undertake global surveys of microbial diversity in a range of free-living and host-associated communities. The importance of the project is that it will provide a comparison of microbial diversity in a range of habitats and provide a platform to underpin many studies of community assembly, diversity, etc. Bacteria thrive on and within the human body. One of the largest human-associated microbial habitats is the skin surface, which harbors large numbers of bacteria that can have important effects on health. We examined the palmar surfaces of the dominant and nondominant hands of 51 healthy young adult volunteers to characterize bacterial diversity on hands and to assess its variability within and between individuals. We used a novel pyrosequencing- based method that allowed us to survey hand surface bacterial communities at an unprecedented level of detail. The diversity of skin-associated bacterial communities was surprisingly high; a typical hand surface harbored >150 unique species-level bacterial phylotypes, and we identified a total of 4,742 unique phylotypes across all of the hands examined. Although there was a core set of bacterial taxa commonly found on the palm surface, we observed pronounced intra- and interpersonal variation in bacterial community composition: hands from the same individual shared only 17% of their phylotypes, with different individuals sharing only 13%. Women had significantly higher diversity than men, and community composition was significantly affected by handedness, time since last hand washing, and an individual's sex. The variation within and between individuals in microbial ecology illustrated by this study emphasizes the challenges inherent in defining what constitutes a ""healthy"" bacterial community; addressing these challenges will be critical for the International Human Microbiome Project.``

STUDY_DESCRIPTION

  Use "Targeted Gene Survey" for 16S or other target gene studies

  *Example*: ``Targeted Gene Survey from Human Skin``

CENTER_NAME

  NCBI-approved name of the center coordinating the overall study,
  e.g. WUGSC. If you don't have a center name, you need to get NCBI to
  define one and then use that - you can't do this as free text. This
  will often be the same as other center names.

  *Example*: ``CCME``

CENTER_PROJECT_NAME

  Name of project as used by the center responsible for the study,
  NULL if none.

  *Example*: ``NULL``

PMID

  PubMed ID of paper describing project, if supplied will write out
  STUDY_LINK block, can be multiple (comma-delimited), can be absent
  if no linked publication yet.

  *Example*: ``19004758``


Submission Input File
---------------------

ACCESSION 

  Accession number for the submission: leave blank if not assigned
  yet, e.g. if new submission rather than replacing XML for an
  existing submission.

  *Example*: ``SRA020077``

SUBMISSION_ID

  Internally unique id for the submission: this is used as an ALIAS
  elsewhere so needs to be short, alphanumeric characters with
  underscores, no special characters.

  *Example*: ``knight_handstudy``

CENTER_NAME

  Name of the center (e.g. sequencing center) preparing the
  submission: assigned by NCBI so you need to get a name for your
  institution rather than making something up here.

  *Example*: ``CCME``

SUBMISSION_COMMENT

  Free-text comments regarding submission

LAB_NAME

  Name of lab preparing submission, can differ from center (usually
  refers to the PI's info, not the sequencing center's)
  
  *Example*: ``Knight``

SUBMISSION_DATE

  Timestamp of submission.
  
  *Example*: ``2009-10-22T06:23:00Z``

CONTACT

  Use semicolon to separate email address from name.  Can include
  multiple contact fields.

  *Example*: ``Jane Smith;janesmith@college.edu``

FILE

  Leave blank if not submitting sequence data, put in filename
  otherwise.

  *Example*: ``fierer_hand_study.seqs.tgz``

Sample Input File
-----------------

SAMPLE_ALIAS

  Unique sample identifier, within each center.  (Note: for dbGAP
  submissions, the name of each specimen is provided from dbGAP). Must
  be formatted as alphanumeric characters plus underscores, with no
  special characters.

  *Example*: ``hand_sample_1``

TITLE

  Title of the sample, for display on the SRA website.

  *Example*: ``human hand microbiome``

TAXON_ID

  NCBI taxon ID of the sample.

  *Example*: ``539655`` (human skin metagenome, species, metagenomes)

COMMON_NAME

  Common name of the sample, should match taxon id's name.

  *Example*: ``human skin metagenome``

ANONYMIZED_NAME

  Anonymized name of the subject, if applicable (e.g. deidentified
  subject IDs from dbGAP, deidentified subject ids from your
  study). Only applies to human studies, leave blank if not
  applicable.

  *Example*: ``Subject 1``

DESCRIPTION

  Free-text description of this specific sample.

  *Example*: ``female right palm``

HOST_TAXID

  Required if there is a host (skip otherwise): taxon id that refers
  to the host. e.g. 9606 = Human.

  *Example*: ``Human``

Experiment Input File
---------------------

EXPERIMENT_TITLE

  Title of the experiment. Must be the same for every member of a
  given SRA Experiment. Free text.

  *Example*: ``Sampling and pyrosequencing methods for quantifying
  bacterial communities in the human gut``

  *Output*: Used as the text of the <TITLE> element in the SRA
  Experiment XML file.

STUDY_REF

  The alias of the SRA Study to which this SRA Experiment belongs.
  The STUDY_REF should indicate the official alias of the study
  registered with SRA.  If the SRA Study was created with QIIME, the
  alias is the value of STUDY_ALIAS defined in the study input file.
  Must be the same for every member of a given SRA Experiment.

  If items from multiple SRA Studies are included on the same
  sequencing run, create separate SRA Experiments for each study.

  *Example*: ``hand_study``

  *Output*: Used as the *refname* attribute of the <STUDY_REF> element
  in the SRA Experiment XML file.  It is also used to derive several
  optional fields.

STUDY_CENTER

  Name of the center that registered the SRA Study.  Study center
  names are assigned by the SRA, so you must contact them if your
  institution does not have an official study center designation.
  Needs to be the same for every member of a given SRA Experiment.

  *Example*: ``UPENNBL``

  *Output*: Used as the *refcenter* attribute of the <EXPERIMENT_REF>
  element in the SRA Run XML file.

SAMPLE_ALIAS

  The alias of the SRA Sample to which this pool member belongs.  The
  SAMPLE_ALIAS should indicate the official alias of a sample registered
  with SRA.  If the SRA Sample was created with QIIME, the alias is the
  value of SAMPLE_ALIAS defined in the sample input file.

  *Example*: ``water_blank``

  *Output*: Used as the *refname* attribute of the <MEMBER> element in
  the SRA Experiment XML file.  It is also used to derive several
  optional fields.

POOL_PROPORTION

  Floating-point number representing the fraction of the pool that was
  intended to come from that pool member.

  *Example*: ``0.05``

  *Output*: Used as the *proportion* attribute of the <MEMBER> element
  in the SRA Experiment XML file.

BARCODE

  Barcode sequence used for each pool member.  Each combination of
  barcode, primer and sequencing region must be unique.

  *Example*: ``ACGTCTGTAGCA``

  *Output*: Used as the text of a <BASECALL> element in the SRA
  Experiment XML file.  It is also used to derive several optional
  fields.

RUN_PREFIX

  The 454 instrument usually produces more than one sff file. This
  should be the prefix of the sff file name that was produced by a
  given run (usually these will have 01, 02, etc. sufixes). This
  allows you to designate a pool as per-library rather than per sff
  file (otherwise you would need to duplicate all the info per run for
  each sff file).

  *Example*: ``GAMA2IO``, for an SFF file named ``GAMA2IO01.sff``

  *Output*: Used as the *name* attribute of the <DATA_BLOCK> element
  in the SRA Run XML file.  It is also used to derive several optional
  fields.

EXPERIMENT_DESIGN_DESCRIPTION

  Free text description of the overall motivation for the SRA
  Experiment: why those samples were mixed together, what it was for,
  etc.  Needs to be the same for all entries with the same
  EXPERIMENT_ALIAS value.

  *Example*: ``Pool of samples from a handwashing study providing longitudinal data about recolonization in a small number of subjects.``

  *Output*: Used as the text of the <DESIGN_DESCRIPTION> element in
  the SRA Experiment XML file.

LIBRARY_CONSTRUCTION_PROTOCOL

  Free-text description of how the library was put together (e.g. from
  the methods section of a paper).  Needs to be the same for all
  entries with the same EXPERIMENT_ALIAS value.

  *Example*: ``Each amplicon library was constructed by amplifying the 16S rRNA gene using the 27f/534r, 357f/926r, 968f/1492r, BSF8/BSR534, BSF343/BSR926, or BSF917/BSR492 primer pair. Primers contained DNA barcode sequences such as those described by Hamady et al. 2008, and the recommended 454 adapter sequences. Amplification conditions are described in McKenna et al 2007, with exception that the polymerase used was AccuPrime (Invitrogen, Carlsbad, CA, USA).``

  *Output*: Used as the text of the <LIBRARY_CONSTRUCTION_PROTOCOL>
  element in the SRA Experiment XML file.

SAMPLE_CENTER *

  Name of the center that provided the sample, can be separate for
  each sample.  

  If sample information is stored in dbGAP, the SAMPLE_CENTER should
  be set to "NCBI".

  *Example*: ``UPENNBL``

  *Output*: Used as the *refcenter* attribute of the <MEMBER> and
  <STUDY_REF> elements in the SRA Experiment XML file.  It is also
  used to derive the DEFAULT_SAMPLE_CENTER field.

PLATFORM *

  The sequencing platform, either ``FLX`` or ``Titanium``.  If the
  platform value is not found in a table of supported platforms, the
  QIIME script will halt with an error.

  *Example*: ``Titanium``

  *Output*: Used to generate the contents of the <PLATFORM> element in
  the SRA Experiment XML file.

KEY_SEQ *

  This is a technical aspect of the 454 platform, is usually TCAG, can
  be obtained from the sff file using the sfftools.

  *Example*: ``TCAG``

  *Output*: Used as the text of the <EXPECTED_BASECALL> element inside
  the <READ_SPEC> element for the Adapter sequence in the SRA
  Experiment XML file.

REGION *

  Region of the plate that was sequenced.  If the plate contained only
  a single region, the SRA requires that this be set to '0'.

  *Example*: ``1``

  *Output*: Used as the *region* attribute of the <DATA_BLOCK> element
  in the SRA Run XML file.

RUN_CENTER *

  Name of the institution that performed the run, assigned by
  NCBI. You can use the center name for your lab for this even if you
  had the sequencing done elsewhere according to SRA.

  *Example*: ``UPENNBL``

  *Output*: Used as the *center_name* and *run_center* attributes of
  the <RUN> element in the SRA Run XML file.

EXPERIMENT_CENTER *

  Official abbreviation for the sequencing center associated with the
  experiment, i.e. who made the pool. Needs to be the same for every
  member of a given pool. This is your center name as assigned by NCBI
  and is often the same as the STUDY center.

  *Example*: ``UPENNBL``

  *Output*: Used as the *center_name* attribute of the <EXPERIMENT>
  element in the SRA Experiment XML file.

EXPERIMENT_ALIAS (optional, derived automatically)

  Unique id (within the submission) for the experiment.  This is the
  decisive element on which separate SRA Experiments are created.  If
  absent, the value will be derived as {STUDY_REF}_{RUN_PREFIX}.

  *Example*: ``hand_study_F0FN7DX``

  *Output*: Used as the *alias* attribute of the <EXPERIMENT> element,
  and as the text of the <LIBRARY_NAME> element in the SRA Experiment
  XML file.  Also used as the *refname* attribute of the
  <EXPERIMENT_REF> element in the SRA Run XML.

RUN_ALIAS (optional, derived automatically)

  Alias for the run.  Presently, this should be different for every
  pool member, since each pool member gets a unique RUN element in the
  run XML.  In the future, we plan to change this behavior, and create
  only a single RUN element of multiple pool members share the same
  RUN_ALIAS. Needs to be a short identifier, alphanumeric and
  underscores only (no special characters).  If absent, this field is
  automatically derived as {STUDY_REF}_{SAMPLE_ALIAS}_{RUN_PREFIX}.

  *Example*: ``hand_study_sample1_F0FN7DX``

  *Output*: Used as the *alias* attribute of the <RUN> element in the
  SRA Run XML.

STUDY_ACCESSION (optional)

  Optional accession number for study.  Leave blank or omit this field
  if not already assigned.

  *Example*: ``SRP003284``

  *Output*: Used as the *accession* attribute of the <STUDY_REF>
  element in the SRA Experiment XML file.

DEFAULT_SAMPLE_CENTER (optional)

  Optional default sample center.  If the field is blank or omitted,
  the value from the SAMPLE_CENTER field is used instead.

  *Example*: ``UPENNBL``

  *Output*: Used as the *refcenter* attribute of the
  <SAMPLE_DESCRIPTOR> element in the SRA Experiment XML file.

DEFAULT_SAMPLE_ACCESSION (optional)

  Optional default sample accession number, if available.  Leave blank
  if this has not already been assigned.

  *Example*: ``SRS107395``

  *Output*: Used as the *accession* attribute of the
  <SAMPLE_DESCRIPTOR> element in the SRA Experiment XML file.

DEFAULT_SAMPLE_NAME (optional, derived automatically)

  Optional reference name for the default sample.  If this field is
  blank/absent, and no DEFAULT_SAMPLE_ACCESSION is provided, the name
  is automatically derived as {STUDY_REF}_default.  Otherwise, the
  default sample is specified by the accession number alone, and this
  attribute is not inserted into the XML output.

  *Example*: ``hand_study_default``

  *Output*: Used as the *refname* attribute of the <SAMPLE_DESCRIPTOR>
  element in the SRA Experiment XML file.

POOL_MEMBER_ACCESSION (optional)

  Optional accession number for pool member. This field should be
  blank or absent if an SRA accession number is not already assigned.

  *Example*: ``SRS107397``

  *Output*: Used as the *accession* attribute of the <MEMBER> element
  in the SRA Experiment XML file.

POOL_MEMBER_NAME (optional, derived automatically)

  Unique (within the pool) id for each pool member.  If the field is
  blank or absent, it is derived automatically.  The derived value of
  this field depends on the primer.  If a primer was provided, it is
  derived as {RUN_PREFIX}_{SAMPLE_ALIAS}_{PRIMER_READ_GROUP_TAG}.
  Otherwise, a value of {RUN_PREFIX}_{SAMPLE_ALIAS} is used.

  *Example*: ``F0FN7DX01_sample1_V1-V2``

  *Output*: Used as the *member_name* attribute of the <MEMBER>
  element in the SRA Experiment XML file.  Also used to derive the
  value of POOL_MEMBER_FILENAME.

POOL_MEMBER_FILENAME (optional, derived automatically)

  Filename for SFF file containing sequences from this pool member.
  The SFF files are searched for in a subdirectory of the sff_dir
  named after the RUN_PREFIX.  If the field is blank or absent, a
  default value of {POOL_MEMBER_NAME}.sff is used.

  *Example*: ``F0FN7DX01_sample1_V1-V2.sff``

  *Output*: Used as the *filename* attribute of the <RUN> element in
  the SRA Run XML.  Also used by the QIIME software to find the SFF
  file and automatically generate a checksum.

BARCODE_READ_GROUP_TAG (optional, derived automatically)

  An identifier for each barcode within a pool.  If this field is
  absent, a value of {RUN_PREFIX}_{BARCODE} is derived automatically.

  *Example*: ``F0FN7DX01_AGCTAGCT``

  *Output*: Used as the *read_group_tag* attribute of a <BASECALL>
  element, and of a <READ_LABEL> element, in the SRA Experiment XML
  file.

LINKER (optional)

  Linker sequence between the primer and the barcode (to reduce
  differences in hybridization based on the barcode).  This field may
  be empty or absent.

  *Example*: ``GGCCAG``

  *Output*: Used as the text of a <BASECALL> element in the SRA
  Experiment XML file.

PRIMER (optional)

  Primer sequence that was used for this particular library member. If
  you used more than one primer for a given pool member (which is
  allowed) you need to duplicate the whole row with the additional
  primer information. This needs to be the actual sequence of the
  primer, not the name of the primer (i.e. not V2).

  *Example*: ``CTGCTGCCTYCCGTA``

  *Output*: Used as the text of a <BASECALL> element in the SRA
  Experiment XML file.  It is also used to derive several optional
  fields.

PRIMER_READ_GROUP_TAG (optional, derived automatically)

  Read group that samples will be assigned to based on the primer,
  e.g. V2 for the V2 primers. By default, multidimensional
  demultiplexing on the barcode and primer is performed.  If it is not
  present, this field will be derived using a table of standard primer
  read group tags.  If the primer is not found in the table, a
  KeyError is raised.

  *Example*: ``V1-V3``

  *Output*: Used as the *read_group_tag* attribute of a <BASECALL>
  element, and of a <READ_LABEL> element, in the SRA Experiment XML
  file.

LIBRARY_STRATEGY (optional, default value provided automatically)

  Sequencing technique intended for this library (optional
  field). This will usually be ``AMPLICON`` (default) or
  ``METAGENOMIC``.

  *Example*: ``AMPLICON``

  *Output*: Used as the text of the <LIBRARY_STRATEGY> element in the
  SRA Experiment XML file.

LIBRARY_SOURCE (optional, default value provided automatically)

  Type of source material that is being sequenced (optional
  field). This will usually be ``GENOMIC`` (default) or
  ``METAGENOMIC``.

  *Example*: ``GENOMIC``

  *Output*: Used as the text of the <LIBRARY_SOURCE> element in the
  SRA Experiment XML file.

LIBRARY_SELECTION (optional, default value provided automatically)

  Whether any method was used to select and/or enrich the material
  being sequenced (optional field). This is used in cases where
  e.g. the cells were sorted, if PCR was used to make a specific
  amplicon, if fractionation for viruses was done, etc.  The default
  value is ``PCR``.

  *Example*: ``PCR``

  *Output*: Used as the text of the <LIBRARY_SELECTION> element in the
  SRA Experiment XML file.

RUN_ACCESSION (optional, currently unused)

  Optional accession number for the run.  Leave blank or omit this
  field if not already assigned.

EXPERIMENT_ACCESSION (optional, currently unused)

  Optional accession number for the experiment.  If you already
  created the Experiment accession in SRA, use it -- otherwise, omit
  this field or leave it blank.

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

