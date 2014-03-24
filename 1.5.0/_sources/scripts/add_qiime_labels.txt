.. _add_qiime_labels:

.. index:: add_qiime_labels.py

*add_qiime_labels.py* -- Takes a directory and a mapping file of SampleIDs to fasta file names, combines all files that have valid fasta extensions into a single fasta file, with valid QIIME fasta labels.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A tab separated text file with SampleIDs 
and fasta file names (just the file name itself, not the full or relative 
filepath) is used to generate a combined fasta file with valid
QIIME labels based upon the SampleIDs specified in the mapping file.

This script is to handle situations where fasta data comes already 
demultiplexed into a one fasta file per sample basis.  Apart from altering
the fasta label to add a QIIME compatible label at the beginning (example:
>FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
could become 
>control.sample_1 FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_

Note that limited checking is done on the mapping file.  The only tests
are that every fasta file name is unique, and that SampleIDs are
MIMARKS compliant (alphanumeric and period characters only).  Duplicate 
SampleIDs are allowed, so care should be taken that there are no typos.

No changes are made to the sequences.



**Usage:** :file:`add_qiime_labels.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		SampleID to fasta file name mapping file filepath
	-i, `-`-fasta_dir
		Directory of fasta files to combine and label.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Required output directory for log file and corrected mapping file, log file, and html file. [default: ./]
	-n, `-`-count_start
		Specify the number to start enumerating sequence labels with. [default: 0]


**Output:**

A combined_seqs.fasta file will be created in the output directory


**Example:**

Specify fasta_dir as the input directory of fasta files, use the SampleID to fasta file mapping file example_mapping.txt, start enumerating with 1000000 following SampleIDs, and output the data to the directory combined_fasta

::

	add_qiime_labels.py -i fasta_dir -m example_mapping.txt -n 1000000 -o combined_fasta


