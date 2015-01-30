.. _add_qiime_labels:

.. index:: add_qiime_labels.py

*add_qiime_labels.py* -- Takes a directory, a metadata mapping file, and a column name that contains the fasta file names that SampleIDs are associated with, combines all files that have valid fasta extensions into a single fasta file, with valid QIIME fasta labels.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A metadata mapping file with SampleIDs
and fasta file names (just the file name itself, not the full or relative
filepath) is used to generate a combined fasta file with valid
QIIME labels based upon the SampleIDs specified in the mapping file.

See: http://qiime.org/documentation/file_formats.html#metadata-mapping-files
for details about the metadata file format.

Example mapping file:
#SampleID	BarcodeSequence	LinkerPrimerSequence	InputFileName	Description
Sample.1	AAAACCCCGGGG	CTACATAATCGGRATT	seqs1.fna	sample.1
Sample.2	TTTTGGGGAAAA	CTACATAATCGGRATT	seqs2.fna	sample.2

This script is to handle situations where fasta data comes already
demultiplexed into a one fasta file per sample basis.  Only alters
the fasta label to add a QIIME compatible label at the beginning.

Example:
With the metadata mapping file above, and an specified directory containing the
files seqs1.fna and seqs2.fna, the first line from the seqs1.fna file might
look like this:
>FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
AACAGATTAGACCAGATTAAGCCGAGATTTACCCGA

and in the output combined fasta file would be written like this
>Sample.1_0 FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
AACAGATTAGACCAGATTAAGCCGAGATTTACCCGA

No changes are made to the sequences.



**Usage:** :file:`add_qiime_labels.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		SampleID to fasta file name mapping file filepath
	-i, `-`-fasta_dir
		Directory of fasta files to combine and label.
	-c, `-`-filename_column
		Specify column used in metadata mapping file for fasta file names.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Required output directory for log file and corrected mapping file, log file, and html file. [default: .]
	-n, `-`-count_start
		Specify the number to start enumerating sequence labels with. [default: 0]


**Output:**

A combined_seqs.fasta file will be created in the output directory, with the sequences assigned to the SampleID given in the metadata mapping file.


**Example:**

Specify fasta_dir as the input directory of fasta files, use the metadata mapping file example_mapping.txt, with the metadata fasta file name column specified as InputFileName, start enumerating with 1000000, and output the data to the directory combined_fasta

::

	add_qiime_labels.py -i fasta_dir -m example_mapping.txt -c InputFileName -n 1000000 -o combined_fasta


