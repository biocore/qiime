.. _split_libraries_lea_seq:

.. index:: split_libraries_lea_seq.py

*split_libraries_lea_seq.py* -- Demultiplexes Low-Error Amplicon Sequencing (LEA-Seq) data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


Implements Low-Error Amplicon Sequencing (LEA-Seq) method, described in:

Faith, Jeremiah J., et al.
The long-term stability of the human gut microbiota.Science 341.6141 (2013).

This method is based on redundant sequencing of a set of linear PCR template
extensions of 16S rRNA genes. The oligonucleotide primer that is used for
PCR template extensions is labeled with a random barcode
5' to the universal 16S rRNA primer sequence. This PCR pool is then
amplified with exponential PCR, using primers that specifically
amplify only the linear PCR molecules. An index primer is added to
the amplicons along with a primer specific for each sample.
This exponential PCR pool is then sequenced redundantly (20x coverage).
The resulting sequences are separated by sample, using the index sequence.
The amplicon sequences within each sample are separated by the random
barcodes. The large number of reads for each barcode helps to create an
error-corrected consensus sequence for the initial template molecule.



**Usage:** :file:`split_libraries_lea_seq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-sequence_read_fps
		The forward and reverse sequence read fastq files (comma-separated)
	-o, `-`-output_dir
		Directory to store output files
	-m, `-`-mapping_fp
		Metadata mapping file
	
	**[OPTIONAL]**
		
	-b, `-`-barcode_type
		The type of barcode used. This can be an integer, e.g. 6 for length 6 barcodes, or golay_12 for golay error-correcting barcodes. Error correction will only be applied for golay_12 barcodes [default: golay_12]
	`-`-max_barcode_errors
		The maximum allowable number of errors in the barcode if passing --barcode_type golay_12 [default: 1.5]
	`-`-min_consensus
		Threshold for consensus score: the minimum score allowable at any position in sequence. where the score is calulated as: occurence of base in consensus sequence/ total sequences[default: 6.6]
	`-`-max_cluster_ratio
		Threshold for cluster ratio: the maximum allowable cluster ratio above which you need to find the consensus sequence for the given sequences.[default: 2.5]
	`-`-min_difference_in_bcs
		Threshold for selecting unique barcodes: Barcodes that are more similar to each other than this value will be discarded.[default: 0.86]
	`-`-fwd_length
		Removes phasing from forward readby truncating it to standard length for the region[default: 64]
	`-`-rev_length
		Removes phasing from reverse readby truncating it to standard length for the region[default: 77]
	`-`-min_difference_in_clusters
		The percent identity threshold while using uclust to cluster sequence reads, which is helpfulin measuring quality of sequencing.[default: 0.98]
	`-`-min_reads_per_random_bc
		Minimum number of reads per randombarcode, attempts to remove random barcodes that are sequencing errors of true barcodesmight be useful in saving memory and time[default: 1]
	`-`-header_barcode_column
		Header of barcode column[default: BarcodeSequence]
	`-`-reverse_primer_column
		Header of reverse primer column[default: ReversePrimer]


**Output:**

The split_libraries_lea_seq.py generates: A fasta file called seqs.fna which contains error corrected consensus sequence for the template DNA


**General Example: Specify forward read and reverse read fasta files, use the metadata mapping file map.txt,and output the data to output_dir:**

output_dir

::

	split_libraries_lea_seq.py -i fwd_read.fq,rev_read.fq -m map.txt -o output --b 7


