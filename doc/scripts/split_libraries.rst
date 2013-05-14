.. _split_libraries:

.. index:: split_libraries.py

*split_libraries.py* -- Split libraries according to barcodes specified in mapping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Since newer sequencing technologies provide many reads per run (e.g. the 454 GS FLX Titanium series can produce 400-600 million base pairs with 400-500 base pair read lengths) researchers are now finding it useful to combine multiple samples into a single 454 run. This multiplexing is achieved through the application of a pyrosequencing-tailored nucleotide barcode design (described in (Parameswaran et al., 2007)). By assigning individual, unique sample specific barcodes, multiple sequencing runs may be performed in parallel and the resulting reads can later be binned according to sample. The script split_libraries.py performs this task, in addition to several quality filtering steps including user defined cut-offs for: sequence lengths; end-trimming; minimum quality score. To summarize, by using the fasta, mapping, and quality files, the program split_libraries.py will parse sequences that meet user defined quality thresholds and then rename each read with the appropriate Sample ID, thus formatting the sequence data for downstream analysis. If a combination of different sequencing technologies are used in any particular study, split_libraries.py can be used to perform the quality-filtering for each library individually and the output may then be combined.

Sequences from samples that are not found in the mapping file (no corresponding barcode) and sequences without the correct primer sequence will be excluded. Additional scripts can be used to exclude sequences that match a given reference sequence (e.g. the human genome; `exclude_seqs_by_blast.py <./exclude_seqs_by_blast.html>`_) and/or sequences that are flagged as chimeras (`identify_chimeric_seqs.py <./identify_chimeric_seqs.html>`_).



**Usage:** :file:`split_libraries.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map
		Name of mapping file. NOTE: Must contain a header line indicating SampleID in the first column and BarcodeSequence in the second, LinkerPrimerSequence in the third.
	-f, `-`-fasta
		Names of fasta files, comma-delimited
	
	**[OPTIONAL]**
		
	-q, `-`-qual
		Names of qual files, comma-delimited [default: None]
	-r, `-`-remove_unassigned
		DEPRECATED: pass --retain_unassigned_reads to keep unassigned reads  [default: None]
	-l, `-`-min-seq-length
		Minimum sequence length, in nucleotides [default: 200]
	-L, `-`-max-seq-length
		Maximum sequence length, in nucleotides [default: 1000]
	-t, `-`-trim-seq-length
		Calculate sequence lengths after trimming primers and barcodes [default: False]
	-s, `-`-min-qual-score
		Min average qual score allowed in read [default: 25]
	-k, `-`-keep-primer
		Do not remove primer from sequences
	-B, `-`-keep-barcode
		Do not remove barcode from sequences
	-a, `-`-max-ambig
		Maximum number of ambiguous bases [default: 6]
	-H, `-`-max-homopolymer
		Maximum length of homopolymer run [default: 6]
	-M, `-`-max-primer-mismatch
		Maximum number of primer mismatches [default: 0]
	-b, `-`-barcode-type
		Barcode type, hamming_8, golay_12, variable_length (will disable any barcode correction if variable_length set), or a number representing the length of the barcode, such as -b 4.  [default: golay_12]
	-o, `-`-dir-prefix
		Directory prefix for output files [default: .]
	-e, `-`-max-barcode-errors
		Maximum number of errors in barcode [default: 1.5]
	-n, `-`-start-numbering-at
		Seq id to use for the first sequence [default: 1]
	`-`-retain_unassigned_reads
		Retain sequences which are Unassigned in the output sequence file[default: False]
	-c, `-`-disable_bc_correction
		Disable attempts to find nearest corrected barcode.  Can improve performance. [default: False]
	-w, `-`-qual_score_window
		Enable sliding window test of quality scores.  If the average score of a continuous set of w nucleotides falls below the threshold (see -s for default), the sequence is discarded. A good value would be 50. 0 (zero) means no filtering. Must pass a .qual file (see -q parameter) if this functionality is enabled.  Default behavior for this function is to truncate the sequence at the beginning of the poor quality window, and test for minimal length (-l parameter) of the resulting sequence. [default: 0]
	-g, `-`-discard_bad_windows
		If the qual_score_window option (-w) is enabled, this will override the default truncation behavior and discard any sequences where a bad window is found.  [default: False]
	-p, `-`-disable_primers
		Disable primer usage when demultiplexing.  Should be enabled for unusual circumstances, such as analyzing Sanger sequence data generated with different primers.  [default: False]
	-z, `-`-reverse_primers
		Enable removal of the reverse primer and any subsequence sequence from the end of each read.  To enable this, there has to be a "ReversePrimer" column in the mapping file. Primers a required to be in IUPAC format and written in the 5' to  3' direction.  Valid options are 'disable', 'truncate_only', and 'truncate_remove'.  'truncate_only' will remove the primer and subsequent sequence data from the output read and will not alter output of sequences where the primer cannot be found. 'truncate_remove' will flag sequences where the primer cannot be found to not be written and will record the quantity of such failed sequences in the log file. [default: disable]
	`-`-reverse_primer_mismatches
		Set number of allowed mismatches for reverse primers (option -z). [default: 0]
	-d, `-`-record_qual_scores
		Enables recording of quality scores for all sequences that are recorded.  If this option is enabled, a file named seqs_filtered.qual will be created in the output directory, and will contain the same sequence IDs in the seqs.fna file and sequence quality scores matching the bases present in the seqs.fna file. [default: False]
	-i, `-`-median_length_filtering
		Disables minimum and maximum sequence length filtering, and instead calculates the median sequence length and filters the sequences based upon the number of median absolute deviations specified by this parameter.  Any sequences with lengths outside the number of deviations will be removed. [default: None]
	-j, `-`-added_demultiplex_field
		Use -j to add a field to use in the mapping file as an additional demultiplexing option to the barcode.  All combinations of barcodes and the values in these fields must be unique. The fields must contain values that can be parsed from the fasta labels such as "plate=R_2008_12_09".  In this case, "plate" would be the column header and "R_2008_12_09" would be the field data (minus quotes) in the mapping file.  To use the run prefix from the fasta label, such as ">FLP3FBN01ELBSX", where "FLP3FBN01" is generated from the run ID, use "-j run_prefix" and set the run prefix to be used as the data under the column headerr "run_prefix".  [default: None]
	-x, `-`-truncate_ambi_bases
		Enable to truncate at the first "N" character encountered in the sequences.  This will disable testing for ambiguous bases (-a option) [default: False]


**Output:**

Three files are generated by split_libraries.py:

1. .fna file (e.g. seqs.fna) - This is a FASTA file containing all sequences which meet the user-defined parameters, where each sequence identifier now contains its corresponding sample id from mapping file.

2. histograms.txt- This contains the counts of sequences with a particular length.

3. split_library_log.txt - This file contains a summary of the split_libraries.py analysis. Specifically, this file includes information regarding the number of sequences that pass quality control (number of seqs written) and how these are distributed across the different samples which, through the use of bar-coding technology, would have been pooled into a single 454 run. The number of sequences that pass quality control will depend on length restrictions, number of ambiguous bases, max homopolymer runs, barcode check, etc. All of these parameters are summarized in this file. If raw sequences do not meet the specified quality thresholds they will be omitted from downstream analysis. Since we never see a perfect 454 sequencing run, the number of sequences written should always be less than the number of raw sequences. The number of sequences that are retained for analysis will depend on the quality of the 454 run itself in addition to the default data filtering thresholds in the split_libraries.py\dir_prefix
dir_prefix
 script. The default parameters (minimum quality score = 25, minimum/maximum length = 200/1000, no ambiguous bases allowed, no mismatches allowed in primer sequence) can be adjusted to meet the user's needs.



**Standard Example:**

Using a single 454 run, which contains a single FASTA, QUAL, and mapping file while using default parameters and outputting the data into the Directory "Split_Library_Output":

::

	split_libraries.py -m Mapping_File.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o Split_Library_Output/

**Multiple FASTA and QUAL Files Example:**

For the case where there are multiple FASTA and QUAL files, the user can run the following comma-separated command as long as there are not duplicate barcodes listed in the mapping file:

::

	split_libraries.py -m Mapping_File.txt -f 1.TCA.454Reads.fna,2.TCA.454Reads.fna -q 1.TCA.454Reads.qual,2.TCA.454Reads.qual -o Split_Library_Output_comma_separated/

**Duplicate Barcode Example:**

An example of this situation would be a study with 1200 samples. You wish to have 400 samples per run, so you split the analysis into three runs and reuse barcoded primers (you only have 600). After initial analysis you determine a small subset is underrepresented (<500 sequences per samples) and you boost the number of sequences per sample for this subset by running a fourth run. Since the same sample IDs are in more than one run, it is likely that some sequences will be assigned the same unique identifier by split_libraries.py when it is run separately on the four different runs, each with their own barcode file. This will cause a problem in file concatenation of the four different runs into a single large file. To avoid this, you can use the '-n' parameter which defines a start index for split_libraries.py. From experience, most FLX runs (when combining both files for a single plate) will have 350,000 to 650,000 sequences. Thus, if Run 1 for split_libraries.py uses '-n 1000000', Run 2 uses '-n 2000000', etc., then you are guaranteed to have unique identifiers after concatenating the results of multiple FLX runs. With newer technologies you will just need to make sure that your start index spacing is greater than the potential number of sequences.

To run split_libraries.py, you will need two or more (depending on the number of times the barcodes were reused) separate mapping files (one for each Run, for example one for Run1 and another one for Run2), then you can run split_libraries.py using the FASTA and mapping file for Run1 and FASTA and mapping file for Run2. Once you have run split libraries on each file independently, you can concatenate (e.g. using the 'cat' command) the sequence files that were generated by split_libraries.py. You can also concatenate the mapping files, since the barcodes are not necessary for downstream analyses, unless the same sample IDs are found in multiple mapping files.

Run split_libraries.py on Run 1:

::

	split_libraries.py -m Mapping_File.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o Split_Library_Run1_Output/ -n 1000000

Run split_libraries.py on Run 2. The resulting FASTA files from Run 1 and Run 2 can then be concatenated using the 'cat' command (e.g. cat Split_Library_Run1_Output/seqs.fna Split_Library_Run2_Output/seqs.fna > Combined_seqs.fna) and used in downstream analyses.

::

	split_libraries.py -m Mapping_File.txt -f 2.TCA.454Reads.fna -q 2.TCA.454Reads.qual -o Split_Library_Run2_Output/ -n 2000000

**Barcode Decoding Example:**

The standard barcode types supported by split_libraries.py are golay (Length: 12 NTs) and hamming (Length: 8 NTs). For situations where the barcodes are of a different length than golay and hamming, the user can define a generic barcode type "-b" as an integer, where the integer is the length of the barcode used in the study.

Note: When analyzing large datasets (>100,000 seqs), users may want to use a generic barcode type, even for length 8 and 12 NTs, since the golay and hamming decoding processes can be computationally intensive, which causes the script to run slow. Barcode correction can be disabled with the -c option if desired.

For the case where the 8 base pair barcodes were used, you can use the following command:

::

	split_libraries.py -m Mapping_File_8bp_barcodes.txt -f 1.TCA.454Reads.fna  -q 1.TCA.454Reads.qual -o split_Library_output_8bp/ -b 8

**Linkers and Primers:**

The linker and primer sequence (or all the degenerate possibilities) are associated with each barcode from the mapping file. If a barcode cannot be identified, all the possible primers in the mapping file are tested to find a matching sequence. Using truncated forms of the same primer can lead to unexpected results for rare circumstances where the barcode cannot be identified and the sequence following the barcode matches multiple primers.

In many cases, sequence reads are long enough to sequence through the reverse primer and sequencing adapter.  To remove these primers and all following sequences, the -z option can be used.  By default, this option is set to 'disable'.  If it is set to 'truncate_only', split_libraries will trim the primer and any sequence following it if the primer is found.  If the 'truncate_remove' option is set, split_libraries.py will trim the primer if found, and will not write the sequence if the primer is not found. The allowed mismatches for the reverse primer are set with the --reverse_primer_mismatches parameter (default 0).  To use reverse primer removal, one must include a 'ReversePrimer' column in the mapping file, with the reverse primer recorded in the 5' to 3' orientation.

Example reverse primer removal, where primers are trimmed if found, and sequence is written unchanged if not found.  Mismatches are increased to 1 from the default 0:

::

	split_libraries.py -m Mapping_File_reverse_primer.txt -f 1.TCA.454Reads.fna -q 1.TCA.454Reads.qual -o split_libraries_output_revprimer/ --reverse_primer_mismatches 1 -z truncate_only


