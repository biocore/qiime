.. _merge_denoiser_output:

.. index:: merge_denoiser_output.py

*merge_denoiser_output.py* -- Merge the output of denoising step back into QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script combines the output of the denoising step with the OTU picker results.


**Usage:** :file:`merge_denoiser_output.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map_file
		Path to denoiser mapping file [default: None]
	-p, `-`-otu_picker_map_file
		Path to OTU picker mapping file [REQUIRED]
	-f, `-`-fasta_fp
		Path to fasta input file, output of `split_libraries.py <./split_libraries.html>`_ [REQUIRED]
	-d, `-`-denoised_fasta_fp
		Path to denoised fasta file [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: Denoiser_out_otu_picked/]


**Output:**




Merge the output of denoising (denoised_seqs.fasta and denoiser_mapping.txt) the OTU picker results on denoised_seqs.fasta (uclust_picked_otus/denoised_seqs_otus.txt) and replace the read IDs with the sampleIDs from the output of `split_libraries.py <./split_libraries.html>`_ (seqs.fna)

::

	merge_denoiser_output.py -f seqs.fna -d denoised_seqs.fasta -m denoiser_mapping.txt -p uclust_picked_otus/denoised_seqs_otus.txt


