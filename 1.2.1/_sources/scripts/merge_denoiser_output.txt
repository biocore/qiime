.. _merge_denoiser_output:

.. index:: merge_denoiser_output.py

*merge_denoiser_output.py* -- Merge the output of denoising step back into QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script combines the output of the denoising step with the OTU picker results.
This script has to be run after denoising and OTU picking to combine the denoiser clusters and OTU clusters.
The input to the script is:

   DENOISER_MAP_FILE: the cluster mapping from the denoiser, usually denoiser_mapping.txt in the denoiser output directory
   
   DENOISED_FASTA_FP: path to the output FASTA file of the denoiser (centroids.fasta and singletons.fasta in the denoiser output dir. Concatenate if you want to include both)

   OTU_PICKER_MAP_FILE: path to OTU mapping file from OTU picking on the

   FASTA_FP:  path to FASTA input file, thus has to be the output of `split_libraries.py <./split_libraries.html>`_
             

See the tutorial on 454 Denoising in the QIIME tutorials on how to use this script properly.



**Usage:** :file:`merge_denoiser_output.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map_file
		Path to denoiser cluster mapping file [REQUIRED]
	-p, `-`-otu_picker_map_file
		Path to OTU mapping file from OTU picker [REQUIRED]
	-f, `-`-fasta_fp
		Path to FASTA file, output of `split_libraries.py <./split_libraries.html>`_ [REQUIRED]
	-d, `-`-denoised_fasta_fp
		Path to denoised fasta file [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: Denoiser_out_otu_picked/]


**Output:**


Two files are generated:

1. A new otu mapping file with the combined clusters from denoising and OTU picking,

2. a FASTA file with one sequence per OTU.

These two files need to be provided to pick_rep_set using -m first


**Example usage:**

Merge the output of denoising (denoised_seqs.fasta and denoiser_mapping.txt) the OTU picker results on denoised_seqs.fasta (uclust_picked_otus/denoised_seqs_otus.txt) and replace the read IDs with the sampleIDs from the output of `split_libraries.py <./split_libraries.html>`_ (seqs.fna)



::

	merge_denoiser_output.py -f seqs.fna -d denoised_seqs.fasta -m denoiser_mapping.txt -p uclust_picked_otus/denoised_seqs_otus.txt


