.. _denoiser_postprocess:

.. index:: denoiser_postprocess.py

*denoiser_postprocess.py* -- Merge denoiser output with OTU picker output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Combine mapping files from denoiser and OTU picker, such that we have a combined mapping file that can be used for the subsequent steps of Qiime. Also replace flowgram identifiers with IDs assigned by `split_libraries.py <./split_libraries.html>`_.


**Usage:** :file:`denoiser_postprocess.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-p, `-`-otu_picker_map_file
		Path to OTU picker mapping file [REQUIRED]
	-f, `-`-fasta_fp
		Path to fasta input file, output of `split_libraries.py <./split_libraries.html>`_ [REQUIRED]
	-d, `-`-denoised_fasta_fp
		Path to denoised fasta file [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: Denoiser_out_otu_picked/]
	-m, `-`-map_file
		Path to denoiser mapping file [default: denoiser_mapping.txt]


**Output:**

 The output of `denoiser_postprocess.py <./denoiser_postprocess.html>`_ are two files:

denoised_otu_map.txt: In this mapping the read/flowgram IDs are
                replaced by sample_id from the `split_libraries.py <./split_libraries.html>`_
                fasta file. Also, the lists for each OTU are sorted
                such that the largest cluster from denosing appears
                first. This will be important for the next step,
                picking representative sequences.

denoised_all.fasta: A fasta sequence where the header lines are
                updated with the sample_ids as assigned by `split_libraries.py <./split_libraries.html>`_.



Combine denoiser output with output of QIIME OTU picker, put results into Outdir:

::

	denoiser_postprocess.py -f seqs.fna -d denoised.fasta -m denoiser_mapping.txt -p cdhit_picked_otus/denoised_otus.txt -v -o Outdir


