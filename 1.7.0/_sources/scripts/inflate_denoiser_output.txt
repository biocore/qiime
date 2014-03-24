.. _inflate_denoiser_output:

.. index:: inflate_denoiser_output.py

*inflate_denoiser_output.py* -- Inflate denoiser results so they can be passed directly to OTU pickers.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Inflate denoiser results so they can be passed directly to `pick_otus.py <./pick_otus.html>`_, `parallel_pick_otus_uclust_ref.py <./parallel_pick_otus_uclust_ref.html>`_, or `pick_de_novo_otus.py <./pick_de_novo_otus.html>`_. Note that the results of this script have not be abundance sorted, so they must be before being passed to the OTU picker. The uclust OTU pickers incorporate this abundance presorting by default.

The inflation process writes each centroid sequence n times, where n is the number of reads that cluster to that centroid, and writes each singleton once. Flowgram identifiers are mapped back to post-split_libraries identifiers in this process (i.e., identifiers in fasta fps).



**Usage:** :file:`inflate_denoiser_output.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-c, `-`-centroid_fps
		The centroid fasta filepaths
	-s, `-`-singleton_fps
		The singleton fasta filepaths
	-f, `-`-fasta_fps
		The input (to denoiser) fasta filepaths
	-d, `-`-denoiser_map_fps
		The denoiser map filepaths
	-o, `-`-output_fasta_fp
		The output fasta filepath


**Output:**




Inflate the results of a single denoiser run.

::

	inflate_denoiser_output.py -c centroids.fasta -s singletons.fasta -f seqs.fna -d denoiser_mapping.txt -o inflated_seqs.fna

Inflate the results of multiple denoiser runs to a single inflated_seqs.fna file.

::

	inflate_denoiser_output.py -c centroids1.fasta,centroids2.fasta -s singletons1.fasta,singletons2.fasta -f seqs1.fna,seqs2.fna -d denoiser_mapping1.txt,denoiser_mapping2.txt -o inflated_seqs_combined.fna


