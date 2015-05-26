.. _denoising_454_data:

Denoising of 454 Data Sets 
--------------------------

The pyrosequencing technology employed by 454 sequencing machines produces characteristic sequencing errors, mostly imprecise signals for longer homopolymers runs. Most of the sequences contain none or only a few errors, but a few sequences contain enough errors to be classified as an additional rare OTU. The goal for the denoising procedure is to reduce the amount of erroneous OTUs and thus increasing the accuracy of the whole QIIME pipeline.

If there are multiple, large 454 runs, follow this tutorial to denoise the data set and analyze it with QIIME. In short, each 454 run needs to be preprocessed with `split_libraries.py <../scripts/split_libraries.html>`_ and denoised separately. Afterwards the output files are combined for OTU picking. We will show an example with two 454 runs (:file:`run1.sff` and :file:`run2.sff`).

**Data preparation:**

From the raw, binary sff file, three files need to be generated for each run with the sffinfo tool from 454. You should have this tool if you have a 454 sequencer. Otherwise ask the sequencing facility for the files::

     sffinfo    run_1.sff > run_1.sff.txt
     sffinfo -s run_1.sff > run_1.fasta
     sffinfo -q run_1.sff > run_1.qual

     sffinfo    run_2.sff > run_2.sff.txt
     sffinfo -s run_2.sff > run_2.fasta
     sffinfo -q run_2.sff > run_2.qual


**Quality filtering and barcode assignment:**

Prior to denoising, each read has to be assigned to one barcode/sample
and low quality reads need to be filtered out. This can be done using
`split_libraries.py <../scripts/split_libraries.html>`_. An example command would be::

	split_libraries.py -o run1 -f run1.fasta -q run1.qual -m run1_mapping.txt -w 50 -r -l 150 -L 350
	split_libraries.py -o run2 -f run2.fasta -q run2.qual -m run2_mapping.txt -w 50 -r -l 150 -L 350 -n 1000000

This step has to be done separately for each 454 pool, following the usual guidelines for running several data sets through `split_libraries.py <../scripts/split_libraries.html>`_.


**Flowgram clustering (aka denoising)**

Each run will be denoised using its quality filtered output of `split_libraries.py <../scripts/split_libraries.html>`_ and the initial :file:`.sff.txt` file. All flowgrams without a match in the provided `split_libraries.py <../scripts/split_libraries.html>`_ FASTA file are removed. The sequencing primer will be extracted from the metadata mapping file::

	denoise_wrapper.py -v -i run1.sff.txt -f run1/seqs.fna -o run1/denoised/ -m run1_mapping.txt 
	denoise_wrapper.py -v -i run2.sff.txt -f run2/seqs.fna -o run2/denoised/ -m run2_mapping.txt


Denoising large data sets is computationally demanding. While smaller data sets (< 50,000 sequences) can be run on one single machine within an hour, a typical 454 run with 400,000 sequences after quality filtering requires up to a day on a 24 core cluster. If the denoiser is set up properly on your cluster or multi-core machine, it can be started in parallel mode using the option -n::

	denoise_wrapper.py -v -i run1.sff.txt -f run1/seqs.fna -o run1/denoised/ -m run1_mapping.txt -n 24



The output files of this step is stored in directory :file:`run1/` and :file:`run2/`, respectively:

	#. :file:`denoiser.log`: Information about the clustering procedure if run in verbose mode (-v). Can be used to monitor the program's progress.
	#. :file:`centroids.fasta`: The centroids of clusters with 2 and more members.
	#. :file:`singletons.fasta`: Reads that could not be clustered. 
	#. :file:`denoiser_mapping.txt`: The cluster to read mapping.

Usually the centroid and singleton files are combined for downstream analysis,
but occasionally it might make sense to remove the low confidence singletons.



**Re-integrating the denoised data into QIIME**

The final step in a denoising run usually is the re-integration of the data into the QIIME pipeline. Since the denoiser uses flowgram similarity for clustering there is no guaranteed sequence (dis)-similarity between cluster centroids. In order to create the usual species-level OTUs at 97% sequence similarity, you must inflate the denoiser results and then run one of QIIME's OTU pickers on the combined denoiser output.

Inflating denoiser results refers to process of creating a new fasta file of denoised sequences where each centroid sequence is written `n` times, where `n` is the cluster size, and each singleton is written once. Flowgram identifiers are mapped to sequence identifiers using the original input file.

To inflate the results of a single denoiser run call::

    inflate_denoiser_output.py -c centroids.fna -s singletons.fna -f seqs.fna -d denoiser_mapping.txt -o denoised_seqs.fna

To inflate the results from independent denoise_wrapper.py runs, pass all of the centroid, singleton, input fasta files, and denoiser maps::

    inflate_denoiser_output.py -c centroids1.fna,centroids2.fna -s singletons1.fna,singletons2.fna -f seqs1.fna,seqs2.fna -d denoiser_mapping1.txt,denoiser_mapping2.txt -o denoised_seqs.fna


Your denoised sequences can now be fed directly into QIIME at the OTU picking stage. The next step will be to run one of the OTU pickers or OTU picking workflow scripts (e.g., `pick_otus.py <../scripts/pick_otus.html>`_, `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_, `pick_reference_otus_through_otu_table.py <../scripts/pick_reference_otus_through_otu_table.html>`_, `core_qiime_analyses.py <../scripts/core_qiime_analyses.html>`_. At the OTU picking stage it is very important that you allow for the abundance presorting, which is currently in place for the uclust OTU picker only. We therefore don't recommend using other OTU pickers, and **do not pass the -D/--suppress_presort_by_abundance_uclust option to pick_otus.py**. If possible, it is worth using uclust with ``--optimal`` to assure the best possible choice of OTUs.::

    pick_otus.py -s 0.97 -i denoised_seqs.fna -m uclust --optimal

Passing ``--optimal`` may be prohibitively compute-intensive for large analyses however (for example, greater than a single 454 FLX run). The default QIIME pick_otus.py parameters are likely to be sufficient.


Notes:

* Denoising very small data sets might be ineffective, since there might not be a good read in the data set that can be used to correct a bad read. If there is a small data set (probably from re-sequencing an under-sampled sample) consider combining it with another, larger data set in your study prior to denoising.

* Currently only one sequencing primer per run is supported. If there is more than one primer the run needs to be split. Simply make per per-primer mapping files and run `split_libraries.py <../scripts/split_libraries.html>`_ with each mapping file, then denoise with each output FASTA file separately.

* Using any other OTU picker than uclust with the exact options as specified above might result in systematic differences between your separately denoised runs. Even small sequence differences in the denoiser output can lead to clustering into different OTUs and an artificial separation of samples. We warned you! 
  
