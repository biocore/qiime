.. _denoising_454_data:

============================
 Denoising of 454 Data Sets
============================

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

Note that the qiime since package v1.2 has a replacement for the sfftools.
It's slower but fully functional.

.. warning:: Warning: Since late 2012, 454 machines have a new feature (flow pattern B) that is supposed to allow for longer reads. Unfortunately, files using this feature can not be denoised, but result in nonsense output. To make sure that your file uses the older, more common flow pattern A, open the .sff.txt file and look for the ``Flow Chars:`` section in the header. If it shows a constant repeat of TACG you are fine. If however the pattern deviates after the third repeat, you are looking at the new flow pattern B that can not be denoised. In any case, all other qiime programs are not affected by this and can be used as usual.

For more details on the available options of each script explained in
the following use the -h option.

**Quality filtering and barcode assignment:**

Prior to denoising, each read has to be assigned to one barcode/sample
and low quality reads need to be filtered out. This can be done using
`split_libraries.py <../scripts/split_libraries.html>`_. An example command would be::

	split_libraries.py -o run1 -f run1.fasta -q run1.qual -m run1_mapping.txt -w 50 -g -r -l 150 -L 350
	split_libraries.py -o run2 -f run2.fasta -q run2.qual -m run2_mapping.txt -w 50 -g -r -l 150 -L 350 -n 1000000

This step has to be done separately for each 454 pool, following the usual guidelines for running several data sets through `split_libraries.py <../scripts/split_libraries.html>`_. Note that all options to `split_libraries.py <../scripts/split_libraries.html>`_ that truncate the sequences on the 3' end should not be used as they do not affect the sff.txt files used for denoising. This includes the ``-x``, ``-z truncate_only``, and ``-w`` without ``-g`` options. We recommend though to use the ``-w 50 -g`` combination to discard reads of bad quality.
Also, do not use the `truncate_fasta_qual_files.py  <../scripts/truncate_fasta_qual_files.html>`_ script if you plan to denoise your data.
If you need to truncate your data, use the ``sfffile`` program from the Roche sfftools package and recreate your fasta and qual files from the truncated sff file.

For a single, non-barcoded sample, `split_libraries.py <../scripts/split_libraries.html>`_
can be provided with a mapping file that has an empty field for the BarcodeSequence.

Example:

.. note::

   * #SampleID   BarcodeSequence	LinkerPrimerSequence	 Description
   * Artificial    			ATTAGATACCCNGGTAG	 ArtificialGSFLX_from_Quince_et_al

Note that fields must be separated by a single tab. For the empty barcode there must be two
tabs between SampleID and the primer sequence. Use QIIME's
`check_id_map.py <../scripts/check_id_map.html>`_ to test for validity. Then, use
`split_libraries.py <../scripts/split_libraries.html>`_ as usual, but with
option -b 0.

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
        #. :file:`denoised_clusters.txt`: A cluster mapping in qiime format. Equivalent to 4.
        #. :file:`denoised_seqs.fasta`: Centroids and singletons combined and sorted by cluster size.
        
Usually the centroid and singleton files are combined for downstream analysis,
but occasionally it might make sense to remove the low confidence singletons.
2, 3, and 4 are used as input to the next step.


**Re-integrating the denoised data into QIIME**

The final step in a denoising run usually is the re-integration of the data into the QIIME pipeline. Since the denoiser uses flowgram similarity for clustering there is no guaranteed sequence (dis)-similarity between cluster centroids. In order to create the usual species-level OTUs at 97% sequence similarity, you must inflate the denoiser results and then run one of QIIME's OTU pickers on the combined denoiser output.

Inflating denoiser results refers to process of creating a new fasta file of denoised sequences where each centroid sequence is written `n` times, where `n` is the cluster size, and each singleton is written once. Flowgram identifiers are mapped to sequence identifiers using the original input file.

To inflate the results of a single denoiser run call::

    inflate_denoiser_output.py -c centroids.fna -s singletons.fna -f seqs.fna -d denoiser_mapping.txt -o denoised_seqs.fna

To inflate the results from independent denoise_wrapper.py runs, pass all of the centroid, singleton, input fasta files, and denoiser maps::

    inflate_denoiser_output.py -c centroids1.fna,centroids2.fna -s singletons1.fna,singletons2.fna -f seqs1.fna,seqs2.fna -d denoiser_mapping1.txt,denoiser_mapping2.txt -o denoised_seqs.fna


Your denoised sequences can now be fed directly into QIIME at the OTU picking stage. The next step will be to run one of the OTU pickers or OTU picking workflow scripts (e.g., `pick_otus.py <../scripts/pick_otus.html>`_, `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, or `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_,. At the OTU picking stage it is very important that you allow for the abundance presorting, which is currently in place for the uclust OTU picker only. We therefore don't recommend using other OTU pickers, and **do not pass the -D/--suppress_presort_by_abundance_uclust option to pick_otus.py**. If possible, it is worth using uclust with ``--optimal`` to assure the best possible choice of OTUs.::

    pick_otus.py -s 0.97 -i denoised_seqs.fna -m uclust --optimal

Passing ``--optimal`` may be prohibitively compute-intensive for large analyses however (for example, greater than a single 454 FLX run). The default QIIME pick_otus.py parameters are likely to be sufficient.


Notes:

* Denoising very small data sets might be ineffective, since there might not be a good read in the data set that can be used to correct a bad read. If there is a small data set (probably from re-sequencing an under-sampled sample) consider combining it with another, larger data set in your study prior to denoising.

* Currently only one sequencing primer per run is supported. If there is more than one primer the run needs to be split. Simply make per per-primer mapping files and run `split_libraries.py <../scripts/split_libraries.html>`_ with each mapping file, then denoise with each output FASTA file separately.

* Using any other OTU picker than uclust with the exact options as specified above might result in systematic differences between your separately denoised runs. Even small sequence differences in the denoiser output can lead to clustering into different OTUs and an artificial separation of samples. We warned you! 
  


**Low-level Interface**

`denoise_wrapper.py <../scripts/denoise_wrapper.html>`_ provides an easy to use interface to the denoiser, which is sufficient in most cases.
For power users, we also provide two low level scripts, that allow for more flexibility.

*Cluster phase 1 - prefix clustering*

All flowgrams corresponding to the sequences that are in :file:`seqs.fna`
(presumed to be the output of `split_libraries.py <../scripts/split_libraries.html>`_)
are pulled from the .sff.txt file and primer, barcodes and
the 454 key sequence are removed. Then, the first clustering phase
groups reads based on common prefixes. For a full FLX run this will
usually take less than an hour on a standard computer and requires
less than 1 GB of memory.
 

Example command::

	denoiser_preprocess.py -i 454Reads.sff.txt -f seqs.fna -o example_pp -s -v -p CATGCTGCCTCCCGTAGGAGT

Several files are stored in the specified output directory. To see the
clustering stastics check the file :file:`preprocess.log` in the output
directory. Basically the less clusters there are (especially small
clusters) the faster the next phase  will run. If there are more than
100.000 sequences remaining, the input set might be split, to achieve
a reasonable run time. The files in the output directory are used in
the next step.



*Cluster phase II - Flowgram clustering or Denoising*

This is the main clustering step and the computationally most expensive one. 
Flowgrams are clustered based on their similarity.

Example command::

	denoiser.py -i 454Reads.sff.txt -p example_pp -v -o example_denoised

The preprocessing information in :file:`example_pp` is used and the output is
stored in a randomly named, new direcory in :file:`example_denoised`. Note, that
when the -p option is not specified here,  the preprocessing is invoked
from `denoiser.py <../scripts/denoiser.html>`_ implicitly.

Because of the potential long runtime, we suggest to distribute the work over
many cpus. If you have a multi-core system or cluster available and
set up the required job submission script (:file:`cluster_jobs_fp` in your qiime config)
the following command will distribute the computation over 24
cpus::

	denoiser.py -i 454Reads.sff.txt -p example_pp -v -o example_denoised -c -n 24

Make sure the output directory is shared by all cluster
nodes. Depending on the complexity of the data this step might take up
to a day even on a 24 core system for a full 454 run with 400-500 k
sequences. Smaller data sets will be finished much faster. The output
will be written to a randomly named directory within the specified
output directory. 
The output files are:

* denoiser.log: Information about the clustering procedure if run in verbose mode (-v).
	      	     Can be used to monitor the program's progress.

* centroids.fasta: The centroids of clusters with 2 and more members

* singletons.fasta: Reads that could not be clustered. 

* denoiser_mapping.txt: The cluster to read mapping.

Usually the centroid and singleton files are combined for downstream analysis,
but occasionally it might make sense to remove the low confidence singletons.




**Notes for running on cluster/multicore system**

We use a very simple setup to farm out the flowgram alignments to a cluster.
A master process (`denoiser.py <../scripts/denoiser.html>`_) sends data to each worker
(`denoiser_worker.py <../scripts/denoiser_worker.html>`_).
A worker sleeps while waiting for the data. Once the file appears it processes it and
sends the result back to the master and goes back to sleep. The master collects all results
and iterates. As such, performance is higly dependent on the actual cluster setup:

* The overall speed is governed by the slowest worker node
* The parallel steps will only start when all worker jobs are established. That means as long
   as one jobs remains queued, the other jobs will block your cluster. Decrease the number of workers
   if you run into this problem.


**FAQ**

Q: How does this denoising procedure differ from PyroNoise?

Q: What is the expected run-time?

Q: Can I denoise Titanium data

Q: How can I speed up the computation?

Q: Why are there so few sequences in my output file after denoising? Did something went wrong with my sequencing run?

Q: So where are all the sequences then?

Q: Can I cluster at different sequence/flowgram similarity thresholds?

Q: Denoising on the clusters "hangs" after a while. What is going on?

Q: How and why can I run the preprocessing step separately?

Q: What about different next-gen sequencing platforms?

---


Q: How does this denoising procedure differ from PyroNoise?

A: PyroNoise uses an expectation maximization (EM) algorithm to figure out the most likely sequence for every read. We, instead, use a greedy scheme that can be seen as an approximation to PyroNoise. According to several test data sets, our approximation gives very similar results in a fraction of the time.


Q: What is the expected run-time?

A: The whole heuristic for our method depends on the actual species distribution in your samples.
An ideal data set has few species and a very skewed abundance distribution with a few, very abundant species.
With more species and a flatter abundance distribution run time increases. You can get a rough estimate of the run time after the preprocessing step by looking at the number of reads printed in the log file in verbose mode. Very, very roughly, compute time increases quadratically with the number of reads after preprocessing:

.. note::
    * ...
    * Prefix matching: removed 242038 out of 339647 seqs
    * Remaining number of sequences: 97609
    * ...

If the number of remaining sequences is smaller than 50.000, you can expect <24 hours on 20 cpus.
With 100k seqs you would need 80 cpus to expect it to finish within a day.

Here are some guidelines from runs with actual data:

- partial GSFLX run with 50.000 reads: ~ 1 hour on a single CPU

- Full GSFLX run (~400.000 reads):   6-24 hours on 24 CPUs

- 1/2 Titanium run (450k reads):   35 hours on 200 CPUs

Titanium data takes longer for two reasons:
 a) Reads are longer, meaning longer alignment times
 b) We observed a higher variability in the Titanium reads, especially towards the 3'end leading to a less efficient greedy clustering.


Q: Can I denoise Titanium data?

A: Yes. The algorithm can process Titanium data and we have done it several times. As of (denoiser) version 0.9/Qiime-1.2 we ship an error profile for the titanium platform with the package. Use the switch --titanium to enable the new profile. Be aware that Titanium still takes considerably longer than FLX.

Q: How can I speed up the computation?

A:
1. Use more CPUs if available.

2. Stop clustering early.
   Clustering phase II processes clusters in decreasing order of their size after cluster phase I. As default, the procedure stops with the first singleton cluster being considered as cluster centroid. Setting -b 3 would stop the clustering with clusters of size 3. Note that setting the -b parameter does not hinder these cluster to be dragged in by another, larger cluster either in phase II or phase III. It just limits their role as cluster centroid.

3. Split your data in smaller pieces. 
   For very large data sets(>1 FLX plate), this is the recommended way to go. While we have observed that splitting into too small pieces (e.g. per sample with 5k sequences/sample) might render the denoising less effective, we expect very little difference when denoising is performed on larger chunks of data (100k+ reads). We recommend pooling similar samples, e.g. time series samples from the same person, but encourage to separate samples from different habitats with expected very different communities.

4. As a rather desperate measure for people who have to limit the compute time we provide a new flag in version 0.9 that controls the maximum number of rounds that the greedy clustering should run for. Note that the lower this number is, the worse the final clustering result can be. 

Q: Why are there so few sequences in my output file after denoising? Did something went wrong with my sequencing run?

A: No, this is expected. The denoising procedure (and this also holds for Chris Quince's Pyronoise) technically do not remove any reads from the input set. This is the task of the initial quality filtering, which we suggest to do using Qiime's split_libraries.py. The denoising is basically a clustering approach on the flowgram level, i.e. all reads that look similar enough on the flowgram level are clustered and only the centroid of each cluster is reported in the output file (either in centroids.fasta if the cluster has more than one member or otherwise in singletons.fasta). You can think of the centroids as OTUs on the flowgram level. Since flowgram similarity does not correlate perfectly with sequence-similarity, we usually don't call them OTUs, but only after an extra OTU picking step with, say, cd-hit or uclust on the denoised sequences.


Q: So where are all the sequences then?

A: If you look at the file denoiser_mapping.txt, e.g. like this::

	wc denoiser_mapping.txt

you should see that the number in the middle of the output (i.e. the number of words) is about the number of sequences in your input set. (Sometimes, the denoiser discards a few additional reads due to quality issues that were not captured by split_libraries.py). All reads that are in this mapping file can and will be used e.g. in the downstream Qiime analysis. The first number in the wc output gives the number of lines on the files, which corresponds to the number of clusters after denoising.



Q: Can I cluster at different sequence/flowgram similarity thresholds?

A: Basically, Yes. The default clustering parameters are set and tested to work well at 0.97% sequence similarity. If you want to cluster at, say, 0.95% you have to increase both cut-offs and decrease the percent_ID:

- low_cut-off=LOW_CUTOFF    low clustering threshold for phase II [default: 3.75]

- high_cut-off=HIGH_CUTOFF  high clustering threshold for phase III [default: 4.5]

- percent_id=PERCENT_ID     sequence similarity clustering threshold [default: 0.97]

The ``low_cut_off`` and the ``percent_id`` are used for clustering in the second, greedy clustering step.
The ``high_cut_off`` is used in the third clustering step, where unclustered reads are mapped according to their best match to any of the clusters of phase II. For good values for the thresholds, we refer to the plot S2 in the supplementary material of the denoiser paper (Reeder and Knight, Nature Methods 2010).



Q: Denoising on the clusters "hangs" after a while. What is going on?

A: If not provided with already preprocessed data via the -p option, the denoiser.py script automatically starts the preprocessing phase (cluster phase I in the paper) on one CPU on the cluster. This preprocessing takes from a few minutes for partial GS FLX runs to an hour or more for large Titanium runs. After this step, the parallel cluster phase II starts. First, all requested workers are started one-by-one. Depending on your queueing system and the number of jobs this might take from few seconds to several minutes. If one or more of the jobs are not started by the queueing system, all submitted jobs will block and wait. This is most likely the state your process is in if nothing seems to happen. We know this is not optimally and already thinking about a better solution for the future. In the meantime, make sure you only request as many jobs as you can safely run in your queue and monitor (qstat) the startup phase to see if all jobs are properly scheduled. If you finf that you requested to many CPUs and need to restart, simply kill the master process (denoiser.py) and it should bring down all but the last submitted jobs. The last job might need to be killed by hand.
Once all workers are succesfully started, you can monitor the progress by following the log file in verbose mode (toggled by the -v option)::
	
	tail -f denoiser.log


Q: How and why can I run the preprocessing step separately?

A: If you call denoiser.py without the -p option (or via its wrapper denoise_wrapper.py in QIIME) the preprocessing step (cluster phase I) is implicitly called. You can explicitly run the preprocessing step via the script preprocess.py and provide the output directory to denoiser.py using the -p option. Reasons for running the steps separately could be:

- run preprocess on a very fast single CPU machine, then transfer the data to a slower multi-cpu cluster

- You want to check the cluster statistics of phase I first, before deciding of wether the data needs to be split or how many CPUs

- something went wrong with the compute cluster in phase II and the program aborted. The results of preprocessing will be in the output dir and can be re-used if you restart the process.
 

Q: What about different next-gen sequencing platforms?

A: Denoising in this form only applies to 454 based pyrosequencing.
