.. _otupipe:

==========================================
Using usearch quality filtering with QIIME
==========================================

Introduction
-------------
**usearch_qf** (usearch quality filter) is a pipeline script built using `USEARCH <http://www.drive5.com/usearch>`_ to perform filtering of noisy sequences, chimera checking, and OTU picking on a set of de-multiplexed (i.e. post `split_libraries.py <../scripts/split_libraries.html>`_) sequences. This tutorial explains how to use usearch_qf through QIIME, with details about each of the steps performed and a brief description of basic parameters and their effect.

Usearch quality filtering is modeled after the OTUPipe series of scripts by Robert Edgar.  For details, please check the website `OTUPIPE <http://www.drive5.com/otupipe>`_ where you can also find some benchmark results using an artificial bacterial community `<http://www.drive5.com/usearch/perf/mock_results.html>`_.

.. _basicuse:

Basic usage
-----------
To use usearch_qf in QIIME, you will need a FASTA file resulting from split_libraries.py. In this tutorial we will use data from the main QIIME tutorial, so our input file will be :file:`seqs.fna`. From the directory where this file is located, type: ::

    pick_otus.py -i seqs.fna -m usearch --db_filepath=/path/to/gold.fa -o usearch_qf_results/ --word_length 64

where :file:`/path/to/gold.fa` specifies the full path to the location of the reference set that will be used when doing chimera checking. A copy of this file can be found `here <http://drive5.com/otupipe/gold.tz>`_ (remember to uncompress the file). After executing this command, several files will be created in the :file:`usearch_qf_results/` directory. The only file that you will need at this point is :file:`usearch_qf_results/seqs_otus.txt` (the OTU map file), which can then be used to pick a set of representative sequences with `pick_rep_set.py <../scripts/pick_rep_set.html>`_ as you would do after running :file:`pick_otus.py` with default options (i.e. using uclust).  The word length is the optimized value for the mock community results `<http://www.drive5.com/usearch/perf/mock_results.html>`_.

How does it work
----------------

usearch_qf performs 10 steps to process the input reads. We will assume log files are created at each step, which is the default setting in QIIME. File names correspond to those you will see if you use the command specified in the section `Basic usage`__. The value of certain parameters might be different depending on what you specified.

__ basicuse_

Step 1. Sort by length
^^^^^^^^^^^^^^^^^^^^^^
Sequences are initially sorted by length, and the result is stored in the file :file:`usearch_qf_results/len_sorted.fasta`.

Step 2. De-replication
^^^^^^^^^^^^^^^^^^^^^^
Sequences are de-replicated, so that the resulting file will contain unique sequences only. Results are stored in the file :file:`usearch_qf_results/dereplicated_seqs.fasta`. In this file each sequence description now includes information on how many sequences match exactly this one, so for instance a sequence with two exact copies would appear as:

.. note::

   * >PC.634_73 FLP3FBN01DX9GS orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0;size=2
   * CTGGTCCGTGTCTCAGTACCAGTGTGGGGGGCCTTCCTCTCAGAACCCCTACGCATCGTCGCCTCGGTGGGCCGTTACCC
   * CGCCGACTAGCTAATGCGCCGCATGCCCATCCGTGGCCGGGATTGCTCCCTTTGGCGGCCCGGGGATGCCCCAAGGCCGC
   * GTTACGCGGTATTAGACGGGGTTTCCCCCGCTTATCCCCCTGCCACGGGCAGGTTGCATACGTGTTACTCACCCGTGCGC
   * CGGTCGCCGGCGG

By default, de-replication is performed using --max_rejects=500, which can be time demanding if your input data set is large. Reducing this value to, for instance, 64, can greatly improve the speed of this step and still produce very similar results.

Step 3. Sort by abundance
^^^^^^^^^^^^^^^^^^^^^^^^^
De-replicated sequences are then sorted by abundance using the information generated in the previous step, the result being stored in the file :file:`usearch_qf_results/abundance_sorted.fasta`.

Step 4. Filtering of noisy sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sequences are clustered at the specified identity (by default, 97%) to filter noisy reads, and the resulting consensus sequences are written to :file:`usearch_qf_results/clustered_error_corrected.fasta`. Each sequence header contains a new identifier (a unique cluster number) and the size of the cluster:

.. note::

   * >Cluster0;size=50
   * TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCC
   * CGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCA
   * TCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGG
   * >Cluster1;size=52
   * CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACAC
   * CGCCAACTACCTAATGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGT
   * TATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCG
   * GTCGCCGG

The identity percentage specified for error correction can be set with the option -j or --percent_id_err, by default 0.97. Higher values of this parameter will result in less reads being merged together at this point; "good" reads that are similar to each other other won't be clustered as a unique read (i.e. you are not artificially reducing diversity), but some "noisy" reads will escape detection, thus artificially inflating diversity estimates. In general we have not found cases where this parameter needs to be modified. Additionally, running time can be affected by larger values of the parameter --max_rejects in this step.

Step 5. Chimera filter, de novo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once the sequences have been corrected for errors, chimera checking is performed using **UCHIME** (1_). In this step "*de novo*" checking is performed, without using any external set of reference sequences. This is particularly useful when are using data for which a good reference set does not exist. However, "*de novo*" chimera checking can be computationally expensive for large datasets, and we suggest to disable it in such cases using the parameter -k or --de_novo_chimera_detection. Results from this step are stored in files :file:`de_novo_non_chimeras.fasta` and :file:`de_novo_chimeras.fasta`.

The parameter -a or --abundance_skew can be used to control the abundance skew for chimera detection.

Step 6. Chimera filter, ref-based
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Reference-based chimera checking is performed against :file:`gold.fa` (or another user-provide set of sequences), and results are stored in files :file:`reference_non_chimeras.fasta` and :file:`reference_novo_chimeras.fasta`.

The parameter -f or --db_filepath can be used to specify the path to the sequence set to be used for ref-based chimera checking. To skip this step altogether, use the option -x or --reference_chimera_detection. 

Step 7. Merging chimera-checked sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sequences tagged as non-chimeric during steps 6 and 7 can be combined either by taking the intersection (only sequences flagged as non-chimeric by both methods) or union (sequences recognized by one of the methods as non-chimeric). Results are stored in :file:`combined_non_chimeras.fasta`.

The parameter -F or --non_chimeras_retention is used to set the merging as the union or intersection of the sets of non-chimeric sequences obtained from "*de novo*" and reference-based chimera checking.

**Example**
Assume there are 4 sequences (A, B, C, D) before chimera checking and "*de novo*" tags sequence A and B as chimeric while ref-based tags sequences B and C. Using --non_chimeras_retention=union will result in sequence B tagged as chimeric and A, C, and D as non-chimeric, while --non_chimeras_retention=intersection will tag A, B, and C as chimeras and only D as a non chimera.

Step 8. Sort by abundance chimera-free sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once sequences tagged as chimeras have been removed, the sequences are again sorted by abundance and clusters with less than MINSIZE reads are discarded. Results are stored in :file:`abundance_sorted_minsize_4.fasta` (this assume MINSIZE is set to the default value of 4). To modify the minimum number of reads that a cluster can have, use the parameter -g or --minsize. A value of 2, for instance, would remove all singletons (clusters of size 1). To skip this step use the parameter -l or --cluster_size_filtering.

Step 9. Cluster chimera-free sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This step corresponds to what is usually known as "*OTU picking*", i.e. sequences are clustered at the desired identity level. Different to regular OTU picking, by using usearch_qf you have also performed error correction and chimera checking, producing a 'cleaner' set of OTUs that will contain less artifacts. Results are stored in :file:`clustered_seqs.fasta`.

The identity percentage to cluster reads can be specified with the parameter -s or --similarity. In general the default of 0.97 works well for most datasets. The parameter --max_rejects can be modified to reduce running time during this step.

Step 10. Assign sequential IDs to OTUs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The OTUs calculated in the previous step get their IDs replaced by a sequential number, and the result is stored in :file:`enumerated_otus.fasta`.

Step 11. Classify reads
^^^^^^^^^^^^^^^^^^^^^^^
Each non-chimeric reads is assigned to the specific OTU identifier it belongs to. This creates the OTU map file (:file:`seqs_otus.txt`), which can be later used by :file:`pick_rep_set.py`.

References
------------
.. [1] Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. UCHIME improves sensitivity and speed of chimera detection. Bioinformatics. 2011 Aug 15;27(16):2194-200. Epub 2011 Jun 23. (`link <http://www.ncbi.nlm.nih.gov/pubmed/21700674>`_)
