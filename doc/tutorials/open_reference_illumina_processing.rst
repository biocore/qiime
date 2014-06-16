.. _open_reference_illumina:

============================================================
Discussion of subsampled open reference OTU picking in QIIME
============================================================

This document describes QIIME's subsampled open-reference OTU picking protocol. These protocol is presented in Rideout et al. 2014, under review (`PeerJ PrePrint available here <https://peerj.com/preprints/411/>`_). The *subsampled open-reference OTU picking protocol* is optimized for large datasets, and yields nearly identical results to *legacy open-reference OTU picking*, so there there is no reason to ever use the *legacy* method anymore. (Note: In QIIME 1.6.0-dev and earlier, we referred to *legacy* open-reference OTU picking as *standard* open-reference OTU picking).

This document very briefly covers legacy open-reference OTU picking. Most of the text covers subsampled open-reference OTU picking, including a description of how to use it, what exactly is happening, and test results from applying it to some well-understood data sets.

 .. note:: In this document we make use of the Greengenes reference OTUs. You can always find a link to the latest version of the Greengenes reference OTUs on the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_.

 .. toctree::

---------------------------------------------------------------
Subsampled open-reference OTU picking
---------------------------------------------------------------

Subsampled reference OTU picking is suitable for any analysis that legacy open-reference OTU picking can be used for, but additionally scales to much larger data sets (such as multiple HiSeq runs, which may require several days on ~100 processors to analyze).

 .. warning:: If processing multiple HiSeq lanes, don't combine the sequence data into a single file. Instead, see :ref:`iterative-mode`.

This is an open-reference OTU picking protocol, meaning that sequences are clustered against a reference database, and reads which fail to hit the reference are subsequently clustered de novo. This differs from legacy open-reference OTU picking as it was optimized at several steps to enable running on massive numbers of sequences (e.g., hundreds of millions). The steps in this workflow are as follows.

Step 0: Optional prefilter (parallel)
-------------------------------------
Prefilter the input sequence collection by searching reads against the reference set with a low percent identity threshold to discard sequences that are not representatives of the targeted marker gene (e.g., sequences that are the result of sequencing error, or host genomic sequences). This filter is disabled by default as it is very slow to run, and the same sequences are often discarded later in the process as they fail to align with PyNAST. If you'd like to use this filter and are working with 16S sequence data, a threshold of 60% is a good place to start (passed as ``--prefilter_percent_id``). The choice of 60% is described :ref:`here <prefilter-threshold>`. All reads which fail to hit the reference set are discarded as likely sequencing error.

 .. warning:: If most or all of your sequences are being filtered at this step, your sequences may be in the reverse orientation with respect to your reference database. To address this, you should add the following to your parameters file (creating one, if necessary) and pass this file as ``-p`` to ``pick_open_reference_otus.py``: ``pick_otus:enable_rev_strand_match True``. This is included in the instructions below, but be aware that this doubles the memory used in this step of the workflow.

Step 1: Closed reference (parallel)
-----------------------------------
Apply closed-reference OTU picking against the reference collection. Generate a fasta file containing all reads that fail to hit the reference collection.

Step 2: De novo clustering of subsampled failures (serial)
----------------------------------------------------------
Randomly subsample the sequences that failed to hit the reference collection, and write these to a new fasta file (default subsampling percentage is 0.1%, modify with ``-s/--percent_subsample``). Cluster these reads de novo, and choose a representative set of sequences as the centroid of each OTU cluster. These are the *new reference* OTUs.

Step 3: Closed reference, round 2 (parallel)
--------------------------------------------
Pick closed-reference OTUs against the representative sequences from the previous step. Write all sequences that fail to hit the reference collection to a fasta file.

Step 4: De novo (serial)
------------------------
Pick de novo OTUs on all reads that failed to hit the reference collection in the previous step. These are the *clean-up* OTUs. This step can be suppress by passing ``--suppress_step4``.

Post-OTU processing (parallel and serial, depending on step)
------------------------------------------------------------

#. Assemble the reference OTUs, the new reference OTUs, and the clean-up OTUs into a new OTU map, and construct an OTU table. At this stage, all OTUs with a sequence count of smaller than 2 (i.e., the singleton OTUs) are discarded. This can be modified with the ``--min_otu_size`` option, and disabled by passing ``--min_otu_size=1``.

#. Construct a new reference collection based on this OTU picking run. This new reference collection will be the combination of the full input reference collection, the new reference OTU representative sequences, and the clean-up OTU representative sequences. Note that this will not include representatives of the singleton OTUs by default. Also note that this differs from the representative set of sequences for this run in that it contains *all* of the input reference sequences, not only the ones that are represented in the current data set (which is what the representative sequences for this run contains).

#. Taxonomy will be assigned to all OTUs (using the uclust consensus taxonomy assigner by default) and representative sequences will be aligned and a tree will be constructed. Finally, an additional OTU table will be constructed that excludes reads that failed to align with PyNAST. It is recommended that this OTU table be used in downstream analysis.

To apply this analysis to ``seqs1.fna``, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command.

You should *always use full paths* which are represented here by ``$PWD``, but will usually look something like ``$HOME/my_data/`` (in other words, they should start with a ``/``). In this example your input sequences (``seqs1.fna``), and your metadata mapping file (``map.txt``) are all in the same directory represented by ``$PWD``. If you work from the directory containing those files, you can leave ``$PWD`` in the commands instead of specifying the full paths::

	pick_open_reference_otus.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss/ -aO 8 -p $PWD/ucrss_params.txt

This command should be run in parallel. Each job will need approximately 4GB of RAM, so if running on EC2 and you want to start 8 parallel jobs (recommended setting for EC2), your instance type should be ``m2.4xlarge``. The ``-aO 8`` specifies that we want to start 8 parallel jobs - adjust this according to the resources you have available.

.. _ucrss-core-diversity-analyses:

When this job completes, you're almost ready to begin running diversity analyses. The last thing you'll need to do before you begin ischoose an even sampling depth for alpha and beta diversity analyses. To do this, you should review the number of reads per sample::

	biom summarize-table -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures_summary.txt

This will create an output file with information on the number of reads per sample. Choose a depth of sampling that maximizes the number of sequences you'll include, and also the number of samples that have at least that many sequences: samples with fewer sequences will be excluded from your beta diversity/PCoA analysis. **Even sampling is absolutely critical to getting meaningful UniFrac distances between your samples.**

After choosing an even sampling depth you can use the ``core_diversity_analyses.py`` script to rarify your OTU table, run alpha and beta diversity analyses, and generate taxonomic summaries of your samples (in this example we have chosen an even sampling depth of 25,000 sequences/sample)::

	core_diversity_analyses.py -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -e 25000 -o $PWD/ucrss/cd_even25000/ -t $PWD/ucrss/rep_set.tre -m $PWD/map.txt -aO8

Again the ``-aO8`` specifies that the job should be run in parallel on 8 processors. Adjust this according to your resources. When this completes you can access your results by opening ``$PWD/ucrss/cd_even25000/index.html`` in a web browser.

--------------------------------------------
 Subsampled OTU picking workflow evaluation
--------------------------------------------

Analyses were performed to confirm that results are highly correlated between the subsampled open-reference OTU picking workflow and the legacy open-reference OTU picking workflow. These results are presented in Rideout et al. 2014, under review (`PeerJ PrePrint available here <https://peerj.com/preprints/411/>`_).

.. _filter_to_closed_ref:

---------------------------------------------------------------
Filtering an open-reference OTU table to reference OTUs only
---------------------------------------------------------------

There are cases where you may be interested in working with the closed reference subset of your open reference OTU table (meaning only those OTUs that hit the reference collection, excluding the new OTUs). Following from the above commands, to do that you can filter the new OTUs from the OTU table with the following command::

	filter_otus_from_otu_table.py -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.reference_only.biom --negate_ids_to_exclude -e $PWD/refseqs.fna

What this does is filter exclude all OTUs with identifiers that are not present in ``$PWD/refseqs.fna``, so all of the new OTUs.

.. _iterative-mode:

----------------------------------------------------------------------------
 Using the subsampled open-reference OTU picking workflow in iterative mode
----------------------------------------------------------------------------

The subsampled open-reference OTU picking workflow can be run in iterative mode to support multiple different sequence collections, such as several HiSeq runs. In iterative mode, the list of sequence files will be processed in order, and the new reference sequences generated at each step will be used as the reference collection for the subsequent step. After all input collections have been processed a single OTU table and tree, covering all of the input collections, will be generated.

To apply this analysis to ``seqs1.fna`` and ``seqs2.fna`` in iterative mode, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command.


To apply this analysis to ``seqs1.fna``, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command.

You should *always use full paths* which are represented here by ``$PWD``, but will usually look something like ``$HOME/my_data/`` (in other words, they should start with a ``/``). In this example your input sequences (``seqs1.fna``), and your metadata mapping file (``map.txt``) are all in the same directory represented by ``$PWD``. If you work from the directory containing those files, you can leave ``$PWD`` in the commands instead of specifying the full paths::

	pick_open_reference_otus.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter/ -aO 8 -p $PWD/ucrss_params.txt

This command should be run in parallel. Each job will need approximately 4GB of RAM, so if running on EC2 and you want to start 8 parallel jobs (recommended setting for EC2), your instance type should be ``m2.4xlarge``. The ``-aO 8`` specifies that we want to start 8 parallel jobs - adjust this according to the resources you have available.

After iterative OTU picking you can continue on with diversity analyses as described :ref:`here <ucrss-core-diversity-analyses>`.
