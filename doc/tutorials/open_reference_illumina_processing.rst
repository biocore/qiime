.. _open_reference_illumina:

============================================================
Discussion of subsampled open reference OTU picking in QIIME
============================================================

This document describes QIIME's subsampled open-reference OTU picking protocol. As of May, 2013 a paper on this workflow is in preparation. The *subsampled open-reference OTU picking protocol* is optimized for large datasets, and yields identical results to *legacy open-reference OTU picking*, so there there is no reason to ever use the *legacy* method anymore. (Note: In QIIME 1.6.0-dev and earlier, we referred to *legacy* open-reference OTU picking as *standard* open-reference OTU picking).

This document very briefly covers legacy open-reference OTU picking. Most of the text covers subsampled open-reference OTU picking, including a description of how to use it, what exactly is happening, and test results from applying it to some well-understood data sets. 

 .. note:: In this document we make use of the Greengenes reference OTUs. You can always find a link to the latest version of the Greengenes reference OTUs on the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_.

 .. toctree::

---------------------------------------------------------------
 Option 1: Legacy open-reference OTU picking
---------------------------------------------------------------

Legacy open-reference OTU picking is suitable for a single HiSeq2000 lane (unless it's very high diversity, in which case runtime may be a limiting factor). You'll use the ``pick_de_novo_otus.py`` workflow script in QIIME with a custom parameters file. Note that the name of the script is misleading here. While this is called ``pick_de_novo_otus.py``, we'll be overriding parameters to have it run open-reference OTU picking. Because we no longer recommend legacy open-reference OTU picking (you should instead use subsampled open reference OTU picking via ``pick_open_reference_otus.py``), we won't modify the name of this script. 

Your parameters file should look like the following::

	pick_otus:otu_picking_method uclust_ref
	pick_otus:refseqs_fp <PATH TO REFERENCE COLLECTION>
	pick_otus:enable_rev_strand_match True

Where ``<PATH TO REFERENCE COLLECTION>`` is replaced with the path to the reference data set you'd like to pick OTUs against. The QIIME development group frequently uses the Greengenes reference OTUs. This might be something like ``$HOME/qiime_software/gg_12_10_otus-release/rep_set/97_otus.fasta``. 

This command should be run in parallel. Each job will need approximately 4GB of RAM, so if running on EC2 and you want to start 8 parallel jobs (recommended setting for EC2), your instance type should be ``m2.4xlarge``.

You can then use the following commands. You should *always use full paths* which are represented here by ``$PWD``, but will usually look something like ``$HOME/my_data/`` (in other words, they should start with a ``/``). In this example your input sequences (``seqs.fna``), your parameters file (``ucr_params.txt``), and your metadata mapping file (``map.txt``) are all in the same directory represented by ``$PWD``. If you work from the directory containing those files, you can leave ``$PWD`` in the commands instead of specifying the full paths.

First, pick otus, choose representative sequences, assign taxonomy to OTUs, and build a phylogenetic tree. The ``-aO 8`` specifies that we want to start 8 parallel jobs - adjust this according to the resources you have available. This is open-reference OTU picking, so reads will be clustered against the reference database (in parallel) and reads which fail to hit the reference data set will subsequently be clustered de novo (serially)::
	
	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/ucr/ -p $PWD/ucr_params.txt -aO 8

When working with Illumina data you typically want to filter singleton OTUs (i.e., OTUs with only one sequence) as these are likely to represent sequencing or PCR errors. In QIIME 1.4.0 (and most earlier versions) you can do that with this command::
	
	filter_otu_table.py -i $PWD/ucr/otu_table.txt -o $PWD/ucr/otu_table_mc2.txt -c 2 -s 1

In QIIME 1.4.0-dev and later, you can filter singleton OTUs with this command::
	
	filter_otus_from_otu_table.py -i $PWD/ucr/otu_table.biom -o $PWD/ucr/otu_table_mc2.biom -n 2

You'll notice that depending on your version of QIIME, the extension on your OTU table file will differ. In QIIME 1.4.0 and earlier, it will be ``.txt``. In QIIME 1.4.0-dev and later it will be ``.biom``. We'll continue this example assuming that your OTU table ends with ``.txt`` (if you're working with QIIME 1.4.0-dev or later, you likely decided to go with option 2 for OTU picking).

As PCoA of UniFrac distances between samples is a frequent result of interest in microbial ecology, we'll cover how to generate PCoA plots next. The first thing you'll want to do is evenly sample your OTU table. To choose an even sampling depth, review the number of reads per sample::
	
	biom summarize-table -i $PWD/ucr97/otu_table_mc2.biom -o $PWD/ucr97/otu_table_mc2_summary.txt

This will create an output file with information on the number of reads per sample. Choose a depth of sampling that maximizes the number of sequences you'll include, and also the number of samples that have at least that many sequences: samples with fewer sequences will be excluded from your beta diversity/PCoA analysis. **Even sampling is absolutely critical to getting meaningful UniFrac distances between your samples.**

After choosing an even sampling depth you can use the ``beta_diversity_through_plots.py`` script to rarify your OTU table, compute UniFrac distances between samples, and run Principal Coordinates Analysis with the following command (in this example we have chosen an even sampling depth of 25,000 sequences/sample)::
	
	beta_diversity_through_plots.py -i $PWD/ucr/otu_table_mc2.biom -e 25000 -o $PWD/ucr/bdiv_even25000/ -t $PWD/ucr/rep_set.tre -m $PWD/map.txt -aO8

Again the ``-aO8`` specifies that the job should be run in parallel on 8 processors. Adjust this according to your resources. When this completes you can open 3D interactive PCoA plots by opening the file ``bdiv_even25000/unweighted_unifrac_3d_continuous/unweighted_unifrac_pc_3D_PCoA_plots.html`` in a web browser. This may take a minute to load for very large sets of samples.


---------------------------------------------------------------
 Option 2: Subsampled open-reference OTU picking
---------------------------------------------------------------

Subsampled reference OTU picking is suitable for any analysis that legacy open-reference OTU picking can be used for, but additionally scales to much larger data sets (such as multiple HiSeq runs, which may require several days on ~100 processors to analyze).

 .. warning:: If processing multiple HiSeq lanes, don't combine the sequence data into a single file. Instead, see :ref:`iterative-mode`.

This is an open-reference OTU picking protocol, meaning that sequences are clustered against a reference database, and reads which fail to hit the reference are subsequently clustered de novo. This differs from legacy open-reference OTU picking as it was optimized at several steps to enable running on massive numbers of sequences (hundreds of millions, which is massive as of this writing). The steps in this workflow are as follows.

Step 0: Prefilter (parallel)
----------------------------
Prefilter the input sequence collection by searching reads against the reference set with a low percent identity threshold (default is 60%, modify with ``--prefilter_percent_id``). The choice of 60% is described :ref:`here <prefilter-threshold>`. All reads which fail to hit the reference set are discarded as likely sequencing error.

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

To apply this analysis to ``seqs1.fna``, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command. Your parameters file should contain the following::

	pick_otus:otu_picking_method uclust_ref
	pick_otus:enable_rev_strand_match True

You should *always use full paths* which are represented here by ``$PWD``, but will usually look something like ``$HOME/my_data/`` (in other words, they should start with a ``/``). In this example your input sequences (``seqs1.fna``), and your metadata mapping file (``map.txt``) are all in the same directory represented by ``$PWD``. If you work from the directory containing those files, you can leave ``$PWD`` in the commands instead of specifying the full paths::

	pick_open_reference_otus.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss/ -aO 8 -p $PWD/ucrss_params.txt

This command should be run in parallel. Each job will need approximately 4GB of RAM, so if running on EC2 and you want to start 8 parallel jobs (recommended setting for EC2), your instance type should be ``m2.4xlarge``. The ``-aO 8`` specifies that we want to start 8 parallel jobs - adjust this according to the resources you have available.

.. _ucrss-beta-diversity:

As PCoA of UniFrac distances between samples is a frequent result of interest in microbial ecology, we'll cover how to generate PCoA plots next. The first thing you'll want to do is evenly sample your OTU table. To choose an even sampling depth, review the number of reads per sample::
	
	biom summarize-table -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures_summary.txt

This will create an output file with information on the number of reads per sample. Choose a depth of sampling that maximizes the number of sequences you'll include, and also the number of samples that have at least that many sequences: samples with fewer sequences will be excluded from your beta diversity/PCoA analysis. **Even sampling is absolutely critical to getting meaningful UniFrac distances between your samples.**

After choosing an even sampling depth you can use the ``beta_diversity_through_plots.py`` script to rarify your OTU table, compute UniFrac distances between samples, and run Principal Coordinates Analysis with the following command (in this example we have chosen an even sampling depth of 25,000 sequences/sample)::
	
	beta_diversity_through_plots.py -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom  -e 25000 -o $PWD/ucrss/bdiv_even25000/ -t $PWD/ucrss/rep_set.tre -m $PWD/map.txt -aO8

Again the ``-aO8`` specifies that the job should be run in parallel on 8 processors. Adjust this according to your resources. When this completes you can open 3D interactive PCoA plots by opening the file ``bdiv_even25000/unweighted_unifrac_3d_continuous/unweighted_unifrac_pc_3D_PCoA_plots.html`` in a web browser. This may take a minute to load for very large sets of samples.

.. _filter_to_closed_ref:

---------------------------------------------------------------
Filtering an open-reference OTU table to reference OTUs only
---------------------------------------------------------------

There are cases where you may be interested in working with the closed reference subset of your open reference OTU table (meaning only those OTUs that hit the reference collection, excluding the new OTUs). Following from the above commands, to do that you can filter the new OTUs from the OTU table with the following command::

	filter_otus_from_otu_table.py -i $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.biom -o $PWD/ucrss/otu_table_mc2_w_tax_no_pynast_failures.reference_only.biom --negate_ids_to_exclude -e $PWD/refseqs.fna

What this does is filter exclude all OTUs with identifiers that are not present in ``$PWD/refseqs.fna``, so all of the new OTUs.

--------------------------------------------
 Subsampled OTU picking workflow evaluation
--------------------------------------------

Several analyses were performed to confirm that results are comparable between the subsampled open-reference OTU picking workflow and the legacy open-reference OTU picking workflow. These include analyses on two different data sets: one host-associated (the `Costello Whole Body <http://www.ncbi.nlm.nih.gov/pubmed/19892944>`_ study) and one free-living (the `Lauber 88 soils <http://www.ncbi.nlm.nih.gov/pubmed/19502440>`_ study). These two were chosen as Greengenes (the reference set being used) is known to be biased toward human-associated microbes, so I wanted to confirm that the method works when few sequences fail to hit the reference set (whole body) and when many sequences fail to hit the reference set (88 soils).

Several tests were performed:
 * beta diversity (procrustes analysis to compare subsampled OTU results to de novo, open-reference, and closed-reference OTU picking)
 * alpha diversity (test for correlation in observed OTU count between subsampled OTU results and de novo, open-reference, and closed-reference OTU picking)
 * otu category significance (reviewed significant OTUs - need a good way to quantitate this)

For all analyses, sequences that fail to align with PyNAST and singleton OTUs were discarded (these are defaults in the subsampled OTU picking workflow).

------------------
88 soils analysis
------------------
This analysis is based on the data presented in the `Lauber 88 soils <http://www.ncbi.nlm.nih.gov/pubmed/19502440>`_ study.


Alpha diversity
---------------
Here I checked whether the subsampled reference OTU alpha diversities for all samples were correlated with the de novo OTU picking, legacy open-reference OTU picking, and closed-reference OTU picking alpha diversities. The *observed species/OTUs* metric was calculated on add data sets (``alpha_diversity.py -m observed_species``), and the Pearson correlations were computed for subsampled reference OTU picking with the three other sets of values.

Results
```````
 * subsampled open-reference OTU picking versus de novo OTU picking: r=0.995 p=4.836e-88
 * subsampled open-reference OTU picking versus legacy open-reference OTU picking: r=1.000 p=0.000
 * subsampled open-reference OTU picking versus closed-reference OTU picking: r=0.8634 p=1.405e-27

Conclusions
```````````
Subsampled open-reference OTU picking alpha diversity values are significantly correlated with de novo, legacy open-reference, and closed-reference OTU picking results. This suggests that we will derive the same biological conclusions between regarding alpha diversity when using the subsampled OTU picking workflow.

Beta diversity
--------------
Here I checked whether Procrustes comparisons of unweighted UniFrac PCoA plots between subsampled open-reference OTU picking and de novo OTU picking, legacy open-reference OTU picking, and closed-reference OTU picking yield significant results. This was calculated using ``transform_coordinate_matrices.py`` which is described in the `Procrustes tutorial <./procrustes_analysis.html>`_. p-values are based on 1000 Monte Carlo iterations.

Results
```````
 * subsampled open-reference OTU picking versus de novo OTU picking: M2=0.009 p<0.001
 * subsampled open-reference OTU picking versus legacy open-reference OTU picking: M2=0.007 p<0.001
 * subsampled open-reference OTU picking versus closed-reference OTU picking: M2=0.039 p<0.001

Conclusions
```````````
Procrustes results are highly significant for the three comparisons, suggesting that we will derive the same biological conclusions regardless of which of these OTU picking workflows is used.


OTU category significance
-------------------------
Here I confirm that the same taxonomy groups are identified as significantly different across the pH gradient in these soils using ANOVA, regardless of which OTU picking workflow is applied. These results were computed with the ``otu_category_significance.py`` script. To define a category for this test I binned the pH values by truncating the values to integers (so 5.0, 5.3, and 5.9 are all binned to pH 5) and using this binned pH as the category. Since I'm just looking for consistent results across the different OTU picking methods I don't think it's important that this isn't the most biologically relevant binning strategy. **Note that OTU ids are not directly comparable across all analyses, so it is best to compare the taxonomies.**

Results
```````


Top 5 OTUs that differ across bins for subsampled open-reference OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
New.CleanUp.ReferenceOTU26927 1.99e-11                      k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Chromatiales; f__Sinobacteraceae
New.CleanUp.ReferenceOTU34053 7.06e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
212596                        9.26e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
112859                        1.22e-08                      k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae
New.CleanUp.ReferenceOTU36189 4.35e-08                      k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Rubrobacterales; f__Rubrobacteraceae
============================= ============================= ==============================================================================================

Top 5 OTUs that differ across bins for de novo OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
26819                         3.19e-11                      k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Chromatiales; f__Sinobacteraceae
28062                         5.92e-10                      k__Bacteria; p__Acidobacteria; c\_\_; o\_\_; f__Koribacteraceae
35264                         2.43e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
45059                         5.48e-09                      k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o\_\_; f\_\_
7687                          2.056e-08	                    k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
============================= ============================= ==============================================================================================


Top 5 OTUs that differ across bins for legacy open-reference OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
DeNovoOTU26928                1.99e-11                      k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Chromatiales; f__Sinobacteraceae
DeNovoOTU34054                7.06e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
212596                        9.26e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae
112859                        1.22e-08                      k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae
DeNovoOTU36190                4.35e-08                      k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Rubrobacterales; f__Rubrobacteraceae
============================= ============================= ==============================================================================================

Top 5 OTUs that differ across bins for closed-reference OTU picking:

============================= ============================= ===================================================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ===================================================================================================================
212596                        4.03e-09                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae; g__CandidatusSolibacter; s\_\_
112859                        4.62e-09                      k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f\_\_; g\_\_; s\_\_
544749                        5.56e-08                      k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Chromatiales; f__Sinobacteraceae; g\_\_; s\_\_
541300                        1.28e-07                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae; g__CandidatusSolibacter; s\_\_
563862                        1.95e-07                      k__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae; g__CandidatusSolibacter; s\_\_
============================= ============================= ===================================================================================================================

Conclusions
```````````
In lieu of a solid statistical approach to compare these results, the results appear consistent across the different OTU picking workflows. The legacy open-reference and subsampled open-reference are remarkably consistent. 

Additional sanity check: is the new reference dataset sane?
-----------------------------------------------------------
To confirm that the new reference data set works as expected, I applied legacy open-reference OTU picking on the original input sequences against the new reference collection generated by the subsampled OTU analysis. The idea here is that most reads should now hit the reference collection. A number of reads still fail, but on close investigation these turn out to all cluster into singleton OTUs. So, this is expected as singletons are not included in the reference collection (possible to adjust this with the ``--min_otu_size`` parameter [default = 2]). The new reference collection that is generated does appear to be sane. The command used for this analysis was::
	
	pick_de_novo_otus.py -i $HOME/data/lauber_88soils/seqs.fna -o $HOME/data/lauber_88soils/subsample_ref_otus_eval/ucr97_v_new_ref/ -p $HOME/data/lauber_88soils/subsample_ref_otus_eval/ucr_v_newref_params.txt -aO 3

The parameters file (``-p``) for this analysis contained the following lines::

	pick_otus:otu_picking_method uclust_ref
	pick_otus:refseqs_fp $HOME/data/lauber_88soils/subsample_ref_otus_eval/prefilter60/new_refseqs.fna
	pick_otus:enable_rev_strand_match True


--------------------
Whole body analysis
--------------------
This analysis is based on the data presented in the `Costello Whole Body <http://www.ncbi.nlm.nih.gov/pubmed/19892944>`_ study.

Alpha diversity
---------------
Here I checked whether the subsampled reference OTU alpha diversities for all samples were correlated with the de novo OTU picking, legacy open-reference OTU picking, and closed-reference OTU picking alpha diversities. The *observed species/OTUs* metric was calculated on add data sets (``alpha_diversity.py -m observed_species``), and the Pearson correlations were computed for subsampled reference OTU picking with the three other sets of values.

Results
```````
 * subsampled open-reference OTU picking versus de novo OTU picking: r=0.99  p=0.0
 * subsampled open-reference OTU picking versus legacy open-reference OTU picking: r=1.0 p=0.0
 * subsampled open-reference OTU picking versus closed-reference OTU picking: r=0.95 p=0.0

Conclusions
```````````
Subsampled open-reference OTU picking alpha diversity values are significantly correlated with de novo, legacy open-reference, and closed-reference OTU picking results. This suggests that we will derive the same biological conclusions between regarding alpha diversity when using the subsampled OTU picking workflow.

Beta diversity
--------------
Here I checked whether Procrustes comparisons of unweighted UniFrac PCoA plots between subsampled open-reference OTU picking and de novo OTU picking, legacy open-reference OTU picking, and closed-reference OTU picking yield significant results. This was calculated using ``transform_coordinate_matrices.py`` which is described in the `Procrustes tutorial <./procrustes_analysis.html>`_. p-values are based on 1000 Monte Carlo iterations.

Results
```````
 * subsampled open-reference OTU picking versus de novo OTU picking: M2=0.056 p<0.001
 * subsampled open-reference OTU picking versus legacy open-reference OTU picking: M2=0.053 p<0.001
 * subsampled open-reference OTU picking versus closed-reference OTU picking: M2=0.072 p<0.001

Conclusions
```````````
Procrustes results are highly significant for the three comparisons, suggesting that we will derive the same biological conclusions regardless of which of these OTU picking workflows is used.


OTU category significance
-------------------------
Here I confirm that the same taxonomy groups are identified as significantly different across the body sites in this study using ANOVA, regardless of which OTU picking workflow is applied. These results were computed with the ``otu_category_significance.py`` script. **Note that OTU ids are not directly comparable across all analyses, so it is best to compare the taxonomies.**

Results
```````


Top 5 OTUs that differ across bins for subsampled open-reference OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
470747                        5.37e-151                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae
471122                        4.35e-143                     k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae
470477                        4.14e-135                     k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Carnobacteriaceae
94166                         1.61e-125                     k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae
442949                        1.33e-124                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae
============================= ============================= ==============================================================================================

Top 5 OTUs that differ across bins for de novo OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
6389                          4.54e-148                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae
17234                         1.16e-141                     k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae
18227                         3.05e-139                     k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Carnobacteriaceae
7262                          1.22e-134                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae
18575                         2.74e-122                     k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae
============================= ============================= ==============================================================================================

Top 5 OTUs that differ across bins for legacy open-reference OTU picking:

============================= ============================= ==============================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ==============================================================================================
470747                        5.37e-151                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae
471122                        4.35e-143                     k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae
470477                        4.14e-135                     k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Carnobacteriaceae
94166                         1.61e-125                     k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae
442949                        1.33e-124                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae
============================= ============================= ==============================================================================================

Top 5 OTUs that differ across bins for closed-reference OTU picking:

============================= ============================= ===================================================================================================================
OTU ID                        Bonferroni-adjusted p-value   Taxonomy
============================= ============================= ===================================================================================================================
470747                        3.25e-150                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae
471122                        2.20e-148                     k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae
470477                        2.51e-138                     k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Carnobacteriaceae
94166                         4.30e-130                     k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae
442949                        2.04e-121                     k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae
============================= ============================= ===================================================================================================================

Conclusions
```````````
In lieu of a solid statistical approach to compare these results, the results are remarkably consistent across the different OTU picking workflows.

.. _prefilter-threshold:

Additional sanity check: what reads are being discarded by the prefilter?
-------------------------------------------------------------------------
To investigate what reads get discarded at the prefilter stage, I evaluated a subset of the reads discarded when the prefilter was set to 80% (``--prefilter_percent_id 0.80``) versus when the prefilter was set to 60% (default).

Sequences filtered at 80% but not at 60%
````````````````````````````````````````

These three have high percent id matches to 16S sequences in NCBI::

	>F12Pinl.140479_129272 FFLHOYS02GCJLO orig_bc=ATACAGAGCTCC new_bc=ATACAGAGCTCC bc_diffs=0
	CTGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCATCCGAAAAGACCAGCTAAGCATCATTGGCTTGGTCAGCCTTTACCTAACCAACTACCTGATACTACGTGGGCTCATCGAACAGCGCGAATTAGCTTGCTTTATGAATTATTCAGGATTTGGAGTGAACTATTCGGCAGATTCCCACGCGTTACGCACCCGTTCGCCACTTTGCTTG
	>F32Indr.140459_1174716 FFO92CG02IYZBA orig_bc=GCTATCACGAGT new_bc=GCTATCACGAGT bc_diffs=0
	CCGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCATCCGAAAAGACCAGCTAAGCATCATTGGCTTGGTCAGCCTTTACCTGACCAACTACCTAATACTACGCAGGCTCATCAAACAGCGCTTTTTAGCTTTCTTCAGGATTTGGCCCGAACTGTTCGGCAGATTCCCACGCGTTACGCACCCGTTCGCCACTTTGTTCTCAACTGTTCCCACCTCCTGGGCGAGA
	>F32Forr.140528_1210712 FFO92CG02IKGYS orig_bc=GCGTTACACACA new_bc=GCGTTACACACA bc_diffs=0
	CCGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCATCCGAAAAGACCAGCTAAGCATCATTGGCTTGGTCAGCCTTTACCTGACCAACTACCTAATACTACGCAGGCTCATCAAACAGCGCTTTTGAGCTTTCTTCAGGATTTGGCCCGAACTGTTCGGCAGATTCCCACGCGTTACGCACCCGTTCGCCACTTTGTTCTCAACTATTCCGATTCTTTTTTCGGTAGGCCGTTA

Sequences filtered 80% and at 60%
`````````````````````````````````
These three reads hit a small fragment, a human sequence, and nothing in NCBI, respectively::

	>M22Pinr.140692_1148864 FFO92CG01EQIWQ orig_bc=CGCACATGTTAT new_bc=CGCACATGTTAT bc_diffs=0
	GGAAAAGGGAAAAACAGATGAGACAAATAGAAAACAAATAGCAAATTAGTAGGTGTTAACATGACTTTATCAATAATTACATCAAATGTAGATGATGTTAACCATGGATTGACAAACTTTTTCTTTATAGGACCAGACAGTCAATATTTTAGGTCTTTGAGGCCATATGGTATCTGTCATAACCACTCAACTGAGCCAGGATCAAACTCTGA
	>F31Nstl.140789_1153834 FFO92CG02FSK33 orig_bc=GCAGCCGAGTAT new_bc=GCAGCCGAGTAT bc_diffs=0
	CANNOT INCLUDE THIS READ DUE TO IRB RESTRICTIONS
	>F32Nstl.140804_1160735 FFO92CG01BRQNZ orig_bc=GCTGCTGCAATA new_bc=GCTGCTGCAATA bc_diffs=0
	CTGAAACCCTGGGTCACCAAAAGGCAGGAGGAGGAGGGACAGGGCAAGGCAGGGGAAGAGAGGGGAGGCTGACTCACATACACACATATGCATGCACACATCACACCCACATTCATGTACACACACACAGATTCACATGCATGCACAGCACAATCGCACACTTGTATACACACACAGGCACA

Conclusions
```````````
Based on this analysis (and currently unpublished data -- will fill in when available), a threshold of 60% was chosen as the default value for discarding sequences that are likely not rRNA.

.. _iterative-mode:

----------------------------------------------------------------------------
 Using the subsampled open-reference OTU picking workflow in iterative mode
----------------------------------------------------------------------------

The subsampled open-reference OTU picking workflow can be run in iterative mode to support multiple different sequence collections, such as several HiSeq runs. In iterative mode, the list of sequence files will be processed in order, and the new reference sequences generated at each step will be used as the reference collection for the subsequent step. After all input collections have been processed a single OTU table and tree, covering all of the input collections, will be generated. 

To apply this analysis to ``seqs1.fna`` and ``seqs2.fna`` in iterative mode, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command. 


To apply this analysis to ``seqs1.fna``, picking OTUs against the reference collection ``refseqs.fna`` you can run the following command. Your parameters file should contain the following::

	pick_otus:otu_picking_method uclust_ref
	pick_otus:enable_rev_strand_match True

You should *always use full paths* which are represented here by ``$PWD``, but will usually look something like ``$HOME/my_data/`` (in other words, they should start with a ``/``). In this example your input sequences (``seqs1.fna``), and your metadata mapping file (``map.txt``) are all in the same directory represented by ``$PWD``. If you work from the directory containing those files, you can leave ``$PWD`` in the commands instead of specifying the full paths::

	pick_open_reference_otus.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter/ -aO 8 -p $PWD/ucrss_params.txt

This command should be run in parallel. Each job will need approximately 4GB of RAM, so if running on EC2 and you want to start 8 parallel jobs (recommended setting for EC2), your instance type should be ``m2.4xlarge``. The ``-aO 8`` specifies that we want to start 8 parallel jobs - adjust this according to the resources you have available. 

After iterative OTU picking you can continue on with beta diversity (and other) analyses as described :ref:`here <ucrss-beta-diversity>`.
