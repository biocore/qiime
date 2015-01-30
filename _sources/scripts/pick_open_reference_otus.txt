.. _pick_open_reference_otus:

.. index:: pick_open_reference_otus.py

*pick_open_reference_otus.py* -- Perform open-reference OTU picking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script is broken down into 4 possible OTU picking steps, and 2 steps
involving the creation of OTU tables and trees. The commands for each step are
described below, including what the input and resulting output files are.
Additionally, the optional specified parameters of this script that can be passed
are referenced.

Step 1) Prefiltering and picking closed reference OTUs
The first step is an optional prefiltering of the input fasta file to remove
sequences that do not hit the reference database with a given sequence
identity (PREFILTER_PERCENT_ID). This step can take a very long time, so is
disabled by default. The prefilter parameters can be changed with the options:
--prefilter_refseqs_fp
--prefilter_percent_id
This filtering is accomplished by picking closed reference OTUs at the specified
prefilter percent id to produce:
prefilter_otus/seqs_otus.log
prefilter_otus/seqs_otus.txt
prefilter_otus/seqs_failures.txt
prefilter_otus/seqs_clusters.uc
Next, the seqs_failures.txt file is used to remove these failed sequences from
the original input fasta file to produce:
prefilter_otus/prefiltered_seqs.fna
This prefiltered_seqs.fna file is then considered to contain the reads
of the marker gene of interest, rather than spurious reads such as host
genomic sequence or sequencing artifacts.

If prefiltering is applied, this step progresses with the prefiltered_seqs.fna.
Otherwise it progresses with the input file. The Step 1 closed reference OTU
picking is done against the supplied reference database. This command produces:
step1_otus/_clusters.uc
step1_otus/_failures.txt
step1_otus/_otus.log
step1_otus/_otus.txt

The representative sequence for each of the Step 1 picked OTUs are selected to
produce:
step1_otus/step1_rep_set.fna

Next, the sequences that failed to hit the reference database in Step 1 are
filtered from the Step 1 input fasta file to produce:
step1_otus/failures.fasta

Then the failures.fasta file is randomly subsampled to PERCENT_SUBSAMPLE of
the sequences to produce:
step1_otus/subsampled_failures.fna.
Modifying PERCENT_SUBSAMPLE can have a big effect on run time for this workflow,
but will not alter the final OTUs.

Step 2) The subsampled_failures.fna are next clustered de novo, and each cluster
centroid is then chosen as a "new reference sequence" for use as the reference
database in Step 3, to produce:
step2_otus/subsampled_seqs_clusters.uc
step2_otus/subsampled_seqs_otus.log
step2_otus/subsampled_seqs_otus.txt
step2_otus/step2_rep_set.fna

Step 3) Pick Closed Reference OTUs against Step 2 de novo OTUs
Closed reference OTU picking is performed using the failures.fasta file created
in Step 1 against the 'reference' de novo database created in Step 2 to produce:
step3_otus/failures_seqs_clusters.uc
step3_otus/failures_seqs_failures.txt
step3_otus/failures_seqs_otus.log
step3_otus/failures_seqs_otus.txt

Assuming the user has NOT passed the --suppress_step4 flag:
The sequences which failed to hit the reference database in Step 3 are removed
from the Step 3 input fasta file to produce:
step3_otus/failures_failures.fasta

Step 4) Additional de novo OTU picking
It is assumed by this point that the majority of sequences have been assigned
to an OTU, and thus the sequence count of failures_failures.fasta is small
enough that de novo OTU picking is computationally feasible. However, depending
on the sequences being used, it might be that the failures_failures.fasta file
is still prohibitively large for de novo clustering, and the jobs might take
too long to finish. In this case it is likely that the user would want to pass
the --suppress_step4 flag to avoid this additional de novo step.

A final round of de novo OTU picking is done on the failures_failures.fasta file
to produce:
step4_otus/failures_failures_cluster.uc
step4_otus/failures_failures_otus.log
step4_otus/failures_failures_otus.txt

A representative sequence for each cluster is chosen to produce:
step4_otus/step4_rep_set.fna

Step 5) Produce the final OTU map and rep set
If Step 4 is completed, the OTU maps from Step 1, Step 3, and Step 4 are
concatenated to produce:
final_otu_map.txt

If Step 4 was not completed, the OTU maps from Steps 1 and Step 3 are
concatenated together to produce:
final_otu_map.txt

Next, the minimum specified OTU size required to keep an OTU is specified with
the --min_otu_size flag. For example, if the user left the --min_otu_size as the
default value of 2, requiring each OTU to contain at least 2 sequences, the any
OTUs which failed to meet this criteria would be removed from the
final_otu_map.txt to produce:
final_otu_map_mc2.txt

If --min_otu_size 10 was passed, it would produce:
final_otu_map_mc10.txt

The final_otu_map_mc2.txt is used to build the final representative set:
rep_set.fna

Step 6) Making the OTU tables and trees
An OTU table is built using the final_otu_map_mc2.txt file to produce:
otu_table_mc2.biom

As long as the --suppress_taxonomy_assignment flag is NOT passed,
then taxonomy will be assigned to each of the representative sequences
in the final rep_set produced in Step 5, producing:
rep_set_tax_assignments.log
rep_set_tax_assignments.txt
This taxonomic metadata is then added to the otu_table_mc2.biom to produce:
otu_table_mc_w_tax.biom

As long as the --suppress_align_and_tree is NOT passed, then the rep_set.fna
file will be used to align the sequences and build the phylogenetic tree,
which includes the de novo OTUs. Any sequences that fail to align are
omitted from the OTU table and tree to produce:
otu_table_mc_no_pynast_failures.biom
rep_set.tre

If both --suppress_taxonomy_assignment and --suppress_align_and_tree are
NOT passed, the script will produce:
otu_table_mc_w_tax_no_pynast_failures.biom

It is important to remember that with a large workflow script like this that
the user can jump into intermediate steps. For example, imagine that for some
reason the script was interrupted on Step 2, and the user did not want to go
through the process of re-picking OTUs as was done in Step 1. They can simply
rerun the script and pass in the:
--step_1_otu_map_fp
--step1_failures_fasta_fp
parameters, and the script will continue with Steps 2 - 4.

**Note:** If most or all of your sequences are failing to hit the reference
during the prefiltering or closed-reference OTU picking steps, your sequences
may be in the reverse orientation with respect to your reference database. To
address this, you should add the following line to your parameters file
(creating one, if necessary) and pass this file as -p:

pick_otus:enable_rev_strand_match True

Be aware that this doubles the amount of memory used in these steps of the
workflow.




**Usage:** :file:`pick_open_reference_otus.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The input sequences filepath or comma-separated list of filepaths
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-m, `-`-otu_picking_method
		The OTU picking method to use for reference and de novo steps. Passing usearch61, for example, means that usearch61 will be used for the de novo steps and usearch61_ref will be used for reference steps. [default: uclust]
	-r, `-`-reference_fp
		The reference sequences [default: /Users/jairideout/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta]
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	`-`-prefilter_refseqs_fp
		The reference sequences to use for the prefilter, if different from the reference sequecnces to use for the OTU picking [default: same as passed for --reference_fp]
	-n, `-`-new_ref_set_id
		Unique identifier for OTUs that get created in this ref set (this is useful to support combining of reference sets) [default:New]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]
	-s, `-`-percent_subsample
		Percent of failure sequences to include in the subsample to cluster de novo, expressed as a fraction between 0 and 1 (larger numbers should give more comprehensive results but will be slower) [default:0.001]
	`-`-prefilter_percent_id
		Sequences are pre-clustered at this percent id (expressed as a fraction between 0 and 1) against the reference and any reads which fail to hit are discarded (a quality filter); pass 0.0 to disable [default:0.0]
	`-`-step1_otu_map_fp
		Reference OTU picking OTU map, to avoid rebuilding if one has already been built. This must be an OTU map generated by this workflow, not (for example) `pick_closed_reference_otus.py <./pick_closed_reference_otus.html>`_.
	`-`-step1_failures_fasta_fp
		Reference OTU picking failures fasta filepath, to avoid rebuilding if one has already been built. This must be a failures file generated by this workflow, not (for example) `pick_closed_reference_otus.py <./pick_closed_reference_otus.html>`_.
	`-`-minimum_failure_threshold
		The minimum number of sequences that must fail to hit the reference for subsampling to be performed. If fewer than this number of sequences fail to hit the reference, the de novo clustering step will run serially rather than invoking the subsampled open reference approach to improve performance. [default: 100000]
	`-`-suppress_step4
		Suppress the final de novo OTU picking step  (may be necessary for extremely large data sets) [default: False]
	`-`-min_otu_size
		The minimum otu size (in number of sequences) to retain the otu [default: 2]
	`-`-suppress_taxonomy_assignment
		Skip the taxonomy assignment step, resulting in an OTU table without taxonomy [default: False]
	`-`-suppress_align_and_tree
		Skip the sequence alignment and tree-building steps [default: False]


**Output:**




Run the subsampled open-reference OTU picking workflow on seqs1.fna using refseqs.fna as the reference collection and using sortmerna and sumaclust as the OTU picking methods. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will genenerally look like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss_sortmerna_sumaclust/ -p $PWD/ucrss_smr_suma_params.txt -m sortmerna_sumaclust

Run the subsampled open-reference OTU picking workflow on seqs1.fna using refseqs.fna as the reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss/ -s 0.1 -p $PWD/ucrss_params.txt

Run the subsampled open-reference OTU picking workflow on seqs1.fna using refseqs.fna as the reference collection and using usearch61 and usearch61_ref as the OTU picking methods. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss_usearch/ -s 0.1 -p $PWD/ucrss_params.txt -m usearch61

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter/ -s 0.1 -p $PWD/ucrss_params.txt

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection. This is useful if you're working with marker genes that do not result in useful alignment (e.g., fungal ITS). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter_no_tree/ -s 0.1 -p $PWD/ucrss_params.txt --suppress_align_and_tree

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection, suppressing assignment of taxonomy. This is useful if you're working with a reference collection without associated taxonomy. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_open_reference_otus.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter_no_tax/ -s 0.1 -p $PWD/ucrss_params.txt --suppress_taxonomy_assignment


