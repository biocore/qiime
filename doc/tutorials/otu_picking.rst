.. _otu_picking:

===============================
OTU picking strategies in QIIME
===============================

QIIME provides three high-level protocols for OTU picking. These can be described as de novo, closed-reference, and open-reference OTU picking, and are accessible through `pick_de_novo_otus.py <>`_, `pick_closed_reference_otus.py <>`_, and `pick_open_reference_otus.py <>`_. Each of these protocol are described in this document, and commands are provided which illustrate how to run each of these with uclust and usearch 6.1 (i.e, usearch61).

Description of OTU picking processes
====================================

TODO: Discussion of pros and cons of each

Running the OTU picking workflows
=================================

The same workflow commands are used for running OTU picking with usearch61 and uclust. To run the methods with usearch, that is either passed in a parameters file or on the command line, depending on the workflow. See :ref:`qiime_parameter_files` for information on parameter files.

Conventions used in these examples
----------------------------------

It's a good idea, particularly for when running these workflows in parallel, to specify absolute paths for your input and output files. That is indicated here with ``$PWD``, but in practice it will often looks something like ``/home/ubuntu/my-analysis/seqs.fna``.

The reference-based OTU picking workflows require that the user provide reference files. Here we define some environment variables to point to those locations. These paths will likely be different on your system. You can download QIIME-compatible reference files from the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_. In this example we're working with the Greengenes 12_10 reference OTU collection. You can set environment variables to point to these as follows::

	export reference_seqs=/home/ubuntu/qiime_software/gg_otus-12_10-release/rep_set/97_otus.fasta
	export reference_tree=/home/ubuntu/qiime_software/gg_otus-12_10-release/trees/97_otus.tree
	export reference_tax=/home/ubuntu/qiime_software/gg_otus-12_10-release/taxonomy/97_otu_taxonomy.txt

De novo OTU picking
-------------------

With uclust::

	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/dn_uc/

With usearch61::
	
	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/dn_us/ -p $PWD/usearch_params.txt

where the following information is in ``usearch_params.txt``::
	
	pick_otus:otu_picking_method usearch61

You can find an additional example using de novo OTU picking in :ref:`tutorial`.


Closed-reference OTU picking
----------------------------

With uclust::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/cr_uc/ -r $reference_seqs -t $reference_tax

With usearch61::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/cr_us/ -r $reference_seqs -t $reference_tax -p usearch_params.txt

where the following information is in ``usearch_params.txt``::
	
	pick_otus:otu_picking_method usearch61

You can find an additional example using closed-reference OTU picking in :ref:`illumina_overview_tutorial`.

Open-reference OTU picking
--------------------------

With uclust::

	pick_open_reference_otus.py -i seqs.fna -o or_uc/ -r $reference_seqs

With usearch61::

	pick_open_reference_otus.py -i seqs.fna -o or_us/ -r $reference_seqs
                        
You can find an additional example using open-reference OTU picking in :ref:`open_reference_illumina`.