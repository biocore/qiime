.. _unweight_fasta:

.. index:: unweight_fasta.py

*unweight_fasta.py* -- Transform fasta files with abundance weighting into unweighted
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

E.g. makes 3 fasta records from a weighted input fasta file containing the following record: 
>goodsample1_12_3 bc_val=20
AATGCTTGTCACATCGATGC



**Usage:** :file:`unweight_fasta.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta
		The input fasta file
	-o, `-`-output_file
		The output fasta filepath
	-l, `-`-label
		Sequence label used for all records. fasta label lines will look like: >label_423


**Output:**

a .fasta file


make 3 fasta records from the following record:
>goodsample1_12_3 bc_val=20
AATGCTTGTCACATCGATGC

resulting in:
>goodsample_0
AATGCTTGTCACATCGATGC
>goodsample_1
AATGCTTGTCACATCGATGC
>goodsample_2
AATGCTTGTCACATCGATGC

::

	unweight_fasta.py -i input.fna -o output.fna -l goodsample


