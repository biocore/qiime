.. _fix_arb_fasta:

.. index:: fix_arb_fasta.py

*fix_arb_fasta.py* -- Reformat ARB FASTA files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script fixes ARB FASTA formatting by repairing incorrect line break chararcters, stripping spaces and replacing "." with "-" characters.


**Usage:** :file:`fix_arb_fasta.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		Path where output will be written [default: print to screen]


**Output:**

The reformatted sequences are written to stdout or to the file path provided with -o.


**Example:**

Fix the input ARB FASTA format file arb.fasta and print the result to stdout:

::

	fix_arb_fasta.py -f arb.fasta

**Example saving to an output file:**

Fix the input ARB FASTA format file arb.fasta and print the result to fixed.fasta:

::

	fix_arb_fasta.py -f arb.fasta -o fixed.fasta


