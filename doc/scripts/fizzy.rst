.. _fizzy:

.. index:: fizzy.py



*fizzy.py* -- Feature Selection.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



**Description:**

Script for that implements several information theoretic feature selection methods given data in biom and the map file 



**Usage:** :file:`fizzy.py [options]`




**Input Arguments:**

.. note::

	**[REQUIRED]**
		
	-i, `-`-input_path
		Path to biom file containing the abundance data [REQUIRED]
	-m, `-`-map_path
		Metadata mapping filepath [REQUIRED]
	-c, `-`-column_label
		Column label in the map file that contains the class 
		labels [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-output_path
		Result file produced after feature selection [default: output.txt]
	-k, `-`-n_select
		Number of features to selection [default: 15]
	-f, `-`-fs_method
		Feature selection method. Available methods are CIFE, CMIM, CondMI, Condred, ICAP, JMI, MIM, MIFS, mRMR. [default: MIM]



**Fizzy example:**

Runs feature selection using the joint mutual information (JMI) objective function with a forwars search algorithm, in abundance data in data.biom and map.txt. The class labels can be found in the "Class" column in the map file. There are 15 features being selected and the result from the feature selection is save in output.txt

::
	
	fizzy.py -i data.biom -m map.txt -c Class -k 15 -f JMI -o output.txt

