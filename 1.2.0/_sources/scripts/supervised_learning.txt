.. _supervised_learning:

.. index:: supervised_learning.py

*supervised_learning.py* -- Run supervised classification using OTUs as predictors and a mapping file category as class labels.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script trains a supervised classifier using OTUs (or other continuous input sample x observation data) as predictors, and a mapping file column containing discrete values as the class labels.

    Outputs:
        predictions.txt: the labels predicted by the classifier for the given
            samples. Each sample is predicted by a model that was trained 
            without it. 
        probabilities.txt: the label probabilities for each of the given 
            samples. (if available)
        summary.txt: a summary of the results, including the expected
            generalization error of the classifier
        features.txt: a list of discriminative OTUs with their associated
            importance scores (if available)
        params.txt: a list of any non-default parameters used in training
            the model.
    
By default, this script removes OTUs that are present in less than 10% of the samples. Run with "--show_params" to see how to set this manually. For an overview of the application of supervised classification to microbiota, see PubMed ID 21039646.

This script requires that R is installed and in the search path. To install R visit: http://www.r-project.org/. Once R is installed, run R and excecute the command "install.packages('randomForest')", then type q() to exit.


**Usage:** :file:`supervised_learning.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_data
		Input data file containing predictors (e.g. otu table)
	-m, `-`-mapping_file
		File containing meta data (response variables)
	-c, `-`-category
		Name of meta data category to predict
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-s, `-`-method
		Comma-separated list of supervised learning methods to apply. Currently one option is available: "random_forest" [default: random_forest].
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-p, `-`-param_file
		File containing parameters for the supervised learning model inference [default: None]
	`-`-show_params
		Show sample parameters file for a given method [default: None]
	-k, `-`-keepfiles
		Keep R-formatted input files [default: None]


**Output:**

Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities according to the classifier.


**Simple example of random forests classifier:**

::

	supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml

**Getting a sample params file for the random forests classifier:**

::

	supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml --show_params

**Running with a user-specified params file for the random forests classifier:**

::

	supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml -p params.txt


