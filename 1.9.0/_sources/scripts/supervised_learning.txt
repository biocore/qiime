.. _supervised_learning:

.. index:: supervised_learning.py

*supervised_learning.py* -- Run supervised classification using OTUs as predictors and a mapping file category as class labels.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script trains a supervised classifier using OTUs (or other continuous input sample x observation data) as predictors, and a mapping file column containing discrete values as the class labels.

Outputs:
    * cv_probabilities.txt: the label probabilities for each of the         given samples. (if available)
    * mislabeling.txt: A convenient presentation of cv_probabilities         for mislabeling detection.
    * confusion_matrix.txt: confusion matrix for hold-out predictions.
    * summary.txt: a summary of the results, including the expected         generalization error of the classifier
    * feature_importance_scores.txt: a list of discriminative OTUs with their         associated importance scores (if available)

It is recommended that you remove low-depth samples and rare OTUs before running this script. This can drastically reduce the run-time, and in many circumstances will not hurt performance. It is also recommended to perform rarefaction to control for sampling effort before running this script. For example, to rarefy at depth 200, then remove OTUs present in < 10 samples run:

`single_rarefaction.py <./single_rarefaction.html>`_ -i otu_table.biom -d 200 -o otu_table_rarefied200.biom
`filter_otus_from_otu_table.py <./filter_otus_from_otu_table.html>`_ -i otu_table_rarefied200.biom -s 10 -o otu_table_rarefied200.present10.biom

For an overview of the application of supervised classification to microbiota, see PubMed ID 21039646.

This script also has the ability to collate the supervised learning results produced on an input directory. For example, in order to reduce any variation introduced through producing a rarefied OTU table, the user can run `multiple_rarefactions_even_depth.py <./multiple_rarefactions_even_depth.html>`_ on the OTU table, and then pass that directory into `supervised_learning.py <./supervised_learning.html>`_. The user can then pass a -w collate_results filepath to produce a single results file that contains the average estimated generalization error of the classified, and the pooled standard deviation (for cv5 and cv10 errortypes).

This script requires that R be installed and in the search path. To install R visit: http://www.r-project.org/. Once R is installed, run R and excecute the command "install.packages("randomForest")", then type q() to exit.


**Usage:** :file:`supervised_learning.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_data
		Input data file containing predictors (e.g. otu table) or a directory of otu tables
	-m, `-`-mapping_file
		File containing meta data (response variables)
	-c, `-`-category
		Name of meta data category to predict
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	`-`-ntree
		Number of trees in forest (more is better but slower) [default: 500]
	-e, `-`-errortype
		Type of error estimation. Valid choices are: oob, loo, cv5, cv10. oob: out-of-bag, fastest, only builds one classifier, use for quick estimates; cv5: 5-fold cross validation, provides mean and standard deviation of error, use for good estimates on very large data sets; cv10: 10-fold cross validation, provides mean and standard deviation of error, use for best estimates; loo: leave-one-out cross validation, use for small data sets (less than ~30-50 samples) [default oob]
	-w, `-`-collate_results_fp
		When passing in a directory of OTU tables that are rarefied at an even depth, this option will collate the results into a single specified output file, averaging the estimated errors and standard deviations. [default: None]


**Output:**

Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities according to the classifier.


**Simple example of random forests classifier:**

::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml

**Running with 10-fold cross-validation for improved estimates of generalization error and feature importances:**

::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml_cv10 -e cv10

**Running with 1,000 trees for improved generalization error:**

::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml_ntree1000 --ntree 1000

**Run 10-fold cross validation on a directory of OTU tables rarefied at an even depth:**

::

	supervised_learning.py -i rarefied_tables/ -m Fasting_Map.txt -c Treatment -o sl_rarefied_tables_cv10 -e cv10

**Run 10-fold cross validation on a directory of OTU tables rarefied at an even depth and collate the results into a single file:**

::

	supervised_learning.py -i rarefied_tables/ -m Fasting_Map.txt -c Treatment -o sl_rarefied_tables_cv10_sweep -e cv10 -w sl_cv10_sweep.txt


