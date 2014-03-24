.. _running_supervised_learning:

============================
Running Supervised Learning
============================

This document describes how to run supervised classification with QIIME. The goal of supervised classification is to classify new, unlabeled communities based on a set of labeled training communities. See (1_) for a general discussion of the application of supervised classification to microbiota. Supervised classification using the Random Forests (2_) classifier is implemented in the QIIME script `supervised_learning.py <../scripts/supervised_learning.html>`_. When you run this script you will get several output files:

* :file:`summary.txt`. This gives the estimated generalization error (an estimate of how much error the classifier would have on a novel data set), a description of the method used to estimate the generalization error, the expected "baseline" error for a model that always guesses the most common class, and the number of trees used in the forest.The estimated generalization error should be compared to the "baseline" error. For your convenience, this file also shows the ratio of the baseline error to the estimated generalization error of the random forests classifier. A reasonable criterion for good classification is that this ratio should be at least 2, or in other words, the random forests classifier does at least twice as well as random guessing. The baseline error will be very low for data sets with unbalanced classes. For example, if 90% of your samples belong to a single class, then your baseline error is 10%, and you should not be impressed if a classifier achieves, say, 8% error. In contrast, if your data is evenly divided into 10 classes, then your baseline error for random guessing is 90%, and even 45% classifier error is an indication that the classifier is doing quite a bit better than random. If you are interested in the performance of the classifier in a specific class, see the file :file:`confusion_matrix.txt` (explained below).
* :file:`cv_probabilities.txt`. Cross-validation estimates of class probabilities for the samples. For each sample, gives the classifier's estimated probability that the sample belonged to each of the possible classes or categories. To avoid overfitting, the estimates for a given sample are always predicted by models that did not contain that sample in their training sets.
* :file:`mislabeling.txt`. For each sample, gives the estimated probability of the alleged class (the class shown in the mapping file), and the probability for the most likely other class. This is just a convenient summary of the information contained in :file:`cv_probabilities.txt`. It also contains columns “mislabeled_probability_above_0.xx” for thresholds 5%, 10%, ..., 95%, 99%, indicating whether the probability that a given sample is mislabeled exceeds the given threshold. For example, if a particular sample has the value “TRUE” in the column “mislabeled_probability_above_0.95”, then the model has estimated that there is a 95% chance that the sample is mislabeled. See also the `tutorial on predicting mislabeled samples <predicting_mislabeled_samples.html>`_.
* :file:`feature_importance_scores.txt`. The classifier's estimates of the importance (i.e., discriminative power) of each of the features (e.g. OTUs). For Random Forests, the importance reported is the expected mean decrease in accuracy when the feature is ignored. There is no accepted threshold for determining which features are significantly discriminative, as this may depend on the overall error rate and other factors. For example, if the overall error rate is only 5%, then a feature whose removal increases error by 1% (i.e. importance score 0.01) might be considered highly discriminative.
* :file:`confusion_matrix.txt`. Entry in row i and column j shows the number of samples whose true class was i that were classified (when held out from model training) in class j.

Input files
------------------
This script requires a QIIME OTU table (or equivalent) and a QIIME metadata mapping file. 

Running the script
--------------------------------------------------------------------------------

To run supervised classification on the `QIIME tutorial <./tutorial.html>`_ data set, where the "Treatment" metadata column gives the class labels::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v
	
All of the result files described above will be contained in the folder :file:`ml`. The `-v` flag causes verbose output, including a trace of the classifier's progress while it is running. This runs Random Forests with the default setting of 500 trees. For larger data sets, the expected generalization error may decrease slightly if more trees are used. You can run random forests with 1,000 trees with the following::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v --ntree 1000

Both of these example build a single random forests classifier, and use "out-of-bag" predictions (that is, each tree in the forest makes predictions about samples that were absent from its bootstrapped set of samples) to estimate the generalization error. If you have a very small data set you may wish to perform leave-one-out cross validation, in which the class label for each sample is predicted using a separate random forests classifier trained on the other n-1 samples::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v -e loo

To obtain more robust estimates of the generalization error and feature importances (including standard deviations), you can run the script with 5-fold or 10-fold cross validation::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v -e cv5

or ::

	supervised_learning.py -i otu_table.biom -m Fasting_Map.txt -c Treatment -o ml -v -e cv10

Cautions
---------
Supervised classification is most useful for larger data sets. When data sets are too small, the estimates of the generalization error, feature importance, and class probabilities may be quite variable. How large a data set needs to be depends on, among other things, how subtle are the differences between classes, and how many noisy features (e.g. OTUs) there are.

Note: we recommend running :file:`single_rarefaction.py` on your OTU table before using it as input to :file:`supervised_learning.py`, to control for variation in sequencing effort.

References
------------
.. [1] Knights D, Costello EK, Knight R (2010). "Supervised Classification of Human Microbiota". FEMS Microbiology Reviews 35, 343-359 (`link <http://www.ncbi.nlm.nih.gov/pubmed/21039646>`_)

.. [2] Breiman L (2001). "Random forests". Maching Learning 45: 5–32. (`link <http://www.springerlink.com/content/u0p06167n6173512/>`_)





