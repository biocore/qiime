.. _compare_categories:

.. index:: compare_categories.py

*compare_categories.py* -- Analyzes distance matrices for statistical significance of sample grouping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script allows for the analysis of distance matrices using several statistical methods. These methods are Adonis, Anosim, BEST, DFA, Moran's I, MRPP, PERMANOVA, PERMDISP, RDA.

Adonis - This method takes a distance matrix and mapping file. It then identifies important points in the data and performs F-tests on the initial data, and random permutations (via shuffling) the category data. Then, it finally returns the information that was identified in the samples. It's stated that it partitions (or seperates the data) for this analysis in order to find underlying relationships.

Anosim - This method takes in a distance matrix and a mapping file. This method tests whether two or more categories are significantly different. You can specify a category in the metadata mapping file to separate samples into groups and then test whether there are significant differences between those groups. For example, you might test whether Control samples are significantly different from Fast samples. Since ANOSIM is non-parametric, significance is determined through permutations.

BEST - This method looks at the numerical environmental variables relating samples in a distance matrix. For instance, the unifrac distance matrix and pH and latitude (or any other number of variables) in soil samples, and ranks them in order of which best explain patterns in the communities. This method will only accept categories that are numerical.

Moran's I - This method takes in a distance matrix and mapping file. Then it uses the geographical data supplied to identify what type of spatial configuration occurs in the samples. Are they dispersed, clustered, or of no distinctly noticeable configuration when compared to each other? This method will only accept a category that is numerical.

MRPP - This method takes in a distance matrix and a mapping file. It then tests whether two or more categories are significantly different. You can specify a category in the metadata mapping file to separate samples into groups and then test whether there are significant differences between those groups. For example, you might test whether Control samples are significantly different from Fast samples. Since MRPP is non-parametric, significance is determined through permutations.

PERMANOVA - This method takes distance matrix and a mapping file. This method is for testing the simultaneous response of one or more variables to one or more factors in an ANOVA experimental design on the basis of any distance metric. It returns a R value and a P value.

PERMDISP - This method takes a distance matrix and a mapping file. This method is a procedure for the analysis of multivariate homogeneity of group dispersions (variances). Permutations can be utilized to measure the dissimilatities between groups.

RDA - This method takes a distance matrix and a mapping file. This method is an ordination method that shows grouping/clustering of samples based on a category in the metadata mapping file and a distance matrix. This category is used to explain the variability between samples. Thus, RDA is similar to PCoA except that it is constrained, while PCoA is unconstrained (you must specify which category should be used to explain the variability in your data).

For more information and examples pertaining to this script, please refer to the accompanying tutorial, which can be found at http://qiime.org/tutorials/category_comparison.html.



**Usage:** :file:`compare_categories.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	`-`-method
		The category analysis method. Valid options: [adonis, anosim, best, morans_i, mrpp, permanova, permdisp, rda]
	-i, `-`-input_dm
		The input distance matrix
	-m, `-`-mapping_file
		The metadata mapping file
	-c, `-`-categories
		A comma-delimited list of categories from the mapping file (NOTE: many methods take just a single category, if multiple are passed only the first will be selected.)
	-o, `-`-output_dir
		Path to the output directory
	
	**[OPTIONAL]**
		
	-n, `-`-num_permutations
		The number of permutations to perform. Only applies to adonis, anosim, mrpp, and permanova [default: 999]


**Output:**


Adonis:
One file is created and outputs the results into it. The results will be: Analysis of variance(AOV) table, degrees of freedom, sequential sums of squares, mean squares, F statistics, partial R-squared and p-values, based on the N permutations.

Anosim:
One file is output to the designated location under the name of anosim_results.txt. The information in the file will be an R-value and p-value.

Best:
This outputs one file 'best_results.txt' It will have teh method name, The number of variables. The list of varibles supplied. And lastly, the RHO values, which are ranked pearson correlations for the best combination of variables that describe the community.

Moran's I:
The output file is placed in a directory specified by -o. The file will be a text file with 4 values: observed, expected, sd, and p.value. The observed value is Morans I index of x. This is computed based on the values passed in to be compared with the weights. The expected value is the value of I under the null hypothesis. The sd is the standard deviation of I under the null hypothesis. P Value is the p-value of the test of the null hypothesis against the alternative hypothesis specified in alternative. Each of these values, except for the p-value, should be between -1 and 1.

MRPP:
The command in the previous section creates a single output file in the directory specified by the -o arguement, if not it will be sent to the directory location it was called from. The file will be named mrpp_results.txt. The file will contain a dissimilarity index, the class mean and counts. It will also conatin information about the chance corrected within-group agreement A, as well as the result Based on observed delta, and expected delta. There will also be the Significance of delta and the amount of permutations performed.

PERMANOVA:
Permanova returns one output file containing the the file passed in, the F-value and the p-value.

PERMDISP:
This method returns one file that outputs an analysis of varaiance table. Responses with the distances will be shown. There will be the strata relationship, then the sample information as well. Lastly you will be able to see the f-value and p-value.

RDA:
RDA outputs a two files. One is calles rda_results.txt, the other file is rda_plot.pdf. rda.txt contains the Inertia Proportion Rank, the Eigenvalues for constrained axes, and the Eigenvalues for unconstrained axes.



**Adonis:**

Performs the Adonis statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category and 999 permutations. Then it outputs the results to the 'adonis' directory. The full file path will be: ./adonis/adonis_results.txt

::

	compare_categories.py --method adonis -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o adonis -n 999

**Anosim:**

Performs the Anosim statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category and 999 perutations. Then it outputs the results to the 'anosim' directory. The full file path will be: ./anosim/anosim_results.txt

::

	compare_categories.py --method anosim -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim -n 999

**BEST:**

Performs the BEST statistical method on a distance matrix and mapping file using the LATITUDE and LONGITUDE categories. Then it outputs the results to the 'best' directory. The full file path will be: ./best/best_results.txt

::

	compare_categories.py --method best -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c LATITUDE,LONGITUDE -o best

**Moran's I:**

Performs the Moran's I statistical method on a distance matrix and mapping file using the PH category. Then it outputs the results to the 'morans_i' directory. The full file path will be: ./morans_i/Morans_I_results.txt

::

	compare_categories.py --method morans_i -i  datasets/88_soils/unweighted_unifrac_dm.txt -m datasets/88_soils/map.txt -c PH -o morans_i

**MRPP:**

Performs the MRPP statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category. Then it outputs the results to the 'mrpp' directory. The full file path will be: ./mrpp/mrpp_results.txt

::

	compare_categories.py --method mrpp -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o mrpp -n 999

**PERMANOVA:**

Performs the PERMANOVA statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category. Then it outputs the results to the 'permanova' directory. The full file path will be: ./permanova/permanova_results.txt

::

	compare_categories.py --method permanova -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova -n 999

**PERMDISP:**

Performs the PERMDISP statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category. Then it outputs the results to the 'permdisp' directory. The full file path will be: ./permdisp/betadisper_results.txt

::

	compare_categories.py --method permdisp -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permdisp

**RDA:**

Performs the RDA statistical method on a distance matrix and mapping file using the HOST_SUBJECT_ID category. Then it outputs the results to the 'rda' directory. The full file path will be: ./RDA/rda_results.txt and ./RDA/rda_plot.txt

::

	compare_categories.py --method rda -i datasets/keyboard/unweighted_unifrac_dm.txt -m datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda


