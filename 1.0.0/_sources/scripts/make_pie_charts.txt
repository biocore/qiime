.. _make_pie_charts:

.. index:: make_pie_charts.py

*make_pie_charts.py* -- Make pie charts based on taxonomy assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script automates the construction of pie charts showing the breakdown of taxonomy by given levels. The script creates an html file for easy visualization of all of the pie charts on the same page. It uses the taxonomy or category counts from `summarize_taxa.py <./summarize_taxa.html>`_ for combined samples by level (-i) and user specified labels for each file passed in (-l). Output will be in a randomly generated folder name within the user specified folder (-o) the default is the current working directory. There is also additional functionality that breaks each taxonomic level up by sample (-s). This will create a pie chart for each sample at each specified level. The user can also specify the number of categories displayed in a single pie charts the rest are grouped together as 'other category' (-n) default is 20.



**Usage:** :file:`make_pie_charts.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_files
		List of files with sample counts by taxonomy [REQUIRED]
	-l, `-`-labels
		List of labels for pie chart(i.e. Phylum,Class)[REQUIRED]
	
	**[OPTIONAL]**
		
	-s, `-`-sample_flag
		If True pie charts will be created for each sample
	-n, `-`-num
		Maximum number of individual categories in each pie chart. All additional categories are grouped into an "other" category. [default: 20]
	-o, `-`-dir-prefix
		Directory prefix for all analyses
	-b, `-`-colorby
		This is the samples to make pie charts for in the counts files from  `summarize_taxa.py <./summarize_taxa.html>`_. The sample name must match the name of a sample id in the header of the counts file exactly and multiple categories can be list by comma separating them without spaces. If you want to see the pie charts broken up by all samples -s is still funtional. If -s is set and -b is used  it will just be broken up by all samples. If neither -s or -b are set the  pie charts will be based on all samples put together, one for each level.  [default: None]
	-p, `-`-prefs_path
		This is the user-generated preferences file. NOTE: This is a file with a dictionary containing preferences for the analysis. The label taxonomy_coloring is used for the coloring, see example prefs file preferences_file. [default: None]
	-k, `-`-background_color
		This is the background color to use in the plots. [default: None]


**Output:**

The script generates an output folder, which contains several files. For each pie chart there is a png and a pdf file. The best way to view all of the pie charts is by opening up the file taxonomy_summary_pie_chart.html.


**Examples:**

If you wish to run the code using default parameters, you must supply a counts file (Class.txt) along with the taxon level label (Class), by using the following command:

::

	make_pie_charts.py -i Class.txt -l Class

If you want to make pie charts for multiple levels at a time (phylum.txt,class.txt,genus.txt) use the following command:

::

	make_pie_charts.py -i phylum.txt,class.txt,genus.txt -l phylum,class,genus

If you want specify an output directory (e.g. "pie_charts/", regardless of whether the directory exists, use the following command:

::

	make_pie_charts.py -i Class.txt -l Class -o pie_charts/

Additionally, if you would like to display on a set number of taxa ("-n 10") and generate pie charts for all samples ("-s"), you can use the following command:

::

	make_pie_charts.py -i Class.txt -l Class -o pie_charts/ -n 10 -s

If you would like to display generate pie charts for samples samples: 'sample1' and 'sample2' that are in the counts file header, you can use the following command:

::

	make_pie_charts.py -i Class.txt -l Class -o pie_charts/ -b sample1,sample2


