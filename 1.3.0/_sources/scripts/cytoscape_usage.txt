.. index:: cytoscape

*cytoscape_usage* -- Loading Results with Cytoscape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


For this visualization, we will use the real edge table (:file:`real_edge_table.txt`) and real node file (:file:`real_node_table`).

The following are the directions to use once Cytoscape is installed and open:

	1. File -> Import -> Network from Table
		a. Select: "real_edge_table.txt" file
		b. Click: Show Text File Import Options
		c. Click: Transfer first line as attribute names
		d. Select: Source Interaction = Column 1
		e. Select: Target Interaction = Column 2
		f. Click the headers for all other columns to import them (they will turn blue: note, might need to clicked more than once)

	* As a result, you should get successful import message.

	2. File -> Import -> Attribute from Table
		a. Select: "real_node_table.txt" file
		b. Click: Show Text File Import Options
		c. Click: Transfer first line as attribute name
		d. Click: Import

	* As a result, you should get successful import message.

	3. Click: VizMapper

	* The user can set different properties, e.g. node color, by double-clicking.
		* Select: as "Discrete Mapper"
		* Select: Attribute to color by

	* The user can also explore layouts by using Layout from the menu bar.
	
	Note: some layouts will take a very long time to run on large datasets.
	
