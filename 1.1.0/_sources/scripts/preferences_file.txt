.. index:: preferences_file

*preferences_file* -- Example of a prefs file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This shows a generic overview of a preferences file (.txt):

.. note::

	* {
	* 'background_color':'color_name',
	* 'sample_coloring':
		* {
			* 'Name_for_color_scheme1':
			* {
				* 'column':'mapping_column1',
				* 'colors':{'Sample1':'red','Sample2':'blue'},
			* }
		* },
	* 'MONTE_CARLO_GROUP_DISTANCES':
		* {
			* 'mapping_column1': distance_to_use1,
			* 'mapping_column2': distance_to_use2,
		* },
	* 'FIELDS':
		* [
			* 'mapping_column1',
			* 'mapping_column2',
		* ],
	* 'taxonomy_coloring':
		* {
			* 'Taxonomy_Level':
			* {
				* 'column':'summarized_otu_table_column_number2',
				* 'colors':
					* {
						* "Lineage":('color_name1',"color1_in_hex"),
						* "Lineage":('color_name2',"color2_in_hsv"),
					* }
			* }
		* }
	* }

This shows an example of a Prefs file (.txt):

.. note::

	* {
	* 'background_color':'black',
	* 'sample_coloring':
		* {
			* 'Samples':
			* {
				* 'column':'SampleID',
				* 'colors':{'Sample1':'red','Sample2':'blue'},
			* },
			* 'TreatmentType':
			* {
				* 'column':'Treatment',
				* 'colors':(('red',(0,100,100)),('blue',(240,100,100)))
			* }
		* },
	* 'MONTE_CARLO_GROUP_DISTANCES':
		* {
			* 'SampleID': 10,
			* 'Treatment': 10
		* },
	* 'FIELDS':
		* [
			* 'SampleID',
			* 'Treatment'
		* ],
	* 'taxonomy_coloring':
		* {
			* 'Level_3':
			* {
				* 'column':'3',
				* 'colors':
					* {
						* 'Root;Bacteria;Bacteroidetes;Flavobacteria':('red',(0,100,100)),
						* 'Root;Bacteria;Bacteroidetes;Sphingobacteria':('blue',(240,100,100))
					* }
			* }
		* }
	* }
