.. index:: preferences_file

*preferences_file* -- Example of a prefs file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This shows a generic overview of a preferences file (.txt):

.. note::

	* {
		* 'Name_for_color_scheme1':
		* {
			* 'column':'mapping_column1',
			* 'colors':{'Sample1':'red','Sample2':'blue'},
		* },
		* 'Name_for_color_scheme2':
		* {
			* 'column':'mapping_column2',
			* 'colors':(('red',(0,100,100)),('blue',(240,100,100)))
		* }
	* }


This shows an example of a Prefs file (.txt):

.. note::

	* {
		* 'Color_By_Hand':
		* {
			* 'column':'Hand',
			* 'colors':{'Left':'red','Right':'blue'},
		* },
		* 'Color_By_pH_gradient':
		* {
			* 'column':'pH',
			* 'colors':(('red',(0,100,100)),('blue',(240,100,100)))
		* }
	* }

