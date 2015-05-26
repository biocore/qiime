.. _make_qiime_py_file:

.. index:: make_qiime_py_file.py

*make_qiime_py_file.py* -- Create python file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This is a script which will add headers and footers to new python files and make them executable.


**Usage:** :file:`make_qiime_py_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-s, `-`-script
		Pass if creating a script to include option parsing framework [default:False].
	-t, `-`-test
		Pass if creating a unit test file to include relevant information [default:False].
	-a, `-`-author_name
		The script author's (probably you) name to be included the header variables. This will typically need to be enclosed in quotes to handle spaces. [default:AUTHOR_NAME]
	-e, `-`-author_email
		The script author's (probably you) e-mail address to be included the header variables. [default:AUTHOR_EMAIL]
	-c, `-`-copyright
		The copyright information to be included in the header variables. [default:Copyright 2014, The QIIME Project]


**Output:**

The results of this script is either a python script, test, or library file, depending on the input parameters.


**Example usage:**

Create a new script:

::

	make_qiime_py_file.py -s -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_script.py

Create a new test file:

::

	make_qiime_py_file.py -t -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_test.py

Create a basic file (e.g., for library code):

::

	make_qiime_py_file.py -a "Greg Caporaso" -e gregcaporaso@gmail.com -o my_lib.py


