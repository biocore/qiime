***************************
QIIME coding guidelines
***************************

QIIME generally follows the `PyCogent Coding Guidelines <http://www.pycogent.org/coding_guidelines.html>`_. You **must** review this document, and adhere to these standards in your QIIME code. 


QIIME command line interface guidelines
=======================================

For many of our users, their primary mode of interaction with QIIME is the command line interfaces to the scripts. It is very important that our command line interfaces behave in a standard way, and are consistent across scripts. 

To get a good start on defining acceptable command line interfaces, you can use the file ``Qiime/scripts/make_qiime_py_file.py`` with the ``-s`` option. This will generate a template file with the option parsing infrastructure already in place.

Defining your command line interfaces
-------------------------------------

QIIME command line interfaces should adhere to the following standards:

 * All options/arguments to QIIME scripts should be passed via options. NO positional arguments are allowed. 

 * Required options must be listed first in the help string, and must end with the text ``[REQUIRED]`` 

 * Required options must be listed in curly braces on the first line followed by the values they should receive, if applicable.

 * One or more usage examples should be provided as part of the usage string.

 * Non-required options must contain their default values in the help string. Typically your help string should end with ``[default: %default]``. optparse will replace ``%default`` with the default value you provide. In some case using ``%default`` is not possible, for example if -o refers to an output file, and the default behavior is to name the output file based on the name of the input file. In that case you should be as descriptive as possible. In this case ``[default: <input_fp>.out]`` is a good way to list the default value.

 * You must test for the presence of required_options, and call parser.error if the user fails to provide a required option

 * if you have been passing dest= to parser.add_option, note that you do not need to do this -- dest will be generated from the long format option (e.g., the value passed to --input_file will be stored in input_file)
 
 * Your usage string should be defined in a variable called ``usage_str`` as a multiline string outside of the function that is using it. This makes the code look neater, and will help automated documentation tools find that text (eventually). 


Naming command line options
---------------------------

We have a minimal set of required option names for QIIME scripts. If your script takes the follow values, they must use these options:

 * input filepath or input directory: ``-i``, ``--input_fp``, ``--input_dir``
 * output filepath or output directory: ``-o``, ``--output_fp``, ``--output_dir``
 * log filepath or log directory: ``-l``, ``--log_fp``, ``--log_dir``
 * verbose output: ``-v``, ``--verbose``

These options cannot be used to pass any other values.

**The uppercase options** ``-N`` **through** ``-Z`` **are reserved for parallel scripts.** No other scripts can define options with uppercase letters ``N`` through ``Z``. 
