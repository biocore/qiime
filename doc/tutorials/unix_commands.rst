.. _unix_commands:

=============================================
Basic Unix/Linux/OS X commands
=============================================

This document covers some very basic unix commands that will be useful working with QIIME. You can find a lot of resources both on the web and as books.

Web resources
-------------

 * Linux/Unix: http://www.linuxcommand.org/learning_the_shell.php
 * OS X: http://acad.coloradocollege.edu/dept/pc/SciCompLab/UnixTutorial/

Books
-----

 * `Unix Power Tools <http://oreilly.com/catalog/9780596003302>`_


cd : change directory
========================

Change to ``data``::

	cd data

Change to your home directory::

	cd

Change to one directory higher (e.g. if you're in ``/home/greg/data/`` change to ``/home/greg/``)::

	cd ..

mv: move/rename file or directory
=================================

Rename ``old_filename`` to ``new_filename``::

	mv old_filename new_filename

wget: pull a file from a remote location
========================================

Pull the greengenes reference OTU file (4feb2011)::

	wget http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/Caporaso_Reference_OTUs/gg_otus_4feb2011.tgz

cp: copy a file
===============

Copy ``my_file.txt`` to ``my_file.txt.backup``::

	cp my_file.txt my_file.txt.backup

Compress a directory to a tgz file
==================================

Compress the ``my_data/`` directory into a single file::

	tar -czf my_data.tgz my_data/

Uncompress a tgz file into the current directory
================================================

Expand the my_data.tgz file::

	tar -xzf my_data.tgz

Uncompress a zip file
=====================

Uncompress the QIIME tutorial data::

	unzip http://bmf.colorado.edu/QIIME/qiime_tutorial-v1.4.0.zip

Working with files
==================

Print the first 10 lines of my_file.txt to screen::

	head -n 10 my_file.txt

Print the last 10 lines of my_file.txt to screen::

	tail -n 10 my_file.txt

Print all of the my_file.txt to screen::

	cat my_file.txt

Open my_file.txt for editing with the pico text editor (`go here for help with pico <http://www.uic.edu/depts/accc/software/pine/pico.html>`_)::

	pico my_file.txt

View the file interactively in read-only mode. The stop viewing the file hit ``q``::

	less my_file.txt



