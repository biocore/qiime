.. _metadata-description:

===================================================
Describing samples based on their metadata
===================================================

Several scripts in QIIME, including `filter_samples_from_otu_table.py <../scripts/filter_samples_from_otu_table.html>`_, `filter_distance_matrix.py <../scripts/filter_distance_matrix.html>`_, and `filter_fasta.py <../scripts/filter_fasta.html>`_, allow you to describe samples based on metadata about those samples in the mapping file. These scripts will require two parameters to achieve this: the mapping file (usually ``-m``) and the state string (usually ``-valid_states``). This section describes how to use the state string to describe a set of samples. When you *describe a set of samples*, you're usually doing that to tell QIIME that you want to retain that set of samples in some filtering operation.

A state string will contain one or more metadata field names (i.e., column header(s) from the mapping file, call this ``FIELD``), and one or more values of that field (call this ``VALUE``). Special characters used in state strings include ``:`` (delimiter between a ``FIELD`` and a ``VALUE``), ``,`` (delimiter between two ``VALUE``s associated with a single ``FIELD``), ``*`` (special ``VALUE`` that indicates any value), ``!`` (negation modifier, which indicates any value other than the one specified), and ``;`` (delimiter between multiple ``FIELD``/``VALUE`` combinations).

To describe the values associated with some field with a state string, that state string should look like: ``FIELD:VALUE``. For example, if you have a ``BodySite`` field in your mapping file, and you want to describe all of the samples that contain the value ``Gut`` in that field in the mapping file, your state string should look like ``BodySite:Gut``. To negate this, to describe all samples that don't have ``Gut`` in the ``BodySite`` field, your state string should look like ``BodySite:*,!Gut``. In this case the ``*`` is saying any value, and the ``!Gut`` is saying *except Gut*.

To describe data based on more than one field, you can separate ``FIELD:VALUE`` pairs with a ``;``, as in ``FIELD1:VALUE1;FIELD2:VALUE2``. For example, ``BodySite:Gut;Age:42`` would describe all samples with the value ``Gut`` in the ``BodySite`` field, and ``42`` in the ``Age`` field. You can of course use the negation operator here as well. For example, ``BodySite:*,!Gut;Age:42`` would describe all samples that don't have the value ``Gut`` in the ``BodySite`` field, and ``42`` in the ``Age`` field.

**IMPORTANT**: When passing state strings on the command line, you must put these in single quotes to avoid the shell interpreting these directly. For example::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_not_control.biom -m map.txt -s 'Treatment:*,!Control'