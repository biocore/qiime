.. _biom_format:

===========================================
Biological Observation Matrix (biom) format
===========================================

As of QIIME 1.4.0-dev, OTU tables are represented in the newly developed biom format: a JSON-derived format rather than the previous tab-separated format. This change was made for several reasons: first, to facilitate efficient handling and storing of very large OTU tables; second, to support encapsulation of core study data (OTU table data and sample/OTU metadata) in a single file; and third, to facilitate the use of these tables between tools that support this format (e.g., passing of data between QIIME and MG-RAST).

Full documentation of the BIOM format can be found at `http://www.biom-format.org <http://www.biom-format.org>`_.

Motivation for changing the OTU Table Format
=============================================

We realize that BIOM formatted files can be less convenient to work with in many cases than classic QIIME OTU table files. We consider this to be an essential switch for several reasons detailed `here <http://biom-format.org/documentation/biom_format.html#motivation-for-the-biom-format>`_.

Converting between BIOM and tab-delimited table representations
===============================================================
Note that you can always use the ``biom convert`` command in `the biom-format package <http://biom-format.org/documentation/biom_conversion.html>`_ (a QIIME dependency, so you already have it installed), to convert BIOM files to tab-delimited text for use in spreadsheet programs.
