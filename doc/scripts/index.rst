.. _index:

=============
QIIME Scripts
=============

All QIIME analyses are performed using python (``.py``) scripts. See `the QIIME install guide <../install/install.html>`_ if you need help getting the QIIME scripts installed.

All QIIME scripts can take the ``-h`` option to provide usage information. You can get this information for the ``align_seqs.py`` script (for example) by running::

    align_seqs.py -h

The same documentation that is presented when calling a script with ``-h`` is available for all QIIME scripts at the links below.

QIIME script index
==================

.. toctree::
   :maxdepth: 1
   :glob:

   ./*

Scripts in related packages
===========================

Other packages provide scripts that are often useful for QIIME users. This isn't intended to be a comprehensive list, but rather to point users to some of the tools most frequently used with QIIME.

biom-format scripts
===================

Scripts in the ``biom-format`` project provide functionality for working with biom tables generated in QIIME.

- `biom convert <http://biom-format.org/documentation/biom_conversion.html>`_
- `biom add-metadata <http://biom-format.org/documentation/adding_metadata.html>`_
- `biom summarize <http://biom-format.org/documentation/summarizing_biom_tables.html>`_

Other documentation of the ``biom-format`` project can be found at `www.biom-format.org <http://www.biom-format.org>`_.

``biom-format`` is a dependency of QIIME, so these scripts should already be available in your QIIME environment.

emperor scripts
===============

Scripts in the ``emperor`` project provide functionality for generating interactive 3D ordination plots from QIIME ordination results.

- `make_emperor.py <http://biocore.github.io/emperor/build/html/scripts/make_emperor.html>`_

Other documentation of the ``emperor`` project can be found at `biocore.github.io/emperor <http://biocore.github.io/emperor>`_.

``emperor`` is a dependency of QIIME, so these scripts should already be available in your QIIME environment.

PICRUSt scripts
===============

Scripts in PICRUSt provide functionality for performing predictive metagenomics, and for normalizing closed-reference OTU tables by 16S copy number.

See the `PICRUSt script index <http://picrust.github.io/picrust/scripts/index.html#index-scripts>`_.

Other documentation of PICRUSt can be found at `picrust.github.io/picrust <http://picrust.github.io/picrust/>`_
