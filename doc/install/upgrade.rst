.. _upgrade:

Upgrading to the latest version of QIIME
========================================

QIIME 1.2.1 to QIIME 1.3.0
--------------------------

 * We have added support for version 2.2 of the RDP classifier in QIIME, and the QIIME wrappers for this are not backward compatible. You can grab version 2.2 of the RDP classifier `here <http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.2.zip/download>`_, and then follow the instructions `here <./install.html#rdp-install>`_ for making it accessible to QIIME.
 * QIIME 1.3.0 depends on PyCogent 1.5.1. You can get this `here <http://sourceforge.net/projects/pycogent/files/PyCogent/1.5.1/PyCogent-1.5.1.tgz/download>`_. You should download this, unzip it, and install either using ``setup.py`` or by putting the ``PyCogent-1.5.1`` directory in your ``$PYTHONPATH``.
 * Denoiser has been integrated into QIIME, so no longer requires an external version. You can remove previous versions of Denoiser from your system and the reference to Denoiser from your ``$PYHTONPATH`` (if applicable).