.. _upgrade:

Upgrading to the latest version of QIIME
========================================
This page documents how to install a new version of QIIME from an installation of the previous version. It should also be possible to follow the steps to upgrade from older versions to the most recent version, but it's best to review all steps in that case to ensure that you won't be performing unnecessary steps. 

For a detailed history of the most notable changes between versions see the `QIIME ChangeLog <http://qiime.svn.sourceforge.net/viewvc/qiime/trunk/ChangeLog?view=markup>`_. This points to the most recent version, so will usually also give a preview of the changes that are coming in the next version of QIIME.


1.2.1 to 1.3.0
--------------------------

 * We have added support for version 2.2 of the RDP classifier in QIIME, and the QIIME wrappers for this are not backward compatible. You can grab version 2.2 of the RDP classifier `here <http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.2.zip/download>`_, and then follow the instructions `here <./install.html#rdp-install>`_ for making it accessible to QIIME. You'll be able to use assign_taxonomy.py with either RDP 2.0.1 or RDP 2.2.
 * QIIME 1.3.0 depends on PyCogent 1.5.1. You can get this `here <http://sourceforge.net/projects/pycogent/files/PyCogent/1.5.1/PyCogent-1.5.1.tgz/download>`_. You should download this, unzip it, and install either using ``setup.py`` or by putting the ``PyCogent-1.5.1`` directory in your ``$PYTHONPATH``.
 * Denoiser has been integrated into QIIME, so no longer requires an external version. You can remove previous versions of Denoiser from your system and the reference to Denoiser from your ``$PYHTONPATH`` (if applicable). The ``denoise.py`` script has been replaced by ``denoise_wrapper.py``.
 * PyroNoise is no longer supported in QIIME, but has been replaced with AmpliconNoise, the updated program by Chris Quince.