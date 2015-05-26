.. _upgrade:

Upgrading to the latest version of QIIME
========================================
This page describes various ways to upgrade your QIIME installation to a new version. For a description of the changes between versions, see the `QIIME ChangeLog <https://github.com/biocore/qiime/blob/master/ChangeLog.md>`_.

If you're using a QIIME virtual machine
---------------------------------------

If you run QIIME from within a virtual machine, you can either download the latest QIIME image or upgrade an existing one (e.g. to keep your settings, files, etc.) by following the :ref:`virtual_box` install instructions.

If you run QIIME using the AWS EC2 service, visit `this page <http://qiime.org/home_static/dataFiles.html>`_ for the latest QIIME AMI.


If you're using MacQIIME
------------------------

See the backup/install instructions on the `MacQIIME website <http://www.wernerlab.org/software/macqiime/macqiime-installation>`_.

If you've done a QIIME native minimal (base) install
----------------------------------------------------

Just:

``pip install --upgrade qiime``

If you've done a QIIME native full install
------------------------------------------

First, upgrade the base install of QIIME with:

``pip install --upgrade qiime``

Then, you'll want to review the dependencies listed in the `QIIME full install instructions <install.html#installing-qiime-natively-with-a-full-install>`_, and add or upgrade any packages whose dependency versions have changed. If you have QIIME installed on a Linux machine, the easiest way to upgrade is by using the `qiime-deploy <https://github.com/qiime/qiime-deploy>`_ tool. The `qiime-deploy <https://github.com/qiime/qiime-deploy>`_ website has instructions for how to upgrade an existing installation.
