.. _upgrade:

Upgrading to the latest version of QIIME
========================================
This page describes various ways to upgrade your QIIME installation to a new version.

For a detailed history of the most notable changes between versions see the `QIIME ChangeLog <https://github.com/biocore/qiime/blob/master/ChangeLog.md>`_. This points to the most recent version, so will usually also give a preview of the changes that are coming in the next version of QIIME.

Linux users
-----------
If you have QIIME installed on a Linux machine, the easiest way to upgrade is by using the `qiime-deploy <https://github.com/qiime/qiime-deploy>`_ tool. The `qiime-deploy <https://github.com/qiime/qiime-deploy>`_ website has instructions for how to upgrade an existing installation.

Virtual Machine (VM) users
--------------------------
If you run QIIME from within a virtual machine, you can either download the latest QIIME image or upgrade an existing one (e.g. to keep your settings, files, etc.) by following the :ref:`virtual_box` install instructions.

Amazon Web Services (AWS) users
-------------------------------
If you run QIIME using the AWS EC2 service, visit `this page <http://qiime.org/home_static/dataFiles.html>`_ for the latest QIIME AMI.

Mac OS X users / Native QIIME installation
------------------------------------------
If you run QIIME on Mac OS X or have a native install (e.g. you are not using `qiime-deploy <https://github.com/qiime/qiime-deploy>`_ to install QIIME on Linux), you will need to manually upgrade the versions of QIIME's dependencies that have changed. To see a complete list of required dependency versions, please see `this page <install.html>`_. For a list of dependency versions that have changed between QIIME releases, please see the `QIIME ChangeLog <https://github.com/biocore/qiime/blob/master/ChangeLog.md>`_.
