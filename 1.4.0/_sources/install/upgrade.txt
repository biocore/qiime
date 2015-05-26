.. _upgrade:

Upgrading to the latest version of QIIME
========================================
This page documents how to install a new version of QIIME from an installation of the previous version. It should also be possible to follow the steps to upgrade from older versions to the most recent version, but it's best to review all steps in that case to ensure that you won't be performing unnecessary steps. 

For a detailed history of the most notable changes between versions see the `QIIME ChangeLog <http://qiime.svn.sourceforge.net/viewvc/qiime/trunk/ChangeLog?view=markup>`_. This points to the most recent version, so will usually also give a preview of the changes that are coming in the next version of QIIME.

1.4.0 to 1.4.0-dev (Linux)
---------------------------
If you want to use the repository version of QIIME in any Linux box, follow these steps:

::
        
        wget http://bmf.colorado.edu/QIIME/app-deploy-qiime-1.4.0.tgz
        tar zxvf app-deploy-qiime-1.4.0.tgz
        cd app-deploy-qiime-1.4.0
        python app-deploy.py /software/ -f etc/qiime_1.4.0_repository.conf
        
1.3.0 to 1.4.0 (QIIME Virtual Machines)
---------------------------------------
If you are upgrading your virtual machine, follow these commands to upgrade. or a previous app-deploy.py installation

::
        
        wget http://bmf.colorado.edu/QIIME/app-deploy-qiime-1.4.0.tgz
        tar zxvf app-deploy-qiime-1.4.0.tgz
        cd app-deploy-qiime-1.4.0
        python app-deploy.py /software/ -f etc/qiime_1.4.0.conf

1.3.0 to 1.4.0 (general)
------------------------
 * QIIME |release| depends on `Python 2.7.1 <http://www.python.org/ftp/python/2.7.1/Python-2.7.1.tgz>`_ or greater; `numpy 1.5.1 <http://sourceforge.net/projects/numpy/files/NumPy/1.5.1/numpy-1.5.1.tar.gz>`_; and `matplotlib 1.1.0 <http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.1.0/matplotlib-1.1.0.tar.gz>`_.
 * AmpliconNoise support has been updated to use its latest version 1.25, which can be found `here <http://ampliconnoise.googlecode.com/files/AmpliconNoiseV1.25.tar.gz>`_.
 * muscle support has been updated to use its latest version 3.8.31, which can be found `link <http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz>`_.
 