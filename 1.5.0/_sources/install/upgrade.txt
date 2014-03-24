.. _upgrade:

Upgrading to the latest version of QIIME
========================================
This page documents how to install a new version of QIIME from an installation of the previous version. It should also be possible to follow the steps to upgrade from older versions to the most recent version, but it's best to review all steps in that case to ensure that you won't be performing unnecessary steps. 

For a detailed history of the most notable changes between versions see the `QIIME ChangeLog <http://qiime.svn.sourceforge.net/viewvc/qiime/trunk/ChangeLog?view=markup>`_. This points to the most recent version, so will usually also give a preview of the changes that are coming in the next version of QIIME.

1.5.0 to 1.5.0-dev (Linux)
---------------------------
If you want to use the repository version of QIIME in any Linux box, follow these steps:

::
        
        wget ftp://thebeast.colorado.edu/pub/QIIME-v1.5.0-dependencies/app-deploy-qiime-1.5.0.tgz
        tar zxvf app-deploy-qiime-1.5.0.tgz
        cd app-deploy-qiime-1.5.0
        python app-deploy.py $HOME/qiime_software/ -f etc/qiime_1.5.0_repository.conf --force-remove-failed-dirs
        
Note that you can replace the installation folder (/software/) for any other path in your system.
        
1.4.0 to 1.5.0 (QIIME Virtual Machines)
---------------------------------------
The new version of the QIIME virtual machine (VM) is based on Ubuntu 12.04 (LTS) due to this upgrade in the operating system we can not support autoupdates. Thus we suggest downloading a new full VM and following the install instructions :ref:`virtual_box`.

1.4.0 to 1.5.0 (general)
------------------------
 * mothur support has been updated to use its latest version 1.25.0, which can be found `link_mothur <http://www.mothur.org/w/images/6/6d/Mothur.1.25.0.zip>`_.
 * rtax support has been added, which can be found `link_rtax <http://dev.davidsoergel.com/trac/rtax/raw-attachment/wiki/Releases/rtax-0.981.tgz>`_.
 * pplacer support has been added, which can be found `link_pplacer <http://matsen.fhcrc.org/pplacer/builds/pplacer-v1.1-Linux.tar.gz>`_.
 * raxml support has been updated to use its latest version 7.3.0 `link_raxml <ftp://thebeast.colorado.edu/pub/QIIME-v1.5.0-dependencies/stamatak-standard-RAxML-5_7_2012.tgz>`_.
 * userch support has been updated to use its latest version 5.2.32 `link_usearch <http://www.drive5.com/usearch/>`_.
 