.. _ubuntu_install:

QIIME Install in Ubuntu Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide was tested on Ubuntu 11.10 & 12.04. Note that you might be able to install in other Linux distributions following these instructions.

Installing QIIME in Ubuntu 11.10 and Ubuntu LTS (12.04)
=======================================================

* Install Ubuntu 11.10 or 12.04 in your preferred way.
* Uncomment the universe and multiverse repositories from /etc/apt/sources.list, you can use your favorite text editor but we suggest pico for simplicity. Note that at the bottom of the screen you will have the commands to use to save, exit, etc.
  ::
  
     sudo pico /etc/apt/sources.list
* Install the building dependencies in your machine. This step requires admin (sudo) access. If you do not have sudo access, you must ask your system administrator to grant you sudo access, or to run these steps for you. In general, all of this software should already be installed but it may not be. It's therefore best to run this command before the following step.
  
  For Ubuntu 11.10:
  ::
  
     sudo apt-get --force-yes -y install libssl-dev libzmq-dev libgsl0-dev openjdk-6-jdk libxml2 libxslt1.1 libxslt1-dev ant subversion build-essential zlib1g-dev libpng12-dev libfreetype6-dev mpich2 libreadline-dev gfortran unzip libmysqlclient16 libmysqlclient-dev ghc
  For Ubuntu 12.04:
  ::
  
     sudo apt-get --force-yes -y install libssl-dev libzmq-dev libgsl0-dev openjdk-6-jdk libxml2 libxslt1.1 libxslt1-dev ant subversion build-essential zlib1g-dev libpng12-dev libfreetype6-dev mpich2 libreadline-dev gfortran unzip libmysqlclient18 libmysqlclient-dev ghc
* Run the following commands. In this example the software will be installed in your home directory under the qiime_software folder, you can replace this with any path you want.
  ::
  
     mkdir $HOME/qiime_software/
     wget ftp://thebeast.colorado.edu/pub/QIIME-v1.5.0-dependencies/app-deploy-qiime-1.5.0.tgz
     tar zxvf app-deploy-qiime-1.5.0.tgz
     cd app-deploy-qiime-1.5.0
     python app-deploy.py $HOME/qiime_software/ -f etc/qiime_1.5.0.conf --force-remove-failed-dirs
* Happy QIIME-ing!
