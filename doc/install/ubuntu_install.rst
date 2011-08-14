.. _ubuntu_install:

QIIME Install in Ubuntu Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide was tested using the Long Term Support (LTS) version of Ubuntu (10.04). Note that you might be able to install in other Linux distributions following these instructions.

Installing the QIIME Virtual Box
================================

1. Install Ubuntu 10.04 in your prefered way.
2. Install the building dependencies in your machine. This step requieres admin (sudo) access, if you do not have them you must ask your system administrator. Note that in general all these software should be installed already in a production machine so you might be able to run the following steps in your home directory.

::
	
	sudo apt-get --force-yes -y install libgsl0-dev openjdk-6-jdk libxml2 libxslt1.1 libxslt1-dev ant subversion build-essential zlib1g-dev libpng12-dev libfreetype6-dev mpich2 libreadline-dev gfortran unzip libmysqlclient16 libmysqlclient16-dev ghc
	
3. Run the following commands. In this example the software will be installed in your home directory under the qiime_software folder, you can replace this to any path you want.

::
	
	mkdir $HOME/qiime_software/
	wget http://bmf.colorado.edu/QIIME/app-deploy-qiime-1.3.0.tgz
	tar zxvf app-deploy-qiime-1.3.0.tgz
	cd app-deploy-qiime-1.3.0
	python app-deploy.py $HOME/qiime_software/ -f etc/qiime_1.3.0.conf

4. Happy QIIME-img!
