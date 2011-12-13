.. _virtual_box:

QIIME Virtual Box
^^^^^^^^^^^^^^^^^

What is the QIIME Virtual Box?
==============================
Because of the 'pipeline' nature of QIIME, there are many external dependencies and installation can therefore be a challenge. The QIIME Virtual Box should get around that problem, and is a fully functional environment for analyzing microbial community surveys and visualizing results.

The QIIME Virtual Box is a virtual machine based on Ubuntu Linux which comes pre-packaged with QIIME's dependencies. This is the fastest way to get up-and-running with QIIME, and is useful for small analyses (approximately up to a full 454 run); and testing QIIME to determine if it meets your needs before investing time in installing it, for example, in your cluster environment.

Installing the QIIME Virtual Box
================================
1. Download and install the `VirtualBox`_ (VB) version for your machine.
2. Download the `64-bit QIIME Virtual Box`_. This file is large so it may take between a few minutes and a few hours depending on your Internet connection speed. You will need to unzip this file, which you can typically do by double-clicking on it.
3. Create a new virtual machine:

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In this screen type QIIME as the name for the virtual machine. Then select Linux as the Operating System, and Ubuntu (64 bit) as the version. Click Next.
  * Select the amount of RAM (memory). You will need at least 1024MB, but the best option is based on your machine. After selecting the amount of RAM, click Next.
  * Select "Use existing hard drive", and click the folder icon next to the selector (it has a green up arrow). In the new window click 'Add', and locate the virtual hard drive that was downloaded in step 2. Click Select and then click Next.
  * In the new window click Finish.


4. Double click on the new virtual machine created -- it will be called QIIME -- to boot it for the first time.
5. Review any messages that are shown, and select whatever options are best for you.
6. When your new virtual machine boots, you will see a folder on the Desktop called 'Before_you_start'. Double click on that folder to open it, and then double click on the 'Welcome' file in that folder. This will get you started with using your QIIME virtual box.

VirtualBox help video
=====================
A video illustrating these steps can be found `here <http://www.youtube.com/watch?v=1jYupkquaME>`_.

Upgrading from previous versions of the QIIME Virtual Box
=========================================================
A QIIME upgrade/deploy script is provided for existing 64-bit Virtual Box users to upgrade to the latest version of QIIME without having to reinstall the QIIME VB. To upgrade your VB follow these steps, launch your existing QIIME VB, open a terminal, and run the following commands::
	
	sudo apt-get install ghc gsl-bin
	wget http://bmf.colorado.edu/QIIME/app-deploy-qiime-1.4.0.tgz
	tar zxvf app-deploy-qiime-1.4.0.tgz
	cd app-deploy-qiime-1.4.0
	python app-deploy.py /software/ -f etc/qiime_1.4.0.conf
	exit
	
Then open a new terminal, and to confirm that you are now running the latest version of Qiime run::
	
	print_qiime_config.py -t
	
QIIME VB and CloVR
==================

As of the QIIME 1.2.0 release, the QIIME VB and EC2 images are built using `CloVR`_.  CloVR provides a platform for building portable virtual machines. The platform automates builds in formats compatible with VirtualBox, VMware, and Clouds, including Amazon EC2.  The `CloVR developer <http://clovr.org/developers>`_ pages have more information on the platform and build process.

Limitations of the QIIME Virtual Box
====================================
Due to licensing restrictions, we cannot package Roche's ``sff tools`` with the QIIME Virtual Box. These must be obtained from Roche/454. When you run QIIME's ``tests/all_tests.py`` script, you will therefore get failures associated with the wrapper for ``sff tools`` which look like:

::
	
	Failed the following unit tests, in part or whole due to missing external applications.
	Depending on the QIIME features you plan to use, this may not be critical.
	/home/qiime/Qiime/tests/test_process_sff.py


You can safely ignore these as QIIME contains built-in tools for parsing sff files. The QIIME sff tools are about 10x slower than Roche's version however, so if you have a copy of Roche's sff tools you may want to install those in the VB.

Support for the 32-bit QIIME Virtual Box is discontinued
========================================================
As of Qiime 1.2.0, we no longer build and distribute the 32 bit Virtual Box. The majority of our users work on 64 bit platforms, and supporting multiple versions has become too time-consuming. We apologize for the inconvenience. If upgrading to a 64-bit system is not an option for you, you might want to look into the `QIIME EC image <./vm_ec2.html>`_.

If you don't know whether your machine is 32 or 64 bit, go through the following steps.

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In the new window, select Linux as the Operating System. You then have the option to choose a version. If you have the option to choose 'Ubuntu (64 bit)' as a version, your system can support a 64-bit Virtual Box. If you do not see the option 'Ubuntu (64 bit)', your system is 32-bit and you will need to upgrade to a new system before you can run QIIME.

Troubleshooting
===============

Error when starting the 64-bit QIIME Virtual Box on Windows
-----------------------------------------------------------
If you see the following error message when starting the 64-bit QIIME Virtual Box on Windows, follow the instructions below.

::
	
	Hardware acceleration has been enabled but is not operational. Your 64-bit guest will fail to detect a 64-bit CPU and will not be able to boot.

By default VT (virtualization technology) comes disabled in Windows. HP has nice instructions on how to fix this, so we point users to those instructions.

Enabling virtualization in the BIOS

1. Shut down the computer and turn it back on. Repeatedly press esc key at startup.
2. Press the F10 key for BIOS Setup.
3. Press the right arrow key to System Configuration, select Virtualization Technology and then press the enter key.
4. Select Enabled and press the enter key.
5. Press the F10 key and select Yes and press the enter key to save changes.

The computer automatically restarts. If you cannot enable Virtualization Technology on in your BIOS, check if your BIOS needs to be updated.

These instructions were extracted directly from HPs website: `link 
<http://h10025.www1.hp.com/ewfrf/wc/document?docname=c01959244&cc=us&lc=en&dlc=en&product=3744198>`_.

If this doesn't work, you might need to instead hit F9 to enter the BIOS. See this `forum discussion <http://forums11.itrc.hp.com/service/forums/questionanswer.do?admit=109447626+1279028363362+28353475&threadId=1120296>`_.

Briefly, on booting the system, hit F9. Enter Advanced Options -> Processor Options -> Intel(R) Virtualization Technology and then hit Enable.

.. _CloVR: http://clovr.org
.. _64-bit QIIME Virtual Box: http://bmf.colorado.edu/QIIME/QIIME-1.4.0-amd64.vdi.gz
.. _VirtualBox: http://www.virtualbox.org/wiki/Downloads