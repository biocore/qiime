.. _virtual_box:

QIIME Virtual Box
^^^^^^^^^^^^^^^^^

What is the QIIME Virtual Box?
==============================
Because of the 'pipeline' nature of QIIME, there are many external dependencies and installation can therefore be challenging for new users. The QIIME Virtual Box should get around that problem, and is a fully functional environment for analyzing microbial community surveys and visualizing results.

The QIIME Virtual Box is a virtual machine based on Ubuntu Linux which comes pre-packaged with QIIME's dependencies. This is the fastest way to get up-and-running with QIIME, and is useful for small analyses (approximately up to a full 454 run); and testing QIIME to determine if it meets your needs before investing time in installing it, for example, in your cluster environment.

Installing the QIIME Virtual Box
================================
1. Download and install the `VirtualBox`_ (VB) version for your machine.
2. Download the `32-bit QIIME Virtual Box`_ or the `64-bit QIIME Virtual Box`_. Review `these notes <#choosing-vb-version>`_ for help deciding whether to install the 32-bit or 64-bit QIIME Virtual Box. This file is large so it may take between a few minutes and a few hours depending on your Internet connection speed. You will need to unzip this file, which you can typically do by double-clicking on it.
3. Create a new virtual machine:

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In this screen type QIIME as the name for the virtual machine. Then select Linux as the Operating System, and Ubuntu (64 bit) as the version. Click Next.
  * Select the amount of RAM (memory). We recommend at least 1024MB, but the best option is based on your machine. After selecting the amount of RAM, click Next.
  * Select "Use existing hard drive", and click the folder icon next to the selector (it has a green up arrow). In the new window click 'Add', and locate the virtual hard drive that was downloaded in step 2. Click Select and then click Next.
  * In the new window click Finish.

4. Double click on the new virtual machine created -- it will be called QIIME -- to boot it for the first time.
5. Review any messages that are shown, and select whatever options are best for you.
6. When your new virtual machine boots, you will see a folder on the Desktop called 'Before_you_start'. Double click on that folder to open it, and then double click on the 'Welcome' file in that folder. This will get you started with using your QIIME virtual box.


Deciding whether to install the 32-bit or 64-bit QIIME Virtual Box
=====================================================================

.. _choosing-vb-version:

You will need to determine whether you want to install the 32-bit or 64-bit QIIME Virtual Box. If you have a 32-bit machine, you must use the 32-bit Virtual Box. If you have a 64-bit machine you can use either, but we recommend using the 64-bit Virtual Box as it will give you access to the 64-bit applications such as ``sffinfo`` and ``sfffile`` (the 454 off-instrument tools). 

If you don't know whether your machine is 32 or 64 bit, go through the following steps.

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In the new window, select Linux as the Operating System. You then have the option to choose a version. If you have the option to choose 'Ubuntu (64 bit)' as a version, your system can support a 64-bit Virtual Box and you should download the 64-bit QIIME Virtual Box. If you do not see the option 'Ubuntu (64 bit)', you should download the 32-bit QIIME Virtual Box. In the latter case you'll only see an option for 'Ubuntu'.
  * Click Cancel to cancel creating a new VirtualBox. (You'll come back to this after you've downloaded the VirtualBox image.)

Limitations of the QIIME Virtual Box
====================================
Due to licensing restrictions, we cannot package the ``sff tools`` software with the QIIME Virtual Box. These must be obtained from Roche/454. When you run QIIME's ``tests/all_tests.py`` script, you will therefore get failures associated with the wrapper for ``sff tools`` which look like:

::
	
	Failed the following unit tests, in part or whole due to missing external applications.
	Depending on the QIIME features you plan to use, this may not be critical.
	/home/qiime/Qiime/tests/test_process_sff.py


You can safely ignore these, unless you are planning to process ``sff`` files directly (rather than beginning with ``fasta`` and ``qual`` files). If you do plan to process the ``sff`` files directly, you will need to install your copy of ``sff tools`` in the virtual box.

.. _32-bit QIIME Virtual Box: http://bmf.colorado.edu/QIIME/QIIME-1.1.0-i386.vdi.gz
.. _64-bit QIIME Virtual Box: http://bmf.colorado.edu/QIIME/QIIME-1.1.0-amd64.vdi.gz
.. _VirtualBox: http://www.virtualbox.org/wiki/Downloads