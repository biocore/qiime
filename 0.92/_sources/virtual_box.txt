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
2. Download the `QIIME Virtual Box`_. This file is large so it may take between a few minutes and a few hours depending on your Internet connection speed. You will need to unzip this file, which you can typically do by double-clicking on it.
3. Create a new virtual machine:

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In this screen type QIIME as the name for the virtual machine. Then select Linux as the Operating System, and Ubuntu as the version. Click Next.
  * Select the amount of RAM (memory). We recommend at least 1024MB, but the best option is based on your machine. After selecting the amount of RAM, click Next.
  * Select "Use existing hard drive", and click the folder icon next to the selector (it has a green up arrow). In the new window click 'Add', and locate the the virtual hard drive that was downloaded in step 2. Click Select and then click Next.
  * In the new window click Finish.

4. Double click on the new virtual machine created -- it will be called QIIME -- to boot it for the first time.
5. Review any messages that are shown, and select whatever options are best for you.
6. When your new virtual machine boots, you will see a folder on the Desktop called 'Before_you_start'. Double click on that folder to open it, and then double click on the 'Welcome' file in that folder. This will get you started with using your QIIME virtual box.

Limitations of the QIIME Virtual Box
====================================
Due to licensing restrictions, we cannot package the ``sff tools`` software with the QIIME Virtual Box. These must be obtained from Roche/454. When you run QIIME's ``tests/all_tests.py`` script, you will therefore get failures associated with the wrapper for ``sff tools`` which look like:

::
	
	Failed the following unit tests, in part or whole due to missing external applications.
	Depending on the QIIME features you plan to use, this may not be critical.
	/Users/caporaso/code/Qiime/tests/test_process_sff.py


You can safely ignore these, unless you are planning to process ``sff`` files directly (rather than beginning with ``fasta`` and ``qual`` files). If you do plan to process the ``sff`` files directly, you will need to install your copy of ``sff tools`` in the virtual box.

The QIIME parallel scripts are not currently supported in the Virtual Box environment. We're not certain that we'll add this capability as the parallel scripts are intended for processing large data sets, and the virtual box is intended for processing smaller data sets, but let us know if you would use parallel scripts in the virtual box environment. If we get enough requests we'll consider adding parallel support.



.. _QIIME Virtual Box: http://bmf.colorado.edu/QIIME/QIIME-0.92.vdi.gz
.. _VirtualBox: http://www.virtualbox.org/wiki/Downloads