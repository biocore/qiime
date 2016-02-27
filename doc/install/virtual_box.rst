.. _virtual_box:

=================
QIIME Virtual Box
=================

As a consequence of QIIME's *pipeline* architecture, **QIIME has a lot of dependencies and can (but doesn't have to) be very challenging to install**. The QIIME Virtual Box gets around the difficulty of installation by providing a functioning QIIME full install inside an Ubuntu Linux virtual machine. You can use the QIIME Virtual Box on Mac OS X, Windows, or Linux.

It is strongly recommended that your system have 8 gigabytes or more of memory to use the QIIME Virtual Box.

Installing the QIIME Virtual Box
================================
1. Download and install the `VirtualBox`_ (VB) version for your machine.
2. Download the 64-bit QIIME Virtual Box, which is linked from the `QIIME Resources page <http://qiime.org/home_static/dataFiles.html>`_. This file is large so it may take between a few minutes and a few hours depending on your Internet connection speed. You will need to unzip this file, which you can typically do by double-clicking on it.
3. Create a new virtual machine:

  * Launch VirtualBox, and create a new machine (press the New button).
  * A new window will show up. Click 'Next'.
  * In this screen type QIIME as the name for the virtual machine. Then select Linux as the Operating System, and Ubuntu (64 bit) as the version. Click Next.
  * Select the amount of RAM (memory). You will need at least 3 GB, but the best option is based on your machine. After selecting the amount of RAM, click Next.
  * Select "Use existing hard drive", and click the folder icon next to the selector (it has a green up arrow). In the new window click 'Add', and locate the virtual hard drive that was downloaded in step 2. Click Select and then click Next.
  * In the new window click Finish.


4. Double click on the new virtual machine created -- it will be called QIIME -- to boot it for the first time.
5. Review any messages that are shown, and select whatever options are best for you.
6. When your new virtual machine boots, you will see a folder on the Desktop called 'Before_you_start'. Double click on that folder to open it, and then double click on the '0.QIIME_version.txt' file in that folder. This will get you started with using your QIIME virtual box.

VirtualBox help video
=====================
A video illustrating these steps can be found `here <http://www.youtube.com/watch?v=1jYupkquaME>`_.

Limitations of the QIIME Virtual Box
====================================
Due to licensing restrictions, we cannot package Roche's ``sff tools`` or ``usearch`` with the QIIME Virtual Box. You can find information on how to obtain these package in the `QIIME Installation Guide <./install.html#installing-qiime-natively-with-a-full-install>`_.

These are not essential for most of QIIME's functionality, so most users can likely ignore this.

Support for the 32-bit QIIME Virtual Box is discontinued
========================================================
As of QIIME 1.2.0, we no longer build and distribute the 32 bit Virtual Box. The majority of our users work on 64-bit platforms, and supporting multiple versions has become too time-consuming. We apologize for the inconvenience. If upgrading to a 64-bit system is not an option for you, you might want to look into the `QIIME EC image <./vm_ec2.html>`_.

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

If this doesn't work, you might need to instead hit F9 to enter the BIOS.

Briefly, on booting the system, hit F9. Enter Advanced Options -> Processor Options -> Intel(R) Virtualization Technology and then hit Enable.

.. _VirtualBox: http://www.virtualbox.org/wiki/Downloads
