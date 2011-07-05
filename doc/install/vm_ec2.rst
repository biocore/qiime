.. _ec2:

QIIME EC2 Image
^^^^^^^^^^^^^^^^^

What is the QIIME EC2 image?
==============================
The QIIME EC2 image is a pre-built image containing QIIME and its dependencies developed for the Amazon Cloud. Any problems should be brought up in the `QIIME Forum <http://groups.google.com/group/qiime-forum>`_.

The QIIME EC2 image will allow users to load a fully-functional QIIME environment in their own Amazon Web Services (AWS) accounts to make use of the compute resources available in Amazon's cloud. 

Loading the QIIME EC2 image in AWS
==================================

To launch the QIIME EC2 image on the Amazon Cloud you will first need an AWS account. Create one the `AWS`_ site if you don't already have one. To get started with using  AWS, check out their documentation on the `AWS console`_ page. You can load the QIIME EC2 image from the `AWS console`_ by searching for ``qiime1.3.0``, which is the name of the QIIME 1.3.0 EC2 image.

Additional information on getting started with the EC image
===========================================================
A video providing detailed instructions on booting and working with the QIIME EC2 image can be found `here <http://www.youtube.com/watch?v=PEcSL_7D-jo>`_.

A very nice tutorial on getting started with EC2 can be found `here <http://ged.msu.edu/angus/tutorials-2011/day1.html>`_. (Note that this is not maintained by the QIIME development group.)

Getting data into and out of your QIIME EC2 instance
====================================================

After you've launched your EC2 instance you'll likely want to get data into it. After completing your analysis, you'll then want to get your data out. It's a good idea to pull your important data back out of your EC2 instance as soon as possible to avoid any issues with losing data if any Amazon systems go down, or if you accidentally terminate your instance. 

You can get data into and out of your running EC2 instance with a graphical tool or via the command line. Which ever you find more convenient. Command line tools such as scp and sftp work fine for this. A good, free graphical tool for Windows and OS X is `Cyberduck <http://cyberduck.ch/>`_ -- this supports interacting with S3 buckets as well as EC2 instances (via sftp).

Working with command line tools
-------------------------------
The primary tool for moving data into and out of your EC2 instance is ``scp``.

To copy data into your EC2 instance you'll need to know that path to your key file, and the Public DNS entry for your EC2 instance (see Figure 1 for where to find the public DNS entry). 

	.. image:: ../images/public_dns_entry.png
	   :width: 700
	Figure 1: Identifying the Public DNS entry.

Your ``scp`` command will look like the following to transfer data into your EC2 instance::

	scp -i PATH_TO_KEY PATH_TO_SOURCE_FILE ubuntu@PUBLIC_DNS_ENTRY:PATH_TO_DESTINATION_FILE

To transfer a file out of your EC2 instance, your command will look like this::

	scp -i PATH_TO_KEY ubuntu@PUBLIC_DNS_ENTRY:PATH_TO_SOURCE_FILE PATH_TO_DESTINATION_FILE

Examples of each are::

	scp -i ~/Documents/keys/caporaso-geobio.pem ubuntu@ec2-50-19-129-122.compute-1.amazonaws.com:~/Fasting_Map.txt ~/code/Qiime/qiime_tutorial/Fasting_Map.txt
	scp -i ~/Documents/keys/caporaso-geobio.pem ~/code/Qiime/qiime_tutorial/Fasting_Map.txt ubuntu@ec2-50-19-129-122.compute-1.amazonaws.com:~/


Working with Cyberduck
----------------------

 1. Download and install `Cyberduck <http://cyberduck.ch/>`_. Installing should be straight-forward on OS X and Windows. 

 2. Launch Cyberduck by using your normal procedure for launching programs. You should see the Cyberduck window open. A new window or sheet will open. 

 3. To open a connection to your EC2 instance, click the 'Open Connection' icon on the top-left. Choose "SFTP (SSH File Transfer Protocol)" in the dropdown box (not Amazon Simple Storage Service!). In the ``Server`` field, enter the text from the ``Public DNS`` field associated with this instance (see Figure 1 for where to find this information). In the ``Username`` field enter ``ubuntu`` (exactly as written here). Do not enter a password. Expand the ``More Options`` section in the bottom of this window, and click the ``Use Public Key Authentication``. A dialog will open - navigate to where you've stored your key and select the key that is associated with this instance. After filling in the field, the window should look like that in Figure 2. Click ``Connect``.


	.. image:: ../images/cyberduck_open_connection.png
	Figure 2: Initiating a new connection.

 4. After connecting, you'll see a file browser on the remote system (Figure 3). If you've just created this instance you won't see any files listed (because you haven't put anything there yet). You can drag-and-drop files from your local system to your EC2 instance now (Figure 4). Those files will now be visible on the EC2 instance (Figure 5).

	.. image:: ../images/cyberduck_opened_connection.png
	Figure 3. Remote file browser.

	.. image:: ../images/cyberduck_drag_and_drop_input_file.png
	   :width: 700
	Figure 4. Drag and drop a local file to the EC2 instance.

	.. image:: ../images/cyberduck_view_file.png
	Figure 5. Transferred file is visible on the EC2 instance.

 5. When you run commands on the remote system, new output files and directories will be created (Figure 6). To drag those back to your local system, first switch to the Cyberduck window and hit the refresh icon. You should then see the newly created files, and you can drag and drop them back to your local system (Figure 7).

	.. image:: ../images/cyberduck_create_output_file.png
	Figure 6. Create an output file on the EC2 instance.

	.. image:: ../images/cyberduck_drag_and_drop_output_file.png
	   :width: 700
	Figure 7. Transfer file from the EC2 instance to the local system.


Get help with Cyberduck `here <http://trac.cyberduck.ch/wiki/help/en>`_.


Other options for getting data in and out of EC2 instances
----------------------------------------------------------
Other options for getting data into and out of your EC2 instance include installing `Dropbox <http://www.dropbox.com>`_ in your instance, and then putting files in your Dropbox folders on your local system and the EC2 instance; or using the ``wget`` command to pull data from URLs. An example of the latter might look like the following::

	wget http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/Caporaso_Reference_OTUs/gg_otus_4feb2011.tgz

This will download the ``gg_otus_4feb2011.tgz`` file (the Greengenes reference OTUs) to your EC2 instance. You can find details on `wget here <http://www.gnu.org/software/wget/>`_.

Acknowledgements
================
We'd like to acknowledge the support of AWS: the QIIME EC2 image was developed using an AWS in Education (for Researchers) grant of compute resource time. 

The QIIME VB and EC2 images are built using `CloVR`_.  CloVR provides a platform for building portable virtual machines. The platform automates builds in formats compatible with VirtualBox, VMware, and Clouds, including Amazon EC2.  The `CloVR developer <http://clovr.org/developers>`_ pages have more information on the platform and build process.

.. _CloVR: http://clovr.org
.. _AWS: http://aws.amazon.com/
.. _AWS console: http://aws.amazon.com/console/