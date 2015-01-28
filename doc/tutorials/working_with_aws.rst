.. _working_with_ec2:

=============================================
Working with QIIME on Amazon Web Services EC2
=============================================

This document covers general steps for working with QIIME on Amazon Web Services (AWS).  If you're completely new to AWS, you should start with the `AWS Getting Started documentation <http://aws.amazon.com/documentation/gettingstarted/>`_. We refer readers there as this will be more up-to-date than what you'll find on the QIIME website. This document covers the details that you'll need to start an AWS instance, but doesn't walk you through the AWS interface as that often changes.

Starting a QIIME EC2 instance
=============================

These details assume that you have logged into your AWS account, and are starting an EC2 instance. The information that you'll need to start a QIIME instance follows. Defaults can be used except as noted below.

AMI
---

You will need the AMI for the QIIME image that you want to run. This tells Amazon what specific virtual machine you want to launch. You can always find the latest QIIME AMI on the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_. You'll need to search for the QIIME AMI under *Community AMIs*. You must be in the *US East* region when logged into your AWS account.

Instance type
-------------

Choosing the instance type depends on what you want to do with QIIME on AWS. The available instance types and pricing are described `here <http://www.ec2instances.info/>`_. The ``m1.large`` instance should be sufficient if you're working with the QIIME tutorials, and for most QIIME analyses **if you do not plan to run parallel jobs**. See :doc:`./parallel_qiime` for more details about setting up and using parallel QIIME. For larger jobs, such as parallel OTU picking on Illumina data, you likely want to use an ``m2.4xlarge`` instance (which has more memory and processors).

Security group
--------------

To log into your QIIME instance, you'll need to have ssh access (i.e., port 22), which is enabled by default.

If you want to use the `IPython Notebook`_ on your instance, you'll need to add another security rule. The *Type* should be ``Custom TCP Rule``, the *Protocol* should be ``TCP``, the *Port Range* should be ``8888``, and the *Source* should be ``Anywhere``.


Connecting to your QIIME EC2 instance using ssh
===============================================

To connect to your EC2 instance with ssh you'll need your key file, the Public IP address for your running EC2 instance, and an ssh client. You should have created and downloaded a key file when creating your EC2 instance (or indicated that you're using one that you've already created). On OS X and Linux, be sure that the permissions are set to ``400`` on this file, or you won't be able to use it. If your key file is called ``$HOME/my_key.pem`` you can set the permissions with the command::

	chmod 400 $HOME/my_key.pem

If you're working on OS X or Linux, you already have an ssh client installed. If you're working on Windows, one free option is `MobaXterm <http://mobaxterm.mobatek.net/>`_.

To connect to your running EC2 instance using ssh from Linux or OS X you'll use a command like the following::

	ssh -i PATH_TO_KEY ubuntu@PUBLIC_IP_ADDRESS


Connecting to your QIIME EC2 instance using the IPython Notebook
================================================================

To connect to your EC2 instance using the `IPython Notebook`_, you will have first need to configured your security group to allow TCP traffic over port 8888. See the instructions above for a description of how to do that. After you've booted your instance, you'll first need to connect using ssh (see instructions in the previous section) and start the IPython Notebook server in a screen session. (See `here <unix_commands.html>`_ if you need help using screen.) You should enter the following commands to start the IPython Notebook server::

	screen
	ipython notebook

You can then connect to the server by going to ``http://PUBLIC_IP_ADDRESS:8888`` in your web browser. The password will be ``qiime``.

.. warning:: This is not a secure IPython Notebook server for several reasons. Anyone who has your public IP address and this password will be able to access your IPython Notebook server, and effectively run anything they want on your instance. For securing your IPython Notebook server, see `here <http://ipython.org/ipython-doc/2/notebook/public_server.html>`_. Note that the QIIME AMI already has an IPython profile. If you'd like to change settings for the notebook (e.g., to change the password, configure SSL, etc.) you can modify the default IPython profile located in ``$HOME/.ipython/profile_default/``.

Setting the QIIME temporary directory
=====================================

By default the QIIME temporary directory is set to ``$HOME/temp/`` which will work for most cases but it can cause issues if your dataset is large as there may not be enough space on the filesystem containing ``$HOME/temp/``. In this case you will need to set your ``temp_dir`` variable within your ``qiime_config`` file to a different directory on a file system that has more space. Note: if you are running in a parallel environment that makes use of multiple AWS instances, for example :ref:`using-qiime-with-starcluster-and-the-ipython-notebook`, the ``temp_dir`` must point to a location that is accessible (usually NFS-shared) across all of the instances.

For details on modifying your ``qiime_config``, see :doc:`../install/qiime_config`. For details on configuring parallel QIIME (which will be helpful if you're using multiple AWS instances), see :doc:`./parallel_qiime`.

Getting data into and out of your QIIME EC2 instance
====================================================

After you've launched your EC2 instance you'll likely want to get data into it. After completing your analysis, you'll then want to get your data out. It's a good idea to pull your important data back out of your EC2 instance as soon as possible to avoid any issues with losing data if any Amazon systems go down, or if you accidentally terminate your instance.

You can get data into and out of your running EC2 instance with a graphical tool or via the command line, which ever you find more convenient. Command line tools such as ``scp`` and ``sftp`` work fine for this. A good, free graphical tool for Windows and OS X is `Cyberduck <http://cyberduck.ch/>`_. This supports interacting with S3 buckets as well as EC2 instances (via sftp).

You will likely want to compress your files for each transfer. You can do this with the tar/gzip commands, or with a tool like WinZip on Windows. See :ref:`unix_commands` for help with these commands.

Working with command line tools
-------------------------------
The primary tool for moving data into and out of your EC2 instance is ``scp``.

To copy data into your EC2 instance you'll need to know that path to your key file, and the Public IP Address for your EC2 instance.

Your ``scp`` command will look like the following to transfer data into your EC2 instance::

	scp -i PATH_TO_KEY PATH_TO_SOURCE_FILE ubuntu@PUBLIC_IP_ADDRESS:PATH_TO_DESTINATION_FILE

To transfer a file out of your EC2 instance, your command will look like this::

	scp -i PATH_TO_KEY ubuntu@PUBLIC_IP_ADDRESS:PATH_TO_SOURCE_FILE PATH_TO_DESTINATION_FILE


Working with Cyberduck
----------------------

.. note:: The screenshots presented in this section may be outdated. Also, while these instructions and screenshots have you connect using the *Public DNS Entry* for your instance, you can connect with either the *Public DNS Entry* or the *Public IP Address* in exactly the same way.

Follow these steps to use Cyberduck to transfer data in and out of your AWS instance:

 1. Download and install `Cyberduck <http://cyberduck.ch/>`_. Installing should be straight-forward on OS X and Windows.

 2. Launch Cyberduck by using your normal procedure for launching programs. You should see the Cyberduck window open. A new window or sheet will open.

 3. To open a connection to your EC2 instance, click the 'Open Connection' icon on the top-left. Choose "SFTP (SSH File Transfer Protocol)" in the dropdown box (not Amazon Simple Storage Service!). In the ``Server`` field, enter the text from the ``Public DNS`` field associated with this instance (see Figure 1 for where to find this information). In the ``Username`` field enter ``ubuntu`` (exactly as written here). Do not enter a password. Expand the ``More Options`` section in the bottom of this window, and click the ``Use Public Key Authentication``. A dialog will open - navigate to where you've stored your key and select the key that is associated with this instance. After filling in the field, the window should look like that in Figure 1. Click ``Connect``.


	.. image:: ../images/cyberduck_open_connection.png

	Figure 1: Initiating a new connection.

 4. After connecting, you'll see a file browser on the remote system (Figure 2). If you've just created this instance you won't see any files listed (because you haven't put anything there yet). You can drag-and-drop files from your local system to your EC2 instance now (Figure 3). Those files will now be visible on the EC2 instance (Figure 4).

	.. image:: ../images/cyberduck_opened_connection.png

	Figure 2. Remote file browser.

	.. image:: ../images/cyberduck_drag_and_drop_input_file.png
	   :width: 700

	Figure 3. Drag and drop a local file to the EC2 instance.

	.. image:: ../images/cyberduck_view_file.png

	Figure 4. Transferred file is visible on the EC2 instance.

 5. When you run commands on the remote system, new output files and directories will be created (Figure 5). To drag those back to your local system, first switch to the Cyberduck window and hit the refresh icon. You should then see the newly created files, and you can drag and drop them back to your local system (Figure 6).

	.. image:: ../images/cyberduck_create_output_file.png

	Figure 5. Create an output file on the EC2 instance.

	.. image:: ../images/cyberduck_drag_and_drop_output_file.png
	   :width: 700

	Figure 6. Transfer file from the EC2 instance to the local system.


Get help with Cyberduck `here <http://trac.cyberduck.ch/wiki/help/en>`_.


Using wget to get data into your instance
-----------------------------------------
If you want to download data that is publicly hosted on the internet into your instance, you can likely use the ``wget`` command to pull data from its URL. An example of this might look like the following::

	wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz

This will download the ``gg_13_8_otus.tar.gz`` file (the Greengenes reference OTUs) to your EC2 instance. You can find details on `wget here <http://www.gnu.org/software/wget/>`_.

Stopping your EC2 instances
===========================
As long as your EC2 instances are running, you're paying for them by the hour. When you're temporarily done using your EC2 instance, you can stop it from the AWS Management Console by selecting the instance and choosing ``Instance Actions`` > ``Stop`` (see Figure 7). Note that you still pay a very small storage fee for stopped instances, so if you're permanently done with an instance you probably want to terminate it. You can restart a stopped EC2 instance by selecting that instance and choosing ``Instance Actions`` > ``Start``.

	.. image:: ../images/stop_instance.png
	   :width: 700

	Figure 7. Stopping or pausing an EC2 instance.

If you're permanently done with an EC2 instance, you can terminate it by selecting the instance and choosing ``Instance Actions`` > ``Terminate``. Once you've terminated an instance you can never get it back: all data in that instance, as well as any configuration changes you've made, etc, is lost forever, so be sure this is what you want to do.

.. _creating-a-volume-for-persistent-storage-across-different-launches-of-an-instance-or-different-instances:

Creating a volume for persistent storage across different launches of an instance (or different instances)
==========================================================================================================
The disk space is fairly limited on the EC2 instances. To get around this you can create a volume (the equivalent of an external hard drive) and mount that on your instance. Data that you store in this volume can be accessed across different launches of an instance, or across different instances, but can only be attached to one instance at a time.

Use the management console to create a volume. To do this, first click the ``EC2`` tab. Next, select ``Volumes`` on the left sidebar. Then click ``Create Volume``. See Figure 8.

	.. image:: ../images/create_an_ebs_volume.png
	   :width: 800

	Figure 8: Create an EBS Volume.

Next you must configure the volume you want to create. You have three options here. First, define the size of the volume. This will be based on the amount of data that you'll need to store. Creating a volume that is around 10x the size of the raw data you want to analyze should leave you plenty of disk space for your analysis. Next, you must define what ``Availability Zone`` you'd like to launch your instance in. This **must** be the same zone that your instance is running in. This information is available under the 'Description' tab associated with your running instance (see ``Zone`` toward the bottom right of Figure 7). Last, you can define an snapshot that you'd like to create your volume from. You typically won't use that here. See Figure 9.

	.. image:: ../images/configure_ebs_volume_creation.png
	   :width: 700

	Figure 9: Configure EBS volume creation.

Finally, you'll attach your volume to your instance: the equivalent of plugging the USB hard drive into the computer. To do this, click the checkbox next to your volume, select ``More`` and then ``Attach Volume``. Select the instance that you'd like to attach your volume to. If you don't see your instance it may not be running, or you may have not selected the correct ``Availability Zone``. Take note of the value associated with ``Device``. You'll need this in the next step (we'll call this the attachment point). See Figure 10.


	.. image:: ../images/configure_ebs_volume_attachment.png
	   :width: 700

	Figure 10: Configure EBS volume attachment.

ssh into your EC2 instance and run the following commands. In this example, I'm assuming that your attachment point is ``/dev/sdf/``. If it's not, replace all occurrences of ``/dev/sdf/`` with your actual attachment point.

The first time you use your volume you'll need to run this command. Do not run this command on an instance that you already have data in - that will erase your data!
::

	sudo mkfs.ext4 /dev/sdf

One your first time attaching a volume to a new instance, you'll need to run this command::

	mkdir $HOME/data

Anytime you attach or re-attach your volume to an instance (so after starting a new or stopped instance) you'll need to run these commands::

	sudo mount /dev/sdf $HOME/data
	sudo chown ubuntu $HOME/data
	sudo chgrp ubuntu $HOME/data

Once you've created your device, you only need to go through the attachment step to attach to future instances. This is the step illustrated in Figure 10. Note that you'll need to create future instances in the same availability zone as this volume if you'd like to attach this volume.

.. _using-qiime-with-starcluster-and-the-ipython-notebook:

Using QIIME with StarCluster
============================

QIIME instances can be loaded using `StarCluster`_, which provides a convenient means for starting and using virtual clusters on AWS.

To start using `StarCluster`_, you should see their install instructions. Your StarCluster config file should look like this::

	[cluster my.qiime.cluster]
	node_image_id = LATEST-QIIME-AMI # see http://qiime.org/home_static/dataFiles.html
	cluster_user = ubuntu
	keyname = YOUR-KEY
	cluster_size = DESIRED-CLUSTER-SIZE
	node_instance_type = DESIRED-INSTANCE-TYPE

Everything in CAPS should be replaced with the corresponding information.

After launching your cluster, we recommend that you connect as the ``ubuntu`` user. You can do this as follows::

	starcluster sshmaster CLUSTER_TAG -u ubuntu

where ``CLUSTER_TAG`` refers to the cluster_tag that was specified when you launched your cluster with ``starcluster start``.

To run jobs in parallel on this system, you'll next need to create the ``$HOME/.qiime_config`` file on the cluster. This file should contain the following line::

	cluster_jobs_fp	start_parallel_jobs_sc.py

Video tutorial
==============

These steps described here are also covered in the `QIIME EC2 video <http://www.youtube.com/watch?v=PEcSL_7D-jo>`_, though that video is now several years old so the AWS interface in the video may be different from what you see when you log in.

.. _AWS: http://aws.amazon.com/
.. _AWS console: http://aws.amazon.com/console/
.. _StarCluster: http://web.mit.edu/star/cluster/
.. _IPython Notebook: http://ipython.org/ipython-doc/stable/interactive/htmlnotebook.html
