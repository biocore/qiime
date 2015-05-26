.. _ec2:

QIIME EC2 Image
^^^^^^^^^^^^^^^^^

What is the QIIME EC2 image?
==============================
The QIIME EC2 image is a pre-built image containing QIIME and its dependencies developed for the Amazon Cloud. Any problems should be brought up in the `QIIME Forum <http://groups.google.com/group/qiime-forum>`_.

The QIIME EC2 image will allow users to load a fully-functional QIIME environment in their own Amazon Web Services (AWS) accounts to make use of the compute resources available in Amazon's cloud. 

Loading the QIIME EC2 image in AWS
==================================

To launch the QIIME EC2 image on the Amazon Cloud you will first need an AWS account. Create one the `AWS`_ site if you don't already have one. To get started with using  AWS, check out their documentation on the `AWS console`_ page. You can load the QIIME EC2 image from the `AWS console`_ by searching for ``qiime1.2.1``, which is the name of the QIIME 1.2.1 EC2 image. 

Acknowledgements
================
We'd like to acknowledge the support of AWS: the QIIME EC2 image was developed using an AWS in Education (for Researchers) grant of compute resource time. 

The QIIME VB and EC2 images are built using `CloVR`_.  CloVR provides a platform for building portable virtual machines. The platform automates builds in formats compatible with VirtualBox, VMware, and Clouds, including Amazon EC2.  The `CloVR developer <http://clovr.org/developers>`_ pages have more information on the platform and build process.

.. _CloVR: http://clovr.org
.. _AWS: http://aws.amazon.com/
.. _AWS console: http://aws.amazon.com/console/