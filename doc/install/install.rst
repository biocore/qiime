.. index:: Installing QIIME

========================
QIIME Installation Guide
========================

Installing QIIME on Mac OS X and Linux using Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QIIME consists of native Python 2 code and additionally wraps many external applications. These instructions describe how to perform a base installation of QIIME using `Miniconda <https://anaconda.org/>`_.

Step 1: Install Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^

Miniconda is a Python distribution, package manager, and virtual environment solution. While QIIME 1 is Python 2 software, we recommend installing
Miniconda with Python 3 (miniconda3), as many bioinformatics packages
are now transitioning to Python 3. You can still install Python 2
software with miniconda3 by passing the ``python=2.7`` flag when you create
a new environment; otherwise the default Python version will be Python
3.

Begin by downloading `Miniconda <http://conda.pydata.org/miniconda.html>`_
and following the associated installation instructions.

You should next work through the tutorials on the Miniconda install page
linked above. We recommend setting up a QIIME 1.9.1 environment for your
Miniconda installation, so pay particular attention to the virtual
environment management in this tutorial. This will eventually allow you
to install QIIME 2 alongside QIIME 1 on your system. It will also
facilitate easy installation of the many other bioinformatics tools
supported by `bioconda <https://bioconda.github.io/>`_.

You'll primarily interact with Miniconda using the ``conda`` command. 

Step 2: Create your qiime1 environment and install QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

Step 3: Activate your qiime1 environment and test your QIIME installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    source activate qiime1
    print_qiime_config.py -t

Using QIIME after installation with Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Anytime you want to use QIIME after installation with Miniconda, you’ll need to reactivate your ``qiime1`` environment using this command:

::

    source activate qiime1

To exit the virtual environment, simply run the deactivate command:

::

    source deactivate

If you decide later that you don’t want the environment or its packages
anymore, deactivate the environment and then run this command:

::

    conda remove --name qiime1 --all

These instructions have been tested on Mac OS X (10.8–10.11) and Linux
(Ubuntu 14.04 and Mint 17.3 Rosa).

Making additional dependencies accessible in your qiime1 environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You now have the QIIME 1 base installation and access to much of QIIME's
functionality. However, if you want to use additional dependencies, your
installation of QIIME will need to have access to the executable
files for each program you'd like to use. If these executables are accessible through your ``PATH`` environment, they will be accessible in your ``qiime1`` conda environment and accessible to QIIME. Some of QIIME's external dependencies may have special installation instructions. These can be found `here <./alternative.html#additional-install-notes-for-some-external-dependencies>`_.

Installing QIIME using alternative approaches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you cannot install QIIME using conda, there are a few other options. These
are described in the `alternative installation methods document <./alternative.html>`_.
