.. index:: Installing QIIME

========================
QIIME Installation Guide
========================

Installing Conda and QIIME on Mac OS X and Linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 1: Install the `Miniconda <http://conda.pydata.org/miniconda.html>`__ version of Conda.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conda is a Python distribution, package manager, and virtual environment
solution. While QIIME 1 is Python 2 software, we recommend installing
Miniconda with Python 3 (miniconda3), as many bioinformatics packages
are now transitioning to Python 3. You can still install Python 2
software with miniconda3 by passing the python=2.7 flag when you create
a new environment; otherwise the default Python version will be Python
3.

You should next work through the tutorials on the Miniconda install page
linked above. We recommend setting up a QIIME 1.9.1 environment for your
Miniconda installation, so pay particular attention to the virtual
environment management in this tutorial. This will eventually allow you
to install QIIME 2 alongside QIIME 1 on your system. It will also
facilitate easy installation of the many other bioinformatics tools
supported by `bioconda <https://bioconda.github.io/>`__.

Step 2: Create your qiime191 environment and install QIIME.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    conda config --add channels bioconda
    conda create -n qiime191 python=2.7 qiime matplotlib=1.4.3 mock nose

Step 3: Activate your qiime191 environment and test your QIIME installation.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    source activate qiime191
    print_qiime_config.py -t

Now you are ready to start exploring Conda and use QIIME.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Anytime you want to use QIIME 1.9.1, you’ll need to reactivate your
qiime191 environment using this command:

::

    source activate qiime191

Alternatively, you can add the above command to your .profile or .bashrc
file in your home directory so that the environment will automatically
activate every time you open a terminal window.

To exit the virtual environment, simply run the deactivate command:

::

    source deactivate

If you decide later that you don’t want the environment or its packages
anymore, deactivate the environment and then run this command:

::

    conda remove --name qiime191 --all

These instructions have been tested on Mac OS X (10.8–10.11) and Linux
(Ubuntu 14.04 and Mint 17.3 Rosa).

Installing QIIME dependencies in your Conda environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You now have the base installation and access to much of QIIME's
functionality. However, if you want to use additional dependencies, your
Conda installation of QIIME will need to have access to the executable
files for each program you'd like to invoke as a dependency.
Commonly-used QIIME dependencies include UCLUST (now included in the
base installation), RDP Classifier, BLAST, USEARCH, and VSEARCH.

Linking to executables you already have
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you already have the appropriate version of these programs installed,
you can ensure that they run the correct executable by making a symbolic
link in the bin folder of your Conda QIIME environment. To allow your
environment to access dependency executables, you can put symbolic links
to them in your environment's bin folder. The path to this folder is
displayed when you load the environment:

::

    source activate qiime191
    discarding /home/username/miniconda3/bin from PATH
    prepending /home/username/miniconda3/envs/qiime191/bin to PATH

To make a symbolic link:

::

    ln -s /path/to/myexecutable /home/username/miniconda3/envs/qiime191/bin/myexecutable

Downloading and linking to new executables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you don’t have the programs installed already, or would like to use
an install specific to your particular environment, you can make new
installations directly into your Conda environment directory.

First, it is helpful to create a variable for your Conda environment
path:

::

    QIIME_ENV=/home/username/miniconda3/envs/qiime191

Next, make a directory for your QIIME dependencies and change to that
directory:

::

    mkdir $QIIME_ENV/qiime_deps
    cd $QIIME_ENV/qiime_deps

Download the binaries for the dependency (uclust, in this case), and
link to the Conda bin (for macOS use curl instead of wget):

::

    wget http://www.drive5.com/uclust/uclustq1.2.22_i86linux64

Make sure the permissions allow you to execute the binary:

::

    chmod 755 uclustq1.2.22_i86linux64

Finally, link the binary (executable) to a properly named symbolic link
in your environment's bin:

::

    ln -s $QIIME_ENV/qiime_deps/uclustq1.2.22_i86linux64 $QIIME_ENV/bin/uclust
