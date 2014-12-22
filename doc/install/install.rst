.. index:: Installing QIIME

=========================
QIIME installation guide
=========================

QIIME consists of native python code and additionally wraps many external applications. This gives the user flexibility to easily build their own analysis pipelines, making use of popular microbial community analysis tools. QIIME handles the processing of input and output of these applications, so the user can spend time analyzing their data rather than parsing, writing, and converting file formats.

As a consequence of this *pipeline* architecture, **QIIME has a lot of dependencies and can (but doesn't have to) be very challenging to install**.

We therefore break this guide down into several parts:

 - First we tell you how you can `avoid installing QIIME natively`__, and instead use pre-built virtual machines or MacQIIME.

 __ vm-or-macqiime_

 - Then we show how to `install QIIME natively with a minimal (base) install.`__ This is easy, and is sufficient for the vast majority of users.

__ native-base_

 - Finally, we cover how to `install QIIME natively doing a complete install.`__ This can be very challenging, and is likely only needed by users looking to develop custom analysis workflows.

__ native-full_

.. _vm-or-macqiime:

How you can avoid installing QIIME by using a virtual machine or MacQIIME
=========================================================================

Since QIIME can be difficult to install, we have attempted to shift this burden away from our users by providing the following options:

 - `MacQIIME <http://www.wernerlab.org/software/macqiime>`_: The easiest way to get a complete, native QIIME base install (and a nearly complete full install) for Mac OS X. MacQIIME is built and maintained by the `Werner Lab <http://www.wernerlab.org/>`_.
 - `QIIME Virtual Box <./virtual_box.html>`_: A virtual machine that you run on your own hardware.
 - `QIIME Amazon Web Services virtual machine <./vm_ec2.html>`_: A virtual machine that you run on rented hardware.

**We highly recommend going with one of these options if you're new to QIIME, or just want to test it out to see if it will do what you want.**

.. _native-base:

Installing QIIME natively with a minimal (base) install
=======================================================

If you want to customize QIIME, work with QIIME in a multi-user environment (e.g., a Linux cluster), are interested in getting involved in QIIME development, or want to use the development version of QIIME, you may need to install QIIME manually.

**For the vast majority of use cases, the QIIME minimal (base) install will be sufficient.** This is much easier than a *QIIME full install*. Additionally, it's possible to add parts of the QIIME full install to an existing QIIME base install as needed, so the QIIME base install is nearly always what you should do to get started.

The QIIME base install allows use of core QIIME functionality, including, but not limited to:

 - `validate_mapping_file.py <../scripts/validate_mapping_file.html>`_
 - `split_libraries.py <../scripts/split_libraries.html>`_
 - `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_
 - `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_
 - `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_
 - `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_
 - `core_diversity_analyses.py <../scripts/core_diversity_analyses.html>`_

Prepare your build environment...
---------------------------------

First, you need to have a build environment on your machine:

 - On OS X, this involves installing the `developer tools <http://developer.apple.com/technologies/xcode.html>`_.
 - On Debian-based Linux (e.g., Ubuntu), this involves installing the ``build-essential`` package with the command: ``sudo apt-get install build-essential``

Then, just pip install qiime!
-----------------------------

The easiest way to install the latest QIIME release and its base dependencies is with pip::

	pip install numpy
	pip install qiime

If you do not have pip, the easiest way to install it is by running::

	easy_install pip

**Note:** You may need to prefix the above commands with ``sudo`` if you do not have permission to place files in the default locations. You will also need an active Internet connection.

Alternatives to pip install qiime
---------------------------------

If you don't have permission to (or don't want to) install QIIME into your system version of python, you might want to look into these alternatives:

 - `virtualenv and virtualenv-wrapper <https://virtualenvwrapper.readthedocs.org/en/latest/>`_ (this is what most of the QIIME developers use instead of system-wide QIIME installations)
 - `Anaconda <https://store.continuum.io/cshop/anaconda/>`_
 - `Canopy <https://www.enthought.com/products/canopy/>`_
 - `pip-install to a non-system directory <https://pip.pypa.io/en/latest/user_guide.html#user-installs>`_

Testing the QIIME base installation
-----------------------------------

After installing the QIIME base packages, you can test the installation::

	print_qiime_config.py -t

If the output doesn't indicate any failures, you're now ready to use QIIME. If you're working in a cluster environment, you may next want to make some customizations by `setting up your qiime config file <./qiime_config.html>`_. You should also read the documentation on `using parallel qiime <../tutorials/parallel_qiime.html>`_.

.. _native-full:

Installing QIIME natively with a full install
=============================================

You should begin by performing the `QIIME minimal (base) install`__. The following steps are not necessary for the vast majority of QIIME's use cases.

__ native-base_

Then, install the following packages. If you're installing on Linux, the easiest way to install all of these packages is using `qiime-deploy <https://github.com/qiime/qiime-deploy>`_.

Alternatively, you can manually install some or all of the following packages. These are grouped by the features that each dependency will provide access to. Installation instructions should be followed for each individual package (e.g., from the project's website or README/INSTALL file).

Alignment, tree-building, taxonomy assignment, OTU picking, and other data generation steps:

* jre1.6.0_05 (`src_jre <http://java.sun.com/javase/downloads/index.jsp>`_) (license: GPL2)
* rdp_classifier-2.2 (`src_rdp <http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.2.zip/download>`_) See :ref:`RDP install notes <rdp-install>`. (license: GPL)
* blast-2.2.22 (legacy BLAST from NCBI, *NOT* BLAST+) (`OS X <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-universal-macosx.tar.gz>`_ or `linux 32-bit <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/blast-2.2.22-ia32-linux.tar.gz>`_) (license: GNU)
* cd-hit 3.1.1 (`src_cdhit <http://www.bioinformatics.org/download/cd-hit/cd-hit-2007-0131.tar.gz>`_) (license: Free access)
* ChimeraSlayer (via microbiomeutil_2010-04-29) (`src_chimeraslayer <http://sourceforge.net/projects/microbiomeutil/files/>`_) See :ref:`ChimeraSlayer install notes <chimeraslayer-install>`.
* mothur 1.25.0 (`src_mothur <http://www.mothur.org/w/images/6/6d/Mothur.1.25.0.zip>`_) (license: GPL)
* clearcut v1.0.9 (`src_clearcut <http://www.mothur.org/w/images/9/91/Clearcut.source.zip>`_)
* raxml 7.3.0 (`src_raxml <ftp://thebeast.colorado.edu/pub/QIIME-v1.5.0-dependencies/stamatak-standard-RAxML-5_7_2012.tgz>`_)
* infernal 1.0.2 (`src_infernal <ftp://selab.janelia.org/pub/software/infernal/infernal.tar.gz>`_) (license: GPL)
* cdbtools (`src_cdbtools <ftp://occams.dfci.harvard.edu/pub/bio/tgi/software/cdbfasta/cdbfasta.tar.gz>`_)
* muscle 3.8.31 (`src_muscle <http://www.drive5.com/muscle/downloads.htm>`_) (Public domain)
* rtax 0.984 (`src_rtax <http://static.davidsoergel.com/rtax-0.984.tgz>`_) (license: BSD)
* usearch v5.2.236 and/or usearch v6.1 (`src_usearch <http://www.drive5.com/usearch/>`_) (license: see http://www.drive5.com/usearch/nonprofit_form.html) **At this stage two different versions of usearch are supported.** usearch v5.2.236 is referred to as ``usearch`` in QIIME, and usearch v6.1 is referred to as ``usearch61``.
* sumaclust v1.0.00 (`src_sumaclust <ftp://ftp.microbio.me/pub/QIIME-v1.9.0-dependencies/suma_package_V_1.0.00.tar.gz>`_) (license: CeCILL FREE SOFTWARE LICENSE AGREEMENT)
* swarm 1.2.19 (`src_swarm <https://github.com/torognes/swarm/releases/tag/1.2.19>`_) (license: GPL)
* sortmerna 2.0 (`src_sortmerna <https://github.com/biocore/sortmerna/releases/tag/2.0>`_) (license: LGPL)

Processing sff files:

* sfffile and sffinfo (optional, QIIME 1.2.0 and later contain built-in tools for processing sff files although they are about 10x slower than the tools from Roche) (license: proprietary - must be obtained from Roche/454)

Denoising 454 data:

* GNU Science Library (required by AmpliconNoise) (`src_gsl <ftp://ftp.gnu.org/gnu/gsl/gsl-1.9.tar.gz>`_)
* AmpliconNoise 1.27 (`src_ampliconnoise <http://ampliconnoise.googlecode.com/files/AmpliconNoiseV1.27.tar.gz>`_) See :ref:`AmpliconNoise install notes <ampliconnoise-install>`.
* ghc 6.8 (required by the QIIME denoiser) (`src_ghc <http://haskell.org/ghc>`_)

Network visualization:

* cytoscape v2.7.0 (`src_cytoscape <http://www.cytoscape.org/>`_) (license: LGPL)

<<<<<<< HEAD
<<<<<<< HEAD
Supervised learning (``supervised_learning.py``) and ``compare_categories.py``:

* R 3.0.2 (`src_r <http://www.r-project.org/>`_) See :ref:`R install notes <R-install>`. (license: GPL2)

Alternative normalization and OTU differential abundance testing techniques (normalize_table.py and differential_abundance.py):

* R 3.1.1

If you plan to build the QIIME documentation locally:

* Sphinx 1.0.4 (`src <http://pypi.python.org/pypi/Sphinx>`_) See :ref:`Building the QIIME documentation <build-qiime-docs>` (license: BSD)
=======
Advanced statistics, including those performed in ``supervised_learning.py``, ``detrend.py`` and ``compare_categories.py``:
>>>>>>> master
=======
Advanced statistics, including those performed in ``supervised_learning.py``, ``detrend.py`` and ``compare_categories.py``:
>>>>>>> master

* R 3.1.2 (`src_r <http://www.r-project.org/>`_) See :ref:`R install notes <R-install>`. (license: GPL2)

If you plan to use SourceTracker with QIIME:

* SourceTracker 0.9.5 (`src <http://downloads.sourceforge.net/project/sourcetracker/sourcetracker-0.9.5.tar.gz>`_) (license: GPL)

For improved performance with large BIOM files, or if you're working with BIOM format 2.x files:

* HDF5 (`web_hdf5 <http://www.hdfgroup.org/HDF5/>`_)
* h5py (`web_h5py <http://www.h5py.org>`_; this can be installed with ``pip install h5py`` after HDF5 has been installed)

Testing the QIIME full installation
-----------------------------------

After installing the packages above, you can test this for sanity by running::

    print_qiime_config.py -tf

If the output doesn't indicate any failures related to packages that you plan to use, you're now ready to use QIIME. If you're working in a cluster environment, you may next want to make some customizations by `setting up your qiime config file <./qiime_config.html>`_. You should also read the documentation on `using parallel qiime <../tutorials/parallel_qiime.html>`_.

Running the test suite
----------------------

Due to size constraints, pip-installing QIIME won't download QIIME's unit tests. If you want to run QIIME's comprehensive unit test suite after installing QIIME, you should `download the QIIME source code from GitHub <https://github.com/biocore/qiime/releases>`_ for the version of QIIME that you've installed. After unpacking the source distribution, you should run the test suite. Execute the following commands::

	cd qiime-<version>/tests/
	python all_tests.py

Where ``qiime-<version>`` will be dependent on the specific version that you've downloaded.

You will see test output on the terminal indicating test successes and failures. Some failures are OK.

The ``all_tests.py`` command will complete with a summary of test failures. Some tests may fail due to missing external applications -- these will be noted separately from other test failures. If these are related to features of QIIME that you are not using, this is acceptable. Otherwise, you'll want to ensure that you have the external applications installed correctly (be sure to check that you have the right versions, QIIME may not always work with the latest version of an external dependency), and re-run the tests.

License information for external dependencies
=============================================
We have attempted to provide accurate licensing information for the above dependencies for the convenience of our users. This information is by no means definitive and may contain errors. Any questions about licenses or the legality of specific uses of these software packages should be directed to the authors of the software. Do not rely solely on the license information presented above!

Additional install notes for some external dependencies
=======================================================

PATH Environment Variable
-------------------------

External applications used by QIIME need to be visible to the shell by existing in the executable search path (i.e., listed in the ``$PATH`` environment variable). For example, if you plan to use cd-hit, and have the cd-hit executables installed in ``$HOME/bin`` you can add this directory to your system path with the commands::

	echo "export PATH=$HOME/bin/:$PATH" >> $HOME/.bashrc
	source $HOME/.bashrc

.. _rdp-install:

RDP_JAR_PATH Environment Variable
---------------------------------

If you plan to use the RDP classifier for taxonomy assignment you must define an ``RDP_JAR_PATH`` environment variable. If you downloaded and unzipped the RDP classifier folder in ``$HOME/app/``, you can do this with the following commands::

	echo "export RDP_JAR_PATH=$HOME/app/rdp_classifier_2.2/rdp_classifier-2.2.jar" >> $HOME/.bashrc
	source $HOME/.bashrc

Note that you will need the contents inside ``rdp_classifier_2.2`` for the program to function properly.

.. _uclust-install:

uclust Install Notes
--------------------

The uclust binary must be called ``uclust``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``uclust1.2.21q_i86linux64`` as ``$HOME/bin/uclust1.2.21q_i86linux64``, we recommend creating a symbolic link to this file::

	ln -s $HOME/bin/uclust1.2.21q_i86linux64 $HOME/bin/uclust

.. _usearch-install:

usearch Install Notes
---------------------

The usearch binary must be called ``usearch``, which differs from the names of the posted binaries, but is the name of the binary if you build from source. If you've installed the binary ``usearch5.2.236_i86linux32`` as ``$HOME/bin/usearch5.2.236_i86linux32``, we recommend creating a symbolic link to this file::

	ln -s $HOME/bin/usearch5.2.236_i86linux32 $HOME/bin/usearch

.. _chimeraslayer-install:

ChimeraSlayer Install Notes
---------------------------

ChimeraSlayer can only be run from the directory where it was unpacked and built as it depends on several of its dependencies being in specific places relative to the executable (``ChimeraSlayer/ChimeraSlayer.pl``). Carefully follow the ChimeraSlayer install instructions. Then add the directory containing ``ChimeraSlayer.pl`` to your ``$PATH`` environment variable. If your ``ChimeraSlayer`` folder is in ``$HOME/app/`` you can set the ``$PATH`` environment variable as follows::

	echo "export PATH=$HOME/app/ChimeraSlayer:$PATH" >> $HOME/.bashrc
	source $HOME/.bashrc

If you're having trouble getting ChimeraSlayer to work via QIIME, you should first check to see if you can run it directly from a directory other than its install directory. For example, try running ``ChimeraSlayer.pl`` from your home directory.

Once you have configured Qiime, you can test your ChimeraSlayer install by running::

	print_qiime_config.py -tf

This includes a check for obvious problems with your ChimeraSlayer install, and should help you determine if you have it installed correctly.

.. _R-install:

R Install Notes
---------------

To install R visit http://www.r-project.org/ and follow the install instructions. Once R is installed, run R and execute the following commands::

<<<<<<< HEAD
<<<<<<< HEAD
	install.packages('randomForest')
	install.packages('optparse')
	install.packages('vegan')
	install.packages('ape')
	install.packages('MASS')
	install.packages('gtools')
	install.packages('klaR')
	install.packages('RColorBrewer')
	install.packages('biom')
	install.packages('metagenomeSeq')
	install.packages('DESeq')
	install.packages('DESeq2')
	q()
=======
=======
>>>>>>> master
    install.packages(c('ape', 'biom', 'optparse', 'RColorBrewer', 'randomForest', 'vegan'))
    source('http://bioconductor.org/biocLite.R')
    biocLite(c('DESeq2', 'metagenomeSeq'))
    q()
<<<<<<< HEAD
>>>>>>> master
=======
>>>>>>> master

.. _ampliconnoise-install:

AmpliconNoise Install Notes
---------------------------

AmpliconNoise requires that several environment variables are set. After you've installed AmpliconNoise, you can set these with the following commands (assuming your AmpliconNoise install directory is ``$HOME/AmpliconNoiseV1.27/``)::

	echo "export PATH=$HOME/AmpliconNoiseV1.27/Scripts:$HOME/AmpliconNoiseV1.27/bin:$PATH" >> $HOME/.bashrc

	echo "export PYRO_LOOKUP_FILE=$HOME/AmpliconNoiseV1.27/Data/LookUp_E123.dat" >> $HOME/.bashrc
	echo "export SEQ_LOOKUP_FILE=$HOME/AmpliconNoiseV1.27/Data/Tran.dat" >> $HOME/.bashrc

QIIME Denoiser Install Notes
----------------------------

If you do not install QIIME using ``setup.py`` and you plan to use the QIIME Denoiser, you'll need to compile the FlowgramAlignment program. To do this you'll need to have ``ghc`` installed. Then from the ``Qiime/qiime/support_files/denoiser/FlowgramAlignment/`` directory, run the following command::

	make ; make install

.. _build-qiime-docs:
