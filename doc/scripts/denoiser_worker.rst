.. _denoiser_worker:

.. index:: denoiser_worker.py

*denoiser_worker.py* -- Start a denoiser worker process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The workers are automatically started by the `denoiser.py <./denoiser.html>`_ script.
You usually never need to use this script yourself.

A worker waits for data and does flowgram alignments once it gets it.


**Usage:** :file:`denoiser_worker.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-file_path
		Path used as prefix for worker data files[REQUIRED]
	-p, `-`-port
		Server port [REQUIRED]
	-s, `-`-server_address
		Server address[REQUIRED]
	
	**[OPTIONAL]**
		
	-e, `-`-error_profile
		Path to error profile [DEFAULT: /Users/jistombaugh/Dropbox/Qiime_work/qiime/support_files/denoiser/Data/FLX_error_profile.dat]
	-c, `-`-counter
		Round counter to start this worker with  [default: 0]


**Output:**

Denoise worker writes a log file if verbose flag is set.


Start worker and connect to server listening on port 12345 on the same machine (localhost)

::

	denoiser_worker.py -f seqs.fna -f denoiser_out/worker99 -p 12345 -s localhost


