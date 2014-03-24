
Report and fixing bugs
^^^^^^^^^^^^^^^^^^^^^^

I found a bug in QIIME. What should I do?
=========================================

If you find a bug in QIIME, you should `report it as a bug in the issue tracker on GitHub <https://github.com/qiime/qiime/issues>`_. It's useful to report even small bugs with the issue tracker, and it allows us to keep records of how frequently bugs are found and fixed. You may first want to review the `Qiime Forum <http://groups.google.com/group/qiime-forum>`_ to ensure that no one has run into your issue before.

In order to report a bug on the GitHub issue tracker, you will need to create a free GitHub account.

This will alert the project managers about the bug, and the issue tracker is actively monitored by the QIIME team.

We aim to be very responsive regarding QIIME bugs, and we greatly appreciate your help in making QIIME a better software package.

Writing useful bug reports
==========================

Regardless of whether you plan to fix the bug yourself, or wait for the developer to fix it, you should write a good bug report. If you do plan to fix it, this is important for the records. If you are not comfortable with fixing the bug because you don't know what's causing it, or because you are not familiar with the code that is failing, a good bug report is likely to speed up the bug fix. Since you are the person who found it, you know the most about it, so you're in a position to help out a lot. 

The developers will need to reproduce the bug into order to fix it, so part of the purpose of a bug report is to help the developers reproduce the bug. To do this, you should include the following information in your bug report:

 #. The exact command or function call that you issue to create the bug.
 #. All necessary input files for reproducing the bug. These files should only be as large as necessary to create the bug. For example, if you have an input file with 10,000 fasta-formatted sequences but the error only arises due to one of the sequences, create a new fasta file with only that sequence, run the command that was giving you problems, and verify that you still get an error. Then post that command and the trimmed fasta file. This is *extremely* useful to the developer, and it is likely that if you don't provide this information, you'll get a response asking for it.

Fixing bugs yourself
====================

After adding a bug report, there are two paths you can follow. If you are comfortable with fixing the bug, you should do that. Otherwise, you can monitor the issue for progress. If you're going to fix the bug yourself, see :doc:`contributing_to_qiime`.