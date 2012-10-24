
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

After adding a bug report, there are two paths you can follow. If you are comfortable with fixing the bug, you should do that. The steps you should take are as follows:

 #. Contact the primary developer of the relevant file(s), to let them know you found a bug and would like to fix it. The primary developer is the person listed in the ``__maintainer__`` variable at the top of the file.  You should include a description of the bug, the issue number of the bug report that you created, and a proposal for how you plan to fix it.
 #. Once you hear back that it is OK to make changes (i.e., they don't have local edits, they agree with the change you'd like to make, and they're comfortable with you editing their code), you should assign the issue to yourself on GitHub.
 #. Fork the QIIME repository on the GitHub website.
 #. Make a clone of your forked repository.
 #. Create a new branch that you will make your changes in.
 #. Ensure that you have the latest version of all files (especially important for the file(s) you plan to edit).
 #. Run Qiime/tests/all_tests.py to confirm that tests pass before you make any changes. You may get some failures, for example if you don't have an external application (e.g., RDP Classifier) installed. It is acceptable to continue if the failing tests are unrelated to the the code your working with. However, if you want to make changes to align_seqs.py and test_align_seqs.py is failing because of missing external applications, you should not proceed until you have installed the external applications and all tests pass.
 #. Make your changes to the code and the tests.
 #. Run Qiime/tests/all_tests.py to ensure that your changes did not cause anything unexpected to break.
 #. Once all of the tests pass, you should commit your changes to your forked repository. Note that some tests may fail again because you do not have external applications installed. This is OK, but this is why it is important that you run all_tests.py prior to making changes: if the same tests fail, then you should be OK.
 #. Issue a pull request on the GitHub website to merge your branch's changes into QIIME's master branch.

For more details on how to make contributions to QIIME through GitHub, please refer to this `document <https://docs.google.com/document/d/1fUpTyCgpYZDK-He62WPWq4P4aXP2GZxmI9tPgzPxnYU/edit>`_.
