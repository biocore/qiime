.. _contributing_to_qiime:

Contributing to QIIME
^^^^^^^^^^^^^^^^^^^^^^

QIIME is an open source software package with an open developer community. You can find the source code and test code for QIIME under public revision control in `the QIIME git repository on GitHub <https://github.com/qiime/qiime>`_.

While we have a `core development group <https://github.com/organizations/qiime/teams/265429>`_, we very much welcome contributions from other users. Some of the types of contributions we're interested in are new features (big or small, but for big ones it's generally a good idea to ask us if we're interested in including it first), bug fixes, and documentation updates, additions, and fixes. We use `pull requests <https://help.github.com/articles/using-pull-requests>`_ to review and accept user contributions. 

Follow the steps below to get started with submitting code to QIIME.

Submitting code to QIIME
========================

You should begin by `creating an issue <https://github.com/qiime/qiime/issues>`_ describing your code. This should include a description of your code (is it a new feature, a bug fix, etc.), and note in the issue description that you're willing to work on it. 

 #. If you'll be modifying existing QIIME file(s), you'll want to get input from the primary developer of the relevant file(s) via a discussion in the issue tracker to let them know you what you'd like to do. The primary developer is the person listed in the ``__maintainer__`` variable at the top of the file. Once you hear back that it is OK to make changes (i.e., they don't have local edits, they agree with the change you'd like to make, and they're comfortable with you editing their code), you should assign the issue to yourself on GitHub.
 #. `Fork <https://help.github.com/articles/fork-a-repo>`_ the QIIME repository on the GitHub website.
 #. Clone your forked repository to the system where you'll be developing with ``git clone``.
 #. Create a new branch that you will make your changes in with ``git checkout -b``.
 #. Ensure that you have the latest version of all files (especially important for the file(s) you plan to edit). You should do this by adding QIIME as a remote repository and then pulling from that repository. For example, to update a branch that you're working on called ``my-topic-branch`` you could do the following (note that you'll only need to run the ``git remote`` step one time)::
	
	git checkout master
	git remote add upstream https://github.com/qiime/qiime.git
	git pull upstream master
	git checkout my-topic-branch
	git rebase master

 #. Run ``Qiime/tests/all_tests.py`` to confirm that tests pass before you make any changes. You may get some failures, for example if you don't have an external application (e.g., RDP Classifier) installed. It is acceptable to continue if the failing tests are unrelated to the the code your working with. However, if you want to make changes to ``align_seqs.py`` and ``test_align_seqs.py`` is failing because of missing external applications, you should not proceed until you have installed the external applications and all tests pass.
 #. Make your changes to the code and the tests. You should make regular commits to your locally cloned repository with ``git add`` and ``git commit``. Write descriptive commit messages to accompany each commit.
 #. Again, ensure that you have the latest version of all files in case some changed while you were working on your edits. 
 #. Run ``Qiime/tests/all_tests.py`` to ensure that your changes did not cause anything unexpected to break. Note that some tests may fail again because you do not have external applications installed. This is OK, but this is why it is important that you run ``all_tests.py`` prior to making changes: if the same tests fail, then you should be OK.
 #. Once the tests pass, you should commit your changes and push them to your forked repository on GitHub using ``git add``, ``git commit``, and ``git push``. 
 #. Issue a `pull request <https://help.github.com/articles/using-pull-requests>`_ on the GitHub website to request that we merge your branch's changes into QIIME's master branch. One of the QIIME developers will get in touch to discuss your code at this stage. If we request changes (which is very common), you don't need to issue a new pull request. You should make changes on your local branch, and commit and push them to GitHub. Your pull request will update automatically.

Getting help with git
=====================

If you're new to ``git`` you'll probably find `gitref.org <http://gitref.org/>`_ helpful. This is what a lot of the QIIME developers used when we transitioned from Sourceforge and ``svn`` to GitHub and ``git``.

Additional notes
================

We keep general notes on using GitHub in a Google Doc `here <https://docs.google.com/document/d/1fUpTyCgpYZDK-He62WPWq4P4aXP2GZxmI9tPgzPxnYU/edit>`_. These are being updated somewhat frequently (at least at the time of this writing).