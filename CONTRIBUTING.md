Contributing to QIIME
=====================

QIIME is an open source software package with an open developer community. You can find the source code and test code for QIIME under public revision control in `the QIIME git repository on [GitHub](https://github.com/qiime/qiime). While we have a [core development group](https://github.com/organizations/qiime/teams/265429), we very much welcome contributions from other users.

This document covers what you should do to get started with contributing to QIIME.

Coding Guidelines
-----------------

We adhere to the [PEP 8](http://www.python.org/dev/peps/pep-0008/) python coding guidelines. Before submitting any code to QIIME, you should read these carefully and apply the guidelines in your code.

On reviewing QIIME, you will notice that all of our code isn't up to PEP 8 standards. This is something that we're striving toward, and as such all new contributions must adhere.

What submissions are we interested in?
--------------------------------------

Some of the types of contributions we're interested in are new features (big or small, but for big ones it's generally a good idea to ask us if we're interested in including it before starting development), bug fixes, and documentation updates, additions, and fixes.

When considering submitting a new feature to QIIME, you should begin by posting an issue to the [QIIME issue tracker](https://github.com/qiime/qiime/issues). The information that you include in that post will differ based on the type of contribution.

For new features, you'll want to describe why the functionality that you are proposing to add is relevant. For it to be relevant, it should be demonstrably useful to QIIME users. This typically means that a new analytic method is implemented (you should describe why it's useful, ideally including a link to a paper that uses this method), or an existing method is enhanced (your implementation matches the performance of the pre-existing method while reducing runtime, memory consumption, etc, or it improves performance over the pre-existing method).

For bug fixes, you should provide a detailed description of the bug so other developers can reproduce it. You should include the following information in your bug report:

 #. The exact command or function call that you issue to create the bug.
 #. A link to all necessary input files for reproducing the bug. These files should only be as large as necessary to create the bug. For example, if you have an input file with 10,000 fasta-formatted sequences but the error only arises due to one of the sequences, create a new fasta file with only that sequence, run the command that was giving you problems, and verify that you still get an error. Then post that command and link to the trimmed fasta file. This is *extremely* useful to other developer, and it is likely that if you don't provide this information you'll get a response asking for it. Often this process helps you to better understand the bug as well.

For documentation additions, you should first post an issue describing what you propose to add, where you'd like to add it in the documentation, and a description of why you think it's an important addition. For documentation improvements and fixes, you should post an issue describing what is currently wrong or missing, and how you propose to address it.

When you post your issue, the QIIME developers will respond to let you know if we agree with the addition or change. It's very important that you go through this step to avoid wasting time working on a feature that we are not interested in including in QIIME.

Code Review
-----------

When you submit code to QIIME, it will be reviewed by one or more QIIME developers. These reviews are intended to confirm a few points:

 *. Your code is sufficiently well-tested.
 *. Your code adheres to our Coding Guidelines.
 *. Your code is sufficiently well-documented.
 *. Your code provides relevant changes or additions to QIIME.

This process is designed to ensure the quality of QIIME, and can be a very useful experience for new developers. If you'd like feedback on your code in the form of a code review as you work, you can request help in the issue that you created and one of the QIIME developers will work with you to perform regular code reviews. This can greatly reduce development time (and frustration) so we highly recommend that new developers take advantage of this rather than submitting a single massive chunk of code. That can lead to frustration when the developer thinks they are done, but the reviewer requests large amounts of changes.

Submitting code to QIIME
------------------------

QIIME is hosted on [GitHub](http://www.github.com), and we use GitHub's [Pull Request](https://help.github.com/articles/using-pull-requests) mechanism for accepting submissions.

You should begin by [creating an issue](https://github.com/qiime/qiime/issues) describing your code. This should include a description of your code (is it a new feature, a bug fix, etc.), and note in the issue description that you're willing to work on it. 

 #. If you'll be modifying existing QIIME file(s), you'll want to get input from the primary developer of the relevant file(s) via a discussion in the issue tracker to let them know you what you'd like to do. The primary developer is the person listed in the ``__maintainer__`` variable at the top of the file. Once you hear back that it is OK to make changes (i.e., they don't have local edits, they agree with the change you'd like to make, and they're comfortable with you editing their code), you should assign the issue to yourself on GitHub.
 #. [Fork](https://help.github.com/articles/fork-a-repo) the QIIME repository on the GitHub website.
 #. Clone your forked repository to the system where you'll be developing with ``git clone``.
 #. Create a new branch that you will make your changes in with ``git checkout -b``.
 #. Ensure that you have the latest version of all files (especially important for the file(s) you plan to edit). You should do this by adding QIIME as a remote repository and then pulling from that repository. For example, to update a branch that you're working on called ``my-topic-branch`` you could do the following (note that you'll only need to run the ``git remote`` step one time)::
	
	git checkout master
	git remote add upstream https://github.com/qiime/qiime.git
	git pull upstream master
	git checkout my-topic-branch
	git merge master

 #. Run ``Qiime/tests/all_tests.py`` to confirm that tests pass before you make any changes. You may get some failures, for example if you don't have an external application (e.g., RDP Classifier) installed. It is acceptable to continue if the failing tests are unrelated to the the code your working with. However, if you want to make changes to ``align_seqs.py`` and ``test_align_seqs.py`` is failing because of missing external applications, you should not proceed until you have installed the external applications and all tests pass.
 #. Make your changes to the code and the tests. You should make regular commits to your locally cloned repository with ``git add`` and ``git commit``. Write descriptive commit messages to accompany each commit.
 #. Again, ensure that you have the latest version of all files in case some changed while you were working on your edits. 
 #. Run ``Qiime/tests/all_tests.py`` to ensure that your changes did not cause anything unexpected to break. Note that some tests may fail again because you do not have external applications installed. This is OK, but this is why it is important that you run ``all_tests.py`` prior to making changes: if the same tests fail, then you should be OK.
 #. Once the tests pass, you should commit your changes and push them to your forked repository on GitHub using ``git add``, ``git commit``, and ``git push``. 
 #. Issue a [pull request](https://help.github.com/articles/using-pull-requests) on the GitHub website to request that we merge your branch's changes into QIIME's master branch. One of the QIIME developers will get in touch to discuss your code at this stage. If we request changes (which is very common), you don't need to issue a new pull request. You should make changes on your local branch, and commit and push them to GitHub. Your pull request will update automatically.

Getting help with git
=====================

If you're new to ``git``, you'll probably find [gitref.org](http://gitref.org/) helpful.
