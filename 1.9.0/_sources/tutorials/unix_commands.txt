.. _unix_commands:

=====================================
Getting started with the command line
=====================================

If you're new to the command line it can seem daunting, but it's not as hard as it seems. We recommend the `Unix Shell lessons from Software Carpentry <http://software-carpentry.org/lessons.html>`_ for getting started. `Software Carpentry <http://software-carpentry.org/>`_ also offers frequent hands-on workshops on scientific computing skills that cover many topics that will make using QIIME easier.

Ensuring that long jobs keep running despite loss of network connection
=======================================================================

When you start a job on a remote server (e.g., where you are connected over ``ssh``), that job will stay running only as long as your ssh connection is active. So, if you lose your network connection (e.g., if you close your laptop or your wifi connection drops) your job will stop running. This commonly happens accidentally, so it's a good idea to follow these steps for all long-running jobs.

To avoid having your jobs stop running if your ssh connection becomes inactive, you can use the ``screen`` command. These steps illustrate basic usage of ``screen``. To start a screen session to safely run commands::

	screen

This will start a *screen session*, which looks identical to your normal command line session. You can run your commands now, and safely close the terminal window, terminate your ssh session, close your computer, etc. When you come back you'll need to ssh into the server again. After doing that, you can then type::

	screen -r

to return to the *screen session* where your command was running.

You can find additional `details on screen here <http://www.rackaid.com/blog/linux-screen-tutorial-and-how-to/>`_.
