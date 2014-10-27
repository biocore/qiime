from IPython.parallel import Client
from subprocess import Popen, PIPE

from qiime.parallel.util import ComputeError
from qiime.util import load_qiime_config


qiime_config = load_qiime_config()


def system_call(cmd):
    """Call cmd and return (stdout, stderr, return_value).

    cmd: can be either a string containing the command to be run, or a
     sequence of strings that are the tokens of the command.

    This function is ported from QIIME (http://www.qiime.org), previously
    named qiime_system_call. QIIME is a GPL project, but we obtained permission
    from the authors of this function to port it to pyqi (and keep it under
    pyqi's BSD license).
    """
    proc = Popen(cmd,
                 universal_newlines=True,
                 shell=True,
                 stdout=PIPE,
                 stderr=PIPE)

    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    if return_value != 0:
        raise ComputeError("Failed to execute: %s\nstdout: %s\nstderr: %s" %
                           (cmd, stdout, stderr))

    return stdout, stderr, return_value


class Context(object):
    """Parallel context

    Methods
    -------
    submit_async
    submit_async_deps
    submit_sync
    sync

    """
    def __init__(self, profile):
        try:
            self._client = Client(profile=profile)
        except IOError:
            raise ComputeError("It looks like the profile %s does not exist "
                               "or the cluster is not actually running."
                               % profile)

        self._stage_imports(self._client)
        self._lview = self._client.load_balanced_view()

    def _stage_imports(self, cluster):
        with cluster[:].sync_imports(quiet=True):
            from qiime.parallel.context import system_call

    def sync(self, data):
        """Sync data to engines

        Parameters
        ----------
        data : dict
            dict of objects and to sync

        """
        self._client[:].update(data)

    def submit_async_deps(self, deps, cmd, *args, **kwargs):
        """Submit as async command to execute making sure that cmd is executed
        after all its dependencies are executed

        Parameters
        ----------
        deps : list of AsyncResult
            The list of job dependencies for cmd
        cmd : {function, str}
            A function to execute or a system call to execute
        args : list
            Arguments to pass to a function (if cmd is function)
        kwargs : dict
            Keyword arguments to pass to a function (if cmd is function)

        Returns
        -------
        IPython.parallel.client.asyncresult.AsyncResult
        """
        with self._lview.temp_flags(after=deps, block=False):
            if isinstance(cmd, str):
                task = self._lview.apply_async(system_call, cmd)
            else:
                task = self._lview.apply_async(cmd, *args, **kwargs)

        return task

    def submit_async(self, cmd, *args, **kwargs):
        """Submit an async command to execute

        Parameters
        ----------
        cmd : {function, str}
            A function to execute or a system call to execute
        args : list
            Arguments to pass to a function (if cmd is function)
        kwargs : dict
            Keyword arguments to pass to a function (if cmd is function)

        Returns
        -------
        IPython.parallel.client.asyncresult.AsyncResult

        """
        if isinstance(cmd, str):
            task = self._lview.apply_async(system_call, cmd)
        else:
            task = self._lview.apply_async(cmd, *args, **kwargs)

        return task

    def submit_sync(self, cmd, *args, **kwargs):
        """Submit an sync command to execute

        Parameters
        ----------
        cmd : {function, str}
            A function to execute or a system call to execute
        args : list
            Arguments to pass to a function (if cmd is function)
        kwargs : dict
            Keyword arguments to pass to a function (if cmd is function)

        Returns
        -------
        Dependent on cmd

        """
        if isinstance(cmd, str):
            result = self._lview.apply_sync(system_call, cmd)
        else:
            result = self._lview.apply_sync(cmd, *args, **kwargs)

        return result

    def wait(self, handlers):
        """Waits until all async jobs in handlers have finished

        Parameters
        ----------
        handlers : list of AsyncResult
            The AsyncResult objects to wait for
        """
        return self._lview.wait(handlers)

    def get_number_of_workers(self):
        """Returns the number of workers present in the cluster

        Returns
        -------
        int
            The number of workers
        """
        return len(self._client.ids)


# likely want this in qiime.parallel.__init__
context = Context(qiime_config['profile'])
