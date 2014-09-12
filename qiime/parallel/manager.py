from IPython.parallel.apps.ipclusterapp import IPClusterStart, IPClusterStop


def start_cluster(profile, n):
    """Start a cluster

    Parameters
    ----------
    profile : str
        The name of the IPython profile to use
    n : int
        The number of engines to start

    Notes
    -----
    This method produces a daemon
    """
    c = IPClusterStart(profile=profile, log_level=0, daemonize=True)
    c.n = n
    c.initialize(argv=[])
    c.start()


def stop_cluster(profile):
    """Stop a cluster

    Parameters
    ----------
    profile : str
        The name of the IPython profile to use

    """
    c = IPClusterStop(profile=profile, log_level=0)
    c.initialize(argv=[])
    c.start()
