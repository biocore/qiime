#!/usr/bin/env python
import click

from qiime.util import load_qiime_config
from qiime.parallel.manager import (start_cluster as _start_cluster,
                                    stop_cluster as _stop_cluster)


qiime_config = load_qiime_config()


@click.group()
def qiime_env():
    pass


@qiime_env.command()
@click.option('--profile', required=False, default=qiime_config['profile'],
              help="The IPython profile to use")
@click.option('-n', type=int, required=False, help="The number of engines",
              default=qiime_config['n_engines'])
def start(profile, n):
    """Start a compute environment"""
    # qiime_config is a defaultdict which results in None when used to 
    # indicate defaults in the above, so we cannot rely right now on click's
    # parameter validation
    if profile is None:
        raise ValueError("No profile specified!")
    if n is None:
        raise ValueError("Must specify the number of engines to use!")

    _start_cluster(profile, n)


@qiime_env.command()
@click.option('--profile', required=False, default=qiime_config['profile'],
              help="The IPython profile to use")
def stop(profile):
    """Stop a compute environment"""
    if profile is None:
        raise ValueError("No profile specified!")
    
    _stop_cluster(profile)

if __name__ == '__main__':
    qiime_env()
