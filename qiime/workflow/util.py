#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso",
               "Kyle Bittinger",
               "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import sys
from os.path import join
from datetime import datetime
from cogent.util.misc import safe_md5
from qiime.util import (qiime_system_call,
                        get_qiime_library_version)


def generate_log_fp(output_dir,
                    basefile_name='log',
                    suffix='txt',
                    timestamp_pattern='%Y%m%d%H%M%S'):
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s_%s.%s' % (basefile_name, timestamp, suffix)
    return join(output_dir, filename)


class WorkflowError(Exception):
    pass


class WorkflowLogger(object):

    def __init__(self, log_fp=None, params=None,
                 qiime_config=None, open_mode='w'):
        if log_fp:
            self._f = open(log_fp, open_mode)
        else:
            self._f = None
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('Logging started at %s\n' % start_time)
        self.write('QIIME version: %s\n\n' % get_qiime_library_version())
        self.writeQiimeConfig(qiime_config)
        self.writeParams(params)

    def write(self, s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of
            # data is being written to the log files.
            self._f.flush()
        else:
            pass

    def writeQiimeConfig(self, qiime_config):
        if qiime_config is None:
            self.write('No qiime config provided.\n')
        else:
            self.write('qiime_config values:\n')
            for k, v in qiime_config.items():
                if v:
                    self.write('%s\t%s\n' % (k, v))
            self.write('\n')

    def writeParams(self, params):
        if params is None:
            self.write('No params provided.\n')
        else:
            self.write('parameter file values:\n')
            for k, v in params.items():
                for inner_k, inner_v in v.items():
                    val = inner_v or 'True'
                    self.write('%s:%s\t%s\n' % (k, inner_k, val))
            self.write('\n')

    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass


def print_commands(commands,
                   status_update_callback,
                   logger,
                   close_logger_on_success=True):
    """Print list of commands to run """
    logger.write("Printing commands only.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('#%s' % e[0])
            print '%s' % e[1]
            logger.write('# %s command\n%s\n\n' % e)


def call_commands_serially(commands,
                           status_update_callback,
                           logger,
                           close_logger_on_success=True):
    """Run list of commands, one after another """
    logger.write("Executing commands.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            logger.write('# %s command \n%s\n\n' % e)
            stdout, stderr, return_value = qiime_system_call(e[1])
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                    "Command run was:\n %s\n" % e[1] +\
                    "Command returned exit status: %d\n" % return_value +\
                    "Stdout:\n%s\nStderr\n%s\n" % (stdout, stderr)
                logger.write(msg)
                logger.close()
                raise WorkflowError(msg)
            # in the no error case, we write commands' output to the log
            # and also echo to this proc's stdout/stderr
            else:
                # write stdout and stderr to log file
                logger.write("Stdout:\n%s\nStderr:\n%s\n" % (stdout, stderr))
                # write stdout to stdout
                if stdout:
                    print stdout
                # write stderr to stderr
                if stderr:
                    sys.stderr.write(stderr)
    if close_logger_on_success:
        logger.close()


def print_to_stdout(s):
    print s


def no_status_updates(s):
    pass


def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append('--%s' % (param_id))
        if param_value is not None:
            result.append(param_value)
    return ' '.join(result)


def validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser):
    if (jobs_to_start != int(default_jobs_to_start)) and \
       not parallel:
        option_parser.error("Passing -O requires that -a is also passed.")
    params['parallel']['jobs_to_start'] = str(jobs_to_start)


def log_input_md5s(logger, fps):
    logger.write("Input file md5 sums:\n")
    for fp in fps:
        if fp is not None:
            logger.write("%s: %s\n" % (fp, safe_md5(open(fp)).hexdigest()))
    logger.write("\n")
