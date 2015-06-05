#!/usr/bin/env python

""" A simple client waiting for data to clean up 454 sequencing data"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option

from qiime.denoiser.utils import get_denoiser_data_dir
from qiime.denoiser.denoise_worker import setup_worker

options_lookup = get_options_lookup()

# denoiser_worker.py
script_info = {}
script_info['brief_description'] = "Start a denoiser worker process"
script_info['script_description'] = """The workers are automatically started by the denoiser.py script.
You usually never need to use this script yourself.

A worker waits for data and does flowgram alignments once it gets it."""

script_info['script_usage'] = [
    ("",
     "Start worker and connect to server listening on port 12345 on the same machine (localhost)",
     "%prog -f seqs.fna -f denoiser_out/worker99 -p 12345 -s localhost")
]

script_info[
    'output_description'] = "Denoise worker writes a log file if verbose flag is set."

script_info['required_options'] = [

    make_option('-f', '--file_path', action='store',
                type='string', dest='file_path',
                help='path used as prefix for worker data files' +
                '[REQUIRED]'),

    make_option('-p', '--port', action='store',
                type='int', dest='port', help='Server port ' +
                '[REQUIRED]'),

    make_option('-s', '--server_address', action='store',
                type='string', dest='server', help='Server address' +
                '[REQUIRED]')
]

script_info['optional_options'] = [

    make_option('-e', '--error_profile', action='store',
                type='string', dest='error_profile',
                help='Path to error profile' +
                ' [DEFAULT: %default]',
                default=get_denoiser_data_dir() +
                'FLX_error_profile.dat'),

    make_option('-c', '--counter', action='store',
                type='int', dest='counter',
                help='Round counter to start this worker with ' +
                ' [default: %default]', default=0)
]

script_info['version'] = __version__


def main(commandline_args=None):
    parser, opts, args = parse_command_line_parameters(**script_info)

    if not opts.file_path:
        parser.error("Required option --file_path not specified.")
    if not opts.port:
        parser.error("Required option --port not specified.")
    if not opts.server:
        parser.error("Required options --server_address not specified.")

    setup_worker(opts.file_path, opts.server, opts.port,
                 opts.counter, opts.verbose, opts.error_profile)

if __name__ == "__main__":
    main()
