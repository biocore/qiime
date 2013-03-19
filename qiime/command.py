#!/usr/bin/env python
# File created on 30 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os import makedirs
from qiime.util import (make_option,
                        load_qiime_config,
                        parse_command_line_parameters)
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (log_input_md5s,
                            WorkflowLogger,
                            generate_log_fp)

qiime_config = load_qiime_config()

def cl_main(cmd, argv):
    script_info = cmd.getScriptInfo()
    
    option_parser, options, arguments =\
       parse_command_line_parameters(**script_info)
    
    try:
        cmd(options = eval(str(options)), 
            arguments = arguments,
            argv = argv)
    except QiimeCommandError, e:
        option_parser.error(e)

class QiimeCommandError(IOError):
    pass

class QiimeCommand(object):
    """ Base class for abstracted QIIME command
    """
    
    _standard_options = [
        make_option('--master_script_log_dir',type='existing_dirpath',
        help='directory where master script log will be stored [default: %default]',
        default='./')]
    _input_file_parameter_ids = []
    
    _brief_description = """ """
    _script_description = """ """
    _script_usage = []
    _script_usage_output_to_remove = []
    _output_description = """ """
    _required_options = []
    _optional_options = []
    _disallow_positional_arguments = True
    _version = __version__
    
    def __init__(self):
        """
        """
        self._optional_options.extend(self._standard_options)
        self._construct_option_lookup()
    
    def __call__(self,options,arguments,argv,logger=None):
        """
        """
        close_logger_on_success = self._start_logging(options,arguments,argv,logger)
        self.run_command(options,arguments)
        self._stop_logging(options,arguments,argv,close_logger_on_success)
    
    def _construct_option_lookup(self):
        """
        """
        _optional_options = {}
        _required_options = {}
        
        for e in self._optional_options:
            _optional_options[e._long_opts[0].lstrip('--')] = e
        for e in self._required_options:
            _required_options[e._long_opts[0].lstrip('--')] = e
        
        self._optional_options = _optional_options
        self._required_options = _required_options

    def _start_logging(self,
                       params,
                       args,
                       argv,
                       logger):
        if logger == None:
            self.logger = WorkflowLogger(generate_log_fp(params['master_script_log_dir']),
                                    params={},
                                    qiime_config=qiime_config)
            close_logger_on_success = True
        else:
            self.logger = logger
            close_logger_on_success = False
        
        self.logger.write('Command:\n')
        self.logger.write(' '.join(argv))
        self.logger.write('\n\n')
    
        log_input_md5s(self.logger,
                       [params[p] for p in self._input_file_parameter_ids])
        
        return close_logger_on_success

    def _stop_logging(self,
                      params,
                      args,
                      argv,
                      close_logger_on_success):
        if close_logger_on_success:
            self.logger.close()

    def getScriptInfo(self):
        result = {}
        result['brief_description'] = self._brief_description
        result['script_description'] = self._script_description
        result['script_usage'] = self._script_usage
        result['script_usage_output_to_remove'] = self._script_usage_output_to_remove
        result['output_description'] = self._output_description
        result['required_options'] = self._required_options.values()
        result['optional_options'] = self._optional_options.values()
        result['version'] = self._version
        result['disallow_positional_arguments'] = self._disallow_positional_arguments
        return result
    
    def run_command(self,params,args):
        raise NotImplementedError, "All subclasses must implement run_command."

class WorkflowCommand(QiimeCommand):

    def _validate_jobs_to_start(self,
                                jobs_to_start,
                                default_jobs_to_start,
                                parallel):
        if (int(jobs_to_start) != int(default_jobs_to_start)) and not parallel:
            raise QiimeCommandError, "Modifying jobs_to_start requires that parallel is True."
        return str(jobs_to_start)
