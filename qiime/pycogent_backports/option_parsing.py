#!/usr/bin/env python
""" Utilities for handle script options
"""

from copy import copy
import types
import sys
from optparse import (OptionParser, OptionGroup, Option, 
                      OptionValueError, OptionError)
from os import popen, remove, makedirs, getenv
from os.path import join, abspath, exists, isdir, isfile


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Greg Caporaso","Daniel McDonald",
               "Gavin Huttley","Rob Knight",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Production"

## Definition of CogentOption option type, a subclass of Option that
## contains specific types for filepaths and directory paths. This 
## will be particularly useful for graphical interfaces that make 
## use of the script_info dictionary as they can then distinguish 
## paths from ordinary strings
def check_existing_filepath(option, opt, value):
    if not exists(value):
        raise OptionValueError(
            "option %s: file does not exist: %r" % (opt, value))
    elif not isfile(value):
        raise OptionValueError(
            "option %s: not a regular file (can't be a directory!): %r" % (opt, value))
    else:
        return value

def check_existing_filepaths(option, opt, value):
    values = value.split(',')
    for v in values:
        check_existing_filepath(option,opt,v)
    return values

def check_existing_dirpath(option, opt, value):
    if not exists(value):
        raise OptionValueError(
            "option %s: directory does not exist: %r" % (opt, value))
    elif not isdir(value):
        raise OptionValueError(
            "option %s: not a directory (can't be a file!): %r" % (opt, value))
    else:
        return value

def check_new_filepath(option, opt, value):
    return value
        
def check_new_dirpath(option, opt, value):
    return value
    
def check_existing_path(option, opt, value):
    if not exists(value):
        raise OptionValueError(
            "option %s: path does not exist: %r" % (opt, value))
    return value
    
def check_new_path(option, opt, value):
    return value

def check_multiple_choice(option, opt, value):
    #split_char = ';' if ';' in value else ','
    values = value.split(option.split_char)
    for v in values:
        if v not in option.mchoices:
            choices = ",".join(map(repr, option.mchoices))
            raise OptionValueError(
                "option %s: invalid choice: %r (choose from %s)"
                % (opt, v, choices))
    return values

class CogentOption(Option):
    ATTRS = Option.ATTRS + ['mchoices','split_char']

    TYPES = Option.TYPES + ("existing_path",
                            "new_path",
                            "existing_filepath",
                            "existing_filepaths",
                            "new_filepath",
                            "existing_dirpath",
                            "new_dirpath",
                            "multiple_choice")
    TYPE_CHECKER = copy(Option.TYPE_CHECKER)
    # for cases where the user specifies an existing file or directory
    # as input, but it can be either a dir or a file
    TYPE_CHECKER["existing_path"] = check_existing_path
    # for cases where the user specifies a new file or directory
    # as output, but it can be either a dir or a file
    TYPE_CHECKER["new_path"] = check_new_path
    # for cases where the user passes a single existing file
    TYPE_CHECKER["existing_filepath"] = check_existing_filepath
    # for cases where the user passes one or more existing files
    # as a comma-separated list - paths are returned as a list
    TYPE_CHECKER["existing_filepaths"] = check_existing_filepaths
    # for cases where the user is passing a new path to be 
    # create (e.g., an output file)
    TYPE_CHECKER["new_filepath"] = check_new_filepath
    # for cases where the user is passing an existing directory
    # (e.g., containing a set of input files)
    TYPE_CHECKER["existing_dirpath"] = check_existing_dirpath
    # for cases where the user is passing a new directory to be 
    # create (e.g., an output dir which will contain many result files)
    TYPE_CHECKER["new_dirpath"] = check_new_dirpath
    # for cases where the user is passing one or more values
    # as comma- or semicolon-separated list
    # choices are returned as a list
    TYPE_CHECKER["multiple_choice"] = check_multiple_choice

    def _check_multiple_choice(self):
        if self.type == "multiple_choice":
            if self.mchoices is None:
                raise OptionError(
                    "must supply a list of mchoices for type '%s'" % self.type, self)
            elif type(self.mchoices) not in (types.TupleType, types.ListType):
                raise OptionError(
                    "choices must be a list of strings ('%s' supplied)"
                    % str(type(self.mchoices)).split("'")[1], self)
            if self.split_char is None:
                self.split_char = ','
        elif self.mchoices is not None:
            raise OptionError(
                "must not supply mchoices for type %r" % self.type, self)

    CHECK_METHODS = Option.CHECK_METHODS + [_check_multiple_choice]

make_option = CogentOption

## End definition of new option type

def build_usage_lines(required_options,
    script_description,
    script_usage,
    optional_input_line,
    required_input_line):
    """ Build the usage string from components 
    """
    line1 = 'usage: %prog [options] ' + '{%s}' %\
     ' '.join(['%s %s' % (str(ro),ro.dest.upper())\
               for ro in required_options])
    usage_examples = []
    for title, description, command in script_usage:
        title = title.strip(':').strip()
        description = description.strip(':').strip()
        command = command.strip()
        if title:
            usage_examples.append('%s: %s\n %s' %\
             (title,description,command))
        else:
            usage_examples.append('%s\n %s' % (description,command))
    usage_examples = '\n\n'.join(usage_examples)
    lines = (line1,
             '', # Blank line
             optional_input_line,
             required_input_line,
             '', # Blank line
             script_description,
             '', # Blank line
             'Example usage: ',\
             'Print help message and exit',
             ' %prog -h\n',
             usage_examples)
    return '\n'.join(lines)

def set_parameter(key,kwargs,default=None):
    try:
        return kwargs[key]
    except KeyError:
        return default
        
def set_required_parameter(key,kwargs):
    try:
        return kwargs[key]
    except KeyError:
        raise KeyError,\
         "parse_command_line_parameters requires value for %s" % key
        
def parse_command_line_parameters(**kwargs):
    """ Constructs the OptionParser object and parses command line arguments
    
        parse_command_line_parameters takes a dict of objects via kwargs which
         it uses to build command line interfaces according to standards 
         developed in the Knight Lab, and enforced in QIIME. The currently 
         supported options are listed below with their default values. If no 
         default is provided, the option is required.
        
        script_description
        script_usage = [("","","")]
        version
        required_options=None
        optional_options=None
        suppress_verbose=False
        disallow_positional_arguments=True
        help_on_no_arguments=True
        optional_input_line = '[] indicates optional input (order unimportant)'
        required_input_line = '{} indicates required input (order unimportant)'
        
       These values can either be passed directly, as:
        parse_command_line_parameters(script_description="My script",\
                                     script_usage=[('Print help','%prog -h','')],\
                                     version=1.0)
                                     
       or they can be passed via a pre-constructed dict, as:
        d = {'script_description':"My script",\
             'script_usage':[('Print help','%prog -h','')],\
             'version':1.0}
        parse_command_line_parameters(**d)
    
    """
    # Get the options, or defaults if none were provided.
    script_description = set_required_parameter('script_description',kwargs)
    version = set_required_parameter('version',kwargs)
    script_usage = set_parameter('script_usage',kwargs,[("","","")])
    required_options = set_parameter('required_options',kwargs,[])
    optional_options = set_parameter('optional_options',kwargs,[])
    suppress_verbose = set_parameter('suppress_verbose',kwargs,False)
    disallow_positional_arguments =\
     set_parameter('disallow_positional_arguments',kwargs,True)
    help_on_no_arguments = set_parameter('help_on_no_arguments',kwargs,True)
    optional_input_line = set_parameter('optional_input_line',kwargs,\
        '[] indicates optional input (order unimportant)')
    required_input_line = set_parameter('required_input_line',kwargs,\
        '{} indicates required input (order unimportant)')
    # command_line_text will usually be nothing, but can be passed for
    # testing purposes
    command_line_args = set_parameter('command_line_args',kwargs,None)
    
    # Build the usage and version strings
    usage = build_usage_lines(required_options,script_description,script_usage,\
                              optional_input_line,required_input_line)
    version = 'Version: %prog ' + version
    
    # Instantiate the command line parser object
    parser = OptionParser(usage=usage, version=version)
    parser.exit = set_parameter('exit_func',kwargs,parser.exit)
    
    # If no arguments were provided, print the help string (unless the
    # caller specified not to)
    if help_on_no_arguments and (not command_line_args) and len(sys.argv) == 1:
        parser.print_usage()
        return parser.exit(-1)

    
    # Process the required options
    if required_options:
        # Define an option group so all required options are
        # grouped together, and under a common header
        required = OptionGroup(parser, "REQUIRED options",
         "The following options must be provided under all circumstances.")
        for ro in required_options:
            # if the option doesn't already end with [REQUIRED], 
            # add it.
            if not ro.help.strip().endswith('[REQUIRED]'):
                ro.help += ' [REQUIRED]'
            required.add_option(ro)
        parser.add_option_group(required)

    # Add a verbose parameter (if the caller didn't specify not to)
    if not suppress_verbose:
        parser.add_option('-v','--verbose',action='store_true',\
           dest='verbose',help='Print information during execution '+\
           '-- useful for debugging [default: %default]',default=False)

    # Add the optional options
    map(parser.add_option,optional_options)
    
    # Parse the command line
    # command_line_text will None except in test cases, in which 
    # case sys.argv[1:] will be parsed
    opts,args = parser.parse_args(command_line_args)
    
    # If positional arguments are not allowed, and any were provided,
    # raise an error.
    if disallow_positional_arguments and len(args) != 0:
        parser.error("Positional argument detected: %s\n" % str(args[0]) +\
         " Be sure all parameters are identified by their option name.\n" +\
         " (e.g.: include the '-i' in '-i INPUT_DIR')")

    # Test that all required options were provided.
    if required_options:
        required_option_ids = [o.dest for o in required.option_list]
        for required_option_id in required_option_ids:
            if getattr(opts,required_option_id) == None:
                return parser.error('Required option --%s omitted.' \
                             % required_option_id)
            
    # Return the parser, the options, and the arguments. The parser is returned
    # so users have access to any additional functionality they may want at 
    # this stage -- most commonly, it will be used for doing custom tests of 
    # parameter values.
    return parser, opts, args

