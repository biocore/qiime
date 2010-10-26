#!/usr/bin/env python
import sys
from optparse import OptionParser, OptionGroup
from os.path import exists, isdir
from os import makedirs
from cogent.util.misc import curry

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Amanda Birmingham",
                    "Sandra Smit", "Zongzhi Liu", "Daniel McDonald",
                    "Kyle Bittinger", "Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

## Begin functions for handling command line parameters

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
        script_usage
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
    script_usage = set_required_parameter('script_usage',kwargs)
    version = set_required_parameter('version',kwargs)
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
    
    # Build the usage and version strings
    usage = build_usage_lines(required_options,script_description,script_usage,\
                              optional_input_line,required_input_line)
    version = 'Version: %prog ' + version
    
    # Instantiate the command line parser object
    parser = OptionParser(usage=usage, version=version)
    
    # If no arguments were provided, print the help string (unless the
    # caller specified not to)
    if help_on_no_arguments and len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()
    
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
    opts,args = parser.parse_args()
    
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
                parser.error('Required option --%s omitted.' \
                             % required_option_id)
            
    # Return the parser, the options, and the arguments. The parser is returned
    # so users have access to any additional functionality they may want at 
    # this stage -- most commonly, it will be used for doing custom tests of 
    # parameter values.
    return parser, opts, args
    
## End functions for handling command line parameters


#some error codes for creating a dir
def get_create_dir_error_codes():
    return {'NO_ERROR':      0,
            'DIR_EXISTS':    1,
            'FILE_EXISTS':   2,
            'OTHER_OS_ERROR':3}

def create_dir(dir_name, fail_on_exist=False, handle_errors_externally=False):
    """Create a dir safely and fail meaningful.

    dir_name: name of directory to create

    fail_on_exist: if true raise an error if dir already exists
    
    handle_errors_externally: if True do not raise Errors, but return
                   failure codes. This allows to handle errors locally and
                   e.g. hint the user at a --force_overwrite options.
                   
    returns values (if no Error raised):
    
         0:  dir could be safely made
         1:  directory already existed
         2:  a file with the same name exists          
         3:  any other unspecified OSError


    See qiime/denoiser.py for an example of how to use this mechanism.

    Note: Depending  of how thorough we want to be we could add tests,
          e.g. for testing actual write permission in an existing dir.
    """
    error_code_lookup = get_create_dir_error_codes()
    #pre-instanciate function with
    ror = curry(handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            #dir is there
            if fail_on_exist:
                return ror(error_code_lookup['DIR_EXISTS'])
            else:
                return error_code_lookup['DIR_EXISTS']
        else:
            #must be file with same name
            return ror(error_code_lookup['FILE_EXISTS'])
    else:
        #no dir there, try making it
        try:
            makedirs(dir_name)
        except OSError:
            return ror(error_code_lookup['OTHER_OS_ERROR'])
    
    return error_code_lookup['NO_ERROR']

def handle_error_codes(dir_name, supress_errors=False,
                       error_code=None):
    """Wrapper function for error_handling.

    dir_name: name of directory that raised the error
    suppress_errors: if True raise Errors, otherwise return error_code
    error_code: the code for the error
    """
    error_code_lookup = get_create_dir_error_codes()
    
    if error_code == None:
        error_code = error_code_lookup['NO_ERROR']
    
    error_strings = \
        {error_code_lookup['DIR_EXISTS'] :
          "Directory already exists: %s" % dir_name,
         error_code_lookup['FILE_EXISTS'] : 
          "File with same name exists: %s" % dir_name,
         error_code_lookup['OTHER_OS_ERROR']: 
          "Could not create output directory: %s. " % dir_name +
          "Check the permissions."}

    if error_code == error_code_lookup['NO_ERROR']:
        return error_code_lookup['NO_ERROR']
    if supress_errors:
        return error_code
    else:
        raise OSError, error_strings[error_code]