#!/usr/bin/env python

"""Application controller for BWA 0.6.2 (release 19 June 2012)"""

from cogent.app.parameters import FlagParameter, ValuedParameter, \
                                  MixedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, \
                            ApplicationError, get_tmp_filename
from os.path import isabs

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Adam Robbins-Pianka", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Production"

# helper functions for argument checking
def is_int(x):
    # return true if it's an int
    return ((type(x) == int) or \
    # or it's a string that is all digits
    (type(x) == str and x.isdigit()) or \
    # otherwise return False
    False)

def is_float(x):
    return (is_int(x) or \
    # or if it's a float
    (type(x) == float) or \
    # or it's a string with exactly one decimal and all digits on both sides of
    # the decimal
    (type(x) == str and '.' in x and all(map(str.isdigit, x.split('.', 1)))) \
    # otherwise return False
    or False)

#Base class
class BWA(CommandLineApplication):
    """BWA generic application controller. Do not instantiate directly.
    
    Instead of instantiating this class, instantiate a subclass for each
    subcommand.  Available subclasses are:
    BWA_index
    BWA_aln
    BWA_samse
    BWA_sampe
    BWA_bwasw
    """

    # all subclasses will accept dictionaries as input that specify input
    # and output files. The required (and optional) types of input and output
    # files differ by subcommand.
    _input_handler = "_input_as_dict"

    # the main command. The program bwa should be in the PATH
    _command = "bwa"

    # holds the values of the dict handled by the input handler
    _input = {}

    # Each subclass can have a dictionary (keys = option names, e.g., -a
    # and values = boolean fucntions) called _valid_arguments
    # that specifies checks to be made on the parameters.
    def check_arguments(self):
        """Sanity check the arguments passed in.

        Uses the boolean functions specified in the subclasses in the
        _valid_arguments dictionary to determine if an argument is valid
        or invalid.
        """
        for k, v in self.Parameters.iteritems():
            if self.Parameters[k].isOn():
                if k in self._valid_arguments:
                    if not self._valid_arguments[k](v.Value):
                        error_message = 'Invalid argument (%s) ' % v.Value
                        error_message += 'for parameter %s\n' % k
                        raise ApplicationError(error_message)

    def _get_base_command(self):
        """ Returns the full command string 

        Overridden here because there are positional arguments (specifically
        the input and output files).
        """
        command_parts = []
        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])
        if self._command is None:
            raise ApplicationError, '_command has not been set.'
        command = self._command
        # also make sure there's a subcommand!
        if self._subcommand is None:
            raise ApplicationError, '_subcommand has not been set.'
        subcommand = self._subcommand
        # sorting makes testing easier, since the options will be written out
        # in alphabetical order. Could of course use option parsing scripts
        # in cogent for this, but this works as well.
        parameters = sorted([str(x) for x in self.Parameters.values() 
                            if str(x)])
        synonyms = self._synonyms
        
        command_parts.append(cd_command)
        command_parts.append(command)
        # add in subcommand
        command_parts.append(subcommand)
        command_parts += parameters
        # add in the positional arguments in the correct order
        for k in self._input_order:
            # this check is necessary to account for optional positional
            # arguments, such as the mate file for bwa bwasw
            # Note that the input handler will ensure that all required
            # parameters have valid values
            if k in self._input:
                command_parts.append(self._input[k])
      
        return self._command_delimiter.join(command_parts).strip()
    
    BaseCommand = property(_get_base_command)

    def _input_as_dict(self, data):
        """Takes dictionary that sets input and output files.

        Valid keys for the dictionary are specified in the subclasses. File
        paths must be absolute.
        """
        # clear self._input; ready to receive new input and output files
        self._input = {}
        # Check that the arguments to the
        # subcommand-specific parameters are valid
        self.check_arguments()

        # Ensure that we have all required input (file I/O)
        for k in self._input_order:
            # N.B.: optional positional arguments begin with underscore (_)!
            # (e.g., see _mate_in for bwa bwasw)
            if k[0] != '_' and k not in data:
                raise ApplicationError, "Missing required input %s" % k

        # Set values for input and output files
        for k in data:
            # check for unexpected keys in the dict
            if k not in self._input_order:
                error_message = "Invalid input arguments (%s)\n" % k
                error_message += "Valid keys are: %s" % repr(self._input_order)
                raise ApplicationError(error_message + '\n')

            # check for absolute paths
            if not isabs(data[k][0]):
                raise ApplicationError, "Only absolute paths allowed.\n%s" %\
                repr(data)
            self._input[k] = data[k]
        
        # if there is a -f option to specify an output file, force the user to
        # use it (otherwise things to to stdout)
        if '-f' in self.Parameters and not self.Parameters['-f'].isOn():
            raise ApplicationError, "Please specify an output file with -f"

        return ''

class BWA_index(BWA):
    """Controls the "index" subcommand of the bwa application.
    
    Valid input keys are: fasta_in
    """

    # the subcommand for bwa index
    _subcommand = "index"

    _parameters = {
        # which algorithm to use.
        # is
        # IS linear-time algorithm for constructing suffix array. It requires
        # 5.37N memory where N is the size of the database. IS is moderately
        # fast, but does not work with database larger than 2GB. IS is the
        # default algorithm due to its simplicity. The current codes for IS
        # algorithm are reimplemented by Yuta Mori.
        #
        # bwtsw
        # Algorithm implemented in BWT-SW. This method works with the whole
        # human genome, but it does not work with database smaller than 10MB
        # and it is usually slower than IS.
        #
        # DEFAULTs to auto-select (based on input fasta file size)
        '-a':ValuedParameter('-', Delimiter=' ', Name='a'),

        # prefix for the output index.
        # DEFAULTs to the base name of the input fasta file
        '-p':ValuedParameter('-', Delimiter=' ', Name='p'),

        # index files named as <in.fasta>.64.* instead of <in.fasta>.*
        '-6':FlagParameter('-', Name='6')
    }


    # The -a command can take on of only two possible values
    # the -p command allows the user to specify a prefix; for our purposes,
    # this prefix should be an abolute path
    _valid_arguments = {
        '-a': lambda x: x in ['is', 'bwtsw'],
        '-p': isabs
    }

    # For the position specific arguments, this is the order that they will
    # be written in the base command
    # input file keys beginning with _ are optional inputs
    _input_order = ['fasta_in']

    def _get_result_paths(self, data):
        """Gets the results for a run of bwa index.

        bwa index outputs 5 files when the index is created. The filename
        prefix will be the same as the input fasta, unless overridden with
        the -p option, and the 5 extensions are listed below:

        .amb
        .ann
        .bwt
        .pac
        .sa

        and these extentions (including the period) are the keys to the
        dictionary that is returned.
        """

        # determine the names of the files. The name will be the same as the
        # input fasta file unless overridden with the -p option
        if self.Parameters['-p'].isOn():
            prefix = self.Parameters['-p'].Value
        else:
            prefix = data['fasta_in']

        # the 5 output file suffixes
        suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        out_files = {}
        for suffix in suffixes:
            out_files[suffix] = ResultPath(prefix+suffix, IsWritten=True)

        return out_files

class BWA_aln(BWA):
    """Controls the "aln" subcommand of the bwa application.
    
    Valid input keys are: prefix, fastq_in 
    """
    _parameters = {
        # max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
        '-n': ValuedParameter('-', Delimiter=' ', Name='n'),
		#maximum number or fraction of gap opens [1]
        '-o': ValuedParameter('-', Delimiter=' ', Name='o'),

		#maximum number of gap extensions, -1 for disabling long gaps [-1]
        '-e': ValuedParameter('-', Delimiter=' ', Name='e'),

		#do not put an indel within bp towards the ends [5]
        '-i': ValuedParameter('-', Delimiter=' ', Name='i'),

		#maximum occurrences for extending a long deletion [10]
        '-d': ValuedParameter('-', Delimiter=' ', Name='d'),

		#seed length [32]
        '-l': ValuedParameter('-', Delimiter=' ', Name='l'),

		#maximum differences in the seed [2]
        '-k': ValuedParameter('-', Delimiter=' ', Name='k'),

		#maximum entries in the queue [2000000]
        '-m': ValuedParameter('-', Delimiter=' ', Name='m'),

		#number of threads [1]
        '-t': ValuedParameter('-', Delimiter=' ', Name='t'),

		#mismatch penalty [3]
        '-M': ValuedParameter('-', Delimiter=' ', Name='M'),

		#gap open penalty [11]
        '-O': ValuedParameter('-', Delimiter=' ', Name='O'),

		#gap extension penalty [4]
        '-E': ValuedParameter('-', Delimiter=' ', Name='E'),

		#stop searching when there are > equally best hits [30]
        '-R': ValuedParameter('-', Delimiter=' ', Name='R'),

		#quality threshold for read trimming down to 35bp [0]
        '-q': ValuedParameter('-', Delimiter=' ', Name='q'),

		#file to write output to instead of stdout
        '-f': ValuedParameter('-', Delimiter=' ', Name='f'),

		#length of barcode
        '-B': ValuedParameter('-', Delimiter=' ', Name='B'),

		#log-scaled gap penalty for long deletions
        '-L': FlagParameter('-', Name='L'),

		#non-iterative mode: search for all n-difference hits (slooow)
        '-N': FlagParameter('-', Name='N'),

		#the input is in the Illumina 1.3+ FASTQ-like format
        '-I': FlagParameter('-', Name='I'),

		#the input read file is in the BAM format
        '-b': FlagParameter('-', Name='b'),

		#use single-end reads only (effective with -b)
        '-0': FlagParameter('-', Name='0'),

		#use the 1st read in a pair (effective with -b)
        '-1': FlagParameter('-', Name='1'),

		#use the 2nd read in a pair (effective with -b)
        '-2': FlagParameter('-', Name='2'),

		#filter Casava-filtered sequences
        '-Y': FlagParameter('-', Name='Y')
    }

    # the subcommand for bwa aln
    _subcommand = 'aln'

    _valid_arguments = {
        # check to see if this is decimal numbers
        '-n': is_float,

        # check to see if these are integers
        '-o': is_int,
        '-e': is_int,
        '-i': is_int,
        '-d': is_int,
        '-l': is_int,
        '-k': is_int,
        '-m': is_int,
        '-t': is_int,
        '-M': is_int,
        '-O': is_int,
        '-E': is_int,
        '-R': is_int,
        '-q': is_int,
        '-B': is_int,

        # check to see if this is an absolute file path
        '-f': isabs
    }

    # input file keys beginning with _ are optional inputs
    _input_order = ['prefix', 'fastq_in']

    def _get_result_paths(self, data):
        """Gets the result file for a bwa aln run.

        There is only one output file of a bwa aln run, a .sai file
        and it can be retrieved with the key 'output'.
        """
        return {'output': ResultPath(self.Parameters['-f'].Value,
                                    IsWritten=True)}

class BWA_samse(BWA):
    """Controls the "samse" subcommand of the bwa application.
    
    Valid input keys are: prefix, sai_in, fastq_in
    """
    _parameters = {
        # Maximum number of alignments to output in the XA tag for reads
        # paired properly. If a read has more than this number of hits, the
        # XA tag will not be written
        '-n': ValuedParameter('-', Delimiter=' ', Name='n'),

		#file to write output to instead of stdout
        '-f': ValuedParameter('-', Delimiter=' ', Name='f'),

        # Specify the read group in a format like '@RG\tID:foo\tSM:bar'
        '-r': ValuedParameter('-', Delimiter=' ', Name='r')
    }

    # the subcommand for samse
    _subcommand = 'samse'

    _valid_arguments = {
        # make sure that this is an int
        '-n': is_int,

        # check to see if this is an absolute file path
        '-f': isabs
    }

    # input file keys beginning with _ are optional inputs
    _input_order = ['prefix', 'sai_in', 'fastq_in']

    def _get_result_paths(self, data):
        """Gets the result file for a bwa samse run.

        There is only one output file of a bwa samse run, a .sam file
        and it can be retrieved with the key 'output'.
        """
        return {'output': ResultPath(self.Parameters['-f'].Value, 
                                    IsWritten=True)}

class BWA_sampe(BWA):
    """Controls the "sampe" subcommand of the bwa application.
    
    Valid input keys are: prefix, sai1_in, sai2_in, fastq1_in,
    fastq2_in
    """
    _parameters = {
        # Maximum insert size for a read pair to be considered being mapped
        # properly
        '-a': ValuedParameter('-', Delimiter=' ', Name='a'),

        # Maximum occurrences of a read for pairing
        '-o': ValuedParameter('-', Delimiter=' ', Name='o'),

        # Load the entire FM-index into memory to reduce disk operations
        '-P': FlagParameter('-', Name='P'),

        # maximum hits to output for paired reads [3]
        '-n': ValuedParameter('-', Delimiter=' ', Name='n'),

        # maximum hits to output for discordant pairs [10]
        '-N': ValuedParameter('-', Delimiter=' ', Name='N'),

		#file to write output to instead of stdout
        '-f': ValuedParameter('-', Delimiter=' ', Name='f'),

        # Specify the read group in a format like '@RG\tID:foo\tSM:bar'
        '-r': ValuedParameter('-', Delimiter=' ', Name='r'),

        # disable Smith-Waterman for the unmapped mate
        '-s': FlagParameter('-', Name='s'),

        # prior of chimeric rate (lower bound) [1.0e-05]
        '-c': ValuedParameter('-', Delimiter= ' ', Name='c'),

        # disable insert size estimate (force -s)
        '-A': FlagParameter('-', Name='A')
    }

    # the subcommand for sampe
    _subcommand = 'sampe'

    _valid_arguments = {
        # make sure this is a float
        '-c': is_float,

        # make sure these are all ints
        '-a': is_int,
        '-o': is_int,
        '-n': is_int,
        '-N': is_int,

        # check to see if this is an absolute file path
        '-f': isabs
    }

    # input file keys beginning with _ are optional inputs
    _input_order = ['prefix', 'sai1_in', 'sai2_in',
    'fastq1_in', 'fastq2_in']

    def _get_result_paths(self, data):
        """Gets the result file for a bwa sampe run.

        There is only one output file of a bwa sampe run, a .sam file,
        and it can be retrieved with the key 'output'.
        """
        return {'output': ResultPath(self.Parameters['-f'].Value, 
                                    IsWritten=True)}

class BWA_bwasw(BWA):
    """Controls the "bwasw" subcommand of the bwa application.
    
    Valid input keys are: prefix, query_fasta, _query_fasta2
    input keys beginning with an underscore are optional.
    """
    _parameters = {
        #Score of a match [1]
        '-a': ValuedParameter('-', Delimiter=' ', Name='a'),
        
        #Mismatch penalty [3]
        '-b': ValuedParameter('-', Delimiter=' ', Name='b'),
        
        #Gap open penalty [5]
        '-q': ValuedParameter('-', Delimiter=' ', Name='q'),
        
        #Gap  extension  penalty.
        '-r': ValuedParameter('-', Delimiter=' ', Name='r'),

        # mask level [0.50]
        '-m': ValuedParameter('-', Delimiter=' ', Name='m'),
        
        #Number of threads in the multi-threading mode [1]
        '-t': ValuedParameter('-', Delimiter=' ', Name='t'),

        # file to output results to instead of stdout
        '-f': ValuedParameter('-', Delimiter=' ', Name='f'),
        
        #Band width in the banded alignment [33]
        '-w': ValuedParameter('-', Delimiter=' ', Name='w'),
        
        #Minimum score threshold divided by a [30]
        '-T': ValuedParameter('-', Delimiter=' ', Name='T'),
        
        #Coefficient  for  threshold  adjustment  according  to query length.
        #Given an l-long query, the threshold for a hit to be retained is
        #a*max{T,c*log(l)}. [5.5]
        '-c': ValuedParameter('-', Delimiter=' ', Name='c'),
        
        #Z-best heuristics. Higher -z increases accuracy at the cost
        #of speed. [1]
        '-z': ValuedParameter('-', Delimiter=' ', Name='z'),
        
        #Maximum SA interval size for initiating a seed. Higher -s increases
        #accuracy at the cost of speed. [3]
        '-s': ValuedParameter('-', Delimiter=' ', Name='s'),
        
        #Minimum  number  of  seeds  supporting  the  resultant alignment to
        #trigger reverse alignment. [5]
        '-N': ValuedParameter('-', Delimiter=' ', Name='N'),

        # in SAM output, use hard clipping instead of soft clipping
        '-H': FlagParameter('-', Name='H'),

        # mark multi-part alignments as secondary
        '-M': FlagParameter('-', Name='M'),

        # skip Smith-Waterman read pariing
        '-S': FlagParameter('-', Name='S'),

        # ignore pairs with insert >= INT for inferring the size of distr
        # [20000]
        '-I': ValuedParameter('-', Delimiter=' ', Name='I')
    }

    # the subcommand fo bwasw
    _subcommand = 'bwasw'

    # input file keys beginning with _ are optional inputs
    _input_order = ['prefix', 'query_fasta', '_query_fasta_2']

    _valid_arguments = {
        # Make sure this is a float
        '-c': is_float,
        '-m': is_float,

        # Make sure these are ints
        '-a': is_int,
        '-b': is_int,
        '-q': is_int,
        '-r': is_int,
        '-t': is_int,
        '-w': is_int,
        '-T': is_int,
        '-z': is_int,
        '-s': is_int,
        '-N': is_int,
        '-I': is_int,

        # make sure this is an absolute path
        '-f': isabs
    }

    def _get_result_paths(self, data):
        """Gets the result file for a bwa bwasw run.

        There is only one output file of a bwa bwasw run, a .sam file,
        and it can be retrieved with the key 'output'.
        """
        return {'output': ResultPath(self.Parameters['-f'].Value, 
                                    IsWritten=True)}

def create_bwa_index_from_fasta_file(fasta_in, params=None):
    """Create a BWA index from an input fasta file.

    fasta_in: the input fasta file from which to create the index
    params: dict of bwa index specific paramters

    This method returns a dictionary where the keys are the various
    output suffixes (.amb, .ann, .bwt, .pac, .sa) and the values
    are open file objects.

    The index prefix will be the same as fasta_in, unless the -p parameter
    is passed in params.
    """
    if params is None:
        params = {}

    # Instantiate the app controller
    index = BWA_index(params)

    # call the application, passing the fasta file in
    results = index({'fasta_in':fasta_in})
    return results

def assign_reads_to_database(query, database_fasta, out_path, params=None):
    """Assign a set of query sequences to a reference database
    
    database_fasta_fp: absolute file path to the reference database
    query_fasta_fp: absolute file path to query sequences
    output_fp: absolute file path of the file to be output
    params: dict of BWA specific parameters.
            * Specify which algorithm to use (bwa-short or bwasw) using the
            dict key "algorithm"
            * if algorithm is bwasw, specify params for the bwa bwasw
            subcommand
            * if algorithm is bwa-short, specify params for the bwa samse
            subcommand
            * if algorithm is bwa-short, must also specify params to use with
            bwa aln, which is used to get the sai file necessary to run samse.
            bwa aln params should be passed in using dict key "aln_params" and
            the associated value should be a dict of params for the bwa aln
            subcommand
            * if a temporary directory is not specified in params using dict
            key "temp_dir", it will be assumed to be /tmp
    
    This method returns an open file object (SAM format).
    """
    if params is None:
        params = {}

    # set the output path
    params['-f'] = out_path

    # if the algorithm is not specified in the params dict, or the algorithm
    # is not recognized, raise an exception
    if 'algorithm' not in params:
        raise ApplicationError("Must specify which algorithm to use " + \
                               "('bwa-short' or 'bwasw')")
    elif params['algorithm'] not in ('bwa-short', 'bwasw'):
        raise ApplicationError('Unknown algorithm "%s". ' % \
                                params['algorithm'] + \
                                "Please enter either 'bwa-short' or 'bwasw'.")

    # if the temp directory is not specified, assume /tmp
    if 'temp_dir' not in params:
        params['temp_dir'] = '/tmp'

    # if the algorithm is bwa-short, we must build use bwa aln to get an sai
    # file before calling bwa samse on that sai file, so we need to know how
    # to run bwa aln. Therefore, we must ensure there's an entry containing
    # those parameters
    if params['algorithm'] == 'bwa-short':
        if 'aln_params' not in params:
            raise ApplicationError("With bwa-short, need to specify a key " + \
                                   "'aln_params' and its value, a " + \
                                   "dictionary to pass to bwa aln, since " + \
                                   "bwa aln is an intermediate step when " + \
                                   "doing bwa-short.")

    # we have this params dict, with "algorithm" and "temp_dir", etc which are
    # not for any of the subcommands, so make a new params dict that is the
    # same as the original minus these addendums
    subcommand_params = {}
    for k, v in params.iteritems():
        if k not in ('algorithm', 'temp_dir', 'aln_params'):
            subcommand_params[k] = v

    # build index from database_fasta
    # get a temporary file name that is not in use
    index_prefix = get_tmp_filename(tmp_dir=params['temp_dir'], suffix='', \
                                    result_constructor=str)

    create_bwa_index_from_fasta_file(database_fasta, {'-p': index_prefix})

    # if the algorithm is bwasw, things are pretty simple. Just instantiate
    # the proper controller and set the files
    if params['algorithm'] == 'bwasw':
        bwa = BWA_bwasw(params = subcommand_params)
        files = {'prefix': index_prefix, 'query_fasta': query}

    # if the algorithm is bwa-short, it's not so simple
    elif params['algorithm'] == 'bwa-short':
        # we have to call bwa_aln to get the sai file needed for samse
        # use the aln_params we ensured we had above
        bwa_aln = BWA_aln(params = params['aln_params'])
        aln_files = {'prefix': index_prefix, 'fastq_in': query}
        # get the path to the sai file
        sai_file_path = bwa_aln(aln_files)['output'].name

        # we will use that sai file to run samse
        bwa = BWA_samse(params = subcommand_params)
        files = {'prefix': index_prefix, 'sai_in': sai_file_path,
                 'fastq_in': query}

    # run which ever app controller we decided was correct on the files
    # we set up
    result = bwa(files)

    # they both return a SAM file, so return that
    return result['output']

def assign_dna_reads_to_dna_database(query_fasta_fp, database_fasta_fp, out_fp,
                                    params = {}):
    """Wraps assign_reads_to_database, setting various parameters.

    The default settings are below, but may be overwritten and/or added to
    using the params dict:

    algorithm:      bwasw
    """
    my_params = {'algorithm': 'bwasw'}
    my_params.update(params)

    result = assign_reads_to_database(query_fasta_fp, database_fasta_fp,
                                        out_fp, my_params)
    
    return result

def assign_dna_reads_to_protein_database(query_fasta_fp, database_fasta_fp, 
                                        out_fp, temp_dir='/tmp', 
                                        params = {}):
    """Wraps assign_reads_to_database, setting various parameters.

    Not yet implemented, as BWA can only align DNA reads to DNA databases.
    """
    raise NotImplementedError, "BWA cannot at this point align DNA to protein"
