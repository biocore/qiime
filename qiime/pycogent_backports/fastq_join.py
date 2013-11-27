#!/usr/bin/env python
# file: fastq_join.py

# Application controller for ea-utils v1.1.2-537 
# fastq processing utilities
# http://code.google.com/p/ea-utils/
# 

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError
import os 
import tempfile

__author__ = "Michael Robeson"
__copyright__ = "Copyright 2007-2013, The Cogent Project"
__credits__ = ["Michael Robeson"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Michael Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

class FastqJoin(CommandLineApplication):
    """fastq-join (v1.1.2) application controller for joining paired-end reads."""
    
    _command = 'fastq-join'
    
    _parameters = {
    # Description copied from 'fastq-join'
    # Usage: fastq-join [options] <read1.fq> <read2.fq> [mate.fq] -o <read.%.fq>
    
    # Output: 
    # You can supply 3 -o arguments, for un1, un2, join files, or one 
    # argument as a file name template.  The suffix 'un1, un2, or join' is 
    # appended to the file, or they replace a %-character if present.
    # If a 'mate' input file is present (barcode read), then the files
    # 'un3' and 'join2' are also created.
    
    # we'll only handle one output base path / file name
    # -o FIL:  See 'Output' above
    '-o':ValuedParameter(Prefix='-', Delimiter=' ', Name='o'),

    # -v C:  Verifies that the 2 files probe id's match up to char C
    # use ' ' (space) for Illumina reads
    '-v':ValuedParameter(Prefix='-', Delimiter=' ', Name='v'),

    # -p N:  N-percent maximum difference (8)
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p'),
    
    # -m N:  N-minimum overlap (6)
    '-m':ValuedParameter(Prefix='-', Delimiter=' ', Name='m'),
   
    # -r FIL:  Verbose stitch length report
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='r')}

    _input_handler = '_input_as_paths'

    def _get_output_path(self):
        """Checks if a base file label / path is set. Returns absolute path."""
        if self.Parameters['-o'].isOn():
            output_path = self._absolute(str(self.Parameters['-o'].Value))
        else:
            raise ValueError, "No output path specified."
        return output_path

    def _get_stitch_report_path(self):
        """Checks if stitch report label / path is set. Returns absolute path."""
        if self.Parameters['-r'].isOn():
            stitch_path = self._absolute(str(self.Parameters['-r'].Value))
            return stitch_path
        elif self.Parameters['-r'].isOff():
            return None

    def _get_result_paths(self, data):
        """Capture fastq-join output.
        
        Three output files are produced, in the form of
            outputjoin : assembled paired reads
            outputun1 : unassembled reads_1
            outputun2 : unassembled reads_2

        If a barcode / mate-pairs file is also provided then the following 
        additional files are output:
            outputjoin2
            outputun3

        If a verbose stitch length report (-r) is chosen to be written by the 
        user then use a user specified filename.
        """
        output_path = self._get_output_path()
        
        result = {}

        # always output:
        result['Assembled'] = ResultPath(Path = output_path + 'join',
                                         IsWritten=True)
        result['UnassembledReads1']  = ResultPath(Path = output_path + 'un1',
                                                  IsWritten=True)
        result['UnassembledReads2']  = ResultPath(Path = output_path + 'un2',
                                                  IsWritten=True)
       
        # check if stitch report is requested:
        stitch_path = self._get_stitch_report_path()
        if stitch_path:
            result['Report'] = ResultPath(Path = stitch_path,
			                              IsWritten=True)

        # Check if mate file / barcode file is present.
        # If not, return result
        # We need to check this way becuase there are no infile parameters.
        mate_path_string = output_path + 'join2'
        mate_unassembled_path_string = output_path + 'un3'
        if os.path.exists(mate_path_string) and \
            os.path.exists(mate_unassembled_path_string):
            result['Mate'] = ResultPath(Path = mate_path_string, 
                                        IsWritten=True)
            result['MateUnassembled'] = ResultPath(Path = 
                                                   mate_unassembled_path_string,
                                                   IsWritten=True)
        else:
            pass
        return result


    def getHelp(self):
        """fastq-join (v1.1.2) help"""
        help_str = """
        For issues with the actual program 'fastq-join', see the following:
    
        For basic help, type the following at the command line:
            'fastq-join'

        Website:
           http://code.google.com/p/ea-utils/

        For questions / comments subit an issue to:
        http://code.google.com/p/ea-utils/issues/list
        """
        return help_str


def join_paired_end_reads_fastqjoin(
    reads1_infile_path,
    reads2_infile_path,
    perc_max_diff=None, # typical default is 8
    min_overlap=None, # typical default is 6
    outfile_label = 'fastqjoin',
    params={},    
    working_dir=tempfile.gettempdir(),
    SuppressStderr=True,
    SuppressStdout=True,
    HALT_EXEC=False): 
    """ Runs fastq-join, with default parameters to assemble paired-end reads.
        Returns file path string.

        -reads1_infile_path : reads1.fastq infile path
        -reads2_infile_path : reads2.fastq infile path
        -perc_max_diff : maximum % diff of overlap differences allowed 
        -min_overlap : minimum allowed overlap required to assemble reads
        -outfile_label : base name for output files.
        -params : dictionary of application controller parameters

    """    
    abs_r1_path = os.path.abspath(reads1_infile_path)
    abs_r2_path = os.path.abspath(reads2_infile_path)
     
    infile_paths = [abs_r1_path, abs_r2_path]

    # check / make absolute infile paths
    for p in infile_paths:
        if not os.path.exists(p):
            raise IOError, 'File not found at: %s' % p
  
    fastq_join_app = FastqJoin(params=params,
                               WorkingDir=working_dir,
                               SuppressStderr=SuppressStderr,
                               SuppressStdout=SuppressStdout,
                               HALT_EXEC=HALT_EXEC)
  
    # set param. Helps with QIIME integration to have these values
    # set to None by default. This way we do not have to worry
    # about changes in default behaviour of the wrapped
    # application
    if perc_max_diff is not None:
        if isinstance(perc_max_diff, int) and 0 <= perc_max_diff <= 100: 
            fastq_join_app.Parameters['-p'].on(perc_max_diff)
        else:
            raise ValueError, "perc_max_diff must be int between 0-100!"

    if min_overlap is not None:
        if isinstance(min_overlap, int) and 0 < min_overlap: 
            fastq_join_app.Parameters['-m'].on(min_overlap)
        else:
            raise ValueError, "min_overlap must be an int >= 0!"

    if outfile_label is not None:
        if isinstance(outfile_label, str): 
            fastq_join_app.Parameters['-o'].on(outfile_label +'.')
        else:
            raise ValueError, "outfile_label must be a string!"
    else:
        pass

  
    # run assembler
    result = fastq_join_app(infile_paths)
    
    # Store output file path data to dict    
    path_dict = {}
    path_dict['Assembled'] = result['Assembled'].name
    path_dict['UnassembledReads1'] = result['UnassembledReads1'].name
    path_dict['UnassembledReads2'] = result['UnassembledReads2'].name
   
    # sanity check that files actually exist in path lcoations
    for path in path_dict.values():
        if not os.path.exists(path):
            raise IOError, 'Output file not found at: %s' % path

    return path_dict



