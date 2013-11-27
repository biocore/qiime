#!/usr/bin/env python
# file: seqprep.py

# Application controller for SeqPrep 
# https://github.com/jstjohn/SeqPrep 
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

# SeqPrep help:
# Usage:
# SeqPrep [Required Args] [Options]
# NOTE 1: The output is always gziped compressed.
# NOTE 2: If the quality strings in the output contain characters less than 
# ascii 33 on an ascii table (they look like lines from a binary file), try 
# running again with or without the -6 option.
#

class SeqPrep(CommandLineApplication):
    """SeqPrep application controller for joining paired-end reads"""
    _command = 'SeqPrep'
    _parameters = {
    # Required Arguments
    # -f <first read input fastq filename>
    # -r <second read input fastq filename>
    # -1 <first read output fastq filename>
    # -2 <second read output fastq filename>
    '-f':ValuedParameter(Prefix='-', Delimiter=' ', Name='f'),
    '-r':ValuedParameter(Prefix='-', Delimiter=' ', Name='r'),
    '-1':ValuedParameter(Prefix='-', Delimiter=' ', Name='1'),
    '-2':ValuedParameter(Prefix='-', Delimiter=' ', Name='2'),
    
    # General Arguments (Optional):
    # -3 <first read discarded fastq filename>
    # -4 <second read discarded fastq filename>
    # -h Display this help message and exit (also works with no args) 
    # -6 Input sequence is in phred+64 rather than phred+33 format, the 
    #    output will still be phred+33 
    # -q <Quality score cutoff for mismatches to be counted in overlap; default = 13>
    # -L <Minimum length of a trimmed or merged read to print it; default = 30>
    '-3':ValuedParameter(Prefix='-', Delimiter=' ', Name='3'),
    '-4':ValuedParameter(Prefix='-', Delimiter=' ', Name='4'),
    '-h':FlagParameter(Prefix='-', Name='h'),
    '-6':FlagParameter(Prefix='-', Name='6'),
    '-q':ValuedParameter(Prefix='-', Delimiter=' ', Name='q'),
    '-L':ValuedParameter(Prefix='-', Delimiter=' ', Name='L'),
  
    # Arguments for Adapter/Primer Trimming (Optional):
    # -A <forward read primer/adapter sequence to trim as it would appear at the 
    #   end of a read (recommend about 20bp of this)
    #	(should validate by grepping a file); 
    #   default (genomic non-multiplexed adapter1) = AGATCGGAAGAGCGGTTCAG>
    # -B <reverse read primer/adapter sequence to trim as it would appear at the 
    #   end of a read (recommend about 20bp of this)
    #	(should validate by grepping a file); 
    #   default (genomic non-multiplexed adapter2) = AGATCGGAAGAGCGTCGTGT>
    # -O <minimum overall base pair overlap with adapter sequence to trim; 
    #   default = 10>
    # -M <maximum fraction of good quality mismatching bases for primer/adapter
    #    overlap; default = 0.020000>
    # -N <minimum fraction of matching bases for primer/adapter overlap; 
    #   default = 0.870000>
    # -b <adapter alignment band-width; default = 50>
    # -Q <adapter alignment gap-open; default = 8>
    # -t <adapter alignment gap-extension; default = 2>
    # -e <adapter alignment gap-end; default = 2>
    # -Z <adapter alignment minimum local alignment score cutoff 
    #   [roughly (2*num_hits) - (num_gaps*gap_open) - (num_gaps*gap_close) - 
    #   (gap_len*gap_extend) - (2*num_mismatches)]; default = 26>
    # -w <read alignment band-width; default = 50>
    # -W <read alignment gap-open; default = 26>
    # -p <read alignment gap-extension; default = 9>
    # -P <read alignment gap-end; default = 5>
    # -X <read alignment maximum fraction gap cutoff; default = 0.125000>
    '-A':ValuedParameter(Prefix='-', Delimiter=' ', Name='A'),
    '-B':ValuedParameter(Prefix='-', Delimiter=' ', Name='B'),
    '-O':ValuedParameter(Prefix='-', Delimiter=' ', Name='O'),
    '-M':ValuedParameter(Prefix='-', Delimiter=' ', Name='M'),
    '-N':ValuedParameter(Prefix='-', Delimiter=' ', Name='N'),
    '-b':ValuedParameter(Prefix='-', Delimiter=' ', Name='b'),
    '-Q':ValuedParameter(Prefix='-', Delimiter=' ', Name='Q'),
    '-t':ValuedParameter(Prefix='-', Delimiter=' ', Name='t'),
    '-e':ValuedParameter(Prefix='-', Delimiter=' ', Name='e'),
    '-Z':ValuedParameter(Prefix='-', Delimiter=' ', Name='Z'),
    '-w':ValuedParameter(Prefix='-', Delimiter=' ', Name='w'),
    '-W':ValuedParameter(Prefix='-', Delimiter=' ', Name='W'),
    '-p':ValuedParameter(Prefix='-', Delimiter=' ', Name='p'),
    '-P':ValuedParameter(Prefix='-', Delimiter=' ', Name='P'),
    '-X':ValuedParameter(Prefix='-', Delimiter=' ', Name='X'),

    # Optional Arguments for Merging:
    # -y <maximum quality score in output ((phred 33) default = ']' )>
    # -g <print overhang when adapters are present and stripped (use this if 
    #   reads are different length)>
    # -s <perform merging and output the merged reads to this file>
    # -E <write pretty alignments to this file for visual Examination>
    # -x <max number of pretty alignments to write (if -E provided);
    #   default = 10000>
    # -o <minimum overall base pair overlap to merge two reads; default = 15>
    # -m <maximum fraction of good quality mismatching bases to overlap reads;
    #   default = 0.020000>
    # -n <minimum fraction of matching bases to overlap reads;
    #   default = 0.900000>
    '-y':ValuedParameter(Prefix='-', Delimiter=' ', Name='y'),
    '-g':FlagParameter(Prefix='-', Name='y'),
    '-s':ValuedParameter(Prefix='-', Delimiter=' ', Name='s'),
    '-E':ValuedParameter(Prefix='-', Delimiter=' ', Name='E'),
    '-x':ValuedParameter(Prefix='-', Delimiter=' ', Name='x'),
    '-o':ValuedParameter(Prefix='-', Delimiter=' ', Name='o'),
    '-m':ValuedParameter(Prefix='-', Delimiter=' ', Name='m'),
    '-n':ValuedParameter(Prefix='-', Delimiter=' ', Name='n')}

    def _unassembled_reads1_out_file_name(self):
        """Checks file name is set for reads1 output. 
           Returns absolute path."""
        if self.Parameters['-1'].isOn():
            unassembled_reads1 = self._absolute(str(self.Parameters['-1'].Value))
        else:
            raise ValueError, "No reads1 (flag: -1) output path specified"
        return unassembled_reads1

    def _unassembled_reads2_out_file_name(self):
        """Checks if file name is set for reads2 output. 
           Returns absolute path."""
        if self.Parameters['-2'].isOn():
            unassembled_reads2 = self._absolute(str(self.Parameters['-2'].Value))
        else:
            raise ValueError, "No reads2 (flag -2) output path specified"
        return unassembled_reads2

    def _discarded_reads1_out_file_name(self):
        """Checks if file name is set for discarded reads1 output. 
           Returns absolute path."""
        if self.Parameters['-3'].isOn():
            discarded_reads1 = self._absolute(str(self.Parameters['-3'].Value))
        else:
            raise ValueError, "No discarded-reads1 (flag -3) output path specified"
        return discarded_reads1

    def _discarded_reads2_out_file_name(self):
        """Checks if file name is set for discarded reads2 output. 
           Returns absolute path."""
        if self.Parameters['-4'].isOn():
            discarded_reads2 = self._absolute(str(self.Parameters['-4'].Value))
        else:
            raise ValueError, "No discarded-reads2 (flag -4) output path specified"
        return discarded_reads2

    def _assembled_out_file_name(self):
        """Checks file name is set for assembled output. 
           Returns absolute path."""
        if self.Parameters['-s'].isOn():
            assembled_reads = self._absolute(str(self.Parameters['-s'].Value))
        else:
            raise ValueError, "No assembled-reads (flag -s) output path specified"
        return assembled_reads

    def _pretty_alignment_out_file_name(self):
        """Checks file name is set for pretty alignment output. 
           Returns absolute path."""
        if self.Parameters['-E'].isOn():
            pretty_alignment = self._absolute(str(self.Parameters['-E'].Value))
        else:
            raise ValueError, "No pretty-=alignment (flag -E) output path specified"
        return pretty_alignment

    def _get_result_paths(self, data):
        """Captures SeqPrep output.
        
        """
        result = {}
        
        # Always output:
        result['UnassembledReads1'] = ResultPath(Path = 
                                                 self._unassembled_reads1_out_file_name(), 
                                                 IsWritten=True)
        result['UnassembledReads2'] = ResultPath(Path = 
                                                 self._unassembled_reads2_out_file_name(),
                                                 IsWritten=True)
        
        # optional output, so we check for each
        # check for assembled reads file
        if self.Parameters['-s'].isOn():
            result['Assembled'] = ResultPath(Path = 
                                             self._assembled_out_file_name(),
                                             IsWritten=True)
        
        # check for discarded (unassembled) reads1 file
        if self.Parameters['-3'].isOn():
            result['Reads1Discarded'] = ResultPath(Path = 
                                                   self._discarded_reads1_out_file_name(), 
                                                   IsWritten=True)

        # check for discarded (unassembled) reads2 file
        if self.Parameters['-4'].isOn():
            result['Reads2Discarded'] = ResultPath(Path = 
                                                   self._discarded_reads2_out_file_name(),
                                                   IsWritten=True)
        
        # check for pretty-alignment file
        if self.Parameters['-E'].isOn():
            result['PrettyAlignments'] = ResultPath(Path = 
                                                    self._pretty_alignment_out_file_name(), 
                                                    IsWritten=True)
        
        return result


    def getHelp(self):
        """seqprep help"""
        help_str = """
        For basic help, type the following at the command line:
            'SeqPrep -h'

        Website:
            https://github.com/jstjohn/SeqPrep
        """
        return help_str


def join_paired_end_reads_seqprep(
    reads1_infile_path,
    reads2_infile_path,
    outfile_label='seqprep',
    max_overlap_ascii_q_score='J',
    min_overlap=15,
    max_mismatch_good_frac=0.02,
    min_frac_matching=0.9,
    phred_64='False',
    params={},
    working_dir=tempfile.gettempdir(),
    SuppressStderr=True,
    SuppressStdout=True,
    HALT_EXEC=False):
    """ Runs SeqPrep parameters to assemble paired-end reads.
        -reads1_infile_path : reads1.fastq infile path
        -reads2_infile_path : reads2.fastq infile path
        -max_overlap_ascii_q_score : 'J' for Illumina 1.8+ phred+33, 
                                    representing a score of 41. See: 
                                    http://en.wikipedia.org/wiki/FASTQ_format
        -min_overlap : minimum overall base pair overlap to merge two reads
        -max_mismatch_good_frac : maximum fraction of good quality mismatching
                                  bases to overlap reads
        -min_frac_matching : minimum fraction of matching bases to overlap 
                             reads
        -phred_64 : if input is in phred+64. Output will always be phred+33.
        -params : other optional SeqPrep parameters
 
         NOTE: SeqPrep always outputs gzipped files
    """

    abs_r1_path = os.path.abspath(reads1_infile_path)
    abs_r2_path = os.path.abspath(reads2_infile_path)
 
    infile_paths = [abs_r1_path, abs_r2_path]

    # check / make absolute infile paths
    for p in infile_paths:
        if not os.path.exists(p):
            raise IOError, 'Infile not found at: %s' % p


    # required by SeqPrep to assemble:
    params['-f'] = abs_r1_path
    params['-r'] = abs_r2_path
    params['-s'] = outfile_label + '_assembled.gz'
    params['-1'] = outfile_label + '_unassembled_R1.gz' 
    params['-2'] = outfile_label + '_unassembled_R2.gz'
    params['-o'] = min_overlap
    params['-m'] = max_mismatch_good_frac
    params['-n'] = min_frac_matching
    params['-y'] = max_overlap_ascii_q_score


    # set up controller
    seqprep_app=SeqPrep(params = params,
                        WorkingDir=working_dir,
                        SuppressStderr=SuppressStderr,
                        SuppressStdout=SuppressStdout,
                        HALT_EXEC=HALT_EXEC)
    
    # if input is phred+64
    if phred_64 == 'True':
        seqprep_app.Parameters['-6'].on()

    # run assembler 
    result = seqprep_app()
  
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



    
