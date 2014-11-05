#!/usr/bin/env python
from __future__ import division
 
__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
 
from unittest import TestCase, main

from qiime.generate_join_paired_ends_commands import (get_pairs,
                                                      create_commands)
 
class GenerateJoinPairedEndsCommands(TestCase):

    def setUp(self):
        # create the temporary input files that will be used
        # Shouldn't need any for this script, as all file IO is in scripts/
        pass

    def tearDown(self):
        # Should not need any file/folder removal
        pass

    def test_get_pairs(self):
        """ Properly returns pairs of matching fastq forward/reverse reads """
        
        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq', '/home/not_a_matched_read/_R1_']
        actual_paired_reads = get_pairs(input_forward_reads, read1_indicator,
            read2_indicator)
        
        expected_paired_reads =\
            {'/home/reads/read1_R1_000.fastq':'/home/reads/read1_R2_000.fastq'}
        
        self.assertEqual(actual_paired_reads, expected_paired_reads)
        
    def test_get_pairs_halts_on_bad_data(self):
        """ Properly raises errors when fastq has no read identifier """
        
        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq',
            '/home/no_matching_indicator_name']
        self.assertRaises(ValueError, get_pairs, input_forward_reads,
            read1_indicator, read2_indicator)
        
    def test_create_commands(self):
        """ Properly creates commands """
        
        pairs = {'f_read1.fastq':'r_read1.fastq',
            'f_read2.fastq':'r_read2.fastq'}
        base_output_dir = "output_dir/"
        optional_params = ""
        leading_text = ""
        trailing_text = ""
        include_input_dir_path=False
        remove_filepath_in_name=False
        
        commands = create_commands(pairs, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)
        
        # Because of dictionary used, must get set of output commands to ensure
        # that order does not cause spurious unit test failures.
        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)
        
        expected_commands = set(['join_paired_ends.py  -f f_read1.fastq -r r_read1.fastq -o output_dir/f_read1 ',
            'join_paired_ends.py  -f f_read2.fastq -r r_read2.fastq -o output_dir/f_read2 '])
            
        self.assertEqual(actual_commands, expected_commands)
        
    def test_create_commands_added_options(self):
        """ Properly creates commands with all optional parameters """
        
        pairs = {'dir1/f_read1.fastq':'dir1/r_read1.fastq',
            'dir2/f_read2.fastq':'dir2/r_read2.fastq'}
        base_output_dir = "output_dir/"
        optional_params = "-m OtherMethod"
        leading_text = "Echo "
        trailing_text = "| qsub COMMANDxxx"
        include_input_dir_path=True
        remove_filepath_in_name=True
        
        commands = create_commands(pairs, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)
        
        # Because of dictionary used, must get set of output commands to ensure
        # that order does not cause spurious unit test failures.
        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)
        
        expected_commands = set(['Echo join_paired_ends.py -m OtherMethod -f dir2/f_read2.fastq -r dir2/r_read2.fastq -o output_dir/dir2 | qsub COMMANDxxx',
            'Echo join_paired_ends.py -m OtherMethod -f dir1/f_read1.fastq -r dir1/r_read1.fastq -o output_dir/dir1 | qsub COMMANDxxx'])
            
        self.assertEqual(actual_commands, expected_commands)

if __name__ == '__main__':
    main()