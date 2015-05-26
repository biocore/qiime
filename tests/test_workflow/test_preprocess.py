#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from unittest import TestCase, main

from qiime.workflow.preprocess import (get_pairs, get_matching_files,
    create_commands_jpe, create_commands_eb, create_commands_slf)

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
        actual_paired_reads, actual_bc_reads =\
            get_pairs(input_forward_reads, read1_indicator, read2_indicator)

        expected_paired_reads =\
            {'/home/reads/read1_R1_000.fastq':'/home/reads/read1_R2_000.fastq'}
        expected_bc_reads = {}

        self.assertEqual(actual_paired_reads, expected_paired_reads)
        self.assertEqual(actual_bc_reads, expected_bc_reads)

    def test_get_pairs_bcs(self):
        """ Properly returns pairs of matching reads plus barcodes """

        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq', '/home/not_a_matched_read/_R1_',
            '/home/reads/read1_I1_000.fastq']
        actual_paired_reads, actual_bc_reads =\
            get_pairs(input_forward_reads, read1_indicator, read2_indicator,
            match_barcodes = True)

        expected_paired_reads =\
            {'/home/reads/read1_R1_000.fastq':'/home/reads/read1_R2_000.fastq'}
        expected_bc_reads = \
            {'/home/reads/read1_R1_000.fastq':'/home/reads/read1_I1_000.fastq'}

        self.assertEqual(actual_paired_reads, expected_paired_reads)
        self.assertEqual(actual_bc_reads, expected_bc_reads)

    def test_get_pairs_halts_on_bad_data(self):
        """ Properly raises errors when fastq has no read identifier """

        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq',
            '/home/no_matching_indicator_name']
        self.assertRaises(ValueError, get_pairs, input_forward_reads,
            read1_indicator, read2_indicator)

    def test_get_pairs_halts_on_bad_data_bcs(self):
        """ Properly raises errors when fastq has no read/bcs identifier """

        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        bc_indicator = "_wrong_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq',
            '/home/reads/read1_notmatching_000.fastq']
        self.assertRaises(ValueError, get_pairs, input_forward_reads,
            read1_indicator, read2_indicator, match_barcodes = True,
            barcode_indicator = bc_indicator)

    def test_get_pairs_halts_on_nonmatched_bcs(self):
        """ Properly raises errors when forward read found without barcode """

        read1_indicator = "_R1_"
        read2_indicator = "_R2_"
        input_forward_reads = ['/home/reads/read1_R1_000.fastq',
            '/home/reads/read1_R2_000.fastq']
        self.assertRaises(ValueError, get_pairs, input_forward_reads,
            read1_indicator, read2_indicator, match_barcodes = True,
            barcode_indicator = '_barcode_')

    def test_get_matching_files(self):
        """ Properly returns matching reads/barcodes/mapping files """

        all_fastq = ['sample1_r1_000.fastq', 'sample2_r1_000.fastq',
            'sample1_bc_000.fastq', 'sample2_bc_000.fastq']
        all_mapping = ['sample1_mapping_000.txt', 'sample2_mapping_000.txt']
        read_indicator = 'r1'
        barcode_indicator = 'bc'
        mapping_indicator = 'mapping'

        actual_matching_files = get_matching_files(all_fastq, all_mapping,
            read_indicator, barcode_indicator, mapping_indicator)

        actual_reads = set(actual_matching_files.keys())
        actual_bcs_mapping = set(actual_matching_files.values())

        expected_matching_reads = set(['sample1_r1_000.fastq',
            'sample2_r1_000.fastq'])
        expected_matching_bcs_reads = set([('sample1_bc_000.fastq',
            'sample1_mapping_000.txt'),
            ('sample2_bc_000.fastq', 'sample2_mapping_000.txt')])

        self.assertEqual(actual_reads, expected_matching_reads)
        self.assertEqual(actual_bcs_mapping, expected_matching_bcs_reads)

    def test_get_matching_files_error_on_invalid_mapping(self):
        """ Raises error when mapping file text not matched"""

        all_fastq = ['sample1_r1_000.fastq', 'sample2_r1_000.fastq',
            'sample1_bc_000.fastq', 'sample2_bc_000.fastq']
        all_mapping = ['sample1_badMap_000.txt', 'sample2_mapping_000.txt']
        read_indicator = 'r1'
        barcode_indicator = 'bc'
        mapping_indicator = 'mapping'

        self.assertRaises(IndexError, get_matching_files, all_fastq,
            all_mapping, read_indicator, barcode_indicator, mapping_indicator)


    def test_get_matching_files_error_on_invalid_reads(self):
        """ Raises error when read file text not matched"""

        all_fastq = ['sample1_BadRead_000.fastq', 'sample2_r1_000.fastq',
            'sample1_bc_000.fastq', 'sample2_bc_000.fastq']
        all_mapping = ['sample1_mapping_000.txt', 'sample2_mapping_000.txt']
        read_indicator = 'r1'
        barcode_indicator = 'bc'
        mapping_indicator = 'mapping'

        self.assertRaises(ValueError, get_matching_files, all_fastq,
            all_mapping, read_indicator, barcode_indicator, mapping_indicator)

    def test_get_matching_files_error_on_invalid_barcodes(self):
        """ Raises error when barcode file text not matched"""

        all_fastq = ['sample1_r1_000.fastq', 'sample2_r1_000.fastq',
            'sample1_BadBarcodes_000.fastq', 'sample2_bc_000.fastq']
        all_mapping = ['sample1_mapping_000.txt', 'sample2_mapping_000.txt']
        read_indicator = 'r1'
        barcode_indicator = 'bc'
        mapping_indicator = 'mapping'

        self.assertRaises(ValueError, get_matching_files, all_fastq,
            all_mapping, read_indicator, barcode_indicator, mapping_indicator)

    def test_get_matching_files_error_on_missing_mapping(self):
        """ Raises error when reads has no matching mapping file"""

        all_fastq = ['sample1_r1_000.fastq', 'sample2_r1_000.fastq',
            'sample1_bc_000.fastq', 'sample2_bc_000.fastq']
        all_mapping = ['sample2_mapping_000.txt']
        read_indicator = 'r1'
        barcode_indicator = 'bc'
        mapping_indicator = 'mapping'

        self.assertRaises(KeyError, get_matching_files, all_fastq,
            all_mapping, read_indicator, barcode_indicator, mapping_indicator)

    def test_create_commands_jpe(self):
        """ Properly creates commands """

        pairs = {'f_read1.fastq':'r_read1.fastq',
            'f_read2.fastq':'r_read2.fastq'}
        base_output_dir = "output_dir/"
        optional_params = ""
        leading_text = ""
        trailing_text = ""
        include_input_dir_path=False
        remove_filepath_in_name=False

        commands = create_commands_jpe(pairs, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)

        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)

        expected_commands = set(['join_paired_ends.py  -f f_read1.fastq -r r_read1.fastq -o output_dir/f_read1 ',
            'join_paired_ends.py  -f f_read2.fastq -r r_read2.fastq -o output_dir/f_read2 '])

        self.assertEqual(actual_commands, expected_commands)

    def test_create_commands_jpe_bcs(self):
        """ Properly creates commands with matched bcs """

        pairs = {'f_read1.fastq':'r_read1.fastq',
            'f_read2.fastq':'r_read2.fastq'}
        base_output_dir = "output_dir/"
        optional_params = ""
        leading_text = ""
        trailing_text = ""
        include_input_dir_path=False
        remove_filepath_in_name=False
        match_barcodes = True
        bc_pairs = {'f_read1.fastq':'b_read1.fastq',
            'f_read2.fastq':'b_read2.fastq'}

        commands = create_commands_jpe(pairs, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name,
            match_barcodes, bc_pairs)

        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)

        expected_commands = set(['join_paired_ends.py  -b b_read1.fastq -f f_read1.fastq -r r_read1.fastq -o output_dir/f_read1 ',
            'join_paired_ends.py  -b b_read2.fastq -f f_read2.fastq -r r_read2.fastq -o output_dir/f_read2 '])

        self.assertEqual(actual_commands, expected_commands)

    def test_create_commands_jpe_added_options(self):
        """ Properly creates commands with all optional parameters """

        pairs = {'dir1/f_read1.fastq':'dir1/r_read1.fastq',
            'dir2/f_read2.fastq':'dir2/r_read2.fastq'}
        base_output_dir = "output_dir/"
        optional_params = "-m OtherMethod"
        leading_text = "Echo"
        trailing_text = "| qsub COMMANDxxx"
        include_input_dir_path=True
        remove_filepath_in_name=True

        commands = create_commands_jpe(pairs, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)

        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)

        expected_commands = set(['Echo join_paired_ends.py -m OtherMethod -f dir2/f_read2.fastq -r dir2/r_read2.fastq -o output_dir/dir2 | qsub COMMANDxxx',
            'Echo join_paired_ends.py -m OtherMethod -f dir1/f_read1.fastq -r dir1/r_read1.fastq -o output_dir/dir1 | qsub COMMANDxxx'])

        self.assertEqual(actual_commands, expected_commands)

    def test_create_commands_eb_unpaired(self):
        """ Properly creates commands """

        pairs = ['read1.fastq','read2.fastq']
        ispaired = False
        base_output_dir = "output_dir/"
        optional_params = ""
        leading_text = ""
        trailing_text = ""
        include_input_dir_path=False
        remove_filepath_in_name=False

        commands = create_commands_eb(pairs, ispaired, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)

        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)

        expected_commands = set(['extract_barcodes.py  -f read1.fastq -o output_dir/read1 ',
            'extract_barcodes.py  -f read2.fastq -o output_dir/read2 '])

        self.assertEqual(actual_commands, expected_commands)

    def test_create_commands_eb_added_options(self):
        """ Properly creates commands with all optional parameters """

        pairs = {'dir1/f_read1.fastq':'dir1/r_read1.fastq',
            'dir2/f_read2.fastq':'dir2/r_read2.fastq'}
        ispaired = True
        base_output_dir = "output_dir/"
        optional_params = "-m OtherMethod"
        leading_text = "Echo"
        trailing_text = "| qsub COMMANDxxx"
        include_input_dir_path=True
        remove_filepath_in_name=True

        commands = create_commands_eb(pairs, ispaired, base_output_dir,
            optional_params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)

        actual_commands = []
        for n in commands:
            actual_commands.append(n[0][1])
        actual_commands = set(actual_commands)

        expected_commands = set(['Echo extract_barcodes.py -m OtherMethod -f dir2/f_read2.fastq -r dir2/r_read2.fastq -o output_dir/dir2 | qsub COMMANDxxx',
            'Echo extract_barcodes.py -m OtherMethod -f dir1/f_read1.fastq -r dir1/r_read1.fastq -o output_dir/dir1 | qsub COMMANDxxx'])

        self.assertEqual(actual_commands, expected_commands)

    def test_create_commands_slf_by_sampleID(self):
        """ Properly creates commands for SampleIDs options """

        all_files = ['sample1_r1.fastq', 'sample2_r1.fastq']
        demultiplexing_method = 'sampleid_by_file'
        output_dir = "sl_out"

        actual_command = create_commands_slf(all_files, demultiplexing_method,
            output_dir)[0][0][1]

        expected_command = "split_libraries_fastq.py  -i sample1_r1.fastq,sample2_r1.fastq --sample_ids sample1,sample2 -o sl_out  --barcode_type 'not-barcoded'"

        self.assertEqual(actual_command, expected_command)

    def test_create_commands_slf_mapping_barcodes(self):
        """ Properly creates commands for barcode/mapping files option """

        all_files = {'sample1_r1.fastq':('sample1_mapping.txt',
            'sample1_bc.fastq'), 'sample2_r1.fastq':('sample2_mapping.txt',
            'sample2_bc.fastq')}
        demultiplexing_method = 'mapping_barcode_files'
        output_dir = "sl_out"

        actual_command = create_commands_slf(all_files, demultiplexing_method,
            output_dir)[0][0][1]

        expected_command = "split_libraries_fastq.py  -i sample1_r1.fastq,sample2_r1.fastq --barcode_read_fps sample1_mapping.txt,sample2_mapping.txt --mapping_fps sample1_bc.fastq,sample2_bc.fastq -o sl_out "
        self.assertEqual(actual_command, expected_command)

    def test_create_commands_slf_added_options(self):
        """ Properly creates slf commands with added parameters """

        all_files = ['sample1/sample_r1.fastq', 'sample2/sample_r1.fastq']
        demultiplexing_method = 'sampleid_by_file'
        output_dir = "sl_out"
        params = "--max_bad_run_length 15"
        leading_text = "echo"
        trailing_text = " | qsub -N BigErn -k oe"
        include_input_dir_path = True
        remove_filepath_in_name = True

        actual_command = create_commands_slf(all_files, demultiplexing_method,
            output_dir, params, leading_text, trailing_text,
            include_input_dir_path, remove_filepath_in_name)[0][0][1]

        expected_command = "echo split_libraries_fastq.py --max_bad_run_length 15 -i sample1/sample_r1.fastq,sample2/sample_r1.fastq --sample_ids sample1,sample2 -o sl_out  | qsub -N BigErn -k oe --barcode_type 'not-barcoded'"

        self.assertEqual(actual_command, expected_command)

if __name__ == '__main__':
    main()
