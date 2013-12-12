#!/usr/bin/env python
# File created February 29, 2012
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import exists, join, basename
from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename, create_dir
from cogent.util.misc import remove_files, get_random_directory_name

from qiime.truncate_reverse_primer import get_rev_primer_seqs,\
    get_output_filepaths, truncate_rev_primers, truncate_reverse_primer
    
class FakeOutFile(object):
    
    def __init__(self):
        self.data = ""
    
    def write(self,s):
        self.data += s

class TruncateRemoveReversePrimerTests(TestCase):
    """ Unit tests for the truncate_reverse_primers.py module """
    
    def setUp(self):
        # create the temporary input files that will be used
        
        self._files_to_remove = []
        
        self.sample_fasta_file1_data = sample_fasta_file1
        self.sample_fasta_file_bad_labels_data =\
         sample_fasta_file_bad_labels
        
        self.sample_mapping_file1_data = sample_mapping_file1
        self.sample_mapping_file_no_revprimer_header =\
         sample_mapping_file_no_revprimer_header
        self.sample_mapping_file_bad_revprimer =\
         sample_mapping_file_bad_revprimer
        self.expected_truncation_default_settings =\
         expected_truncation_default_settings
        self.expected_truncation_zero_mismatches =\
         expected_truncation_zero_mismatches
        self.expected_truncation_zero_mismatches_truncate_remove =\
         expected_truncation_zero_mismatches_truncate_remove
         
        self.fasta_fp = get_tmp_filename(\
         prefix = 'fasta_seqs_',
         suffix = '.fna')
        seq_file = open(self.fasta_fp, 'w')
        seq_file.write(self.sample_fasta_file1_data)
        seq_file.close()
        
        self.fasta_badlabels_fp = get_tmp_filename(\
         prefix = "fasta_seqs_badlabels_",
         suffix = ".fna")
        seq_file = open(self.fasta_badlabels_fp, "w")
        seq_file.write(self.sample_fasta_file_bad_labels_data)
        seq_file.close()
        
        self.mapping_fp = get_tmp_filename(\
         prefix = 'sample_mapping_',
         suffix = '.txt')
        mapping_file = open(self.mapping_fp, "w")
        mapping_file.write(self.sample_mapping_file1_data)
        mapping_file.close()
        
        self.mapping_bad_header_fp = get_tmp_filename(\
         prefix = 'sample_mapping_badheader_',
         suffix = ".txt")
        mapping_file = open(self.mapping_bad_header_fp, "w")
        mapping_file.write(self.sample_mapping_file_no_revprimer_header)
        mapping_file.close()
        
        self.mapping_bad_primer_fp = get_tmp_filename(\
         prefix = 'sample_mapping_badprimer_',
         suffix = ".txt")
        mapping_file = open(self.mapping_bad_primer_fp, "w")
        mapping_file.write(self.sample_mapping_file_bad_revprimer)
        mapping_file.close()
        
        self.output_dir = get_random_directory_name(prefix = '/tmp/')
        self.output_dir += '/'
        
        create_dir(self.output_dir)
        
        self._files_to_remove =\
         [self.fasta_fp, self.mapping_fp, self.mapping_bad_header_fp,
         self.mapping_bad_primer_fp, self.fasta_badlabels_fp]
        
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir ):
            rmtree(self.output_dir )
            
    def test_get_rev_primer_seqs(self):
        """ Builds dict of reverse primers from mapping file """
        
        actual_rev_primers = get_rev_primer_seqs(open(self.mapping_fp, "U"))
        
        expected_rev_primers = {'PC.481': ['CTCTCCG'], 'PC.634': ['CTCTCAG'],
        'PC.635': ['CTCTCAG'], 'PC.636': ['CTCTCAG'], 'PC.354': ['CTCTCAG',
        'TTCTCRG']}
	    
        self.assertEqual(actual_rev_primers, expected_rev_primers)
	    
    def test_get_rev_primer_seqs_errors(self):
        """Raises errors with invalid mapping file """
        
        # Raises error if missing ReversePrimer column
        self.assertRaises(KeyError, get_rev_primer_seqs,
         open(self.mapping_bad_header_fp, "U"))
         
        # Raises error if invalid characters in primer
        self.assertRaises(ValueError, get_rev_primer_seqs,
         open(self.mapping_bad_primer_fp, "U"))
	    
    def test_get_output_filepaths(self):
        """ Properly returns output filepaths """
        
        
        actual_fna_fp, actual_log_fp = get_output_filepaths(".", 
         '/home/tests/seqs.fna')
         
        expected_fna_fp = "./seqs_rev_primer_truncated.fna"
        expected_log_fp = "./rev_primer_truncation.log"
        
        self.assertEqual(actual_fna_fp, expected_fna_fp)
        self.assertEqual(actual_log_fp, expected_log_fp)
	    
    def test_truncate_rev_primers(self):
        """ Properly truncates remove primers """
        
        out_f = FakeOutFile()
        
        rev_primers = {'PC.481': ['CTCTCCG'], 'PC.634': ['CTCTCAG'],
        'PC.635': ['CTCTCAG'], 'PC.636': ['CTCTCAG'], 'PC.354': ['CTCTCAG']}
        
        # Use default options, all sequences should get truncated and written
        log_data = truncate_rev_primers(open(self.fasta_fp, "U"),
         out_f, rev_primers)

        expected_log_data = {'seqs_written': 5, 'total_seqs': 6,
        'sample_id_not_found': 0, 'reverse_primer_not_found': 0}
        
        self.assertEqual(log_data, expected_log_data)
        
        # Note that because these are short sequences, two mismatches allows
        # for a very short truncation of one of the sequences
        self.assertEqual(out_f.data, self.expected_truncation_default_settings)
        
        # With zero mismatches will not truncate all seqs
        out_f = FakeOutFile()
        
        # Use default options, all sequences should get truncated and written
        log_data = truncate_rev_primers(open(self.fasta_fp, "U"),
         out_f, rev_primers, primer_mismatches=0)

        expected_log_data = {'seqs_written': 5, 'total_seqs': 6,
        'sample_id_not_found': 0, 'reverse_primer_not_found': 2}
        
        self.assertEqual(log_data, expected_log_data)
        
        # With zero mismatches allowed, 2 seqs should not be truncated
        self.assertEqual(out_f.data, self.expected_truncation_zero_mismatches)
        
        # With zero mismatches and truncate_remove option, should only write
        # 3 of the sequences
        out_f = FakeOutFile()
        
        log_data = truncate_rev_primers(open(self.fasta_fp, "U"),
         out_f, rev_primers, truncate_option="truncate_remove",
         primer_mismatches=0)

        expected_log_data = {'seqs_written': 3, 'total_seqs': 6,
        'sample_id_not_found': 0, 'reverse_primer_not_found': 2}
        
        self.assertEqual(log_data, expected_log_data)
        
        # With zero mismatches allowed, 3 seqs total should be written
        self.assertEqual(out_f.data,
         self.expected_truncation_zero_mismatches_truncate_remove)
         
        # Should count sample ids not found in log
        out_f = FakeOutFile()
        
        rev_primers = {'PC.481': ['CTCTCCG'], 'PC.634': ['CTCTCAG'],
        'PC.635': ['CTCTCAG'], 'PC.636': ['CTCTCAG'], 'PC.354': ['CTCTCAG']}
        
        # Use default options, all sequences should get truncated and written
        log_data = truncate_rev_primers(open(self.fasta_badlabels_fp, "U"),
         out_f, rev_primers)

        expected_log_data = {'seqs_written': 5, 'total_seqs': 5,
        'sample_id_not_found': 5, 'reverse_primer_not_found': 0}
        
        self.assertEqual(log_data, expected_log_data)
        
        # No matches to sample IDs, so sequences are written unmodified
        self.assertEqual(out_f.data, self.sample_fasta_file_bad_labels_data)
        
	    
    def test_truncate_reverse_primer(self):
        """ Overall module functionality test """
        
        
        truncate_reverse_primer(self.fasta_fp, self.mapping_fp,
         self.output_dir, truncate_option='truncate_only', primer_mismatches=2)
         
        output_log_fp = open(join(self.output_dir,
        "rev_primer_truncation.log"), "U")
        
        actual_log_lines = [line.strip() for line in output_log_fp]
        
        # Because filepaths used are recorded, need the tmp filepaths
        expected_log_lines = ['Details for removal of reverse primers',
        'Original fasta filepath: %s' % self.fasta_fp,
        'Total seqs in fasta: 6',
        'Mapping filepath: %s' % self.mapping_fp,
        'Truncation option: truncate_only',
        'Mismatches allowed: 2',
        'Total seqs written: 5',
        'SampleIDs not found: 0',
        'Reverse primers not found: 0']
        
        self.assertEqual(actual_log_lines, expected_log_lines)
        
        expected_fna_fp = basename(self.fasta_fp.replace('.fna', '')) +\
        "_rev_primer_truncated.fna"
        
        fasta_fp = open(join(self.output_dir, expected_fna_fp), "U")
        
        actual_fasta_output = [line.strip() for line in fasta_fp]
        
        self.assertEqual(actual_fasta_output,
         self.expected_truncation_default_settings.strip().split('\n'))


            
sample_fasta_file1 = """>PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTAC
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCCTCTCAGAACCCCTATC
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTAT
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTAC
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTAC
>PC.634_6 RevPrimerAtStart FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTCTCAGCACCGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACC
"""

sample_fasta_file_bad_labels = """>CaptHammer_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTAC
>JohnnySnow_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCCTCTCAGAACCCCTATC
>BadHorse_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTAT
>Moist_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTAC
>DrHorrible_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTAC
"""

expected_truncation_default_settings = """>PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACC
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCC
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGT
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAAC
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACC
"""

expected_truncation_zero_mismatches = """>PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACC
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCC
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTAT
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTAC
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACC
"""

expected_truncation_zero_mismatches_truncate_remove = """>PC.634_1 FLP3FBN01ELBSX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACC
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCC
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACC
"""

sample_mapping_file1 = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	ReversePrimer	Description
"#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009)."						
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	CTGAGAG,CYGAGAA	Control_mouse_I.D._354
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	CGGAGAG	Control_mouse_I.D._481
PC.634	ACAGAGTCGGCT	&&&YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._636"""

sample_mapping_file_no_revprimer_header = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	XXXX	Description
"#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009)."						
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	CTGAGAG	Control_mouse_I.D._354
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	CGGAGAG	Control_mouse_I.D._481
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._636"""

sample_mapping_file_bad_revprimer = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	ReversePrimer	Description
"#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009)."						
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	CTGAGAG	Control_mouse_I.D._354
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	CGGAGAG	Control_mouse_I.D._481
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGXG	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	CTGAGAG	Fasting_mouse_I.D._636"""


if __name__ =='__main__':
    main()	