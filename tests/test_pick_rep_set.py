#!/usr/bin/env python

"""Tests of code for representative set picking"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Kyle Bittinger", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"

from os import remove
from cogent import LoadSeqs
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename
from qiime.pick_rep_set import (RepSetPicker, GenericRepSetPicker, first_id,
                                first, random_id, longest_id, unique_id_map, label_to_name,
                                make_most_abundant, fasta_parse, ReferenceRepSetPicker)


class RepSetPickerTests(TestCase):

    """Tests of the abstract RepSetPicker class"""

    def test_init(self):
        """Abstract RepSetPicker __init__ should store name, params"""
        p = RepSetPicker({})
        self.assertEqual(p.Name, 'RepSetPicker')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract RepSetPicker __call__ should raise NotImplementedError"""
        p = RepSetPicker({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs',
                          '/path/to/otus')


class SharedSetupTestCase(TestCase):

    """Wrapper for shared setup stuff"""

    def setUp(self):
        # create the temporary input files
        self.tmp_seq_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath, 'w')
        seq_file.write(dna_seqs)
        seq_file.close()

        self.tmp_otu_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest_',
            suffix='.otu')
        otu_file = open(self.tmp_otu_filepath, 'w')
        otu_file.write(otus)
        otu_file.close()

        self.files_to_remove = [self.tmp_seq_filepath, self.tmp_otu_filepath]

        self.params = {'Algorithm': 'first', 'ChoiceF': first_id}

    def tearDown(self):
        remove_files(self.files_to_remove)


class GenericRepSetPickerTests(SharedSetupTestCase):

    """ Tests of the generic RepSet picker """

    def test_call_default_params(self):
        """GenericRepSetPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp = {'0': 'R27DLI_4812',
               '1': 'U1PLI_7889',
               '2': 'W3Cecum_4858',
               '3': 'R27DLI_3243',
               }
        app = GenericRepSetPicker(params={'Algorithm': 'first',
                                          'ChoiceF': first_id})
        obs = app(self.tmp_seq_filepath, self.tmp_otu_filepath)
        self.assertEqual(obs, exp)

    def test_call_wrapped_function(self):
        """GenericRepSetPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp = {'0': 'R27DLI_4812',
               '1': 'U1PLI_7889',
               '2': 'W3Cecum_4858',
               '3': 'R27DLI_3243',
               }
        app = GenericRepSetPicker(params={'Algorithm': 'most_abundant',
                                          'ChoiceF': make_most_abundant, 'ChoiceFRequiresSeqs': True})
        obs = app(self.tmp_seq_filepath, self.tmp_otu_filepath)
        self.assertEqual(obs, exp)

    def test_call_output_to_file(self):
        """GenericRepSetPicker.__call__ output to file functions as expected
        """

        tmp_result_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest.test_call_output_to_file_',
            suffix='.txt')

        app = GenericRepSetPicker(params=self.params)
        obs = app(self.tmp_seq_filepath, self.tmp_otu_filepath,
                  result_path=tmp_result_filepath)

        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)

        # compare data in result file to fake expected file
        self.assertEqual(result_file_str, rep_seqs_result_file_exp)
        # confirm that nothing is returned when result_path is specified
        self.assertEqual(obs, None)

    def test_call_output_to_file_sorted(self):
        """GenericRepSetPicker.__call__ output to file sorts when requested
        """

        tmp_result_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest.test_call_output_to_file_',
            suffix='.txt')

        app = GenericRepSetPicker(params=self.params)
        obs = app(self.tmp_seq_filepath, self.tmp_otu_filepath,
                  result_path=tmp_result_filepath, sort_by='seq_id')

        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)

        # compare data in result file to fake expected file
        self.assertEqual(result_file_str, rep_seqs_result_file_sorted_exp)
        # confirm that nothing is returned when result_path is specified
        self.assertEqual(obs, None)

    def test_call_log_file(self):
        """GenericRepSetPicker.__call__ writes log when expected
        """

        tmp_log_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest.test_call_output_to_file_l_',
            suffix='.txt')
        tmp_result_filepath = get_tmp_filename(
            prefix='GenericRepSetPickerTest.test_call_output_to_file_r_',
            suffix='.txt')

        app = GenericRepSetPicker(params=self.params)
        obs = app(self.tmp_seq_filepath, self.tmp_otu_filepath,
                  result_path=tmp_result_filepath, log_path=tmp_log_filepath)

        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        # remove the temp files before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)

        log_file_exp = ["GenericRepSetPicker parameters:",
                        'Algorithm:first',
                        "Application:None",
                        'ChoiceF:first',
                        'ChoiceFRequiresSeqs:False',
                        "Result path: %s" % tmp_result_filepath, ]
        # compare data in log file to fake expected log file
        for i, j in zip(log_file_str.splitlines(), log_file_exp):
            if not i.startswith('ChoiceF:'):  # can't test, different each time
                self.assertEqual(i, j)


class ReferenceRepSetPickerTests(SharedSetupTestCase):

    """Tests of the ReferenceRepSetPickerclass """

    def setUp(self):
        # create the temporary input files
        self.tmp_seq_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath, 'w')
        seq_file.write(dna_seqs)
        seq_file.close()

        self.ref_seq_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.fasta')
        seq_file = open(self.ref_seq_filepath, 'w')
        seq_file.write(reference_seqs)
        seq_file.close()

        self.tmp_otu_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.otu')
        otu_file = open(self.tmp_otu_filepath, 'w')
        otu_file.write(otus_w_ref)
        otu_file.close()

        self.result_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.fasta')
        otu_file = open(self.result_filepath, 'w')
        otu_file.write(otus_w_ref)
        otu_file.close()

        self.files_to_remove = [self.tmp_seq_filepath,
                                self.tmp_otu_filepath,
                                self.ref_seq_filepath,
                                self.result_filepath]

        self.params = {'Algorithm': 'first', 'ChoiceF': first_id}

    def test_call_default_params(self):
        """ReferenceRepSetPicker.__call__ expected clusters default params"""

        exp = {'0': ('R27DLI_4812', 'CTGGGCCGTATCTC'),
               'ref1': ('ref1', 'GGGGGGGAAAAAAAAAAAAA'),
               '2': ('W3Cecum_4858', 'TTGGGCCGTGTCTCAGT'),
               'ref0': ('ref0', 'CCCAAAAAAATTTTTT'),
               }
        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})
        obs = app(self.tmp_seq_filepath,
                  self.tmp_otu_filepath,
                  self.ref_seq_filepath)
        self.assertEqual(obs, exp)

    def test_call_write_to_file(self):
        """ReferenceRepSetPicker.__call__ otu map correctly written to file"""
        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})
        app(self.tmp_seq_filepath,
            self.tmp_otu_filepath,
            self.ref_seq_filepath,
            result_path=self.result_filepath)
        exp = rep_seqs_reference_result_file_exp
        self.assertEqual(LoadSeqs(self.result_filepath, aligned=False),
                         LoadSeqs(data=exp, aligned=False))

    def test_non_ref_otus(self):
        """ReferenceRepSetPicker.__call__ same result as Generic when no ref otus
        """
        exp = {'0': ('R27DLI_4812', 'CTGGGCCGTATCTC'),
               '1': ('U1PLI_7889', 'TTGGACCGTG'),
               '2': ('W3Cecum_4858', 'TTGGGCCGTGTCTCAGT'),
               '3': ('R27DLI_3243', 'CTGGACCGTGTCT')}
        tmp_otu_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.otu')
        otu_file = open(tmp_otu_filepath, 'w')
        otu_file.write(otus)
        otu_file.close()

        self.files_to_remove.append(tmp_otu_filepath)

        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})
        obs = app(self.tmp_seq_filepath,
                  tmp_otu_filepath,
                  self.ref_seq_filepath)
        self.assertEqual(obs, exp)

    def test_call_invalid_id(self):
        """ReferenceRepSetPicker.__call__ expected clusters default params"""
        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})

        tmp_otu_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.otu')
        otu_file = open(tmp_otu_filepath, 'w')
        # replace a valid sequence identifier with an invalid
        # sequence identifier (i.e., one that we don't have a sequence for)
        otu_file.write(otus_w_ref.replace('R27DLI_4812', 'bad_seq_identifier'))
        otu_file.close()
        self.files_to_remove.append(tmp_otu_filepath)

        # returning in dict
        self.assertRaises(KeyError,
                          app,
                          self.tmp_seq_filepath,
                          tmp_otu_filepath,
                          self.ref_seq_filepath)
        # writing to file
        self.assertRaises(KeyError,
                          app,
                          self.tmp_seq_filepath,
                          tmp_otu_filepath,
                          self.ref_seq_filepath,
                          result_path=self.result_filepath)

    def test_call_ref_only(self):
        """ReferenceRepSetPicker.__call__ functions with no non-refseqs"""

        tmp_otu_filepath = get_tmp_filename(
            prefix='ReferenceRepSetPickerTest_',
            suffix='.otu')
        otu_file = open(tmp_otu_filepath, 'w')
        otu_file.write(otus_all_ref)
        otu_file.close()
        self.files_to_remove.append(tmp_otu_filepath)

        exp = {'ref1': ('ref1', 'GGGGGGGAAAAAAAAAAAAA'),
               'ref0': ('ref0', 'CCCAAAAAAATTTTTT')}

        # passing only reference (not input seqs)
        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})
        obs = app(None,
                  tmp_otu_filepath,
                  self.ref_seq_filepath)
        self.assertEqual(obs, exp)

        # passing reference and input seqs
        app = ReferenceRepSetPicker(params={'Algorithm': 'first',
                                            'ChoiceF': first_id})
        obs = app(self.tmp_seq_filepath,
                  tmp_otu_filepath,
                  self.ref_seq_filepath)
        self.assertEqual(obs, exp)

    def test_call_alt_non_ref_picker(self):
        """ReferenceRepSetPicker.__call__ handles alt non-ref picking method"""

        exp = {'0': ('U1PLI_9526', 'CTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCG'
                     'GTCTCTCAACCCGGCTACCCATCGCGGGCTAGGTGGGCCGTT'
                     'ACCCCGCCTACTACCTAATGGGCCGCGACCCCATCCCTTGCCGTCTGGGC'
                     'TTTCCCGGGCCCCCCAGGAGGGGGGCGAGGAGTATCCGGTATTAGCCTCGGTT'
                     'TCCCAAGGTTGTCCCGGAGCAAGGGGCAGGTTGGTCACGTGTTACTCACCCGT'
                     'TCGCCACTTCATGTCCGCCCGAGGGCGGTTTCATCG'),
               'ref1': ('ref1', 'GGGGGGGAAAAAAAAAAAAA'),
               '2': ('W3Cecum_4858', 'TTGGGCCGTGTCTCAGT'),
               'ref0': ('ref0', 'CCCAAAAAAATTTTTT'),
               }
        app = ReferenceRepSetPicker(params={'Algorithm': 'longest',
                                            'ChoiceF': longest_id})
        obs = app(self.tmp_seq_filepath,
                  self.tmp_otu_filepath,
                  self.ref_seq_filepath)
        self.assertEqual(obs, exp)


class TopLevelTests(SharedSetupTestCase):

    """Tests of top-level functions"""

    def test_first(self):
        """first should always return first item"""
        vals = [3, 4, 2]
        self.assertEqual(first(vals), 3)
        vals.reverse()
        self.assertEqual(first(vals), 2)

    def test_first_id(self):
        """first_id should return first id from list"""
        ids = \
            "R27DLI_4812 R27DLI_600  R27DLI_727  U1PLI_403   U1PLI_8969".split(
            )
        self.assertEqual(first_id(ids, {}), 'R27DLI_4812')

    def test_random_id(self):
        """random_id should return random id from list"""
        ids = \
            "R27DLI_4812 R27DLI_600  R27DLI_727  U1PLI_403   U1PLI_8969".split(
            )
        assert random_id(ids, {}) in ids
        # just test we got something from the list, don't add stochastic test

    def test_longest_id(self):
        """longest_id should return id associated with longest seq"""
        ids = \
            "R27DLI_4812 R27DLI_600  R27DLI_727  U1PLI_403   U1PLI_8969".split(
            )
        seqs = dict(fasta_parse(dna_seqs.splitlines(),
                                       label_to_name=label_to_name))
        self.assertEqual(longest_id(ids, seqs), 'U1PLI_403')

    def test_unique_id_map(self):
        """unique_id_map should return map of seqs:unique representatives"""
        seqs = {'a': 'AG', 'b': 'AG', 'c': 'CC', 'd': 'CT'}
        obs = unique_id_map(seqs)
        exp = {'c': ['c'], 'd': ['d'], 'a': ['a', 'b'], 'b': ['a', 'b']}
        # can't predict if a or b
        for k in obs:
            assert obs[k] in exp[k]

    def test_make_most_abundant(self):
        """make_most_abundant should return function with correct behavior"""
        ids = \
            "R27DLI_4812 R27DLI_600  R27DLI_727  U1PLI_403   U1PLI_8969".split(
            )
        seqs = dict(fasta_parse(dna_seqs.splitlines(),
                                       label_to_name=label_to_name))
        f = make_most_abundant(seqs)
        result = f(ids, seqs)
        assert result in ['R27DLI_4812', 'R27DLI_727', 'U1PLI_8969']


dna_seqs = """>R27DLI_4812 FMSX0OV01EIYV5 orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGGCCGTATCTC
>R27DLI_600 FMSX0OV01D110Y orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGGCCGTATCTCA
>R27DLI_727 FMSX0OV01D5X55 orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGGCCGTATCTC
>U1PLI_403 FMSX0OV01DVG99 orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGGCCGTATCTCAGTCCCAA
>U1PLI_8969 FMSX0OV01ARWY7 orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGGCCGTATCTC
>U1PLI_9080 FMSX0OV01C9JUX orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGGCCG
>U1PLI_9526 FMSX0OV01EUN7B orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGGCCGTATCTCAGTCCCAATGTGGCCGGTCGGTCTCTCAACCCGGCTACCCATCGCGGGCTAGGTGGGCCGTTACCCCGCCTACTACCTAATGGGCCGCGACCCCATCCCTTGCCGTCTGGGCTTTCCCGGGCCCCCCAGGAGGGGGGCGAGGAGTATCCGGTATTAGCCTCGGTTTCCCAAGGTTGTCCCGGAGCAAGGGGCAGGTTGGTCACGTGTTACTCACCCGTTCGCCACTTCATGTCCGCCCGAGGGCGGTTTCATCG
>W3Cecum_6642 FMSX0OV01CW7FI orig_bc=GATACGTCCTGA new_bc=GATACGTCCTGA bc_diffs=0
CTGGGCCGTATCTCAGT
>W3Cecum_8992 FMSX0OV01C3YXK orig_bc=GATACGTCCTGA new_bc=GATACGTCCTGA bc_diffs=0
CTGGGCCGTGTCTC
>U1PLI_7889 FMSX0OV01C6HRL orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
TTGGACCGTG
>W3Cecum_4858 FMSX0OV01BX4KM orig_bc=GATACGTCCTGA new_bc=GATACGTCCTGA bc_diffs=0
TTGGGCCGTGTCTCAGT
>R27DLI_3243 FMSX0OV01DH41R orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGACCGTGTCT
>R27DLI_4562 FMSX0OV01EJKLT orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGACCGTGTCT
>R27DLI_6828 FMSX0OV01BCWTL orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGACCGTGTCT
>R27DLI_9097 FMSX0OV01APUV6 orig_bc=CTTGATGCGTAT new_bc=CTTGATGCGTAT bc_diffs=0
CTGGACCGTGTCT
>U1PLI_2780 FMSX0OV01E2K1S orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGACCGTGTCTC
>U1PLI_67 FMSX0OV01DO1NS orig_bc=TACAGATGGCTC new_bc=TACAGATGGCTC bc_diffs=0
CTGGACCGTGT
>U9PSI_10475 FMSX0OV01BB4Q3 orig_bc=GATAGCTGTCTT new_bc=GATAGCTGTCTT bc_diffs=0
CTGGACCGTGTCTC
>U9PSI_4341 FMSX0OV01B8SXV orig_bc=GATAGCTGTCTT new_bc=GATAGCTGTCTT bc_diffs=0
CTGGACCGTGTCT
>W3Cecum_5191 FMSX0OV01BMU6R orig_bc=GATACGTCCTGA new_bc=GATACGTCCTGA bc_diffs=0
CTGGACCGTGTCT
"""

otus = """0	R27DLI_4812	R27DLI_600	R27DLI_727	U1PLI_403	U1PLI_8969	U1PLI_9080	U1PLI_9526	W3Cecum_6642	W3Cecum_8992
1	U1PLI_7889
2	W3Cecum_4858
3	R27DLI_3243	R27DLI_4562	R27DLI_6828	R27DLI_9097	U1PLI_2780	U1PLI_67	U9PSI_10475	U9PSI_4341	W3Cecum_5191
"""

rep_seqs_result_file_exp = """>0 R27DLI_4812
CTGGGCCGTATCTC
>1 U1PLI_7889
TTGGACCGTG
>2 W3Cecum_4858
TTGGGCCGTGTCTCAGT
>3 R27DLI_3243
CTGGACCGTGTCT
"""

rep_seqs_result_file_sorted_exp = """>3 R27DLI_3243
CTGGACCGTGTCT
>0 R27DLI_4812
CTGGGCCGTATCTC
>2 W3Cecum_4858
TTGGGCCGTGTCTCAGT
>1 U1PLI_7889
TTGGACCGTG
"""

otus_w_ref = """0	R27DLI_4812	R27DLI_600	R27DLI_727	U1PLI_403	U1PLI_8969	U1PLI_9080	U1PLI_9526	W3Cecum_6642	W3Cecum_8992
ref1	U1PLI_7889
2	W3Cecum_4858
ref0	R27DLI_3243	R27DLI_4562	R27DLI_6828	R27DLI_9097	U1PLI_2780	U1PLI_67	U9PSI_10475	U9PSI_4341	W3Cecum_5191
"""

otus_all_ref = """ref1	U1PLI_7889
ref0	R27DLI_3243	R27DLI_4562	R27DLI_6828	R27DLI_9097	U1PLI_2780	U1PLI_67	U9PSI_10475	U9PSI_4341	W3Cecum_5191
"""

reference_seqs = """>ref0
CCCAAAAAAATTTTTT
>ref1 some comment
GGGGGGGAAAAAAAAAAAAA
>ref2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAA
"""

rep_seqs_reference_result_file_exp = """>0 R27DLI_4812
CTGGGCCGTATCTC
>ref1 ref1
GGGGGGGAAAAAAAAAAAAA
>2 W3Cecum_4858
TTGGGCCGTGTCTCAGT
>ref0 ref0
CCCAAAAAAATTTTTT
"""


# run unit tests if run from command-line
if __name__ == '__main__':
    main()
