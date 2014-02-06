#!/usr/bin/env python

"""Tests of code for aligning 16S sequences"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import remove
from os.path import getsize
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment, Alignment
from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename
from qiime.align_seqs import (compute_min_alignment_length,
                              Aligner, CogentAligner, PyNastAligner, InfernalAligner,
                              alignment_module_names,
                              )


def remove_files(list_of_filepaths, error_on_missing=True):
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError(
            "Some filepaths were not accessible: %s" %
            '\t'.join(missing))


class AlignerTests(TestCase):

    """Tests of the abstract Aligner class"""

    def test_init(self):
        """Abstract Aligner __init__ should store name, params"""
        p = Aligner({})
        self.assertEqual(p.Name, 'Aligner')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract Aligner __call__ should raise NotImplementedError"""
        p = Aligner({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs')


class SharedSetupTestCase(TestCase):

    """Shared setup for aligner tests"""

    def tearDown(self):
        remove_files(self._paths_to_clean_up)


class CogentAlignerTests(SharedSetupTestCase):

    """Tests of the CogentAligner class"""

    def setUp(self):
        self.input_fp = get_tmp_filename(
            prefix='CogentAlignerTests_', suffix='.fasta')
        open(self.input_fp, 'w').write(seqs_for_muscle)

        self._paths_to_clean_up =\
            [self.input_fp]
        self.muscle_module = alignment_module_names['muscle']

    def test_call_correct_alignment(self):
        """CogentAligner: output expected alignment file
        """
        p = CogentAligner({'Module': self.muscle_module})
        log_fp = get_tmp_filename(
            prefix='CogentAlignerTests_', suffix='.log')
        self._paths_to_clean_up.append(log_fp)

        actual = p(result_path=None, seq_path=self.input_fp,
                   log_path=log_fp)
        expected = expected_muscle_alignment
        # note: lines in diff order w/ diff versions
        self.assertEqualItems(str(actual).splitlines(), expected.splitlines())

    def test_muscle_max_memory(self):
        """CogentAligner: muscle_max_memory should be passed to alignment fcn
        """
        p = CogentAligner({
            'Module': self.muscle_module,
            '-maxmb': '200',
        })
        self.assertEqual(p.Params["-maxmb"], "200")

        log_fp = get_tmp_filename(
            prefix='CogentAlignerTests_', suffix='.log')
        self._paths_to_clean_up.append(log_fp)

        actual = p(result_path=None, seq_path=self.input_fp,
                   log_path=log_fp)
        expected = expected_muscle_alignment
        # note: lines in diff order w/ diff versions
        self.assertEqualItems(str(actual).splitlines(), expected.splitlines())


class InfernalAlignerTests(SharedSetupTestCase):

    """Tests of the InfernalAligner class"""

    def setUp(self):
        self.infernal_test1_input_fp = get_tmp_filename(
            prefix='InfernalAlignerTests_', suffix='.fasta')
        open(
            self.infernal_test1_input_fp,
            'w').write(
            infernal_test1_input_fasta)

        self.infernal_test1_template_fp = get_tmp_filename(
            prefix='InfernalAlignerTests_', suffix='template.sto')
        open(self.infernal_test1_template_fp, 'w').\
            write(infernal_test1_template_stockholm)

        # create temp file names (and touch them so we can reliably
        # clean them up)
        self.result_fp = get_tmp_filename(
            prefix='InfernalAlignerTests_', suffix='.fasta')
        open(self.result_fp, 'w').close()

        self.log_fp = get_tmp_filename(
            prefix='InfernalAlignerTests_', suffix='.log')
        open(self.log_fp, 'w').close()

        self._paths_to_clean_up = [
            self.infernal_test1_input_fp,
            self.result_fp,
            self.log_fp,
            self.infernal_test1_template_fp,
        ]

        self.infernal_test1_aligner = InfernalAligner({
            'template_filepath': self.infernal_test1_template_fp,
        })
        self.infernal_test1_expected_aln = \
            LoadSeqs(data=infernal_test1_expected_alignment, aligned=Alignment,
                     moltype=DNA)

    def test_call_infernal_test1_file_output(self):
        """InfernalAligner writes correct output files for infernal_test1 seqs
        """
        # do not collect results; check output files instead
        actual = self.infernal_test1_aligner(
            self.infernal_test1_input_fp, result_path=self.result_fp,
            log_path=self.log_fp)

        self.assertTrue(actual is None,
                        "Result should be None when result path provided.")

        expected_aln = self.infernal_test1_expected_aln
        actual_aln = LoadSeqs(self.result_fp, aligned=Alignment)
        self.assertEqual(actual_aln, expected_aln)

    def test_call_infernal_test1(self):
        """InfernalAligner: functions as expected when returing objects
        """
        actual_aln = self.infernal_test1_aligner(self.infernal_test1_input_fp)
        expected_aln = self.infernal_test1_expected_aln

        expected_names = ['seq_1', 'seq_2', 'seq_3']
        self.assertEqual(sorted(actual_aln.Names), expected_names)
        self.assertEqual(actual_aln, expected_aln)


class PyNastAlignerTests(SharedSetupTestCase):

    """Tests of the PyNastAligner class"""

    def setUp(self):
        self.pynast_test1_input_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='.fasta')
        open(self.pynast_test1_input_fp, 'w').write(pynast_test1_input_fasta)

        self.pynast_test1_template_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='template.fasta')
        open(self.pynast_test1_template_fp, 'w').\
            write(pynast_test1_template_fasta)

        self.pynast_test_template_w_dots_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='template.fasta')
        open(self.pynast_test_template_w_dots_fp, 'w').\
            write(pynast_test1_template_fasta.replace('-', '.'))

        self.pynast_test_template_w_u_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='template.fasta')
        open(self.pynast_test_template_w_u_fp, 'w').\
            write(pynast_test1_template_fasta.replace('T', 'U'))

        self.pynast_test_template_w_lower_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='template.fasta')
        open(self.pynast_test_template_w_lower_fp, 'w').\
            write(pynast_test1_template_fasta.lower())

        # create temp file names (and touch them so we can reliably
        # clean them up)
        self.result_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='.fasta')
        open(self.result_fp, 'w').close()
        self.failure_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='.fasta')
        open(self.failure_fp, 'w').close()
        self.log_fp = get_tmp_filename(
            prefix='PyNastAlignerTests_', suffix='.log')
        open(self.log_fp, 'w').close()

        self._paths_to_clean_up = [
            self.pynast_test1_input_fp,
            self.result_fp,
            self.failure_fp,
            self.log_fp,
            self.pynast_test1_template_fp,
            self.pynast_test_template_w_dots_fp,
            self.pynast_test_template_w_u_fp,
            self.pynast_test_template_w_lower_fp
        ]

        self.pynast_test1_aligner = PyNastAligner({
            'template_filepath': self.pynast_test1_template_fp,
            'min_len': 15,
        })

        self.pynast_test1_expected_aln = \
            LoadSeqs(
                data=pynast_test1_expected_alignment,
                aligned=DenseAlignment)
        self.pynast_test1_expected_fail = \
            LoadSeqs(data=pynast_test1_expected_failure, aligned=False)

    def test_call_pynast_test1_file_output(self):
        """PyNastAligner writes correct output files for pynast_test1 seqs
        """
        # do not collect results; check output files instead
        actual = self.pynast_test1_aligner(
            self.pynast_test1_input_fp, result_path=self.result_fp,
            log_path=self.log_fp, failure_path=self.failure_fp)

        self.assertTrue(actual is None,
                        "Result should be None when result path provided.")

        expected_aln = self.pynast_test1_expected_aln
        actual_aln = LoadSeqs(self.result_fp, aligned=DenseAlignment)
        self.assertEqual(actual_aln, expected_aln)

        actual_fail = LoadSeqs(self.failure_fp, aligned=False)
        self.assertEqual(actual_fail.toFasta(),
                         self.pynast_test1_expected_fail.toFasta())

    def test_call_pynast_test1_file_output_alt_params(self):
        """PyNastAligner writes correct output files when no seqs align
        """
        aligner = PyNastAligner({
            'template_filepath': self.pynast_test1_template_fp,
            'min_len': 1000})

        actual = aligner(
            self.pynast_test1_input_fp, result_path=self.result_fp,
            log_path=self.log_fp, failure_path=self.failure_fp)

        self.assertTrue(actual is None,
                        "Result should be None when result path provided.")

        self.assertEqual(getsize(self.result_fp), 0,
                         "No alignable seqs should result in an empty file.")

        # all seqs reported to fail
        actual_fail = LoadSeqs(self.failure_fp, aligned=False)
        self.assertEqual(actual_fail.getNumSeqs(), 3)

    def test_call_pynast_test1(self):
        """PyNastAligner: functions as expected when returing objects
        """
        actual_aln = self.pynast_test1_aligner(self.pynast_test1_input_fp)
        expected_aln = self.pynast_test1_expected_aln

        expected_names = ['1 description field 1..23', '2 1..23']
        self.assertEqual(actual_aln.Names, expected_names)
        self.assertEqual(actual_aln, expected_aln)

    def test_call_pynast_template_aln_with_dots(self):
        """PyNastAligner: functions when template alignment contains dots
        """
        pynast_aligner = PyNastAligner({
            'template_filepath': self.pynast_test_template_w_dots_fp,
            'min_len': 15,
        })
        actual_aln = pynast_aligner(self.pynast_test1_input_fp)
        expected_aln = self.pynast_test1_expected_aln

        expected_names = ['1 description field 1..23', '2 1..23']
        self.assertEqual(actual_aln.Names, expected_names)
        self.assertEqual(actual_aln, expected_aln)

    def test_call_pynast_template_aln_with_lower(self):
        """PyNastAligner: functions when template alignment contains lower case
        """
        pynast_aligner = PyNastAligner({
            'template_filepath': self.pynast_test_template_w_lower_fp,
            'min_len': 15,
        })
        actual_aln = pynast_aligner(self.pynast_test1_input_fp)
        expected_aln = self.pynast_test1_expected_aln

        expected_names = ['1 description field 1..23', '2 1..23']
        self.assertEqual(actual_aln.Names, expected_names)
        self.assertEqual(actual_aln, expected_aln)

    def test_call_pynast_template_aln_with_U(self):
        """PyNastAligner: error message when template contains bad char
        """
        pynast_aligner = PyNastAligner({
            'template_filepath': self.pynast_test_template_w_u_fp,
            'min_len': 15,
        })
        self.assertRaises(KeyError, pynast_aligner, self.pynast_test1_input_fp)

    def test_call_pynast_alt_pairwise_method(self):
        """PyNastAligner: alternate pairwise alignment method produces correct alignment
        """
        aligner = PyNastAligner({
            'pairwise_alignment_method': 'muscle',
            'template_filepath': self.pynast_test1_template_fp,
            'min_len': 15,
        })
        actual_aln = aligner(self.pynast_test1_input_fp)
        expected_aln = self.pynast_test1_expected_aln
        self.assertEqual(actual_aln, expected_aln)

    def test_call_pynast_test1_alt_min_len(self):
        """PyNastAligner: returns no result when min_len too high
        """
        aligner = PyNastAligner({
            'template_filepath': self.pynast_test1_template_fp,
            'min_len': 1000})

        actual_aln = aligner(
            self.pynast_test1_input_fp)
        expected_aln = {}

        self.assertEqual(actual_aln, expected_aln)

    def test_call_pynast_test1_alt_min_pct(self):
        """PyNastAligner: returns no result when min_pct too high
        """
        aligner = PyNastAligner({
            'template_filepath': self.pynast_test1_template_fp,
            'min_len': 15,
            'min_pct': 100.0})

        actual_aln = aligner(self.pynast_test1_input_fp)
        expected_aln = {}

        self.assertEqual(actual_aln, expected_aln)

    def tearDown(self):
        """
        """
        remove_files(self._paths_to_clean_up)


class TopLevelTests(TestCase):

    """ tests of top-level functions """

    def setUp(self):
        """ """
        self.min_length_computation_seqs = min_length_computation_seqs.split(
            '\n')

    def test_compute_min_alignment_length(self):
        """compute_min_alignment_length: returns n std devs below mean seq len
        """
        self.assertEqual(compute_min_alignment_length(
                         self.min_length_computation_seqs), 16)
        self.assertEqual(compute_min_alignment_length(
                         self.min_length_computation_seqs, 0.60), 13)


seqs_for_muscle = \
    """>abc
ACACACAC
>def
ACAGACAC
>ghi
ACAGACACTT
>jkl
TTACAC"""

expected_muscle_alignment = """>jkl\n--TTACAC--\n>abc\nACACACAC--\n>ghi\nACAGACACTT\n>def\nACAGACAC--\n"""

infernal_test1_input_fasta = """>seq_1
ACTGCTAGCTAGTAGCGTACGTA
>seq_2
GCTACGTAGCTAC
>seq_3
GCGGCTATTAGATCGTA"""

infernal_test1_template_stockholm = """# STOCKHOLM 1.0
seq_a           TAGGCTCTGATATAATAGC-TCTC---------
seq_b           ----TATCGCTTCGACGAT-TCTCTGATAGAGA
seq_c           ------------TGACTAC-GCAT---------
#=GC SS_cons    ............((.(....)))..........
//"""

infernal_test1_expected_alignment = """>seq_1
-----ACTGCTA-GCTAGTAGCGTACGTA----
>seq_2
--------GCTACG-TAGCTAC-----------
>seq_3
-----GCGGCTATTAGATC-GTA----------
"""

pynast_test1_template_fasta = """>1
ACGT--ACGTAC-ATA-C-----CC-T-G-GTA-G-T---
>2
AGGTTTACGTAG-ATA-C-----CC-T-G-GTA-G-T---
>3
AGGTACT-CCAC-ATA-C-----CC-T-G-GTA-G-T---
>4
TCGTTCGT-----ATA-C-----CC-T-G-GTA-G-T---
>5
ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---
"""

pynast_test1_input_fasta = """>1 description field
ACCTACGTTAATACCCTGGTAGT
>2
ACCTACGTTAATACCCTGGTAGT
>3
AA
"""

min_length_computation_seqs = """>1 description field
ACCTACGTTAATACCCTGGTA
>2
ACCTACGTTAATACCCTGGTAGT
>3
ACCTACGTTAATACCCTGGTAA
"""

pynast_test1_expected_alignment = """>1 description field 1..23
ACCTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---
>2 1..23
ACCTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---
"""

pynast_test1_expected_failure = """>3
AA
"""

# run unit tests if run from command-line
if __name__ == '__main__':
    main()
