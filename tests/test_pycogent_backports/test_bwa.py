#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from os.path import join, exists
from os import remove, getcwd
from cogent.app.bwa import BWA, BWA_index, BWA_aln, BWA_samse, \
    BWA_sampe, BWA_bwasw, create_bwa_index_from_fasta_file, \
    assign_reads_to_database
from cogent.app.util import get_tmp_filename, ApplicationError

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Adam Robbins-Pianka", "Daniel McDonald", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Production"


class BWAtests(TestCase):

    """Tests for the BWA app controller
    """

    # keeps track of which files are created during the tests so that they
    # can be removed during tearDown
    files_to_remove = []

    def setUp(self):
        """Performs setup for the tests.

        Nothing to set up for these tests.
        """

        pass

    def tearDown(self):
        """Properly and politely terminates the test.

        Removes files created during the tests.
        """

        for f in self.files_to_remove:
            if exists(f):
                remove(f)

    def test_check_arguments(self):
        """Tests the "check_arguments" method of the BWA base class.

        Arguments passed to certain parameters of the various subcommands can
        take only certain values. The check_arguments function enforces these
        constraints. This function ensures that the rules are being enforced
        as expected.
        """

        # set up test parameters
        # should pass
        index_params_is = {'-a': 'is'}
        # should pass
        index_params_bwtsw = {'-a': 'bwtsw'}
        # should fail, -a must be one of "is" or "bwtsw"
        index_params_invalid = {'-a': 'invalid'}
        # should fail, -p must specify a prefix that is an absolute path
        index_params_invalid_prefix = {'-p': 'invalid'}
        # should pass
        index_params_valid_prefix = {'-p': '/prefix'}

        # instantiate objects built from the above parameters
        index_is = BWA_index(params=index_params_is, HALT_EXEC=True)
        index_bwtsw = BWA_index(params=index_params_bwtsw, HALT_EXEC=True)
        index_invalid = BWA_index(params=index_params_invalid, HALT_EXEC=True)
        index_invalid_prefix = BWA_index(params=index_params_invalid_prefix,
                                         HALT_EXEC=True)
        index_valid_prefix = BWA_index(params=index_params_valid_prefix,
                                       HALT_EXEC=True)

        # Should not be allowed
        self.assertRaisesRegexp(ApplicationError, "Invalid argument",
                                index_invalid.check_arguments)
        self.assertRaisesRegexp(ApplicationError, "Invalid argument",
                                index_invalid_prefix.check_arguments)

        # Should execute and not raise any exceptions
        index_is.check_arguments()
        index_bwtsw.check_arguments()
        index_valid_prefix.check_arguments()

        # The rest of the _valid_arguments are for checking is_int and is_float
        # and they all use the same function from the base-class, so testing
        # just one of the subcommands should suffice

        # -n must be a float (expressed either as a float or as a string)
        # -o must be an int (expressed either as an int or as a string)
        # pass, both valid
        aln_params_valid = {'-n': 3.0, '-o': 5, '-f': '/sai_out'}
        # fail, second invalid
        aln_params_invalid1 = {'-n': 3.0, '-o': 'nope', '-f': '/sai_out'}
        # fail, first invalid
        aln_params_invalid2 = {'-n': '3.5.1', '-o': 4, '-f': '/sai_out'}
        # fail, did not specify -f
        aln_params_invalid3 = {'-n': 3.0, '-o': 5}

        # instantiate objects
        aln_valid = BWA_aln(params=aln_params_valid, HALT_EXEC=True)
        aln_invalid1 = BWA_aln(params=aln_params_invalid1, HALT_EXEC=True)
        aln_invalid2 = BWA_aln(params=aln_params_invalid2, HALT_EXEC=True)
        aln_invalid3 = BWA_aln(params=aln_params_invalid3, HALT_EXEC=True)

        test_paths = {'prefix': '/fa_in', 'fastq_in': '/fq_in'}

        # Should Halt Exec (AssertionError) right before execution
        self.assertRaisesRegexp(AssertionError, 'Halted exec', aln_valid,
                                test_paths)
        # also need to make sure the base command is correct
        self.assertIn('; bwa aln -f /sai_out -n 3.0 -o 5 /fa_in /fq_in',
                      aln_valid.BaseCommand)

        # Should fail
        self.assertRaisesRegexp(ApplicationError,
                                "Invalid argument", aln_invalid1,
                                test_paths)

        self.assertRaisesRegexp(ApplicationError,
                                "Invalid argument", aln_invalid2,
                                test_paths)

        self.assertRaisesRegexp(ApplicationError,
                                "Please specify an output file",
                                aln_invalid3, test_paths)

    def test_input_as_dict(self):
        """Tests the input handler (_input_as_dict)

        The input handler should throw exceptions if there are not enough
        arguments, or if there are unrecognized arguments, or if a file path
        appears to be a relative filepath.
        """

        # Arguments for BWA_bwasw, which was chosen since it is the only one
        # that also has an optional argument (optional arguments are denoted
        # by a leading underscore)
        missing = {'prefix': '/fa_in', '_query_fasta_2': '/mate'}
        extra = {'prefix': '/fa_in', 'query_fasta': '/query_fasta',
                 'extra': '/param'}
        rel_fp = {'prefix': 'fa_in', 'query_fasta': '/query_fasta'}
        valid = {'prefix': '/fa_in', 'query_fasta': '/query_fasta'}
        valid_with_mate = {'prefix': '/fa_in', 'query_fasta': '/query_fasta',
                           '_query_fasta_2': '/mate'}

        # instantiate the object
        bwasw = BWA_bwasw(params={'-f': '/sam_out'}, HALT_EXEC=True)

        # should raise ApplicationError for wrong I/O files; failure
        self.assertRaisesRegexp(ApplicationError, "Missing required input",
                                bwasw, missing)
        self.assertRaisesRegexp(ApplicationError, "Invalid input arguments",
                                bwasw, extra)
        self.assertRaisesRegexp(ApplicationError, "Only absolute paths",
                                bwasw, rel_fp)

        # should raise AssertionError (Halt Exec); success
        # tests valid arguments with and without the optional
        # _query_fasta_2 argument
        self.assertRaisesRegexp(AssertionError, 'Halted exec', bwasw, valid)
        self.assertRaisesRegexp(AssertionError, 'Halted exec', bwasw,
                                valid_with_mate)

    def test_get_base_command(self):
        """Tests the function that generates the command string.

        Tests whether an object can be instantiated and then called using
        one set of files, and then another set of files.

        Since the structure of the various sublcasses is consistent, testing
        that the correct command is generated by one of the subclasses should
        suffice here.
        """

        # instantiate one instance
        aln = BWA_aln(params={'-n': 1.0, '-f': '/sai_out'}, HALT_EXEC=True)

        # set up two different sets of files
        first_files = {'prefix': '/fa_in1', 'fastq_in': '/fq_in1'}
        second_files = {'prefix': '/fa_in2', 'fastq_in': '/fq_in2'}

        # make sure both sets run, and that the command appears to be correct
        self.assertRaisesRegexp(AssertionError,
                                'Halted exec', aln, first_files)
        self.assertIn('; bwa aln -f /sai_out -n 1.0 /fa_in1 /fq_in1',
                      aln.BaseCommand)

        self.assertRaisesRegexp(AssertionError, 'Halted exec', aln,
                                second_files)
        self.assertIn('; bwa aln -f /sai_out -n 1.0 /fa_in2 /fq_in2',
                      aln.BaseCommand)

        # instantiate another object, to test that there is no cross-talk
        # between instances with the same baseclass
        aln2 = BWA_aln(params={'-n': 2.5, '-o': 7, '-f': '/sai_out'},
                       HALT_EXEC=True)

        self.assertRaisesRegexp(AssertionError, 'Halted exec', aln2,
                                first_files)
        self.assertIn('; bwa aln -f /sai_out -n 2.5 -o 7 /fa_in1 /fq_in1',
                      aln2.BaseCommand)

    def test_get_result_paths(self):
        """Tests the function that retrieves the result paths.

        aln, sampe, samse, bwasw return only one file.
        BWA_index returns 5 files, and the name depends on whether or not the
        -p option is on or not
        """

        # instantiate objects
        index = BWA_index(params={}, HALT_EXEC=True)
        index2 = BWA_index(params={'-p': '/prefix'}, HALT_EXEC=True)
        aln = BWA_aln(params={'-f': '/sai_out'}, HALT_EXEC=True)
        samse = BWA_samse(params={'-f': '/sam_out'}, HALT_EXEC=True)
        sampe = BWA_sampe(params={'-f': '/sam_out'}, HALT_EXEC=True)
        bwasw = BWA_bwasw(params={'-f': '/sam_out'}, HALT_EXEC=True)

        # pass in the data, and make sure the output paths are as expected.
        # -p is off here
        index_data = {'fasta_in': '/fa_in'}
        results = index._get_result_paths(index_data)
        self.assertEqual(results['.amb'].Path, '/fa_in.amb')
        self.assertEqual(results['.ann'].Path, '/fa_in.ann')
        self.assertEqual(results['.bwt'].Path, '/fa_in.bwt')
        self.assertEqual(results['.pac'].Path, '/fa_in.pac')
        self.assertEqual(results['.sa'].Path, '/fa_in.sa')

        # pass in the data, and make sure the output paths are as expected.
        # -p is on here
        results = index2._get_result_paths(index_data)
        self.assertEqual(results['.amb'].Path, '/prefix.amb')
        self.assertEqual(results['.ann'].Path, '/prefix.ann')
        self.assertEqual(results['.bwt'].Path, '/prefix.bwt')
        self.assertEqual(results['.pac'].Path, '/prefix.pac')
        self.assertEqual(results['.sa'].Path, '/prefix.sa')

        # pass in the data, and make sure the output path is as expected
        aln_data = {'prefix': '/fa_in', 'fastq_in': '/fq_in'}
        results = aln._get_result_paths(aln_data)
        self.assertEqual(results['output'].Path, '/sai_out')

        samse_data = {'prefix': '/fa_in', 'sai_in': '/sai_in',
                      'fastq_in': '/fq_in'}
        results = samse._get_result_paths(samse_data)
        self.assertEqual(results['output'].Path, '/sam_out')

        sampe_data = {'prefix': '/fa_in', 'sai1_in': '/sai1_in',
                      'sai2_in': '/sai2_in', 'fastq1_in': '/fq1_in',
                      'fastq2_in': '/fq2_in'}
        results = sampe._get_result_paths(sampe_data)
        self.assertEqual(results['output'].Path, '/sam_out')

    def test_create_bwa_index_from_fasta_file(self):
        """Test create_bwa_index_from_fasta_file

        Makes sure that the file paths are as expected.
        """

        # get a new temp file for the input fasta
        fasta_in = get_tmp_filename(suffix=".fna")
        # write the test fasta (see end of this file) to the temp file
        fasta = open(fasta_in, 'w')
        fasta.write(test_fasta)
        fasta.close()

        # make sure to remove this fasta file upon tearDown
        self.files_to_remove.append(fasta_in)

        # run the function
        results = create_bwa_index_from_fasta_file(fasta_in, {})

        # for each of the 5 output files (not counting stdout, stderr, and
        # the exitStatus), make sure the file paths are as expcted.
        for filetype, result in results.iteritems():
            if filetype not in ('ExitStatus'):
                # be sure to remove these 5 files
                self.files_to_remove.append(result.name)
            if filetype not in ('StdOut', 'ExitStatus', 'StdErr'):
                self.assertEqual(fasta_in + filetype, result.name)

    def test_assign_reads_to_database(self):
        """Tests for proper failure in assign_reads_to_database
        """

        # sets of params that should cause failure
        no_alg = {}
        wrong_alg = {'algorithm': 'not_an_algorithm'}
        no_aln_params = {'algorithm': 'bwa-short'}

        # dummy files -- checking for failure as expected, so the function
        # won't get as far as actually running the program
        database = '/db'
        query = '/query'
        out = '/sam'

        self.assertRaisesRegexp(ApplicationError,
                                "Must specify which algorithm",
                                assign_reads_to_database, query, database, out,
                                no_alg)

        self.assertRaisesRegexp(ApplicationError, "Unknown algorithm",
                                assign_reads_to_database, query, database, out,
                                wrong_alg)

        self.assertRaisesRegexp(ApplicationError,
                                "aln is an intermediate step",
                                assign_reads_to_database, query, database, out,
                                no_aln_params)

test_fasta = '''>NZ_GG770509_647533119
UACUUGGAGUUUGAUCCUGGCUCAGAACGAACGCUGGCGGCAGGCUUAACACAUGCAAGUCGAGCGAGCGGCAGACGGGUGAGUAACGCGUGGGAACGUACCAUUUGCUACGGAAUAACUCAGGGAAACUUGUGCUAAUACCGUAUGUGGAAAGUCGGCAAAUGAUCGGCCCGCGUUGGAUUAGCUAGUUGGUGGGGUAAAGGCUCACCAAGGCGACGAUCCAUAGCUGGUCUGAGAGGAUGAUCAGCCACACUGGGACUGAGACACGGCCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGGACAAUGGGCGCAAGCCUGAUCCAGCCAUGCCGCGUGAGUGAUGAAGGCCCUAGGGUUGUAAAGCUCUUUCACCGGUGAAGAUGACGGUAACCGGAGAAGAAGCCCCGGCUAACUUCGUGCCAGCAGCCGCGGUAAUACGAAGGGGGCUAGCGUUGUUCGGAUUUACUGGGCGUAAAGCGCACGUAGGCGGACUUUUAAGUCAGGGGUGAAAUCCCGGGGCUCAACCCCGGAACUGCCUUUGAUACUGGAAGUCUUGAGUAUGGUAGAGGUGAGUGGAAUUCCGAGUGUAGAGGUGAAAUUCGUAGAUAUUCGGAGGAACACCAGUGGCGAAGGCGGCUCACUGGACCAACUGACGCUGAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGAAUGUUAGCCGUCGGGGCUUCGGUGGCGCAGCUAACGCAUUAAACAUUCCGCCUGGGGAGUGCGGUCGCAAGAUUAAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAAGCAACGCGCAGAACCUUACCAGCCCUUGACAUCGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUCGCCCUUAGUUGCCAGCAUGGGCACUCUAAGGGGACUGCCGGUGAUAAGCCGGAGGAAGGUGGGGAUGACGUCAAGUCCUCAUGGCCCUUACGGGCUGGGCUACACACGUGCUACAAUGGUGGUCAGUGGGCAGCGAGCACGCGAGUGUGAGCUAAUCUCCGCCAUCUCAGUUCGGAUGCACUCUGCAACUCGAGUGCAGAAGUUGGAAUCGCUAGUAAUCGCGGAUCAGCAUGCCGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUUGGUUUUACCCGAAGGCGCUUGCUAGGCAGGCGACCACGGUAGGGUCAGCGACUGGGGUGAAGUCGUAACAAGGUAGCCGUAGGGGAACCUGCGGCUGGAUCACCUCCUUUCU
>NZ_GG739926_647533195
UAAUGGGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCUACAGGCUUAACACAUGCAAGUCGAGGGACCGGCGCACGGGUGAGUAACGCGUAUCCAACCUUCCCGCGACCAAGGGAUAACCUGCCGAAAGGCAGACUAAUACCUUAUGUCCAAAGUCGGUCACGGAUGGGGAUGCGUCCGAUUAGCUUGUUGGCGGGGCAACGGCCCACCAAGGCAUCGAUCGGUAGGGGUUCUGAGAGGAAGGCCCCCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGAGGAAUAUUGGUCAAUGGGCGGAAGCCUGAACCAGCCAAGUAGCGUGCAGGACGACGGCCUACGGGUUGUAAACUGCUUUUAUGCGGGGAUAUGCAGGUACCGCAUGAAUAAGGACCGGCUAAUUCCGUGCCAGCAGCCGCGGUAAUACGGAAGGUCCGGGCGUUAUCCGGAUUUAUUGGGUUUAAAGGGAGCGCAGGCCGCCGUGCAAGCGUGCCGUGAAAAGCAGCGGCCCAACCGCUGCCCUGCGGCGCGAACUGCUUGGCUUGAGUGCGCCGGAAGCGGGCGGAAUUCGUGGUGUAGCGGUGAAAUGCUUAGAUAUCACGAAGAACCCCGAUUGCGAAGGCAGCCCGCUGUGGCGACUGACGCUGAGGCUCGAAGGUGCGGGUAUCGAACAGGAUUAGAUACCCUGGUAGUCCGCACGGUAAACGAUGGAUACCCGCUGUCCGGCUCUGGGCGGCCAAGCGAAAGCGUUAAGUAUCCCACCUGGGGAGUACGCCGGCAACGGUGAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGAGGAACAUGUGGUUUAAUUCGAUGAUACGCGAGGAACCUUACCCGGGCUUGAAUUGUGAAGGUGCUGCAUGGUUGUCGUCAGCUCGUGCCGUGAGGUGUCGGCUCAAGUGCCAUAACGAGCGCAACCCCUCUCCGCAGUUGCCAUCGGCCGGGCACUCUGCGGACACUGCCGCCGCAAGGUGGAGGAAGGUGGGGAUGACGUCAAAUCAGCACGGCCCUUACGUCCGGGGCCACACACGUGUUACAAUGGCCGGCAGAGGGCUGUCCGCGCGCAAGUGCGGGUGAAUCCCCUCCGGUCCCAGUUCGGAUGGGGUCUGCAACCCGACCCCAGAAGCUGGAUUCGCUAGUAAUCGCGCAUCAGCCAUGGCGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCAAGCCAUGAAAGCCGGGGGUGCCUGAAGUCCGUGUCGGCCUAGGGCAAAACCGGUGAUUGGGGCUAAGUCGUAACAAGGUAGCCGUACCGGAAGGUGCGGCUGGAACACCUCCUUUCU
>NZ_ACIZ01000148_643886127
AAUAUGGAGUUUGAUCCUGGCUCAGGAUGAACGCUGGCGGCGUGCCUAAUACAUGCAAGUCGAACGAGUGGCGGACGGGUGAGUAACACGUGGGUAACCUGCCCUUAAGUGGGGGAUAACAUUUGGAAACAGAUGCUAAUACCGCAUAAAGAAAGUCGCUUUUGGAUGGACCCGCGGCGUAUUAGCUAGUUGGUGAGGUAACGGCUCACCAAGGCAAUGAUACGUAGCCGAACUGAGAGGUUGAUCGGCCACAUUGGGACUGAGACACGGCCCAAACUCCUACGGGAGGCAGCAGUAGGGAAUCUUCCACAAUGGACGCAAGUCUGAUGGAGCAACGCCGCGUGAGUGAAGAAGGCUUUCGGGUCGUAAAACUCUGUUGUUGGAGAAGAUGACGGUAUCCAACCAGAAAGCCACGGCUAACUACGUGCCAGCAGCCGCGGUAAUACGUAGGUGGCAAGCGUUAUCCGGAUUUAUUGGGCGUAAAGCGAGCGCAGGCGGUUUUUUAAGUCUGAUGUGAAAGCCCUCGGCUUAACCGAGGAAGUGCAUCGGAAACUGGGAAACUUGAGUGCAGAAGAGGACAGUGGAACUCCAUGUGUAGCGGUGAAAUGCGUAGAUAUAUGGAAGAACACCAGUGGCGAAGGCGGCUGUCUGGUCUGACUGACGCUGAGGCUCGAAAGCAUGGGUAGCGAACAGGAUUAGAUACCCUGGUAGUCCAUGCCGUAAACGAUGAAUGCUAGGUGUUGGAGCUUCAGUGCCGCAGCUAACGCAUUAAGCAUUCCGCCUGGGGAGUACGACCGCAAGGUUGAAACUCAAAGGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAAGCAACGCGAAGAACCUUACCAGGUCUUGACAUCGACAGGUGGUGCAUGGUUGUCGUCAGCUCGUGUCGUGAGAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUGACUAGUUGCCAGCAUGGGCACUCUAGUAAGACUGCCGGUGACAAACCGGAGGAAGGUGGGGAUGACGUCAAAUCAUCAUGCCCCUUAUGACCUGGGCUACACACGUGCUACAAUGGAUGGCAACGAGUUGCGAGACCGCGAGGUCAAGCUAAUCUCUUCCAUUCUCAGUUCGGAUGUAGGCUGCAACUCGCCUACAGAAGUCGGAAUCGCUAGUAAUCGCGGAUCAGCACGCCGCGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGAGAGUUUGUAACACCCGAAGCCGGUGCGUAGCGAGCCGUCUAAGGUGGGACAAAUGAUUAGGGUGAAGUCGUAACAAGGUAGCCGUAGGAGAACCUGCGGCUGGAUCACCUCCUUUCU'''

if __name__ == "__main__":
    main()
