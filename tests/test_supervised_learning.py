#!/usr/bin/env python

"""Tests of code for assigning taxonomy"""

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Dan Knights"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Release"


from os import remove, system, mkdir
from shutil import rmtree
from os.path import join, exists
from tempfile import NamedTemporaryFile, mkdtemp
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationError
from qiime.util import get_tmp_filename

from cogent.util.misc import remove_files
from qiime.supervised_learning import RSupervisedLearner,\
    RSupervisedLearnerFilter, R_format_table
from numpy import array

def is_float(input_string):
    """True if string can be cast as a float"""
    try:
        float(input_string)
        return True
    except ValueError:
        return False

class RSupervisedLearnerTests(TestCase):
    """Tests of the RSupervisedLearner class"""
    
    def setUp(self):
        
        # Temporary input file
        self.tmp_otu_filepath = get_tmp_filename(
            prefix='R_test_otu_table_',
            suffix='.txt'
            )
        seq_file = open(self.tmp_otu_filepath, 'w')
        seq_file.write(test_otu_table)
        seq_file.close()

        self.tmp_map_filepath = get_tmp_filename(
            prefix='R_test_map_',
            suffix='.txt'
            )
        seq_file = open(self.tmp_map_filepath, 'w')
        seq_file.write(test_map)
        seq_file.close()


        self.files_to_remove = \
         [self.tmp_otu_filepath, self.tmp_map_filepath]
        self.dirs_to_remove = []
   
        
    def tearDown(self):
        remove_files(set(self.files_to_remove))
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_R_format_table(self):
        """Correctly formats otu table and mapping file for R
        """
        # expected value has comment lines, but no comment char in header
        exp = test_otu_table.split('\n')
        exp[1] = exp[1][1:] # remove '#' from header
        converted = R_format_table(self.tmp_otu_filepath,
            write_to_file=False)
        self.assertEqual(converted, exp)

        # mapping file
        exp = test_map.split('\n')
        exp[0] = exp[0][1:] # remove '#' from header
        converted = R_format_table(self.tmp_map_filepath,
            write_to_file=False)
        self.assertEqual(converted, exp)
    
    def test_RSupervisedLearner(self):
        """Verify that classification algorithms work through R app controller.
        """

        # Prep input files in R format
        otu_fp, map_fp, output_dir = self.prep_learning_input_files()

        # test random forests
        # Temporary param file
        params = ['params$ntree=100','params$seed=0']
        tmp_param_fp = get_tmp_filename(prefix='params_', suffix='.txt')
        param_file = open(tmp_param_fp, 'w')
        param_file.write('\n'.join(params))
        param_file.close()

        self.files_to_remove.append(tmp_param_fp)
        mkdir(join(output_dir, 'random_forest'))
        results = RSupervisedLearner()(
            otu_fp, map_fp,
            'Individual', ['random_forest'], output_dir,
            param_file=tmp_param_fp)
        self.verify_random_forests(results)

    def test_RSupervisedLearnerFilter(self):
        """Verify that classification with pre-filter works through R app controller.
           
           Note: the elastic net code below is commented-out. This classifier
           is likely to work on all data sets, but is not yet enabled
           pending further validation. See doc string for "verify_elastic_net()
           for more information.
        """

        # Prep input files in R format
        otu_fp, map_fp, output_dir = self.prep_learning_input_files()

        # test random forests
        # Temporary param file
        params = ['params$ntree=100','params$seed=0']
        tmp_param_fp = get_tmp_filename(prefix='params_', suffix='.txt')
        param_file = open(tmp_param_fp, 'w')
        param_file.write('\n'.join(params))
        param_file.close()

        self.files_to_remove.append(tmp_param_fp)
        mkdir(join(output_dir, 'random_forest'))
        results = RSupervisedLearnerFilter()(
            otu_fp, map_fp,
            'Individual', ['random_forest'], output_dir,
            param_file=tmp_param_fp, filter='BSSWSS', filter_min=5,
            filter_max=15,filter_step=5,filter_reps=2)
        self.verify_random_forests_filter(results)

    def test_param_file(self):
        """Verify that R correctly reads in params file.        
        """
        otu_fp, map_fp, output_dir = self.prep_learning_input_files()
        mkdir(join(output_dir, 'random_forest'))

        # Temporary param file
        params = ['params$ntree=100','params$seed=0']
        tmp_param_fp = get_tmp_filename(
                prefix='params_',
                suffix='.txt'
            )
        param_file = open(tmp_param_fp, 'w')
        param_file.write('\n'.join(params))
        param_file.close()

        self.files_to_remove.append(tmp_param_fp)

        # test random forests
        results = RSupervisedLearner()(
            otu_fp, map_fp,
            'Individual', ['random_forest'], output_dir, param_file=tmp_param_fp)
        obs_params = results['random_forest']['params'].readlines()
        obs_param_dict = dict([l.strip().split('\t') for l in obs_params[2:]])
        self.assertEqual(obs_param_dict['ntree:'],'100')

    
    def prep_learning_input_files(self):
        # temporarily reformat otu_table, map_file, put in tmp files
        otu_fp = R_format_table(self.tmp_otu_filepath)
        map_fp = R_format_table(self.tmp_map_filepath)

        # make tmp output dir
        output_dir = mkdtemp()
        self.files_to_remove += \
            [otu_fp, map_fp]
        self.dirs_to_remove.append(output_dir)
        return otu_fp, map_fp, output_dir

    def verify_features_format(self, results):
        features_output = results['features'].readlines()

        # ensure that at least one feature is listed (skip header and comment)
        num_features_returned = len(features_output) - 2
        self.assertGreaterThan(num_features_returned, 0)

        # ensure that each line has two elements, and that the first one
        # is the name of one of the OTUs, the second is a float
        for line in features_output[2:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 2)
            self.assertEqual(words[0] in test_OTU_IDs, True)
            self.assertEqual(is_float(words[1]), True)

    def verify_cv_probabilities_format(self,results):
        # verify FORMAT of cross-validation probabilities file
        probabilities_output = results['cv_probabilities'].readlines()

        # ensure that all input samples were predicted
        num_samples_returned = len(probabilities_output)-1
        self.assertEqual(num_samples_returned, len(test_sample_IDs))

        # ensure that each line has five elements, and that the first one
        # is the name of one of the samples, the others are floats
        for line in probabilities_output[1:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 4)
            self.assertEqual(words[0] in test_sample_IDs, True)            
            for word in words[1:3]:
                self.assertEqual(is_float(word),True)

    def verify_mislabeling_format(self,results):
        # verify FORMAT of mislabeling predictions file
        mislabeling_output = results['mislabeling'].readlines()

        # ensure that all input samples were predicted
        num_samples_returned = len(mislabeling_output)-1
        self.assertEqual(num_samples_returned, len(test_sample_IDs))

        # ensure that each line has five elements, and that the first one
        # is the name of one of the samples, the others are floats
        exp = "SampleID\tP(alleged label)\tP(second best)\tP(alleged label)-P(second best)"
        self.assertEqual(mislabeling_output[0].strip(), exp)
        for line in mislabeling_output[1:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 4)
            self.assertEqual(words[0] in test_sample_IDs, True)            
            for word in words[1:3]:
                self.assertEqual(is_float(word),True)

    def verify_filter_features_format(self, results):
        features_output = results['filter_features'].readlines()

        # ensure that at least one feature is listed (skip header and comment)
        num_features_returned = len(features_output) - 1
        self.assertGreaterThan(num_features_returned, 0)

        self.assertEqual(features_output[0].strip(), 'OTU ID')
        # ensure that each line contains an OTU
        for line in features_output[1:]:
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 1)
            self.assertEqual(words[0] in test_OTU_IDs, True)

    def verify_filter_error_format(self, errors_output):

        # ensure that at least one feature is listed (skip header and comment)
        self.assertEqual(len(errors_output)-1, 3)
        
        # check first line
        exp = "Number of features\tMean error\tStandard error"
        self.assertEqual(errors_output[0].strip(), exp)
        
        # ensure that each line has three elements,
        # the number of features, the error, and the std. error
        nfeatures = ['5','10','15']
        for i,line in enumerate(errors_output[1:]):
            words = line.strip().split('\t')
            line_length = len(words)
            self.assertEqual(line_length, 3)
            self.assertEqual(words[0], nfeatures[i])
            self.assertEqual(is_float(words[1]), True)
            self.assertEqual(is_float(words[2]), True)

    def verify_random_forests(self, results):
        # verify FORMAT of features file, since method is stochastic and 
        # external package is subject to change without notice
        self.verify_features_format(results['random_forest'])
        self.verify_mislabeling_format(results['random_forest'])
        self.verify_cv_probabilities_format(results['random_forest'])

        # verify parameters file
        parameters_output = results['random_forest']['params'].readlines()
        exp = ['# values of all non-default parameters\n',
                '# method was "random_forest"\n',
                'ntree:\t100\n',
                'seed:\t0\n']
        self.assertEqual(parameters_output,exp)
        
        # verify summary file (except don't explicitly test error value)
        summary_output = results['random_forest']['summary'].readlines()
        # check generalization error as a float
        exp_error_message = 'Estimated generalization error'
        result_error = summary_output[0].strip().split(' = ')
        self.assertEqual(result_error[0],exp_error_message)
        self.assertEqual(is_float(result_error[1]), True)
        # make sure error is between 0 and 1
        assert(float(result_error[1]) >= 0)
        assert(float(result_error[1]) <= 1)        
        exp = 'Error estimation method = Out-of-bag prediction of training data\n'
        self.assertEqual(summary_output[1],exp)
        exp = 'Number of features used'
        words = summary_output[2].strip().split(' = ')
        self.assertEqual(words[0], exp)
        self.assertEqual(is_float(words[1]), True)       

    def verify_otu_subset_format(self, results, errors_output):
        output = results['otu_subset'].readlines()
        nfeatures = array([int(line.split('\t')[0]) for line in errors_output[1:]])
        errs = array([float(line.split('\t')[1]) for line in errors_output[1:]])
        best_n = nfeatures[errs.argmin()]
        
        # ensure that the correct number of features is listed
        num_features_returned = len(output) - 2
        self.assertEqual(num_features_returned, nfeatures[errs.argmin()])
        
        # check header format (and initial comment line)
        self.assertEqual(output[0][0], '#')
        self.assertEqual(output[1].strip().split('\t')[0], 'OTU ID')
        self.assertEqual(output[1].strip().split('\t')[-1], 'Consensus Lineage')
        self.assertEqual(output[1].strip().split('\t')[1:-1],
            test_sample_IDs)
        
        # ensure that each line is extracted from the OTU table
        # make dict of test_otu_table
        test_otu_lines = test_otu_table.split('\n')[2:]
        test_otu_lineages = {}
        for line in test_otu_lines:
            words = line.strip().split('\t')
            otuid = words[0]
            lineage = words[-1]
            test_otu_lineages[otuid] = lineage
            
        for line in output[2:]:
            words = line.strip().split('\t')            
            otuid = words[0]
            
            self.assertEqual(otuid in test_OTU_IDs, True)
            self.assertEqual(words[-1], test_otu_lineages[otuid])
            for count in words[1:-1]:
                self.assertFloatEqual(is_float(count), True)


    def verify_random_forests_filter(self, results):
        # verify FORMAT of filter_features file, since method is stochastic and 
        # external package is subject to change without notice
        self.verify_filter_features_format(results['random_forest'])
        errors_output = results['random_forest']['filter_errors'].readlines()
        self.verify_filter_error_format(errors_output)
        self.verify_otu_subset_format(results['random_forest'], errors_output)

        # verify parameters file
        parameters_output = results['random_forest']['params'].readlines()
        exp = ['# values of all non-default parameters\n',
                '# method was "random_forest"\n',
                '# filter was "BSSWSS"\n',
                'ntree:\t100\n',
                'seed:\t0\n']
        self.assertEqual(parameters_output,exp)

        # verify summary file (except don't explicitly test error value)
        summary_output = results['random_forest']['filter_summary'].readlines()
        # check generalization error as a float
        exp_error_message = 'Estimated generalization error'
        result_error = summary_output[0].strip().split(' = ')
        self.assertEqual(result_error[0],exp_error_message)
        self.assertEqual(is_float(result_error[1]), True)
        # make sure error is between 0 and 1        
        assert(float(result_error[1]) >= 0)
        assert(float(result_error[1]) <= 1)       
        exp = 'Error estimation method = Out-of-bag prediction of training data\n'
        self.assertEqual(summary_output[1],exp)
        exp = 'Optimal feature subset size'
        words = summary_output[2].strip().split(' = ')
        self.assertEqual(words[0], exp)
        self.assertEqual(is_float(words[1]), True)       


test_sample_IDs = ['S1RingL', 'S1keyM', 'S1keySpace', 'S1IndexL', 'S1keyK', 'S1ThumbR', 'S1keyV', 'S1IndexR', 'S1keyA', 'S1RingR', 'S1MiddleR', 'S1keyD', 'S2keySpace', 'S2keyJ', 'S2keyLeftShift', 'S2keyN', 'S2keyZ', 'S2IndexL', 'S2keyA', 'S2PinkyL', 'S2keyK', 'S2keyRightShift', 'S2keyM', 'S2keyI', 'S2PinkyR', 'S3keySpace', 'S3keyEnter', 'S3keyS', 'S3IndexR', 'S3ThumbR', 'S3MiddleR', 'S3keyY', 'S3ThumbL', 'S3keyF', 'S3IndexL', 'S3keyW', 'S3keyQ', 'S3keyL']
test_OTU_IDs = ['88', '131', '144', '158', '193', '225', '260', '588', '634', '721', '821', '843', '883', '891', '976', '979', '983', '1035', '1088', '1156', '1287', '1314', '1351', '1373', '1487', '1582', '1591', '1784', '1848', '1886', '2007', '2059', '2096', '2187', '2218', '2270', '2328', '2360', '2366', '2407', '2519', '2526', '2810', '2915', '2932', '2956', '3006', '3060', '3108', '3127']

test_otu_table = """#Full OTU Counts
#OTU ID	S1RingL	S1keyM	S1keySpace	S1IndexL	S1keyK	S1ThumbR	S1keyV	S1IndexR	S1keyA	S1RingR	S1MiddleR	S1keyD	S2keySpace	S2keyJ	S2keyLeftShift	S2keyN	S2keyZ	S2IndexL	S2keyA	S2PinkyL	S2keyK	S2keyRightShift	S2keyM	S2keyI	S2PinkyR	S3keySpace	S3keyEnter	S3keyS	S3IndexR	S3ThumbR	S3MiddleR	S3keyY	S3ThumbL	S3keyF	S3IndexL	S3keyW	S3keyQ	S3keyL	Consensus Lineage
88	0	0	0	0	0	0	0	0	0	0	0	0	0	1	4	2	0	3	4	3	0	0	1	5	3	3	1	1	1	1	2	1	0	0	0	1	0	1	Aa;other
131	0	0	0	0	0	0	0	0	0	1	0	0	1	1	1	0	0	1	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Ab;other
144	0	0	0	0	0	0	0	0	0	0	0	0	2	1	4	0	2	4	6	7	0	0	0	2	4	0	0	0	0	0	0	0	0	0	0	0	0	1	Ac;other
158	0	0	0	0	0	0	0	0	0	0	0	0	2	1	3	5	0	1	7	6	1	1	0	1	3	0	1	1	0	0	0	0	0	0	0	1	0	0	Ad;other
193	2	0	0	0	2	0	0	0	0	0	0	0	3	8	0	1	2	2	1	0	1	3	3	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	Ae;other
225	0	0	0	0	1	0	0	0	3	0	0	2	1	1	0	0	0	0	0	0	0	1	0	1	0	0	2	2	0	3	0	2	0	0	1	0	3	5	Af;other
260	1	0	0	1	0	2	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	Ag;other
588	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	1	0	0	1	0	1	1	0	0	0	0	1	Ah;other
634	0	0	0	0	0	0	0	1	0	0	0	0	3	0	0	2	0	1	1	0	2	2	0	2	0	0	0	0	1	1	0	0	1	0	1	0	1	0	Ai;other
721	0	0	0	0	0	0	0	0	0	0	0	0	0	1	2	0	1	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Aj;other
821	5	0	0	1	1	1	1	2	0	3	1	1	6	10	4	5	2	1	6	2	3	3	10	8	3	0	0	1	0	0	0	2	1	1	2	0	1	0	Ak;other
843	1	0	1	3	0	0	0	0	0	0	2	0	2	4	0	2	1	3	4	0	1	2	0	2	3	0	0	0	0	2	0	0	0	0	0	0	0	0	Al;other
883	4	7	3	1	1	3	5	1	1	2	0	5	0	2	3	1	0	1	1	3	14	10	6	3	2	0	0	0	0	0	1	0	0	0	3	0	0	1	Am;other
891	0	1	5	0	0	0	3	0	1	0	0	1	5	3	10	3	5	12	10	8	2	0	3	3	9	0	2	3	8	2	8	2	1	1	1	3	4	2	An;other
976	0	4	0	0	8	0	0	0	0	1	0	1	4	16	2	4	1	0	1	0	1	0	0	3	0	7	2	9	2	3	5	17	6	10	4	7	12	6	Ao;other
979	0	0	0	0	0	0	0	0	0	0	0	0	2	2	2	2	6	1	0	1	1	0	0	2	1	0	0	0	0	0	0	0	0	0	0	0	0	0	Ap;other
983	2	1	0	0	0	2	1	0	0	1	1	0	7	0	3	0	0	0	0	0	1	1	0	3	0	0	0	0	0	1	0	0	0	0	0	0	0	0	Aq;other
1035	7	10	4	3	3	19	9	11	10	5	4	13	7	18	7	16	6	21	26	16	8	23	17	9	30	4	7	7	6	6	3	4	2	0	2	10	2	3	Ar;other
1088	0	0	0	0	0	0	0	0	0	0	0	0	2	1	4	1	9	12	7	9	3	2	2	2	9	0	0	0	1	1	2	0	0	0	2	0	0	0	As;other
1156	0	2	2	1	0	1	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	At;other
1287	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	1	1	4	2	2	3	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	Au;other
1314	1	0	0	0	0	1	0	0	2	3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Av;other
1351	0	0	0	0	0	0	0	0	0	0	0	0	2	0	1	2	3	0	0	0	0	0	2	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Aw;other
1373	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	0	1	2	0	1	0	0	1	3	0	Ax;other
1487	0	0	0	0	1	0	1	0	1	0	0	0	9	3	2	3	18	3	3	2	3	5	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	Ay;other
1582	0	0	1	0	0	1	1	1	1	1	1	0	4	3	1	1	1	0	2	0	6	3	3	1	3	0	0	0	1	1	1	0	1	0	1	0	0	1	Az;other
1591	0	0	0	0	0	0	0	0	0	0	0	0	2	3	4	1	0	0	3	1	1	4	1	4	1	0	0	0	1	0	0	0	0	0	0	0	0	0	Ba;other
1784	0	0	0	0	0	0	0	0	0	0	0	0	2	0	1	1	0	2	2	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bb;other
1848	150	149	154	58	156	137	153	163	151	157	176	155	16	17	9	38	6	8	3	14	54	37	52	17	29	137	153	140	110	132	97	128	161	177	143	141	126	123	Bc;other
1886	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	1	1	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bd;other
2007	0	0	0	2	0	0	1	0	0	0	0	0	5	8	5	3	1	0	4	0	8	9	5	2	3	0	1	2	0	0	0	0	1	0	0	0	0	1	Be;other
2059	0	0	1	1	0	2	1	0	0	3	1	0	12	2	2	2	1	0	4	0	5	6	2	4	4	14	0	3	23	18	12	11	4	1	4	1	1	4	Bf;other
2096	0	0	0	0	0	1	0	0	2	1	0	0	5	0	8	3	5	43	13	30	1	1	2	3	20	0	0	0	0	0	0	0	0	0	0	0	1	0	Bg;other
2187	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	1	2	0	0	0	2	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bh;other
2218	0	0	0	1	0	0	0	0	0	0	0	0	2	1	0	3	0	0	1	1	2	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bi;other
2270	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	2	1	0	0	0	0	0	0	0	Bj;other
2328	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bk;other
2360	0	0	0	0	0	0	0	0	0	0	0	0	2	0	1	0	0	0	1	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bl;other
2366	0	1	2	0	0	0	0	0	0	0	0	0	2	8	34	2	2	2	6	35	1	4	4	3	3	1	0	0	0	0	0	1	0	0	0	0	0	0	Bm;other
2407	0	0	0	0	0	0	0	0	0	0	0	0	2	4	1	0	0	7	6	3	2	0	0	1	2	0	0	0	0	0	0	0	0	0	0	0	0	0	Bn;other
2519	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	3	2	1	3	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bo;other
2526	0	0	1	0	10	1	0	0	0	2	0	2	5	3	1	2	0	0	2	1	4	2	5	14	7	1	1	2	4	0	14	0	2	2	0	2	2	5	Bp;other
2810	0	0	2	0	0	1	0	0	0	2	0	0	0	0	6	2	0	9	0	3	0	0	0	2	7	3	1	0	2	0	1	0	2	1	3	1	1	0	Bq;other
2915	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	4	1	2	0	0	0	1	3	0	0	0	0	0	0	0	0	0	0	0	1	0	Br;other
2932	2	3	10	96	2	4	7	8	2	5	1	0	18	13	12	39	13	4	16	2	32	21	32	13	8	1	0	2	0	1	5	0	3	0	3	3	3	0	Bs;other
2956	1	0	1	1	0	0	1	1	0	2	1	0	1	1	0	3	0	0	1	0	0	0	2	0	2	1	0	0	0	0	2	0	0	0	0	1	0	0	Bt;other
3006	0	0	0	0	1	0	0	0	0	0	0	0	2	2	3	3	3	19	3	9	0	1	2	5	4	0	1	0	0	0	0	3	0	0	0	0	0	0	Bu;other
3060	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	7	1	1	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bv;other
3108	0	0	0	0	0	0	0	0	0	0	0	0	3	0	3	0	3	0	0	1	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	Bw;other
3127	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	2	1	2	0	0	0	0	1	0	0	Bx;other"""

test_map = \
"""#SampleID	Digit	Hand	Keysize	Key	Source	Individual	Spacebar	SpacebarSex	Typingdigit	Handkey	HandIndiv	TypingHand	TypingHandSource	TypingHandIndiv	TypingHandSourceIndiv	UseFreq	Description
S1RingL	Ring	Left	NA	NA	Fing	S1	NA	NA	Ri	Left	S1Left	RiLeft	RiLeftFing	RiLeftS1	RiLeftFingS1	NA	S1.ring-left
S1keyM	NA	NA	Small	M	Key	S1	NA	NA	In	Right	S1Right	InRight	InRightKey	InRightS1	InRightKeyS1	2to6	S1.key, m
S1keySpace	NA	NA	Large	Space	Key	S1	Private	Male	Th	Mixed	S1Mixed	ThMixed	ThMixedKey	ThMixedS1	ThMixedKeyS1	NA	S1.key, spacebar
S1IndexL	Index	Left	NA	NA	Fing	S1	NA	NA	In	Left	S1Left	InLeft	InLeftFing	InLeftS1	InLeftFingS1	NA	S1.index-left
S1keyK	NA	NA	Small	K	Key	S1	NA	NA	Mi	Right	S1Right	MiRight	MiRightKey	MiRightS1	MiRightKeyS1	1minus	S1.key, k
S1ThumbR	Thumb	Right	NA	NA	Fing	S1	NA	NA	Th	Right	S1Right	ThRight	ThRightFing	ThRightS1	ThRightFingS1	NA	S1.thumb-right
S1keyV	NA	NA	Small	V	Key	S1	NA	NA	In	Left	S1Left	InLeft	InLeftKey	InLeftS1	InLeftKeyS1	1to2	S1.key, v
S1IndexR	Index	Right	NA	NA	Fing	S1	NA	NA	In	Right	S1Right	InRight	InRightFing	InRightS1	InRightFingS1	NA	S1.index-right
S1keyA	NA	NA	Small	A	Key	S1	NA	NA	Pi	Left	S1Left	PiLeft	PiLeftKey	PiLeftS1	PiLeftKeyS1	7plus	S1.key, a
S1RingR	Ring	Right	NA	NA	Fing	S1	NA	NA	Ri	Right	S1Right	RiRight	RiRightFing	RiRightS1	RiRightFingS1	NA	S1.ring-right
S1MiddleR	Middle	Right	NA	NA	Fing	S1	NA	NA	Mi	Right	S1Right	MiRight	MiRightFing	MiRightS1	MiRightFingS1	NA	S1.middle-right
S1keyD	NA	NA	Small	D	Key	S1	NA	NA	Mi	Left	S1Left	MiLeft	MiLeftKey	MiLeftS1	MiLeftKeyS1	2to6	S1.key, d
S2keySpace	NA	NA	Large	Space	Key	S2	Private	Male	Th	Mixed	S2Mixed	ThMixed	ThMixedKey	ThMixedS2	ThMixedKeyS2	NA	S2.key, spacebar
S2keyJ	NA	NA	Small	J	Key	S2	NA	NA	In	Right	S2Right	InRight	InRightKey	InRightS2	InRightKeyS2	1minus	S2.key, j
S2keyLeftShift	NA	NA	Medium	LeftShift	Key	S2	NA	NA	Pi	Left	S2Left	PiLeft	PiLeftKey	PiLeftS2	PiLeftKeyS2	NA	S2.key, shift-left
S2keyN	NA	NA	Small	N	Key	S2	NA	NA	In	Right	S2Right	InRight	InRightKey	InRightS2	InRightKeyS2	6to7	S2.key, n
S2keyZ	NA	NA	Small	Z	Key	S2	NA	NA	Pi	Left	S2Left	PiLeft	PiLeftKey	PiLeftS2	PiLeftKeyS2	1minus	S2.key, z
S2IndexL	Index	Left	NA	NA	Fing	S2	NA	NA	In	Left	S2Left	InLeft	InLeftFing	InLeftS2	InLeftFingS2	NA	S2.index-left
S2keyA	NA	NA	Small	A	Key	S2	NA	NA	Pi	Left	S2Left	PiLeft	PiLeftKey	PiLeftS2	PiLeftKeyS2	7plus	S2.key, a
S2PinkyL	Pinky	Left	NA	NA	Fing	S2	NA	NA	Pi	Left	S2Left	PiLeft	PiLeftFing	PiLeftS2	PiLeftFingS2	NA	S2.pinky-left
S2keyK	NA	NA	Small	K	Key	S2	NA	NA	Mi	Right	S2Right	MiRight	MiRightKey	MiRightS2	MiRightKeyS2	1minus	S2.key, k
S2keyRightShift	NA	NA	Medium	RightShift	Key	S2	NA	NA	Pi	Right	S2Right	PiRight	PiRightKey	PiRightS2	PiRightKeyS2	NA	S2.key, shift-right
S2keyM	NA	NA	Small	M	Key	S2	NA	NA	In	Right	S2Right	InRight	InRightKey	InRightS2	InRightKeyS2	2to6	S2.key, m
S2keyI	NA	NA	Small	I	Key	S2	NA	NA	Mi	Right	S2Right	MiRight	MiRightKey	MiRightS2	MiRightKeyS2	6to7	S2.key, i
S2PinkyR	Pinky	Right	NA	NA	Fing	S2	NA	NA	Pi	Right	S2Right	PiRight	PiRightFing	PiRightS2	PiRightFingS2	NA	S2.pinky-right
S3keySpace	NA	NA	Large	Space	Key	S3	Private	Male	Th	Mixed	S3Mixed	ThMixed	ThMixedKey	ThMixedS3	ThMixedKeyS3	NA	S3.key, spacebar
S3keyEnter	NA	NA	Medium	Enter	Key	S3	NA	NA	Unk	Right	S3Right	UnkRight	UnkRightKey	UnkRightS3	UnkRightKeyS3	NA	S3.key, enter
S3keyS	NA	NA	Small	S	Key	S3	NA	NA	Ri	Left	S3Left	RiLeft	RiLeftKey	RiLeftS3	RiLeftKeyS3	6to7	S3.key, s
S3IndexR	Index	Right	NA	NA	Fing	S3	NA	NA	In	Right	S3Right	InRight	InRightFing	InRightS3	InRightFingS3	NA	S3.index-right
S3ThumbR	Thumb	Right	NA	NA	Fing	S3	NA	NA	Th	Right	S3Right	ThRight	ThRightFing	ThRightS3	ThRightFingS3	NA	S3.thumb-right
S3MiddleR	Middle	Right	NA	NA	Fing	S3	NA	NA	Mi	Right	S3Right	MiRight	MiRightFing	MiRightS3	MiRightFingS3	NA	S3.middle-right
S3keyY	NA	NA	Small	Y	Key	S3	NA	NA	In	Right	S3Right	InRight	InRightKey	InRightS3	InRightKeyS3	1to2	S3.key, y
S3ThumbL	Thumb	Left	NA	NA	Fing	S3	NA	NA	Th	Left	S3Left	ThLeft	ThLeftFing	ThLeftS3	ThLeftFingS3	NA	S3.thumb-left
S3keyF	NA	NA	Small	F	Key	S3	NA	NA	In	Left	S3Left	InLeft	InLeftKey	InLeftS3	InLeftKeyS3	2to6	S3.key, f
S3IndexL	Index	Left	NA	NA	Fing	S3	NA	NA	In	Left	S3Left	InLeft	InLeftFing	InLeftS3	InLeftFingS3	NA	S3.index-left
S3keyW	NA	NA	Small	W	Key	S3	NA	NA	Ri	Left	S3Left	RiLeft	RiLeftKey	RiLeftS3	RiLeftKeyS3	2to6	S3.key, w
S3keyQ	NA	NA	Small	Q	Key	S3	NA	NA	Pi	Left	S3Left	PiLeft	PiLeftKey	PiLeftS3	PiLeftKeyS3	1minus	S3.key, q
S3keyL	NA	NA	Small	L	Key	S3	NA	NA	Ri	Right	S3Right	RiRight	RiRightKey	RiRightS3	RiRightKeyS3	2to6	S3.key, l"""

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
