#!/usr/bin/env python
# File created on 19 May 2011

from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase,main
from qiime.parse import parse_mapping_file_to_dict, parse_rarefaction
from qiime.compare_alpha_diversity import (sampleId_pairs,
    compare_alpha_diversities, _correct_compare_alpha_results)
from numpy.random import seed
from numpy import nan

class TopLevelTests(TestCase):
    """Tests of top level functions"""
    
    def setUp(self):
        """define data for tests"""
        # small amount of redundancy here since setUp called at each test, but
        # limited tests means little concern
        self.rarefaction_file = \
            ['\tsequences per sample\titeration\tSam1\tSam2\tSam3\tSam4\tSam5\tSam6',
            'rare480.txt\t480\t0\t2.52800404052\t2.3614611247\t2.59867416108\t3.56970811181\t3.44800265895\t1.9433560517',
            'rare480.txt\t480\t1\t2.06375457238\t3.32293450758\t3.4189896645\t3.35312890712\t3.10763472113\t2.78155253726',
            'rare480.txt\t480\t2\t2.44788730109\t3.42464996459\t2.24541787295\t2.491419231\t2.60106690099\t5.40828403581',
            'rare480.txt\t480\t3\t5.1846120153\t3.67022675065\t1.54879964908\t2.8055801405\t4.3086171269\t3.87761898868',
            'rare910.txt\t910\t0\t2.67580703282\t1.72405794627\t2.15312863498\t2.4300954476\t3.7753658185\t3.36198860355',
            'rare910.txt\t910\t1\t4.10226466956\t2.24587945345\t3.02932964779\t2.98218513619\t3.73316846484\t1.85879566537',
            'rare910.txt\t910\t2\t1.65800670063\t2.42281993323\t3.02400997565\t3.271608097\t2.99265263795\t3.68802382515',
            'rare910.txt\t910\t3\t2.50976021964\t2.43976761056\t3.32119905587\t2.47487750248\t1.901408525\t3.42883742207',
            'rare500.txt\t500\t0\t3.42225118215\tn/a\t4.03758268426\t2.35344629448\t2.26690085385\t1.80164570104',
            'rare850.txt\t850\t0\t4.2389858006\t4.97464230229\t1.53451087057\t3.35785261181\t1.91658777533\t2.32583475424',
            'rare850.txt\t850\t1\t2.81445883827\tn/a\t2.54767461948\t1.38835207925\t3.70018890199\t1.57359105209',
            'rare850.txt\t850\t2\t2.9340493412\t3.95897035158\tn/a\t2.07761860166\t3.42393336685\t2.6927305603']
        self.rarefaction_data = parse_rarefaction(self.rarefaction_file)
        self.mapping_file = \
            ['#SampleID\tDose\tLinkerPrimerSequence\tWeight\tTTD\tDescription',
            '#Comment Line',
            'Sam1\t1xDose\tATCG\tHigh\t31\ts1_desc',
            'Sam2\t1xDose\tACCG\tLow\t67\ts2_desc',
            'Sam3\t2xDose\tACGT\tMed\t21\ts3_desc',
            'Sam4\t2xDose\tAACG\tLow\t55\ts4_desc',
            'Sam5\tControl\tCGTC\tLow\t67\ts5_desc',
            'Sam6\t1xDose\tACCT\tLow\t55\ts6_desc']
        self.mapping_data = parse_mapping_file_to_dict(self.mapping_file)[0]
    
    def test_sampleId_pairs(self):
        """Test that sampleId_pairs returns the correct combos/sampleId's."""
        # expected values
        dose_vps = \
            [('1xDose', '2xDose'), ('1xDose', 'Control'), ('2xDose', 'Control')]
        ttd_vps = \
            [('31', '21'), ('31', '55'), ('31', '67'), ('21', '55'),
             ('21', '67'), ('55', '67')]
        dose_sids = \
            [(['Sam1','Sam2','Sam6'], ['Sam3','Sam4']),
             (['Sam1','Sam2','Sam6'], ['Sam5']),
             (['Sam3','Sam4'], ['Sam5'])]
        ttd_sids = \
            [(['Sam1'], ['Sam3']), 
             (['Sam1'], ['Sam4','Sam6']), 
             (['Sam1'], ['Sam2','Sam5']), 
             (['Sam3'], ['Sam4','Sam6']),
             (['Sam3'], ['Sam2','Sam5']), 
             (['Sam4','Sam6'], ['Sam2','Sam5'])]
        
        # observed values
        obs_dose_sids, obs_dose_vps = sampleId_pairs(self.mapping_data, 
            self.rarefaction_data, 'Dose')
        obs_ttd_sids, obs_ttd_vps = sampleId_pairs(self.mapping_data,
            self.rarefaction_data, 'TTD')
        
        # sort -- order is unimportant and depends on way presented in mf
        self.assertEqual(dose_vps.sort(),obs_dose_vps.sort())
        self.assertEqual(dose_sids.sort(),obs_dose_sids.sort())
        self.assertEqual(ttd_vps.sort(),obs_ttd_vps.sort())
        self.assertEqual(ttd_sids.sort(),obs_ttd_sids.sort())
        
        # check errors when no samples had this category
        self.assertRaises(ValueError, sampleId_pairs, self.mapping_data,
            self.rarefaction_data, 'DNE')

        # check no error if map file has more sampleids than rarefaction data
        superset_mf = \
            ['#SampleID\tDose\tLinkerPrimerSequence\tWeight\tTTD\tDescription',
            '#Comment Line',
            'Sam1\t1xDose\tATCG\tHigh\t31\ts1_desc',
            'Sam2\t1xDose\tACCG\tLow\t67\ts2_desc',
            'Sam3\t2xDose\tACGT\tMed\t21\ts3_desc',
            'Sam4\t2xDose\tAACG\tLow\t55\ts4_desc',
            'Sam5\tControl\tCGTC\tLow\t67\ts5_desc',
            'Sam6\t1xDose\tACCT\tLow\t55\ts6_desc',
            'Sam7\t4xDose\tACCT\tLow\t55\ts7_desc',
            'Sam8\t3xDose\tACCT\tLow\t55\ts8_desc',
            'Sam9\t1xDose\tACCT\tLow\t55\ts9_desc']
        superset_mf = parse_mapping_file_to_dict(superset_mf)[0] #(mf, comments)
        obs_dose_sids, obs_dose_vps = sampleId_pairs(superset_mf, 
            self.rarefaction_data, 'Dose')

        self.assertEqual(dose_vps.sort(),obs_dose_vps.sort())
        self.assertEqual(dose_sids.sort(),obs_dose_sids.sort())


    def test_correct_compare_alpha_results(self):
        """Test that FDR and Bonferroni are applied correctly."""

        input_results = \
            {'1xDose,2xDose': (-1.8939787722170394, 0.14),
             'Control,1xDose': (3.365078231689424, 0.34),
             'Control,2xDose': (0.43262479194397335, 1.0)}
            
        method = 'fdr'
        expected_results = \
            {'1xDose,2xDose': (-1.8939787722170394, 0.42),
             'Control,1xDose': (3.365078231689424, 0.51),
             'Control,2xDose': (0.43262479194397335, 1.0)}
        
        observed_results = _correct_compare_alpha_results(input_results, method)
        # test each key in expected results -- this won't catch if 
        # observed_results has extra entries, but test that via the next call
        for k in expected_results:
            for val0, val1 in zip(expected_results[k],observed_results[k]):
                self.assertAlmostEqual(val0,val1)
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))

        method = 'bonferroni'
        expected_results = \
            {'1xDose,2xDose': (-1.8939787722170394, 0.14*3),
             'Control,1xDose': (3.365078231689424, 1.0), #because maxes at 1
             'Control,2xDose': (0.43262479194397335, 1.0)} #becuase maxes at 1

        observed_results = _correct_compare_alpha_results(input_results, method)
        # test each key in expected results -- this won't catch if 
        # observed_results has extra entries, but test that via the next call
        for k in expected_results:
            for val0, val1 in zip(expected_results[k],observed_results[k]):
                self.assertAlmostEqual(val0,val1)
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))

        method = 'none'
        expected_results = \
            {'1xDose,2xDose': (-1.8939787722170394, 0.14),
             'Control,1xDose': (3.365078231689424, 0.34),
             'Control,2xDose': (0.43262479194397335, 1.0)}
        observed_results = _correct_compare_alpha_results(input_results, method)
        # test each key in expected results -- this won't catch if 
        # observed_results has extra entries, but test that via the next call
        for k in expected_results:
            for val0, val1 in zip(expected_results[k],observed_results[k]):
                self.assertAlmostEqual(val0,val1)
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))

        # check errors if wrong method
        self.assertRaises(ValueError, _correct_compare_alpha_results,
            input_results, 'DNE')
        
    def test_compare_alpha_diversities(self):
        """Tests alpha diversities are correctly calculated."""
        # test 'Dose' at 480 inputs
        category = 'Dose'
        depth = 480
        test_type = 'parametric'
        observed_results = compare_alpha_diversities(self.rarefaction_file,
            self.mapping_file, category=category, depth=depth, 
            test_type=test_type)
        
        # hardcoded order of the terms in the keys otherwise would comps fail
        expected_results = \
            {'Control,2xDose': (1.1746048668554037, 0.44899351189030801),
             '1xDose,2xDose': (1.7650193854830403, 0.17574514418562981),
             'Control,1xDose': (0.43618805086434992, 0.7052689260099092)}
             
        # test each key in expected results -- this won't catch if 
        # observed_results has extra entries, but test that via the next call
        for k in expected_results:
            self.assertEqual(expected_results[k],observed_results[k])
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))

        # test 'Dose' at 480 inputs with nonparametric test
        seed(0) # set the seed to reproduce random MC pvals
        category = 'Dose'
        depth = 480
        test_type = 'nonparametric'
        num_permutations = 100
        observed_results = compare_alpha_diversities(self.rarefaction_file,
            self.mapping_file, category=category, depth=depth, 
            test_type=test_type, num_permutations=num_permutations)

        expected_results = \
            {'Control,2xDose': (1.1746048668554037, 0.63),
             '1xDose,2xDose': (1.7650193854830403, 0.09),
             'Control,1xDose': (0.43618805086434992, 0.76)}
 
        # test each key in expected results -- this won't catch if 
        # observed_results has extra entries, but test that via the next call
        for k in expected_results:
            self.assertEqual(expected_results[k],observed_results[k])
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))

        # test it works with NA values
        # test 'Dose' at 500 inputs with paramteric test
        category = 'Dose'
        depth = 500
        test_type = 'parametric'
        observed_results = compare_alpha_diversities(self.rarefaction_file,
            self.mapping_file, category=category, depth=depth, 
            test_type=test_type)

        expected_results = \
            {'Control,2xDose': (-0.63668873339963239, 0.63906168713487699), 
             '1xDose,2xDose': (None,None), 
             'Control,1xDose': (nan,nan)}

        # the comparison will fail on (nan,nan)==(nan,nan) because nan's don't
        # compare equal. we avoid this check knowing that its producing nans.
        for k in expected_results:
            if k is not 'Control,1xDose':
                for val0, val1 in zip(expected_results[k],observed_results[k]):
                    self.assertEqual(val0,val1)
        self.assertEqual(set(expected_results.keys()),set(observed_results.keys()))


if __name__ == "__main__":
    main()
