#!/usr/bin/env python
#file test_alpha_diversity.py
from __future__ import division
from numpy import array, log, sqrt, exp
from math import e
from cogent.util.unit_test import TestCase, main
from qiime.pycogent_backports.alpha_diversity import (expand_counts, counts, 
    observed_species, singles, doubles, osd, margalef, menhinick, dominance, 
    simpson, simpson_reciprocal, reciprocal_simpson, shannon, equitability, 
    berger_parker_d, mcintosh_d, brillouin_d, strong, kempton_taylor_q, 
    fisher_alpha, mcintosh_e, heip_e, simpson_e, robbins, robbins_confidence, 
    chao1_uncorrected, chao1_bias_corrected, chao1, chao1_var, 
    chao1_confidence, ACE, michaelis_menten_fit, diversity, lorenz_curve, 
    lorenz_curve_integrator, gini_index, goods_coverage, 
    lladser_ci, lladser_pe,
    starr, ul_confidence_bounds, upper_confidence_bound, 
    lower_confidence_bound, lladser_ci_from_r, lladser_ci_series, starr_est,
    esty_ci, get_interval_for_r_new_species, lladser_point_estimates)

from qiime.alpha_diversity import (AlphaDiversityCalc, AlphaDiversityCalcs, 
    single_file_cup, get_cup_metric, list_known_metrics)

from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files
from itertools import izip
import qiime.pycogent_backports.alpha_diversity as alph

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight","Justin Kuczynski, William Van Treuren",
                "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class diversity_tests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Set up shared variables"""
        self.TestData = array([0,1,1,4,2,5,2,4,1,2])
        self.NoSingles = array([0,2,2,4,5,0,0,0,0,0])
        self.NoDoubles = array([0,1,1,4,5,0,0,0,0,0])
        # for Gini index
        self.gini_data = array([4.5, 6.7, 3.4, 15., 18., 3.5, 6.7, 14.1])
        self.gini_lorenz_curve_points = \
            [(0.125, 0.047287899860917935),
             (0.25, 0.095966620305980521),
             (0.375, 0.15855354659248957),
             (0.5, 0.2517385257301808),
             (0.625, 0.34492350486787204),
             (0.75, 0.541029207232267),
             (0.875, 0.74965229485396379),
             (1.0, 1.0)]
        # for Manuel's coverage estimations
        self.tmp_file    = get_tmp_filename(tmp_dir = "./", suffix="test_single_file_cup.biom")
        self.tmp_outfile = get_tmp_filename(tmp_dir = "./", suffix="test_single_file_cup.txt")
        self.files_to_remove = []
    def tearDown(self):
        """ Remove tmp files """
        remove_files(self.files_to_remove)
    def test_expand_counts(self):
        """expand_counts should return correct expanded array"""
        c = array([2,0,1,2])
        self.assertEqual(expand_counts(c), array([0,0,2,3,3]))

    def test_counts(self):
        """counts should return correct array"""
        c = array([5,0,1,1,5,5])
        obs = counts(c)
        exp = array([1,2,0,0,0,3])
        self.assertEqual(obs, exp)
        d = array([2,2,1,0])
        obs = counts(d, obs)
        exp = array([2,3,2,0,0,3])
        self.assertEqual(obs, exp)

    def test_singles(self):
        """singles should return correct # of singles"""
        self.assertEqual(singles(self.TestData), 3)
        self.assertEqual(singles(array([0,3,4])), 0)
        self.assertEqual(singles(array([1])), 1)

    def test_doubles(self):
        """doubles should return correct # of doubles"""
        self.assertEqual(doubles(self.TestData), 3)
        self.assertEqual(doubles(array([0,3,4])), 0)
        self.assertEqual(doubles(array([2])), 1)

    def test_osd(self):
        """osd should return correct # of observeds, singles, doubles"""
        self.assertEqual(osd(self.TestData), (9,3,3))

    def test_margalef(self):
        """margalef should match hand-calculated values"""
        self.assertEqual(margalef(self.TestData), 8/log(22))

    def test_menhinick(self):
        """menhinick should match hand-calculated values"""
        self.assertEqual(menhinick(self.TestData), 9/sqrt(22))

    def test_dominance(self):
        """dominance should match hand-calculated values"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(dominance(c), .34)
        d = array([5])
        self.assertEqual(dominance(d), 1)

    def test_simpson(self):
        """simpson should match hand-calculated values"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(simpson(c), .66)
        d = array([5])
        self.assertFloatEqual(simpson(d), 0)

    def test_reciprocal_simpson(self):
        """reciprocal_simpson should match hand-calculated results"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(reciprocal_simpson(c), 1/.66)

    def test_simpson_reciprocal(self):
        """simpson_reciprocal should match 1/D  results"""
        c = array([1,0,2,5,2])
        self.assertFloatEqual(simpson_reciprocal(c), 1./dominance(c))

    def test_shannon(self):
        """shannon should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(shannon(c), 0)
        c = array([5,5])
        self.assertFloatEqual(shannon(c), 1)
        c = array([1,1,1,1,0])
        self.assertEqual(shannon(c), 2)

    def test_equitability(self):
        """equitability should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(equitability(c), 0)
        c = array([5,5])
        self.assertFloatEqual(equitability(c), 1)
        c = array([1,1,1,1,0])
        self.assertEqual(equitability(c), 1)

    def test_berger_parker_d(self):
        """berger-parker_d should match hand-calculated values"""
        c = array([5])
        self.assertFloatEqual(berger_parker_d(c), 1)
        c = array([5,5])
        self.assertFloatEqual(berger_parker_d(c), 0.5)
        c = array([1,1,1,1,0])
        self.assertEqual(berger_parker_d(c), 0.25)

    def test_mcintosh_d(self):
        """mcintosh_d should match hand-calculated values"""
        c = array([1,2,3])
        self.assertFloatEqual(mcintosh_d(c), 0.636061424871458)

    def test_brillouin_d(self):
        """brillouin_d should match hand-calculated values"""
        c = array([1,2,3,1])
        self.assertFloatEqual(brillouin_d(c), 0.86289353018248782)

    def test_strong(self):
        """strong's dominance index should match hand-calculated values"""
        c = array([1,2,3,1])
        self.assertFloatEqual(strong(c), 0.214285714)

    def test_kempton_taylor_q(self):
        """kempton_taylor_q should approximate Magurran 1998 calculation p143"""
        c = array([2,3,3,3,3,3,4,4,4,6,6,7,7,9,9,11,14,15,15,20,29,33,34,
            36,37,53,57,138,146,170])
        self.assertFloatEqual(kempton_taylor_q(c), 14/log(34/4))

    def test_fisher_alpha(self):
        """fisher alpha should match hand-calculated value."""
        c = array([4,3,4,0,1,0,2])
        obs = fisher_alpha(c)
        self.assertFloatEqual(obs, 2.7823795367398798)

    def test_mcintosh_e(self):
        """mcintosh e should match hand-calculated value."""
        c = array([1,2,3,1])
        num = sqrt(15)
        den = sqrt(19)
        exp = num/den
        self.assertEqual(mcintosh_e(c), exp)

    def test_heip_e(self):
        """heip e should match hand-calculated value"""
        c = array([1,2,3,1])
        h = shannon(c, base=e)
        expected = exp(h-1)/3
        self.assertEqual(heip_e(c), expected)

    def test_simpson_e(self):
        """simpson e should match hand-calculated value"""
        # example1 - a totally even community should have simpson_e=1
        c = array([1,1,1,1,1,1,1])
        self.assertEqual(simpson_e(c), 1)
        # example2
        c = array([0,30,25,40,0,0,5]).astype(float)
        freq_c = c/c.sum()
        D = (freq_c**2).sum()
        exp = 1./(D*4)
        obs = simpson_e(c)
        self.assertEqual(obs, exp)
        # example3 - from https://groups.nceas.ucsb.edu/sun/meetings/calculating-evenness-of-habitat-distributions
        c = array([500, 400, 600, 500])
        freq_c = array([.25, .2, .3, .25])
        D = .0625+.04+.09+.0625
        exp = 1./(D*4)
        self.assertEqual(simpson_e(c), exp)

    def test_robbins(self):
        """robbins metric should match hand-calculated value"""
        c = array([1,2,3,0,1])
        r = robbins(c)
        self.assertEqual(r,2./7) 

    def test_robbins_confidence(self):
        """robbins CI should match hand-calculated value"""
        c = array([1,2,3,0,1])
        r = robbins_confidence(c, 0.05)
        n = 7
        s = 2
        k = sqrt(8/0.05)
        self.assertEqual(r, ((s-k)/(n+1), (s+k)/(n+1))) 


    def test_observed_species(self):
        """observed_species should return # observed species"""
        c = array([4,3,4,0,1,0,2])
        obs = observed_species(c)
        exp = 5
        self.assertEqual(obs, exp)
        c = array([0,0,0])
        obs = observed_species(c)
        exp = 0
        self.assertEqual(obs, exp)
        self.assertEqual(observed_species(self.TestData), 9)

    def test_chao1_bias_corrected(self):
        """chao1_bias_corrected should return same result as EstimateS"""
        obs = chao1_bias_corrected(*osd(self.TestData))
        self.assertEqual(obs, 9.75)

    def test_chao1_uncorrected(self):
        """chao1_uncorrected should return same result as EstimateS"""
        obs = chao1_uncorrected(*osd(self.TestData))
        self.assertEqual(obs, 10.5)

    def test_chao1(self):
        """chao1 should use right decision rules"""
        self.assertEqual(chao1(self.TestData), 9.75)
        self.assertEqual(chao1(self.TestData,bias_corrected=False),10.5)
        self.assertEqual(chao1(self.NoSingles), 4)
        self.assertEqual(chao1(self.NoSingles,bias_corrected=False),4)
        self.assertEqual(chao1(self.NoDoubles), 5)
        self.assertEqual(chao1(self.NoDoubles,bias_corrected=False),5)

    def test_chao1_var(self):
        """chao1_var should match observed results from EstimateS"""
        #NOTE: EstimateS reports sd, not var, and rounds to 2 dp
        self.assertFloatEqual(chao1_var(self.TestData), 1.42**2, eps=0.01)
        self.assertFloatEqual(chao1_var(self.TestData,bias_corrected=False),\
            2.29**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoSingles), 0.39**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoSingles, \
            bias_corrected=False), 0.39**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoDoubles), 2.17**2, eps=0.01)
        self.assertFloatEqualAbs(chao1_var(self.NoDoubles, \
            bias_corrected=False), 2.17**2, eps=0.01)

    def test_chao1_confidence(self):
        """chao1_confidence should match observed results from EstimateS""" 
        #NOTE: EstimateS rounds to 2 dp
        self.assertFloatEqual(chao1_confidence(self.TestData), (9.07,17.45), \
            eps=0.01)
        self.assertFloatEqual(chao1_confidence(self.TestData, \
            bias_corrected=False), (9.17,21.89), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoSingles),\
            (4, 4.95), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoSingles, \
            bias_corrected=False), (4,4.95), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoDoubles), \
            (4.08,17.27), eps=0.01)
        self.assertFloatEqualAbs(chao1_confidence(self.NoDoubles, \
            bias_corrected=False), (4.08,17.27), eps=0.01)
    
    def test_ACE(self):
        """ACE should match values calculated by hand""" 
        self.assertFloatEqual(ACE(array([2,0])), 1.0, eps=0.001)
        # next: just returns the number of species when all are abundant
        self.assertFloatEqual(ACE(array([12,0,9])), 2.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,2,8])), 3.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,2,1])), 4.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,1,2,1])), 7.0, eps=0.001)
        self.assertFloatEqual(ACE(array([12,3,2,1])), 4.6, eps=0.001)
        self.assertFloatEqual(ACE(array([12,3,6,1,10])), 5.62749672, eps=0.001)

    def test_michaelis_menten_fit(self):
        """ michaelis_menten_fit should match hand values in limiting cases"""
        res = michaelis_menten_fit([22])
        self.assertFloatEqual(res,1.0,eps=.01)
        res =  michaelis_menten_fit([42])
        self.assertFloatEqual(res,1.0,eps=.01)
        res =  michaelis_menten_fit([34],num_repeats=3,params_guess=[13,13])
        self.assertFloatEqual(res,1.0,eps=.01)
        res =  michaelis_menten_fit([70,70],num_repeats=5)
        self.assertFloatEqual(res,2.0,eps=.01)

    def test_lorenz_curve(self):
        """Tests that Lorenz curve points are correctly calculated from data."""
        self.assertEqual(self.gini_lorenz_curve_points, lorenz_curve(self.gini_data))

        # check errors on negative data
        self.assertRaises(ValueError, lorenz_curve, [1.0, -3.1, 4.5])

    def test_lorenz_curve_integrator(self):
        """Tests both methods of integration work correctly."""
        expected_trapezoids = 0.33614394993045893
        expected_rectangles = 0.39864394993045893

        self.assertAlmostEqual(expected_trapezoids, 
            lorenz_curve_integrator(self.gini_lorenz_curve_points, 'trapezoids'))

        self.assertAlmostEqual(expected_rectangles, 
            lorenz_curve_integrator(self.gini_lorenz_curve_points, 'rectangles'))

    def test_gini_index(self):
        """Test Gini index is correctly calculated."""
        expected_trapezoids = 0.32771210013908214
        expected_rectangles = 0.20271210013908214

        self.assertAlmostEqual(expected_trapezoids, 
            gini_index(self.gini_data, 'trapezoids'))

        self.assertAlmostEqual(expected_rectangles, 
            gini_index(self.gini_data, 'rectangles'))

    def test_single_file_cup(self):
        """single_file_cup returns matrix with estimates"""
        # Test using a string as metrics
        # convert_biom using otu_table w/o leading #
        bt_string = '{"rows": [{"id": "1", "metadata": null}, {"id": "2",\
 "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata":\
 null}, {"id": "5", "metadata": null}], "format": "Biological Observation\
 Matrix 0.9.1-dev", "data": [[0, 0, 3.0], [0, 1, 4.0], [1, 0, 2.0],\
 [1, 1, 5.0], [2, 0, 1.0], [2, 1, 2.0], [3, 1, 4.0], [4, 0, 1.0]], "columns":\
 [{"id": "S1", "metadata": null}, {"id": "S2", "metadata": null}],\
 "generated_by": "BIOM-Format 0.9.1-dev", "matrix_type": "sparse", "shape":\
 [5, 2], "format_url": "http://biom-format.org", "date":\
 "2012-05-04T09:28:28.247809", "type": "OTU table", "id": null,\
 "matrix_element_type": "float"}'

        fh = open(self.tmp_file,"w")
        fh.write(bt_string)
        fh.close()
        self.files_to_remove.append(self.tmp_file)
        self.files_to_remove.append(self.tmp_outfile)

        # Not much testing here, just make sure we get back a (formatted)
        # matrix with the right dimensions
        single_file_cup(self.tmp_file, 'lladser_pe,lladser_ci',
                        self.tmp_outfile, r=4, alpha=0.95, f=10, ci_type="ULCL")
        observed = open(self.tmp_outfile,"U").readlines()
        self.assertEqual(len(observed), 3)
        self.assertEqual(len(observed[1].split('\t')), 4)

        # Test using a list as metrics
        # convert_biom using otu_table w/o leading #
        bt_string = '{"rows": [{"id": "1", "metadata": null}], "format":\
 "Biological Observation Matrix 0.9.1-dev", "data": [[0, 0, 3.0]], "columns":\
 [{"id": "S1", "metadata": null}], "generated_by": "BIOM-Format 0.9.1-dev",\
 "matrix_type": "sparse", "shape": [1, 1], "format_url":\
 "http://biom-format.org", "date": "2012-05-04T09:36:57.500673", "type":\
 "OTU table", "id": null, "matrix_element_type": "float"}'

        fh = open(self.tmp_file,"w")
        fh.write(bt_string)
        fh.close()

        single_file_cup(self.tmp_file, ['lladser_pe','lladser_ci'],
            self.tmp_outfile, r=4, alpha=0.95, f=10, ci_type="ULCL")
        observed = open(self.tmp_outfile,"U").readlines()
        expected=["\tlladser_pe\tlladser_lower_bound\tlladser_upper_bound\n",
                  "S1\tNaN\tNaN\tNaN"]
        self.assertEqual(observed, expected)

    def test_lladser_point_estimates(self):
        """lladser_point_estimates calculates correct estimates"""
        
        s = [5,1,5,1,2,3,1,5,3,2,5,3]
        r=3
        observed = list(lladser_point_estimates(s,r))
        self.assertEqual(len(observed), 3)
        
        for k in (range(3)):
            x = observed[k]
            t = x[2]
            self.assertEqual(x[0],  (r-1)/t)

        
        #Estimator has variance of (1-p)^2/(r-2),
        # which for r=7 and p=0.5 is 0.05     
        seq="WBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBW"
        reps=1000
        sum=0
        for i in range(reps):
            p,_,_ = list(lladser_point_estimates(seq, r=7))[0]
            sum +=p
        self.assertTrue(0.45 < sum/reps and sum/reps <0.55)
        
    def test_get_interval_for_r_new_species(self):
        """get_interval_for_r_new_species should return the right intervals"""
        
        s = [5,1,5,1,2,3,1,5,3,2,5]
        expected = [(3,set([5]),4,0),
                    (4,set([5,1]),6,1),
                    (4,set([5,1,2]),9,4)]
        for x,y in izip (get_interval_for_r_new_species(s,2), expected):
            self.assertEqual(y,x)

        s = [5,5,5,5,5]
        # never saw new one 
        self.assertEqual(list(get_interval_for_r_new_species(s,2)), [])

    def test_lladser_ci_series_exact(self):
        """lladser_ci_series returns series of predictions"""

        #Values are from Manuel's email of 9/11/09
        #have seen RWB
        urn_1 = 'RWBWWBWRRWRYWRPPZ'
        results = list(lladser_ci_series(urn_1, r=4))
        self.assertEqual(len(results), 3)

    def test_lladser_ci_series_random(self):
        """lladser_ci_series' interval contain true prob with expected alpha."""
                
        seq="WBWBWBWBWBWB"
        observations=[]
        alpha=0.95
        reps = 1000
        for i in range(reps):
            obs = list(lladser_ci_series(seq, r=4, alpha=alpha))[0]
            observations.append(obs)
        tps = filter (lambda (a,b): a < 0.5 and 0.5 < b, observations)
        self.assertTrue(len(tps) >= alpha*reps ) #100%-95%

    def test_lladser_ci_from_r(self):
        """lladser_ci_from_r returns correct conf interval."""

        f=10
        t=10
        r=4
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f)
        self.assertFloatEqual(obs_low, 0.0806026244)
        self.assertFloatEqual(obs_high, 0.806026244)

        r=20
        t=100
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f)
        self.assertFloatEqual(obs_low, 0.02787923964)
        self.assertFloatEqual(obs_high, 0.2787923964)


        # make sure we test with each possible alpha
        alpha=0.99
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        self.assertFloatEqual(obs_low, 0.03184536992)
        self.assertFloatEqual(obs_high, 0.3184536992)

        alpha=0.9
        r=3
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        self.assertFloatEqual(obs_low, 0.005635941995)
        self.assertFloatEqual(obs_high, 0.05635941995)


        # test other ci_types
        ci_type='ULCU'
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0.01095834700)
        self.assertFloatEqual(obs_high, 0.1095834700)

        alpha=0.95
        t=10
        ci_type='U'
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0)
        self.assertFloatEqual(obs_high, 0.6295793622)

        ci_type='L'
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0.0817691447)
        self.assertFloatEqual(obs_high, 1)

        #Requesting CI for not precomputed values raises Exception
        r=500
        self.assertRaises(ValueError, lladser_ci_from_r, r=r, t=t, f=f,
                          alpha=alpha, ci_type=ci_type)

    def test_esty_ci(self):
        """esty_ci computes correct confidence intervals."""

        data = [1,1,2,1,1,3,2,1,3,4]

        (observed_upper, observed_lower) = \
                      zip(*diversity(data, f=esty_ci, step=1))

        expected_upper = [1, 1.38590382, 1.40020259, 0.67434465, 0.55060902,
                          0.71052858, 0.61613483, 0.54041008, 0.43554755,
                          0.53385652]
        expected_lower= [1, -1.38590382, -0.73353593, -0.17434465, -0.15060902,
                         -0.04386191, -0.33042054, -0.29041008, -0.43554755,
                         -0.33385652]

        self.assertFloatEqual(observed_upper, expected_upper)
        self.assertFloatEqual(observed_lower, expected_lower)

    def test_goods_coverage(self):
        counts = [1]*75 + [2,2,2,2,2,2,3,4,4]
        res = goods_coverage(counts)
        self.assertFloatEqual(res, 0.23469387755)

if __name__ == '__main__':
    main()
