#!/usr/bin/env python

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main

from itertools import izip
from qiime.conditional_uncovered_probability import *
from cogent.maths.stats.alpha_diversity import diversity

class CUPtests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Set up shared variables"""

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

    def test_lladser_ci(self):
        """lladser_ci returns correct conf interval."""

        f=10
        t=10
        r=4
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f)
        self.assertFloatEqual(obs_low, 0.0806026244)
        self.assertFloatEqual(obs_high, 0.806026244)

        r=20
        t=100
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f)
        self.assertFloatEqual(obs_low, 0.02787923964)
        self.assertFloatEqual(obs_high, 0.2787923964)


        # make sure we test with each possible alpha
        alpha=0.99
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f, alpha=alpha)
        self.assertFloatEqual(obs_low, 0.03184536992)
        self.assertFloatEqual(obs_high, 0.3184536992)

        alpha=0.9
        r=3
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f, alpha=alpha)
        self.assertFloatEqual(obs_low, 0.005635941995)
        self.assertFloatEqual(obs_high, 0.05635941995)


        # test other ci_types
        ci_type='ULCU'
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0.01095834700)
        self.assertFloatEqual(obs_high, 0.1095834700)

        alpha=0.95
        t=10
        ci_type='U'
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0)
        self.assertFloatEqual(obs_high, 0.6295793622)

        ci_type='L'
        obs_low, obs_high = lladser_ci(r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertFloatEqual(obs_low, 0.0817691447)
        self.assertFloatEqual(obs_high, 1)

        #Requesting CI for not precomputed values raises Exception
        r=500
        self.assertRaises(ValueError, lladser_ci, r=r, t=t, f=f,
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
       
if __name__ == '__main__':
    main()

