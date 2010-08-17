#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.plot_rank_abundance_graph import make_sorted_frequencies
from numpy import array

class PlotRankAbundanace(TestCase):
    
    def test_make_sorted_frequencies(self):
        """make_sorted_frequencies transforms and sorts correctly."""
        
        #works on empty
        counts = array([])
        self.assertEqual(make_sorted_frequencies(counts),[])

         #works on zeros
        counts = array([0,0,0,0,0,0])
        self.assertEqual(make_sorted_frequencies(counts),[])
 

        #works on flat data
        counts = array([3,3,3,3,3])
        expected_freqs = [0.2, 0.2, 0.2, 0.2, 0.2] 
        observed_freqs = make_sorted_frequencies(counts)
        self.assertEqual(observed_freqs, expected_freqs)

        #works on real data
        counts = array([1,2,0,1,0,2,4])
        expected_freqs = [0.4, 0.2, 0.2, 0.1, 0.1] 
        observed_freqs = make_sorted_frequencies(counts)
        self.assertEqual(observed_freqs, expected_freqs)
   
    def test_make_sorted_frequencies_abolute(self):
        """make_sorted_frequencies returns correct absolute values."""
        
        #works on empty
        counts = array([])
        self.assertEqual(make_sorted_frequencies(counts,True),[])
  
        #works on zeros
        counts = array([0,0,0,0,0,0])
        self.assertEqual(make_sorted_frequencies(counts,True),[])
        
        #works on flat data
        counts = array([3,3,3,3,3])
        expected_freqs = [3,3,3,3,3] 
        observed_freqs = make_sorted_frequencies(counts,True)
        self.assertEqual(observed_freqs, expected_freqs)

        #works o real data
        counts = array([1,2,0,1,0,2,4])
        expected_freqs = [4,2,2,1,1] 
        observed_freqs = make_sorted_frequencies(counts,True)
        self.assertEqual(observed_freqs, expected_freqs)

if __name__ == "__main__":
    main()
