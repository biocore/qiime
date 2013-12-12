#!/usr/bin/env python
# File created on 19 Jan 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat_to_dict, parse_mapping_file
from qiime.categorized_dist_scatterplot import get_avg_dists, get_sam_ids
import StringIO
class FunctionTests(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_get_avg_dists(self):
        """get_avg_dists functions as expected """
        dmtx_str = StringIO.StringIO("""\ts1\ts2\ts3
s1\t0\t.5\t.6
s2\t.5\t0\t.7
s3\t.6\t.7\t0.0
""")
        distdict1 = parse_distmat_to_dict(dmtx_str)
        distdict2 = {'s1':{'s2':.5,'s3':.6},'s2':{'s1':.5,'s3':.7},
            's3':{'s2':.7,'s1':.6}}
        state1_samids = ['s1','s2']
        state2_samids = ['s3','s2']
        # note s2 in both
        exp_avgs = [.55, .7]
        obs_avgs = get_avg_dists(state1_samids,state2_samids,distdict1)
        self.assertFloatEqual(exp_avgs, obs_avgs)

    def test_get_sam_ids(self):
        """set of sample ids in get_sam_ids should be correct"""
        map_file = StringIO.StringIO("""#SampleID	Country	AgeYears	Family	AgeCat
    h208A.1	Malawi	0.032854209	h208	Child
    h301A.1	Malawi	0.05	h301	Child
    h301B.1	Malawi	0.05	h301	Child
    USinfTw20.1	USA	0.083333333	USinfTw20	Child
    USinfTw20.2	USA	0.083333333	USinfTw20	Child
    USinfTw1.1	USA	0.083333333	USinfTw1	Child
    h10M	Malawi	26	h10	Adult
    h68M	Malawi	26	h68	Adult
    TS25	USA	26	USts9	Adult
    TS26	USA	26	USts9	Adult""")

        map_data, map_header, comments = parse_mapping_file(map_file)
        colorby = 'Country'
        cat = 'USA'
        primary_state = 'AgeCat:Child'
        ids1, ids2 = get_sam_ids(map_data, map_header, colorby, cat, 
                primary_state, secondary_state=None)
        self.assertEqual(set(ids1),
            set(['USinfTw20.1','USinfTw20.2','USinfTw1.1']))
        self.assertEqual(set(ids2), set(['TS25','TS26']))

if __name__ == "__main__":
    main()