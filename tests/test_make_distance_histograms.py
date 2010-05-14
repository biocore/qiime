#!/usr/bin/env python
# file test_make_distance_histograms.py

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
import shutil
from os import mkdir, listdir
from qiime.parse import parse_mapping_file, parse_distmat, group_by_field,\
    group_by_fields
from numpy import array, arange
from qiime.make_distance_histograms import between_sample_distances, \
    within_category_distances_grouped, between_category_distances_grouped, \
    within_category_distances, within_and_between_fields, \
    all_category_distances, draw_all_histograms, get_histogram_scale, \
    draw_histogram, make_nav_html, make_main_html, get_valid_indices, \
    distances_by_groups, write_distance_files, group_distances, \
    monte_carlo_group_distances, permute_for_monte_carlo, \
    _make_histogram_filenames, _make_path, _make_relative_paths, \
    _make_random_filename, _get_script_dir, matplotlib_rgb_color, \
    average_colors, average_all_colors, assign_unassigned_colors,\
    assign_mapped_colors
    
from qiime.colors import data_colors
from collections import defaultdict

class DistanceHistogramsTests(TestCase):
    """Tests for make_distance_histograms.py
    """
    
    def setUp(self):
        """setup data function for DistanceHistogramsTests."""
        self.working_dir = '/tmp/distance_histogram_tests/'
        try:
            mkdir(self.working_dir)
        except OSError: #except already exisits
            pass
        
        self.histogram_dir = self.working_dir+'histograms/'
        try:
            mkdir(self.histogram_dir)
        except OSError: #except already exisits remove it and make a new one
            pass
            
        #Create distance matrix file
        self.dmat_file = self.working_dir+'dmat.txt'
        dmat_out = open(self.dmat_file,'w')
        dmat_out.write(DISTANCE_MATRIX_STRING)
        dmat_out.close()

        self.distance_header, self.dmat = \
            parse_distmat(open(self.dmat_file,'U'))
        
        #Create mapping file
        self.map_file = self.working_dir+'map.txt'
        map_out = open(self.map_file,'w')
        map_out.write(MAPPING_STRING)
        map_out.close()
        
        mapping, header, comments = parse_mapping_file(open(self.map_file,'U'))
        header[0] = '#'+header[0]
        header = [header]
        header.extend(mapping)
        self.mapping=header

        
        #Create prefs file
        self.prefs_file = self.working_dir+'prefs.txt'
        prefs_out = open(self.prefs_file,'w')
        prefs_out.write(str(PREFS))
        prefs_out.close()
        
        #Build single field dict for 'Treatment' field.
        self.single_field_treatment = defaultdict(dict)
        self.treatment_groups = group_by_field(self.mapping, 'Treatment')
        self.single_field_treatment['Treatment'] = \
            distances_by_groups(self.distance_header,self.dmat,\
                self.treatment_groups)
        self.paired_field_treatment = {'Treatment_to_Treatment':[\
            [('Control','Control'),('Fast','Fast'),\
                             array([[0.729,  0.8  ,  0.721, 0.765],
                                    [0.776,  0.744,  0.749, 0.677],
                                    [0.734,  0.777,  0.733, 0.724],
                                    [0.696,  0.675,  0.654, 0.696],
                                    [0.731,  0.758,  0.738, 0.737]])],\
            [('Control','Control'),('Control','Control'),\
                                 array([0.625,  0.623,  0.61 ,  0.577, 0.615,
                                        0.642,  0.673, 0.682,  0.737, 0.704])],\
            [('Fast','Fast'),('Fast','Fast'),\
                             array([0.718,  0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}
        
        self.distances_file = self.working_dir+'distances_out.txt'
        dist_out = open(self.distances_file,'w')
        dist_out.write(DISTANCES_OUT)
        dist_out.close()

    def tearDown(self):
        """clean up after running all tests.
        """
        shutil.rmtree(self.working_dir)
    
    def test_matplotlib_rgb_color(self):
        """matplotlib_rgb_color should correctly convert RGB to decimal.
        """
        test_data = [(255,255,255),(0,0,0),(255,0,100),(100,253,18)]
        exp = [(1.0,1.0,1.0),(0.,0.,0.),(1.,0.,0.39215686274509803),\
            (0.39215686274509803,0.99215686274509807,0.070588235294117646)]
        for t, e in zip(test_data,exp):
            self.assertEqual(matplotlib_rgb_color(t),e)
        
    def test_average_colors(self):
        """average_colors should properly average two RGB colors.
        """
        to_average = [((255,255,255),(0,0,0)),((255,0,100),(100,253,18))]
        exp = [(127.5,127.5,127.5),(177.5,126.5,59.)]
        for t,e in zip(to_average,exp):
            self.assertEqual(average_colors(t[0],t[1]),e)
    
    def test_average_all_colors(self):
        """average_all_colors should properly average all colors.
        """
        to_average = ['Treatment_Control_to_Fast', 'DOB_20070314_to_20061126']
        exp = {'Treatment_Control_to_Fast':(0.,.5,.5),\
            'DOB_20070314_to_20061126':(.5,0.,.5)}
        self.assertEqual(average_all_colors(to_average,FIELD_TO_COLOR_PREFS),\
            exp)
    
    def test_assign_unassigned_colors(self):
        """assign_unassigned_colors should correctly assign unassigned.
        """
        unassigned = ['first','second','third']
        exp = {'first':(0., 0., 1.),\
            'second':(0.50196078431372548, 0.50196078431372548,\
                0.50196078431372548),\
            'third':(0.50196078431372548, 0., 0.50196078431372548)}
        self.assertEqual(assign_unassigned_colors(unassigned),exp)
    
    def test_assign_mapped_colors(self):
        """assign_mapped_colors should correctly assign mapped colors.
        """
        assigned = ['Treatment_Within_Control_Distances','DOB_Within_20070314']
        exp = {'Treatment_Within_Control_Distances':(0.,0.,1.),\
            'DOB_Within_20070314':(1.,0.,0.)}
        self.assertEqual(assign_mapped_colors(assigned,FIELD_TO_COLOR_PREFS),\
            exp)
        
    def test_between_sample_distances(self):
        """between_sample_distances should return correct result.
        """
        exp = {'All_Between_Sample_Distances':\
            [0.625, 0.623, 0.61, 0.577, 0.729,  0.8, 0.721, 0.765, 0.615,\
             0.642, 0.673, 0.776, 0.744, 0.749, 0.677,0.682, 0.737, 0.734,\
             0.777, 0.733, 0.724, 0.704, 0.696, 0.675, 0.654, 0.696, 0.731,\
             0.758, 0.738, 0.737, 0.718, 0.666, 0.727, 0.6, 0.578,0.623]}
        self.assertEqual(between_sample_distances(self.dmat),exp)
    
    def test_within_category_distances_grouped(self):
        """within_category_distances_grouped should return correct result.
        """
        exp = {'Treatment_All_Within_Category_Distances':\
            [0.625, 0.623, 0.61 , 0.577, 0.615, 0.642, 0.673, 0.682, 0.737, \
             0.704 ,0.718,  0.666, 0.727, 0.6, 0.578, 0.623]}
        self.assertEqual(within_category_distances_grouped(\
            self.single_field_treatment),exp)

    
    def test_between_category_distances_grouped(self):
        """between_category_distances_grouped should return correct result.
        """
        exp = {'Treatment_All_Between_Category_Distances':\
            [0.729,  0.8  ,  0.721, 0.765, 0.776,  0.744,  0.749, 0.677,\
             0.734,  0.777,  0.733, 0.724, 0.696,  0.675,  0.654, 0.696,\
             0.731,  0.758,  0.738, 0.737]}
        self.assertEqual(between_category_distances_grouped(\
            self.single_field_treatment),exp)

    def test_within_category_distances(self):
        """within_category_distances should return correct result.
        """
        exp = {'Treatment_Within_Control_Distances':\
                [0.625,  0.623,  0.61 ,  0.577, 0.615,\
                 0.642,  0.673, 0.682,  0.737, 0.704],\
               'Treatment_Within_Fast_Distances':\
                [0.718,  0.666, 0.727, 0.6, 0.578, 0.623]}
        self.assertEqual(within_category_distances(\
            self.single_field_treatment),exp)
    
    def test_within_and_between_fields(self):
        """within_and_between_fields should return correct result.
        """
        exp = {'Within_All_Fields':\
            [0.625, 0.623, 0.61 , 0.577, 0.615, 0.642, 0.673, 0.682, 0.737, \
             0.704 ,0.718,  0.666, 0.727, 0.6, 0.578, 0.623],\
               }
        
        self.assertEqual(within_and_between_fields(\
            self.paired_field_treatment),exp)
            
    
    def test_all_category_distances(self):
        """all_category_distances should return correct result.
        """
        exp = {'Treatment_Control_to_Control':\
                                       [0.625, 0.623, 0.61, 0.577, 0.615,\
                                        0.642, 0.673, 0.682, 0.737, 0.704],\
               'Treatment_Control_to_Fast':[0.729,  0.8  ,  0.721, 0.765,\
                                            0.776,  0.744,  0.749, 0.677,\
                                            0.734,  0.777,  0.733, 0.724,\
                                            0.696,  0.675,  0.654, 0.696,\
                                            0.731,  0.758,  0.738, 0.737],\
               'Treatment_Fast_to_Fast':[0.718,  0.666, 0.727, 0.6, 0.578,\
                                         0.623]
                }
        self.assertEqual(all_category_distances(\
            self.single_field_treatment),exp)
    
    def test_draw_all_histograms(self):
        """draw_all_histograms should return correct result.
        """
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field = self.single_field_treatment, \
                                paired_field = self.paired_field_treatment, \
                                dmat=self.dmat,\
                                histogram_dir = self.histogram_dir,\
                                field_to_color_prefs = FIELD_TO_COLOR_PREFS,\
                                background_color='white')
        
        #Iterate through each histogram file and assure it exisits.
        for k,v in label_to_histogram_filename.items():
            obs_file = open(v,'U').read()
            self.assertGreaterThan(len(obs_file),0)
            

    def test_get_histogram_scale(self):
        """get_histogram_scale should return correct result.
        """
        distances_dict = {'All_Category_Pairs':{'Treatment_Control_to_Control':\
                               [0.625, 0.623, 0.61, 0.577, 0.615,\
                                0.642, 0.673, 0.682, 0.737, 0.704],\
       'Treatment_Control_to_Fast':[0.729,  0.8  ,  0.721, 0.765,\
                                    0.776,  0.744,  0.749, 0.677,\
                                    0.734,  0.777,  0.733, 0.724,\
                                    0.696,  0.675,  0.654, 0.696,\
                                    0.731,  0.758,  0.738, 0.737],\
       'Treatment_Fast_to_Fast':[0.718,  0.666, 0.727, 0.6, 0.578,\
                                 0.623]
        }}
        bins = arange(0,1.01,0.05)
        xscale_exp = (0.55,0.85)
        yscale_exp = (0.0,0.5)
        xscale_obs, yscale_obs = get_histogram_scale(distances_dict,bins)
        self.assertFloatEqual(xscale_obs,xscale_exp)
        self.assertFloatEqual(yscale_obs,yscale_exp)

    def test_draw_histogram(self):
        """draw_histogram should return correct result.
        """
        hist_outfile = self.working_dir+'test_hist.png'
        distances = [0.625, 0.623, 0.61 , 0.577, 0.615, 0.642, 0.673, 0.682, \
            0.737, 0.704 ,0.718,  0.666, 0.727, 0.6, 0.578, 0.623]
        color = 'blue'
        draw_histogram(distances,color,10,hist_outfile)
        obs_file = open(hist_outfile,'U').read()
        self.assertGreaterThan(len(obs_file),0)

    def test_make_nav_html(self):
        """make_nav_html should return correct result.
        """
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field = self.single_field_treatment, \
                                paired_field = self.paired_field_treatment, \
                                dmat=self.dmat,\
                                histogram_dir = self.histogram_dir,\
                                field_to_color_prefs = FIELD_TO_COLOR_PREFS,\
                                background_color = 'white')
        nav_html_obs = \
            make_nav_html(distances_dict, label_to_histogram_filename)
        self.assertEqual(nav_html_obs, NAV_HTML)
        
    def test_make_main_html(self):
        """make_main_html should return correct result.
        """
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field = self.single_field_treatment, \
                                paired_field = self.paired_field_treatment, \
                                dmat=self.dmat,\
                                histogram_dir = self.histogram_dir,\
                                field_to_color_prefs = FIELD_TO_COLOR_PREFS,\
                                background_color = 'white')
                                
        make_main_html(distances_dict=distances_dict, \
                       label_to_histogram_filename=label_to_histogram_filename,\
                       root_outdir=self.working_dir, \
                       outfile_name='QIIME_Distance_Histograms.html')
        main_html_obs = \
            open(self.working_dir+\
                'QIIME_Distance_Histograms.html','U').readlines()
        
        self.assertEqual(len(main_html_obs),len(MAIN_HTML.split('\n')))

    def test_get_valid_indices(self):
        """get_valid_indices should return correct result.
        """
        control_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']
        exp_control_indices = [0,1,2,3,4]
        
        fast_ids = ['PC.607', 'PC.634', 'PC.635', 'PC.636']
        exp_fast_indices = [5,6,7,8]
        
        obs_control = get_valid_indices(self.distance_header, control_ids)
        self.assertEqual(obs_control, exp_control_indices)
        
        obs_fast = get_valid_indices(self.distance_header, fast_ids)
        self.assertEqual(obs_fast, exp_fast_indices)

    def test_distances_by_groups(self):
        """distances_by_groups should return correct result.
        """
        exp = [\
            ['Control','Fast',array([[0.729,  0.8  ,  0.721, 0.765],
                                    [0.776,  0.744,  0.749, 0.677],
                                    [0.734,  0.777,  0.733, 0.724],
                                    [0.696,  0.675,  0.654, 0.696],
                                    [0.731,  0.758,  0.738, 0.737]])],\
            ['Control','Control', array([0.625,  0.623,  0.61 ,  0.577, 0.615,
                                        0.642,  0.673, 0.682,  0.737, 0.704])],\
            ['Fast','Fast', array([0.718,  0.666, 0.727, 0.6, 0.578, 0.623])]
        ]
        
        obs = distances_by_groups(self.distance_header,self.dmat,\
            self.treatment_groups)
        self.assertEqual(obs,exp)

    def test_write_distance_files(self):
        """write_distance_files should return correct result.
        """
        exp_path = self.working_dir+'distances/dist_Treatment.xls'
        write_distance_files(self.single_field_treatment,\
            dir_prefix=self.working_dir)
        self.assertEqual(open(exp_path,'U').read(),
            open(self.distances_file,'U').read())

    def test_group_distances(self):
        """group_distances should return correct result.
        """
        single_field_exp = {'Treatment':[\
            ['Control','Fast',array([[0.729,  0.8  ,  0.721, 0.765],
                                    [0.776,  0.744,  0.749, 0.677],
                                    [0.734,  0.777,  0.733, 0.724],
                                    [0.696,  0.675,  0.654, 0.696],
                                    [0.731,  0.758,  0.738, 0.737]])],\
            ['Control','Control', array([0.625,  0.623,  0.61 ,  0.577, 0.615,
                                        0.642,  0.673, 0.682,  0.737, 0.704])],\
            ['Fast','Fast', array([0.718,  0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}
        
        paired_field_exp = {'Treatment_to_Treatment':[\
            [('Control','Control'),('Fast','Fast'),\
                             array([[0.729,  0.8  ,  0.721, 0.765],
                                    [0.776,  0.744,  0.749, 0.677],
                                    [0.734,  0.777,  0.733, 0.724],
                                    [0.696,  0.675,  0.654, 0.696],
                                    [0.731,  0.758,  0.738, 0.737]])],\
            [('Control','Control'),('Control','Control'),\
                                 array([0.625,  0.623,  0.61 ,  0.577, 0.615,
                                        0.642,  0.673, 0.682,  0.737, 0.704])],\
            [('Fast','Fast'),('Fast','Fast'),\
                             array([0.718,  0.666, 0.727, 0.6, 0.578, 0.623])]
        ]}
        
        single_field_obs, paired_field_obs, dmat_obs = \
            group_distances(mapping_file=self.map_file,\
                dmatrix_file=self.dmat_file,\
                fields=['Treatment'],\
                dir_prefix=self.working_dir)
        
        self.assertEqual(single_field_exp,single_field_obs)
        self.assertEqual(paired_field_exp,paired_field_obs)        
        self.assertEqual(self.dmat,dmat_obs)
        
    def test_monte_carlo_group_distances(self):
        """monte_carlo_group_distances should return correct result.
        """
        mc_group_dist_path = self.working_dir+\
            'monte_carlo_group_distances/group_distances_Treatment.xls'
        monte_carlo_group_distances(mapping_file=self.map_file,\
                                    dmatrix_file=self.dmat_file,\
                                    prefs=PREFS,\
                                    dir_prefix=self.working_dir)
        obs_res = open(mc_group_dist_path,'U').readlines()
        exp_res = MONTE_CARLO_DISTANCES.split('\n')
        for obs, exp in zip(obs_res, exp_res):
            obs_fields = obs.split('\t')
            exp_fields = exp.split('\t')
            #Check first 10 fields should be identical from run to run.
            for i in [0,1,2,3,5,6,7,8,9,11,13]:
                self.assertEqual(obs_fields[i],exp_fields[i])

    def test_permute_for_monte_carlo(self):
        """permute_for_monte_carlo should return correct result.
        """
        obs = permute_for_monte_carlo(self.dmat)
        self.assertNotEqual(obs,self.dmat)
        self.assertEqual(len(obs),len(self.dmat))
        self.assertEqual(sorted(obs.flat),sorted(self.dmat.flat))

    def test__make_histogram_filenames(self):
        """_make_histogram_filenames should return correct result.
        """
        distances = {'Distances1':[0,.5,.6,.7],\
                     'Distances2':[.3,.7,.9,.2,1.]}
        obs = _make_histogram_filenames(distances,self.histogram_dir)
        self.assertEqual(obs.keys(),distances.keys())
        for k,v in obs.items():
            exp_prefix = self.histogram_dir+k
            self.assertEquals(v[:len(exp_prefix)],exp_prefix)
            self.assertEquals(v[-4:],'.png')

    def test__make_path(self):
        """_make_path should return correct result.
        """
        exp = 'path/to/files/'
        paths = ['path','to','files']
        self.assertEqual(_make_path(paths),exp)

    def test__make_relative_paths(self):
        """_make_relative_paths should return correct result.
        """
        label_to_path = {'first_path':'/path/to/file1.txt',\
                         'second_path':'/path/to/file2.txt'}
        prefix = '/path/'
        exp = {'first_path':'./to/file1.txt',\
               'second_path':'./to/file2.txt'}
        self.assertEqual(_make_relative_paths(label_to_path,prefix),exp)
            
    def test__make_random_filename(self):
        """_make_random_filename should return correct result.
        """
        prefix = 'test_prefix'
        suffix = 'test_suffix'
        num_chars = 10
        obs = _make_random_filename(prefix=prefix,suffix=suffix,\
            num_chars=num_chars)
        self.assertEqual(len(obs),sum([len(prefix),len(suffix),num_chars]))
        self.assertEqual(obs[:len(prefix)],prefix)
        self.assertEqual(obs[len(prefix)+num_chars:],suffix)

    def test__get_script_dir(self):
        """_get_script_dir should return correct result.
        """
        script_path = '/Qiime/qiime/make_distance_histograms.py'
        exp = '/Qiime/qiime/'
        self.assertEqual(_get_script_dir(script_path),exp)


DISTANCE_MATRIX_STRING = \
"""\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\tPC.607\tPC.634\tPC.635\tPC.636
PC.354\t0.0\t0.625\t0.623\t0.61\t0.577\t0.729\t0.8\t0.721\t0.765
PC.355\t0.625\t0.0\t0.615\t0.642\t0.673\t0.776\t0.744\t0.749\t0.677
PC.356\t0.623\t0.615\t0.0\t0.682\t0.737\t0.734\t0.777\t0.733\t0.724
PC.481\t0.61\t0.642\t0.682\t0.0\t0.704\t0.696\t0.675\t0.654\t0.696
PC.593\t0.577\t0.673\t0.737\t0.704\t0.0\t0.731\t0.758\t0.738\t0.737
PC.607\t0.729\t0.776\t0.734\t0.696\t0.731\t0.0\t0.718\t0.666\t0.727
PC.634\t0.8\t0.744\t0.777\t0.675\t0.758\t0.718\t0.0\t0.6\t0.578
PC.635\t0.721\t0.749\t0.733\t0.654\t0.738\t0.666\t0.6\t0.0\t0.623
PC.636\t0.765\t0.677\t0.724\t0.696\t0.737\t0.727\t0.578\t0.623\t0.0"""

MAPPING_STRING = \
"""#SampleID\tBarcodeSequence\tTreatment\tDOB
PC.354\tAGCACGAGCCTA\tControl\t20061218
PC.355\tAACTCGTCGATG\tControl\t20061218
PC.356\tACAGACCACTCA\tControl\t20061126
PC.481\tACCAGCGACTAG\tControl\t20070314
PC.593\tAGCAGCACTTGT\tControl\t20071210
PC.607\tAACTGTGCGTAC\tFast\t20071112
PC.634\tACAGAGTCGGCT\tFast\t20080116
PC.635\tACCGCAGAGTCA\tFast\t20080116
PC.636\tACGGTGAGTGTC\tFast\t20080116
"""

PREFS = {
'MONTE_CARLO_GROUP_DISTANCES': {'Treatment':10},
'FIELDS':['Treatment'],
}

DISTANCES_OUT = \
"""Control_to_Fast\t0.729\t0.8\t0.721\t0.765\t0.776\t0.744\t0.749\t0.677\t0.734\t0.777\t0.733\t0.724\t0.696\t0.675\t0.654\t0.696\t0.731\t0.758\t0.738\t0.737
Control_to_Control\t0.625\t0.623\t0.61\t0.577\t0.615\t0.642\t0.673\t0.682\t0.737\t0.704
Fast_to_Fast\t0.718\t0.666\t0.727\t0.6\t0.578\t0.623
"""

MONTE_CARLO_DISTANCES = \
"""Control\tto\tFast\tavg\t0.730721926055\tcompared with\tControl\tto\tControl\tavg\t0.648710844983\t: t=\t3.2800428288\tp=\t0.00597272430874\tp_greater:\t0.0\tp_less:\t1.0\tnum_iters:\t10
Control\tto\tFast\tavg\t0.730721926055\tcompared with\tFast\tto\tFast\tavg\t0.651918984745\t: t=\t2.48083030062\tp=\t0.0349438048591\tp_greater:\t0.0\tp_less:\t1.0\tnum_iters:\t10
Control\tto\tControl\tavg\t0.648710844983\tcompared with\tFast\tto\tFast\tavg\t0.651918984745\t: t=\t-0.114952651644\tp=\t0.910115083574\tp_greater:\t0.6\tp_less:\t0.4\tnum_iters:\t10"""

NAV_HTML = \
"""<td>
    <div style="overflow:scroll; width: 300px; height: 400px;">
    
    <p>
<span class="normal">All_Within_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Within_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Within_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Within_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Within_Category_Distances')">Treatment_All_Within_Category_Distances</a></span><br />
<span class="normal">All_Between_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Between_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Between_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Between_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Between_Category_Distances')">Treatment_All_Between_Category_Distances</a></span><br />
<span class="normal">All_Between_Sample_Distances</span><br />
<span class="smnorm"><input type="checkbox" id="check_All_Between_Sample_Distances" checked onclick="visibilityAndOpacity(this, 'All_Between_Sample_Distances')" />
<a onmouseover="mouseoverVisible('All_Between_Sample_Distances')"; onmouseout="mouseoverHidden('All_Between_Sample_Distances')">All_Between_Sample_Distances</a></span><br />
<span class="normal">All_Within_And_Between_Fields</span><br />
<span class="smnorm"><input type="checkbox" id="check_Within_All_Fields"  onclick="visibilityAndOpacity(this, 'Within_All_Fields')" />
<a onmouseover="mouseoverVisible('Within_All_Fields')"; onmouseout="mouseoverHidden('Within_All_Fields')">Within_All_Fields</a></span><br />
<span class="normal">All_Category_Pairs</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Control_to_Control"  onclick="visibilityAndOpacity(this, 'Treatment_Control_to_Control')" />
<a onmouseover="mouseoverVisible('Treatment_Control_to_Control')"; onmouseout="mouseoverHidden('Treatment_Control_to_Control')">Treatment_Control_to_Control</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Control_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment_Control_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment_Control_to_Fast')"; onmouseout="mouseoverHidden('Treatment_Control_to_Fast')">Treatment_Control_to_Fast</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Fast_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment_Fast_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment_Fast_to_Fast')"; onmouseout="mouseoverHidden('Treatment_Fast_to_Fast')">Treatment_Fast_to_Fast</a></span><br />
<span class="normal">All_Within_Categories</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Control_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Control_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Control_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Control_Distances')">Treatment_Within_Control_Distances</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Fast_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Fast_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Fast_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Fast_Distances')">Treatment_Within_Fast_Distances</a></span><br />
    </p>
    </div>
    </td>
"""

MAIN_HTML = \
"""<html><head> <title>
Distance Histograms
</title>
<script type="text/javascript" src="./js/histograms.js"></script>
 <style type="text/css">
.smnorm {color: blue; font-family:Arial,Verdana; font-size:10; font-weight: bold;}
.normal {color: black; font-family:Arial,Verdana; font-size:11; font-weight: bold;}

</style>
</head>
<body> 
<div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;"></div>
<table width="200" border="0" cellspacing="2" cellpadding="2"> <table width="200" border="0" cellspacing="2" cellpadding="2"> <tr><td colspan="2" class="header_qiime" align="center">
    <table width=800 cellpadding=0 cellspacing=0 border=0>
    <tr valign=middle><td class=ntitle width=200 valign="middle"><img src="./web_resources/qiime_header.png" border="0" /></td>
        <td width=300 align=center >
            &nbsp; 
        </td>
    </tr> 
    </table>
</td></tr>
<tr><td colspan="2" align="left" valign="top">&nbsp;</td></tr> 
 </table>  <tr><td colspan="2" align="left" valign="top" class="normal"> <table border="0" cellspacing="1" cellpadding="0" width="800">

</table>  </td> </tr> </table>
<div>
    <table>
        <tr>
            
<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:visible;" id="All_Between_Sample_Distances" name="visible" src=".//histograms/All_Between_Sample_Distances_VDGtFjE1p2qRCvzsBt0w.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Within_Fast_Distances" name="hidden" src=".//histograms/Treatment_Within_Fast_Distances_oMHomERhcoGVy0taqoXy.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Within_All_Fields" name="hidden" src=".//histograms/Within_All_Fields_Ne8iQ6jzEs0QKECgkaOJ.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_All_Between_Category_Distances" name="hidden" src=".//histograms/Treatment_All_Between_Category_Distances_O8U91gtDsRDPfkc1B5Bg.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Control_to_Fast" name="hidden" src=".//histograms/Treatment_Control_to_Fast_BguPgXBBPQHd3JGZiWWy.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Control_to_Control" name="hidden" src=".//histograms/Treatment_Control_to_Control_VwziTye01o8i5masUoN2.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="All_Between_Sample_Distances" name="hidden" src=".//histograms/All_Between_Sample_Distances_VDGtFjE1p2qRCvzsBt0w.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_All_Within_Category_Distances" name="hidden" src=".//histograms/Treatment_All_Within_Category_Distances_y2jFfSgOfxmpRJAbc5LP.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Fast_to_Fast" name="hidden" src=".//histograms/Treatment_Fast_to_Fast_czBTLqB8WcaJyGREPG12.png" border="0"></td>

<td style="position:absolute; top:100; left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:hidden;" id="Treatment_Within_Control_Distances" name="hidden" src=".//histograms/Treatment_Within_Control_Distances_tRN8cTb1QXc0qzknWKQQ.png" border="0"></td>
        </tr>
    </table>
</div>


<table style="position:absolute; top:100; left:600">
    <tr>
    <td>
    <div style="overflow:scroll; width: 300px; height: 400px;">
    
    <p>
<span class="normal">All_Within_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Within_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Within_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Within_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Within_Category_Distances')">Treatment_All_Within_Category_Distances</a></span><br />
<span class="normal">All_Between_Category_Grouped</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_All_Between_Category_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_All_Between_Category_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_All_Between_Category_Distances')"; onmouseout="mouseoverHidden('Treatment_All_Between_Category_Distances')">Treatment_All_Between_Category_Distances</a></span><br />
<span class="normal">All_Between_Sample_Distances</span><br />
<span class="smnorm"><input type="checkbox" id="check_All_Between_Sample_Distances" checked onclick="visibilityAndOpacity(this, 'All_Between_Sample_Distances')" />
<a onmouseover="mouseoverVisible('All_Between_Sample_Distances')"; onmouseout="mouseoverHidden('All_Between_Sample_Distances')">All_Between_Sample_Distances</a></span><br />
<span class="normal">All_Within_And_Between_Fields</span><br />
<span class="smnorm"><input type="checkbox" id="check_Within_All_Fields"  onclick="visibilityAndOpacity(this, 'Within_All_Fields')" />
<a onmouseover="mouseoverVisible('Within_All_Fields')"; onmouseout="mouseoverHidden('Within_All_Fields')">Within_All_Fields</a></span><br />
<span class="normal">All_Category_Pairs</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Control_to_Control"  onclick="visibilityAndOpacity(this, 'Treatment_Control_to_Control')" />
<a onmouseover="mouseoverVisible('Treatment_Control_to_Control')"; onmouseout="mouseoverHidden('Treatment_Control_to_Control')">Treatment_Control_to_Control</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Control_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment_Control_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment_Control_to_Fast')"; onmouseout="mouseoverHidden('Treatment_Control_to_Fast')">Treatment_Control_to_Fast</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Fast_to_Fast"  onclick="visibilityAndOpacity(this, 'Treatment_Fast_to_Fast')" />
<a onmouseover="mouseoverVisible('Treatment_Fast_to_Fast')"; onmouseout="mouseoverHidden('Treatment_Fast_to_Fast')">Treatment_Fast_to_Fast</a></span><br />
<span class="normal">All_Within_Categories</span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Control_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Control_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Control_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Control_Distances')">Treatment_Within_Control_Distances</a></span><br />
<span class="smnorm"><input type="checkbox" id="check_Treatment_Within_Fast_Distances"  onclick="visibilityAndOpacity(this, 'Treatment_Within_Fast_Distances')" />
<a onmouseover="mouseoverVisible('Treatment_Within_Fast_Distances')"; onmouseout="mouseoverHidden('Treatment_Within_Fast_Distances')">Treatment_Within_Fast_Distances</a></span><br />
    </p>
    </div>
    </td>

    </tr>
</table>
"""

DATA_COLOR_ORDER = ['blue', 'lime', 'red', 'aqua', 'fuchsia', 'yellow', \
    'green', 'maroon', 'teal', 'purple', 'olive', 'silver', 'gray']

FIELD_TO_COLOR_PREFS = {'DOB': ({'20070314': ['PC.481'], '20071112': ['PC.607'], '20080116': ['PC.634', 'PC.635', 'PC.636'], '20061126': ['PC.356'], '20061218': ['PC.354', 'PC.355'], '20071210': ['PC.593']}, {'20070314': 'red', '20071112': 'aqua', '20080116': 'yellow', '20061126': 'blue', '20061218': 'lime', '20071210': 'fuchsia'}, data_colors, DATA_COLOR_ORDER),'Treatment': ({'Control': ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593'], 'Fast': ['PC.607', 'PC.634', 'PC.635', 'PC.636']}, {'Control': 'blue', 'Fast': 'lime'}, data_colors, DATA_COLOR_ORDER), 'BarcodeSequence': ({'ACCAGCGACTAG': ['PC.481'], 'ACCGCAGAGTCA': ['PC.635'], 'AACTGTGCGTAC': ['PC.607'], 'AGCAGCACTTGT': ['PC.593'], 'ACAGAGTCGGCT': ['PC.634'], 'AACTCGTCGATG': ['PC.355'], 'ACGGTGAGTGTC': ['PC.636'], 'AGCACGAGCCTA': ['PC.354'], 'ACAGACCACTCA': ['PC.356']}, {'ACCAGCGACTAG': 'fuchsia', 'ACCGCAGAGTCA': 'yellow', 'AACTGTGCGTAC': 'lime', 'ACAGACCACTCA': 'red', 'ACAGAGTCGGCT': 'aqua', 'AACTCGTCGATG': 'blue', 'ACGGTGAGTGTC': 'green', 'AGCAGCACTTGT': 'teal', 'AGCACGAGCCTA': 'maroon'}, data_colors, DATA_COLOR_ORDER),'Description': ({'Fasting mouse, I.D. 607': ['PC.607'], 'Control mouse, I.D. 481': ['PC.481'], 'Control mouse, I.D. 593': ['PC.593'], 'Control mouse, I.D. 356': ['PC.356'], 'Control mouse, I.D. 354': ['PC.354'], 'Control mouse, I.D. 355': ['PC.355'], 'Fasting mouse, I.D. 634': ['PC.634'], 'Fasting mouse, I.D. 635': ['PC.635'], 'Fasting mouse, I.D. 636': ['PC.636']}, {'Fasting mouse, I.D. 607': 'yellow', 'Control mouse, I.D. 481': 'aqua', 'Control mouse, I.D. 593': 'fuchsia', 'Control mouse, I.D. 356': 'red', 'Control mouse, I.D. 354': 'blue', 'Control mouse, I.D. 355': 'lime', 'Fasting mouse, I.D. 634': 'green', 'Fasting mouse, I.D. 635': 'maroon', 'Fasting mouse, I.D. 636': 'teal'}, data_colors, DATA_COLOR_ORDER), 'SampleID': ({'PC.636': ['PC.636'], 'PC.355': ['PC.355'], 'PC.607': ['PC.607'], 'PC.634': ['PC.634'], 'PC.635': ['PC.635'], 'PC.593': ['PC.593'], 'PC.356': ['PC.356'], 'PC.481': ['PC.481'], 'PC.354': ['PC.354']}, {'PC.636': 'teal', 'PC.355': 'lime', 'PC.607': 'yellow', 'PC.634': 'green', 'PC.635': 'maroon', 'PC.593': 'fuchsia', 'PC.356': 'red', 'PC.481': 'aqua', 'PC.354': 'blue'}, data_colors, DATA_COLOR_ORDER)}

#run tests if called from command line
if __name__ == "__main__":
    main()