#!/usr/bin/env python

"""Tests of code for assigning taxonomy"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Greg Caporaso", "Kyle Bittinger"] 
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


from cStringIO import StringIO
from os import remove, system
from glob import glob
from tempfile import NamedTemporaryFile
from shutil import copy as copy_file
from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.util.misc import remove_files
from cogent.parse.fasta import MinimalFastaParser
from qiime.assign_taxonomy import TaxonAssigner, BlastTaxonAssigner,\
 RdpTaxonAssigner


class TaxonAssignerTests(TestCase):
    """Tests of the abstract TaxonAssigner class"""

    def test_init(self):
        """Abstract TaxonAssigner __init__ should store name, params"""
        p = TaxonAssigner({})
        self.assertEqual(p.Name, 'TaxonAssigner')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract TaxonAssigner __call__ should raise NotImplementedError"""
        p = TaxonAssigner({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs')
        

class BlastTaxonAssignerTests(TestCase):
    """Tests of the BlastTaxonAssigner class"""
    
    def setUp(self):
        self.id_to_taxonomy_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.txt') 
        self.input_seqs_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self.reference_seqs_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        
        self._paths_to_clean_up =\
         [self.id_to_taxonomy_fp,self.input_seqs_fp,self.reference_seqs_fp] 
        
        open(self.id_to_taxonomy_fp,'w').write(id_to_taxonomy_string)
        open(self.input_seqs_fp,'w').write(test_seq_coll.toFasta())
        self.test_seqs = test_seq_coll.items()
        open(self.reference_seqs_fp,'w').write(test_refseq_coll.toFasta())
        
        self.expected1 = {
            's1': ('Archaea;Euryarchaeota;Halobacteriales;uncultured', 0.0, "AY800210"),
            's2': ('Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.', 0.0, "EU883771"),
            's3': ('Archaea;Crenarchaeota;uncultured;uncultured', 0.0, "EF503699"),
            's4': ('Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium', 0.0, "DQ260310"),
            's5': ('Archaea;Crenarchaeota;uncultured;uncultured', 0.0, "EF503697"),
            's6': ('No blast hit', None, None),
            }
        
    def tearDown(self):
        remove_files(set(self._paths_to_clean_up))

    def test_init(self):
        """BlastTaxonAssigner __init__ should store name, params"""
        p = BlastTaxonAssigner({})
        self.assertEqual(p.Name, 'BlastTaxonAssigner')
        # default parameters correctly initialized
        default_params = {'Min percent identity':0.90,\
         'Max E value':1e-30,\
         'Application':'blastn/megablast'}
        self.assertEqual(p.Params, default_params)
        
    def test_parse_id_to_taxonomy_file(self):
        """Parsing taxonomy files functions as expected
        """
        lines = id_to_taxonomy_string.splitlines()
        p = BlastTaxonAssigner({})
        expected = {\
         "AY800210":"Archaea;Euryarchaeota;Halobacteriales;uncultured",\
         "EU883771":"Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",\
         "EF503699":"Archaea;Crenarchaeota;uncultured;uncultured",\
         "DQ260310":"Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",\
         "EF503697":"Archaea;Crenarchaeota;uncultured;uncultured"}
        self.assertEqual(p._parse_id_to_taxonomy_file(lines),expected)
        
    def test_map_ids_to_taxonomy(self):
        """Mapping sequence ids to taxonomy functions as expected
        """
        p = BlastTaxonAssigner({})
        id_to_taxonomy_map = {
            "AY800210": "Archaea;Euryarchaeota;Halobacteriales;uncultured",
            "EU883771": "Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",
            "EF503699": "Archaea;Crenarchaeota;uncultured;uncultured",
            "DQ260310": "Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",
            "EF503697": "Archaea;Crenarchaeota;uncultured;uncultured",
            }
        hits = {
            's1': ("AY800210", 1e-99),
            's5': ("EU883771", 'weird confidence value'),
            's3': ("DQ260310", 42.),
            's4': None,
            }
        expected = {
            's1': ("Archaea;Euryarchaeota;Halobacteriales;uncultured", 1e-99, "AY800210"),
            's5': ('Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.',
                   'weird confidence value',"EU883771"),
            's3': ("Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium", 42.,"DQ260310"),
            's4': ('No blast hit', None, None),
            }
        actual = p._map_ids_to_taxonomy(hits,id_to_taxonomy_map)
        self.assertEqual(actual,expected)
        
    def test_get_first_blast_hit_per_seq(self):
        """Extracting the first blast hit for each seq functions as expected
        """
        p = BlastTaxonAssigner({})
        blast_hits = {'s1':[('blah',0.0)],\
                      's3':[('dsasd',1e-42),('rrr',1e-12),('qqq',0.001)],\
                      's2':[]}
        expected = {'s1':('blah',0.0),\
                      's3':('dsasd',1e-42),\
                      's2':None}
        actual = p._get_first_blast_hit_per_seq(blast_hits)
        self.assertEqual(actual,expected)
        
    def test_get_blast_hits(self):
        """BlastTaxonAssigner._get_blast_hits functions w existing db
        
        """
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove
        
        p = BlastTaxonAssigner({})
        seq_coll_blast_results = p._get_blast_hits(blast_db,self.test_seqs)
        # mapping from identifier in test_seq_coll to the id of the sequence
        # in the refseq collection (a silva derivative)
        expected_matches = {\
         's1':'AY800210',
         's2':'EU883771',\
         's3':'EF503699',\
         's4':'DQ260310',\
         's5':'EF503697'}
        
        # no results for s6 (which is a randomly-generated sequence) 
        s6_blast_results = seq_coll_blast_results['s6']
        self.assertEqual(s6_blast_results,[])
        
        # expected results for all other query sequences
        for seq_id in expected_matches:
            blast_results = seq_coll_blast_results[seq_id]
            blast_results_d = dict(blast_results)
            # explicitly checks that the result is in the data before 
            # pulling it out (this is redundant, but allows for a useful
            # error message if the data wasn't in there b/c e.g. there 
            # were no blast results returned)
            self.assertTrue(expected_matches[seq_id] in blast_results_d)
            # now check that the perfect match got a 0.0 e-value as it should
            # on this data
            self.assertEqual(blast_results_d[expected_matches[seq_id]],0.0)
            
    def test_call_existing_blast_db(self):
        """BlastTaxonAssigner.__call__ functions w existing db
        """
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove
        
        p = BlastTaxonAssigner({'blast_db':blast_db,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        actual = p(self.input_seqs_fp)
        
        self.assertEqual(actual,self.expected1)
        
    def test_call_alt_input_types(self):
        """BlastTaxonAssigner.__call__ functions w alt input types """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
         
        # neither seqs or seq_fp passed results in AssertionError
        self.assertRaises(AssertionError,p)
        
        # Functions with a list of (seq_id, seq) pairs
        seqs = list(MinimalFastaParser(open(self.input_seqs_fp)))
        actual = p(seqs=seqs)
        self.assertEqual(actual,self.expected1)
        
        # Functions with input path
        actual = p(self.input_seqs_fp)
        self.assertEqual(actual,self.expected1)
        
        # same result when passing fp or seqs
        self.assertEqual(p(seqs=seqs),p(self.input_seqs_fp))
        
    def test_seqs_to_taxonomy(self):
        """BlastTaxonAssigner._seqs_to_taxonomy: functions as expected
        """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
         
        # build the id_to_taxonomy_map as this test doesn't execute __call__
        id_to_taxonomy_map = {
            "AY800210": \
             "Archaea;Euryarchaeota;Halobacteriales;uncultured",
            "EU883771": \
             "Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",
            "EF503699": \
             "Archaea;Crenarchaeota;uncultured;uncultured",
            "DQ260310": \
             "Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",
            "EF503697": \
             "Archaea;Crenarchaeota;uncultured;uncultured",
            }
        
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove
        
        # read the input file into (seq_id, seq) pairs
        seqs = list(MinimalFastaParser(open(self.input_seqs_fp)))
        
        actual = p._seqs_to_taxonomy(seqs,blast_db,id_to_taxonomy_map)
        self.assertEqual(actual,self.expected1)
        
        # passing empty list of seqs functions as expected
        actual = p._seqs_to_taxonomy([],blast_db,id_to_taxonomy_map)
        self.assertEqual(actual,{})
        
        
    def test_call_on_the_fly_blast_db(self):
        """BlastTaxonAssigner.__call__ functions w creating blast db
        """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        actual = p(self.input_seqs_fp)
        
        self.assertEqual(actual,self.expected1)
            
    def test_call_output_to_file(self):
        """BlastTaxonAssigner.__call__ functions w output to file
        """
        result_path = get_tmp_filename(
            prefix='BlastTaxonAssignerTests_', suffix='.fasta')
        self._paths_to_clean_up.append(result_path)

        p = BlastTaxonAssigner({
            'reference_seqs_filepath': self.reference_seqs_fp,
            'id_to_taxonomy_filepath': self.id_to_taxonomy_fp,
            })
        actual = p(self.input_seqs_fp, result_path=result_path)

        expected_lines = set([
            's1\tArchaea;Euryarchaeota;Halobacteriales;uncultured\t0.0\tAY800210\n',
            's2\tArchaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.\t0.0\tEU883771\n',
            's3\tArchaea;Crenarchaeota;uncultured;uncultured\t0.0\tEF503699\n',
            's4\tArchaea;Euryarchaeota;Methanobacteriales;Methanobacterium\t0.0\tDQ260310\n',
            's5\tArchaea;Crenarchaeota;uncultured;uncultured\t0.0\tEF503697\n',
            's6\tNo blast hit\tNone\tNone\n',
            ])
        f = open(result_path)
        observed_lines = set(f.readlines())
        f.close()
        self.assertEqual(observed_lines, expected_lines)
        
        # Return value is None when result_path is provided (Not sure
        # if this is what we want yet, or if we would want both so 
        # results could be logged to file...)
        self.assertEqual(actual, None)
        
    def test_call_logs_run(self):
        """BlastTaxonAssigner.__call__ logs the run when expected
        """
        log_path = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self._paths_to_clean_up.append(log_path) 
        
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove
        
        p = BlastTaxonAssigner({\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp,\
         'blast_db':blast_db})
        actual = p(self.input_seqs_fp,log_path=log_path)
        
        log_file = open(log_path)
        log_file_str = log_file.read()
        log_file.close()

        log_file_exp = [
            "BlastTaxonAssigner parameters:",
            'Min percent identity:0.9',
            'Application:blastn/megablast',
            'Max E value:1e-30',
            'Result path: None, returned as dict.',
            'blast_db:%s' % str(self.reference_seqs_fp)[1:-1],
            'id_to_taxonomy_filepath:%s' % self.id_to_taxonomy_fp,
            'Number of sequences inspected: 6',
            'Number with no blast hits: 1',
            '',
         ]
        # compare data in log file to fake expected log file
        # NOTE: Since p.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that 
        # the equal unordered lists of lines is present in actual and expected
        self.assertEqualItems(log_file_str.split('\n'), log_file_exp)
        
        
class RdpTaxonAssignerTests(TestCase):
    """Tests for the Rdp-based taxonomy assigner.

    Slow tests are a bug, and currently these tests take about 38s on
    a Dual 2.3GHz Mac.  Presumably most of the time is spent
    initializing the Java VM on each run.  If so, this problem should
    be fixed upstream in PyCogent.
    """

    def setUp(self):
        # Temporary input file
        self.tmp_seq_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTest_',
            suffix='.fasta'
            )
        seq_file = open(self.tmp_seq_filepath, 'w')
        seq_file.write(rdp_test1_fasta)
        seq_file.close()

        # Temporary results filename
        self.tmp_res_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTestResult_',
            suffix='.tsv',
            )
        # touch the file so we don't get an error trying to close it
        open(self.tmp_res_filepath,'w').close()

        # Temporary log filename
        self.tmp_log_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTestLog_',
            suffix='.txt',
            )
        # touch the file so we don't get an error trying to close it
        open(self.tmp_log_filepath,'w').close()

        self._paths_to_clean_up = \
         [self.tmp_seq_filepath, self.tmp_res_filepath, self.tmp_log_filepath]
         
        self.id_to_taxonomy_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.txt')
        self.id_to_taxonomy_file.write(rdp_id_to_taxonomy)
        self.id_to_taxonomy_file.seek(0)

        self.reference_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.fasta')
        self.reference_seqs_file.write(rdp_reference_seqs)
        self.reference_seqs_file.seek(0)
        


        self.default_app = RdpTaxonAssigner({})
        
        # Why is this giving weird results (strings get increasingly longer)
        # print str(self.default_app)

    def tearDown(self):
        remove_files(self._paths_to_clean_up)
 
    def test_init(self):
        """RdpTaxonAssigner.__init__ should set default attributes and params
        """
        a = RdpTaxonAssigner({})
        self.assertEqual(a.Name, 'RdpTaxonAssigner')
        
    def test_train_on_the_fly(self):
        """Training on-the-fly classifies reference sequence correctly with 100% certainty
        """
        input_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.fasta')
        input_seqs_file.write(test_seq_coll.toFasta())
        input_seqs_file.seek(0)

        expected = rdp_trained_test1_expected_dict
        
        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': self.id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        actual = app(self.tmp_seq_filepath)
        
        key = 'X67228 some description'
        self.assertEqual(actual[key], expected[key])

    def test_parse_lineage(self):
        """Lineage in csv format is correctly parsed to a list
        """
        str = 'Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.;a;b'
        actual = RdpTaxonAssigner._parse_lineage(str)
        expected = ['Archaea', 'Euryarchaeota', 'Methanomicrobiales', 
                    'Methanomicrobium et rel.', 'a', 'b']
        self.assertEqual(actual, expected)

    def test_build_tree(self):
        """RdpTaxonAssigner._build_tree() should return a tree with correct Rdp-format taxonomy
        """
        tree = RdpTaxonAssigner._build_tree(self.id_to_taxonomy_file)
        actual = tree.rdp_taxonomy()
        # The order of the lines in this file depends on python's
        # dict() implementation, so we should ideally build two sets
        # of lines and check that their contents match.
        expected = rdp_expected_taxonomy
        self.assertEqual(actual, expected)

    def test_generate_training_seqs(self):
        seqs = RdpTaxonAssigner._generate_training_seqs(
            self.reference_seqs_file, self.id_to_taxonomy_file)
        actual = LoadSeqs(data=seqs, aligned=False).toFasta()
        self.assertEqual(actual, rdp_expected_training_seqs)
        
    def test_generate_training_files(self):
        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': self.id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        actual_taxonomy_file, actual_training_seqs_file = \
            app._generate_training_files()

        # see note in test_build_tree()
        self.assertEqual(actual_taxonomy_file.read(), rdp_expected_taxonomy)
        

    def test_call_result_as_dict(self):
        """RdpTaxonAssigner should return correct taxonomic assignment
        
           This test may periodically fail, but should be rare.
           
        """
        expected = rdp_test1_expected_dict
        min_confidence = self.default_app.Params['Confidence']
        
        # Since there is some variation in the assignments, run
        # 10 trials and make sure we get the expected result at least once
        num_trials = 10
        num_seqs = len(expected)
        seq_ids = expected.keys()
        assignment_comp_results = [False] * num_seqs
        expected_assignment_comp_results = [True] * num_seqs
        
        for i in range(num_trials):
            actual = self.default_app(self.tmp_seq_filepath)
            # seq ids are the same, and all input sequences get a result
            self.assertEqual(actual.keys(),expected.keys())
            for j,seq_id in enumerate(seq_ids):
                # confidence is above threshold
                self.assertTrue(actual[seq_id][1] >= min_confidence)
                # confidence roughly matches expected
                self.assertFloatEqual(\
                 actual[seq_id][1],expected[seq_id][1],0.1)
                # check if the assignment is correct -- this must happen
                # at least once per seq_id for the test to pass
                if actual[seq_id][0] == expected[seq_id][0]:
                    assignment_comp_results[j] = True
            if assignment_comp_results == expected_assignment_comp_results:
                # break once we've seen a correct assignment for each seq
                break
                    
        self.assertEqual(\
         assignment_comp_results,\
         expected_assignment_comp_results,\
         "Taxonomic assignments never correct in %d trials." % num_trials)
    
    def test_call_result_to_file(self):
        """RdpTaxonAssigner should save results to file
        
           This test may periodically fail, but should be rare.
           
        """
        expected_lines = rdp_test1_expected_lines
        
        # Since there is some variation in the assignments, run
        # 10 trials and make sure we get the expected result at least once
        # for each sequence
        num_trials = 10
        num_seqs = len(expected_lines)
        assignment_comp_results = [False] * num_seqs
        expected_assignment_comp_results = [True] * num_seqs
        
        for i in range(num_trials):
            retval = self.default_app(
             seq_path=self.tmp_seq_filepath,
             result_path=self.tmp_res_filepath,
             log_path=None)
            actual = [l.strip() for l in open(self.tmp_res_filepath, 'r')]
            message = "Expected return value of None but observed %s" % retval
            self.assertTrue(retval is None, message)
            for j in range(num_seqs):
                a = actual[j]
                e = expected_lines[j]
                # note we're testing using startswith here to allow
                # for some variability in confidence
                if a.startswith(e):
                    assignment_comp_results[j] = True
            if assignment_comp_results == expected_assignment_comp_results:
                break
                    
        self.assertEqual(\
         assignment_comp_results,\
         expected_assignment_comp_results,\
         "Taxonomic assignments never correct in %d trials." % num_trials)


    def test_log(self):
        """RdpTaxonAssigner should write correct message to log file"""
        # expected result when no result_path is provided
        a = RdpTaxonAssigner({})
        a(seq_path=self.tmp_seq_filepath, 
          result_path=None, 
          log_path=self.tmp_log_filepath)
            
        # open the actual log file and the expected file, and pass into lists
        obs = [l.strip() for l in list(open(self.tmp_log_filepath, 'r'))]
        exp = rdp_test1_log_file_contents.split('\n')
        # sort the lists as the entries are written from a dict,
        # so order may vary
        obs.sort()
        exp.sort()
        self.assertEqual(obs, exp)
        

rdp_test1_fasta = \
""">X67228 some description
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>EF503697
TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC
"""

rdp_test1_log_file_contents = \
"""RdpTaxonAssigner parameters:
Application:RDP classfier
Citation:Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.
Taxonomy:RDP
Confidence:0.8
id_to_taxonomy_fp:None
reference_sequences_fp:None"""

rdp_test1_expected_dict = {\
 'X67228 some description': ('Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium',0.95),\
 'EF503697': ('Root;Archaea;Crenarchaeota;Thermoprotei',0.88)
}


rdp_trained_test1_expected_dict = {
    'X67228 some description': ('Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium', 1.0),
    'EF503697': ('Bacteria;Proteobacteria;Gammaproteobacteria', 0.83999999999999997),
    }

rdp_test1_expected_lines = [\
 "\t".join(["X67228 some description",\
  "Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium",\
  "0.9"]),
 "\t".join(['EF503697','Root;Archaea;Crenarchaeota;Thermoprotei','0.8'])]

rdp_id_to_taxonomy = \
"""X67228	Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium
X73443	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium
AB004750	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
xxxxxx	Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
AB004748	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
AB000278	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium
AB000390	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio
"""

rdp_reference_seqs = \
""">X67228
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>X73443
nnnnnnngagatttgatcctggctcaggatgaacgctggccggccgtgcttacacatgcagtcgaacgaagcgcttaaactggatttcttcggattgaagtttttgctgactgagtggcggacgggtgagtaacgcgtgggtaacctgcctcatacagggggataacagttagaaatgactgctaataccnnataagcgcacagtgctgcatggcacagtgtaaaaactccggtggtatgagatggacccgcgtctgattagctagttggtggggt
>AB004750
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgca
>xxxxxx
ttgaacgctggcggcaggcctaacacatgcaagtcgagcggcagcannnncttcgggaggctggcgagcggcggacgggtgagtaacgcatgggaacttacccagtagtgggggatagcccggggaaacccggattaataccgcatacgccctgagggggaaagcgggctccggtcgcgctattggatgggcccatgtcggattagttagttggtggggtaatggcctaccaaggcgacgatccgtagctggtctgagaggatgatcagccacaccgggactgagacacggcccggactcctacgggaggcagcagtggggaatattggacaatgggggcaaccctgatccagccatgccg
>AB004748
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggccttcgggttg
>AB000278
caggcctaacacatgcaagtcgaacggtaanagattgatagcttgctatcaatgctgacgancggcggacgggtgagtaatgcctgggaatataccctgatgtgggggataactattggaaacgatagctaataccgcataatctcttcggagcaaagagggggaccttcgggcctctcgcgtcaggattagcccaggtgggattagctagttggtggggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggggaaaccctgatgcagccatgccgcgtgta
>AB000390
tggctcagattgaacgctggcggcaggcctaacacatgcaagtcgagcggaaacgantnntntgaaccttcggggnacgatnacggcgtcgagcggcggacgggtgagtaatgcctgggaaattgccctgatgtgggggataactattggaaacgatagctaataccgcataatgtctacggaccaaagagggggaccttcgggcctctcgcttcaggatatgcccaggtgggattagctagttggtgaggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgag
"""


rdp_expected_taxonomy = \
"""1*Bacteria*0*0*domain
2*Firmicutes*1*1*phylum
3*Clostridia*2*2*class
4*Clostridiales*3*3*order
5*Clostridiaceae*4*4*family
6*Clostridium*5*5*genus
7*Proteobacteria*1*1*phylum
8*Alphaproteobacteria*7*2*class
9*Rhizobiales*8*3*order
10*Rhizobiaceae*9*4*family
11*Rhizobium*10*5*genus
12*Gammaproteobacteria*7*2*class
13*Enterobacteriales*12*3*order
14*Enterobacteriaceae*13*4*family
15*Enterobacter*14*5*genus
16*Pseudomonadales*12*3*order
17*Pseudomonadaceae*16*4*family
18*Pseudomonas*17*5*genus
19*Vibrionales*12*3*order
20*Vibrionaceae*19*4*family
21*Photobacterium*20*5*genus
22*Vibrio*20*5*genus
"""

# newline at the end makes a difference
rdp_expected_training_seqs = \
""">AB000278 Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium
caggcctaacacatgcaagtcgaacggtaanagattgatagcttgctatcaatgctgacgancggcggacgggtgagtaatgcctgggaatataccctgatgtgggggataactattggaaacgatagctaataccgcataatctcttcggagcaaagagggggaccttcgggcctctcgcgtcaggattagcccaggtgggattagctagttggtggggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggggaaaccctgatgcagccatgccgcgtgta
>AB000390 Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio
tggctcagattgaacgctggcggcaggcctaacacatgcaagtcgagcggaaacgantnntntgaaccttcggggnacgatnacggcgtcgagcggcggacgggtgagtaatgcctgggaaattgccctgatgtgggggataactattggaaacgatagctaataccgcataatgtctacggaccaaagagggggaccttcgggcctctcgcttcaggatatgcccaggtgggattagctagttggtgaggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgag
>AB004748 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggccttcgggttg
>AB004750 Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgca
>X67228 Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>X73443 Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium
nnnnnnngagatttgatcctggctcaggatgaacgctggccggccgtgcttacacatgcagtcgaacgaagcgcttaaactggatttcttcggattgaagtttttgctgactgagtggcggacgggtgagtaacgcgtgggtaacctgcctcatacagggggataacagttagaaatgactgctaataccnnataagcgcacagtgctgcatggcacagtgtaaaaactccggtggtatgagatggacccgcgtctgattagctagttggtggggt
>xxxxxx Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
ttgaacgctggcggcaggcctaacacatgcaagtcgagcggcagcannnncttcgggaggctggcgagcggcggacgggtgagtaacgcatgggaacttacccagtagtgggggatagcccggggaaacccggattaataccgcatacgccctgagggggaaagcgggctccggtcgcgctattggatgggcccatgtcggattagttagttggtggggtaatggcctaccaaggcgacgatccgtagctggtctgagaggatgatcagccacaccgggactgagacacggcccggactcctacgggaggcagcagtggggaatattggacaatgggggcaaccctgatccagccatgccg"""

id_to_taxonomy_string = \
"""AY800210\tArchaea;Euryarchaeota;Halobacteriales;uncultured
EU883771\tArchaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.
EF503699\tArchaea;Crenarchaeota;uncultured;uncultured
DQ260310\tArchaea;Euryarchaeota;Methanobacteriales;Methanobacterium
EF503697\tArchaea;Crenarchaeota;uncultured;uncultured"""

test_seq_coll = LoadSeqs(data=[\
 ('s1','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('s2','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('s3','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('s4','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('s5','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC'),\
 ('s6','ATAGTAGGTGATTGCGAAGACCGCGGAACCGGGACCTAGCACCCAGCCTGTACCGAGGGATGGGGAGCTGTGGCGGTCCACCGACGACCCTTTGTGACAGCCGATTCCTACAATCCCAGCAACTGCAATGATCCACTCTAGTCGGCATAACCGGGAATCGTTAACCTGGTAGGGTTCTCTACGTCTGAGTCTACAGCCCAGAGCAGTCAGGCTACTATACGGTTTGCTGCATTGCATAGGCATCGGTCGCGGGCACTCCTCGCGGTTTCAGCTAGGGTTTAAATGGAGGGTCGCTGCATGAGTATGCAAATAGTGCCACTGCTCTGATACAGAGAAGTGTTGATATGACACCTAAGACCTGGTCACAGTTTTAACCTGCCTACGCACACCAGTGTGCTATTGATTAACGATATCGGTAGACACGACCTTGGTAACCTGACTAACCTCATGGAAAGTGACTAGATAAATGGACCGGAGCCAACTTTCACCCGGAAAACGGACCGACGAATCGTCGTAGACTACCGATCTGACAAAATAAGCACGAGGGAGCATGTTTTGCGCAGGCTAGCCTATTCCCACCTCAAGCCTCGAGAACCAAGACGCCTGATCCGGTGCTGCACGAAGGGTCGCCTCTAGGTAAGGAGAGCTGGCATCTCCAGATCCGATATTTTACCCAACCTTTGCGCGCTCAGATTGTTATAGTGAAACGATTTAAGCCTGAACGGAGTTCCGCTCCATATGTGGGTTATATATGTGAGATGTATTAACTTCCGCAGTTGTCTCTTTCGGTGCAGTACGCTTGGTATGTGTCTCAAATAATCGGTATTATAGTGATCTGAGAGGTTTTAAG')],aligned=False)
 
test_refseq_coll = LoadSeqs(data=[\
 ('AY800210','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('EU883771','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('EF503699','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('DQ260310','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('EF503697','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC')],aligned=False)


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
