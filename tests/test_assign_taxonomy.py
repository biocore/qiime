#!/usr/bin/env python

"""Tests of code for assigning taxonomy"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project" 
#remember to add yourself if you make changes
__credits__ = ["Greg Caporaso", "Kyle Bittinger"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs
from cogent.app.util import get_tmp_filename
from os import remove
from pipe454.assign_taxonomy import TaxonAssigner, BlastTaxonAssigner,\
 RdpTaxonAssigner
 
def remove_files(list_of_filepaths,error_on_missing=True):
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError,\
         "Some filepaths were not accessible: %s" % '\t'.join(missing) 

 
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
        """ """        
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
        open(self.reference_seqs_fp,'w').write(test_refseq_coll.toFasta())
        
        self.expected1 = {\
         's1':('Archaea,Euryarchaeota,Halobacteriales,uncultured',None),
         's2':('Archaea,Euryarchaeota,Methanomicrobiales,Methanomicrobium et rel.',None),\
         's3':('Archaea,Crenarchaeota,uncultured,uncultured',None),\
         's4':('Archaea,Euryarchaeota,Methanobacteriales,Methanobacterium',None),\
         's5':('Archaea,Crenarchaeota,uncultured,uncultured',None),\
         's6':None}
        
    def tearDown(self):
        remove_files(self._paths_to_clean_up)

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
         "AY800210":"Archaea,Euryarchaeota,Halobacteriales,uncultured",\
         "EU883771":"Archaea,Euryarchaeota,Methanomicrobiales,Methanomicrobium et rel.",\
         "EF503699":"Archaea,Crenarchaeota,uncultured,uncultured",\
         "DQ260310":"Archaea,Euryarchaeota,Methanobacteriales,Methanobacterium",\
         "EF503697":"Archaea,Crenarchaeota,uncultured,uncultured"}
        self.assertEqual(p._parse_id_to_taxonomy_file(lines),expected)
        
    def test_map_ids_to_taxonomy(self):
        """Mapping sequence ids to taxonomy functions as expected
        """
        p = BlastTaxonAssigner({})
        id_to_taxonomy_map = {\
         "AY800210":"Archaea,Euryarchaeota,Halobacteriales,uncultured",\
         "EU883771":"Archaea,Euryarchaeota,Methanomicrobiales,Methanomicrobium et rel.",\
         "EF503699":"Archaea,Crenarchaeota,uncultured,uncultured",\
         "DQ260310":"Archaea,Euryarchaeota,Methanobacteriales,Methanobacterium",\
         "EF503697":"Archaea,Crenarchaeota,uncultured,uncultured"}
        
        hits = {'s1':("AY800210",1e-99),\
         's5':("EU883771",'weird confidence value'),
         's3':("DQ260310",42.),\
         's4':None}
        expected = {\
         's1':("Archaea,Euryarchaeota,Halobacteriales,uncultured",None),\
         's5':(
          'Archaea,Euryarchaeota,Methanomicrobiales,Methanomicrobium et rel.',\
          None),
         's3':("Archaea,Euryarchaeota,Methanobacteriales,Methanobacterium",None),\
         's4':None}
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
        
            This test is kind of a hack until I figure out a better way. The
             sequences in test_seq_coll are pulled directly from the sequences
             in the blast database that is being tested against. So, this 
             test requires that for each of the three query sequences, the 
             exact match in the database is returned with an E-value of 0.0. 
             All other blast results for each query sequence are ignored. 
        
        
        """
        p = BlastTaxonAssigner({})
        blast_db = 'ribosomal_v11'
        seq_coll_blast_results = p._get_blast_hits(blast_db,test_seq_coll)
        # mapping from identifier in test_seq_coll to the id of the sequence
        # in my ribosomal_v11 database (a silva derivative)
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
            self.assertTrue(expected_matches[seq_id] in blast_results_d, \
             msg='** May be missing blast or necessary blast database: %s'\
             % blast_db)
            # now check that the perfect match got a 0.0 e-value as it should
            # on this data
            self.assertEqual(blast_results_d[expected_matches[seq_id]],0.0,\
             msg='** May be missing blast or necessary blast database: %s'\
             % blast_db)
            
    def test_call_existing_blast_db(self):
        """BlastTaxonAssigner.__call__ functions w existing db
        """
        blast_db = 'ribosomal_v11'
        
        p = BlastTaxonAssigner({'blast_db':blast_db,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        actual = p(self.input_seqs_fp)
        
        self.assertEqual(actual,self.expected1,\
             msg='** May be missing blast or necessary blast database: %s'\
             % blast_db)
        
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
         
        result_path = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self._paths_to_clean_up.append(result_path) 

        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        
        actual = p(self.input_seqs_fp,result_path=result_path)
        
        # because the order of the lines is not guaranteed, check
        # that on parsing each line we end up with the expected data
        of = open(result_path)
        for line in of:
            fields = line.strip().split('\t')
            self.assertEqual(\
             (fields[1],fields[2]),(self.expected1[fields[0]][0],'None'))
        of.close()
        
        # Return value is None when result_path is provided (Not sure
        # if this is what we want yet, or if we would want both so 
        # results could be logged to file...)
        self.assertEqual(actual,None)
        
    def test_call_logs_run(self):
        """BlastTaxonAssigner.__call__ logs the run when expected
        """
        log_path = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self._paths_to_clean_up.append(log_path) 
        
        p = BlastTaxonAssigner({\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp,\
         'blast_db':'ribosomal_v11'})
        actual = p(self.input_seqs_fp,log_path=log_path)
        
        log_file = open(log_path)
        log_file_str = log_file.read()
        log_file.close()
        
        log_file_exp = ["BlastTaxonAssigner parameters:",\
         'Min percent identity:0.9','Application:blastn/megablast',\
         'Max E value:1e-30',\
         'Result path: None, returned as dict.',
         'blast_db:ribosomal_v11',\
         'id_to_taxonomy_filepath:%s' % self.id_to_taxonomy_fp,\
         '']
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
        """ Attmepting to train on-the-fly raises error
        """
        app = RdpTaxonAssigner({\
         'id_to_taxonomy_fp':'/some/path',
         'reference_sequences_fp':'/some/other/path'})
        self.assertRaises(NotImplementedError,app,self.tmp_seq_filepath)

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
 'X67228 some description': ('Root,Bacteria,Proteobacteria,Alphaproteobacteria,Rhizobiales,Rhizobiaceae,Rhizobium',0.95),\
 'EF503697': ('Root,Archaea,Crenarchaeota,Thermoprotei',0.88)
}

rdp_test1_expected_lines = [\
 "\t".join(["X67228 some description",\
  "Root,Bacteria,Proteobacteria,Alphaproteobacteria,Rhizobiales,Rhizobiaceae,Rhizobium",\
  "0.9"]),
 "\t".join(['EF503697','Root,Archaea,Crenarchaeota,Thermoprotei','0.8'])]

id_to_taxonomy_string = \
"""AY800210\tArchaea,Euryarchaeota,Halobacteriales,uncultured
EU883771\tArchaea,Euryarchaeota,Methanomicrobiales,Methanomicrobium et rel.
EF503699\tArchaea,Crenarchaeota,uncultured,uncultured
DQ260310\tArchaea,Euryarchaeota,Methanobacteriales,Methanobacterium
EF503697\tArchaea,Crenarchaeota,uncultured,uncultured"""

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
