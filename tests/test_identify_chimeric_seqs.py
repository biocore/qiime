#!/usr/bin/env python
# File created on 06 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso","Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from os.path import exists, split, splitext
from cogent import LoadSeqs, DNA
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_file
from qiime.identify_chimeric_seqs import BlastFragmentsChimeraChecker,\
    chimeraSlayer_identify_chimeras, parse_CPS_file,\
    get_chimeras_from_Nast_aligned

class BlastFragmentsChimeraCheckerTests(TestCase):
    """ """
    
    def setUp(self):
        """ """    
        self.id_to_taxonomy_fp = get_tmp_filename(\
         prefix='BlastFragmentsChimeraCheckerTests_',suffix='.txt') 
        self.input_seqs_fp = get_tmp_filename(\
         prefix='BlastFragmentsChimeraCheckerTests_',suffix='.fasta')
        self.reference_seqs_fp = get_tmp_filename(\
         prefix='BlastFragmentsChimeraCheckerTests_',suffix='.fasta')
        
        self._paths_to_clean_up =\
         [self.id_to_taxonomy_fp,self.input_seqs_fp,self.reference_seqs_fp] 
        
        open(self.id_to_taxonomy_fp,'w').write(id_to_taxonomy_string)
        open(self.input_seqs_fp,'w').write(test_seq_coll.toFasta())
        open(self.reference_seqs_fp,'w').write(test_refseq_coll.toFasta())
        
        self.bcc = None
        
        self.expected1 = [\
            ('c1',\
                ['Archaea;Euryarchaeota;Halobacteriales;uncultured',\
                 'Archaea;Crenarchaeota;uncultured;uncultured'])]
                 
    def tearDown(self):
        remove_files(self._paths_to_clean_up)
        if self.bcc != None:
            self.bcc.cleanUp()

    def test_init_creates_db(self):
        """BlastFragmentsChimeraCheckerTests: blast database created/cleaned
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        self.assertEqual(self.bcc.num_fragments,3)
        # blast db created
        db_name = split(self.reference_seqs_fp)[1]
        db_fp = '/tmp/%s.nsq' % db_name
        self.assertTrue(exists(db_fp))
        self.bcc.cleanUp()
        # db is cleaned up
        self.assertFalse(exists(db_fp))
        
        
        
    def test_init_creates_db(self):
        """BlastFragmentsChimeraChecker: parameters initialized correctly
        """        
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp,\
                  'num_fragments':2}
        self.bcc = BlastFragmentsChimeraChecker(params)
        self.assertEqual(self.bcc.Params['num_fragments'],2)
        
    def test_function_w_preexisting_blastdb(self):
        blast_db, db_files_to_remove = \
         build_blast_db_from_fasta_file(test_refseq_coll.toFasta().split('\n'))
        self._paths_to_clean_up += db_files_to_remove
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':None,\
                  'blast_db':blast_db,\
                  'num_fragments':2}
        self.bcc = BlastFragmentsChimeraChecker(params)
        actual = list(self.bcc(self.input_seqs_fp))
        self.assertEqual(actual,self.expected1)
        
    def test_fragmentSeq_even_len_frags(self):
        """BlastFragmentsChimeraChecker:_fragment_seq fns when frags are evenly divisible
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        
        # default: 3 frags
        self.assertEqual(self.bcc._fragment_seq('ACCGTTATATTT'),\
            ['ACCG','TTAT','ATTT'])
            
        self.bcc.Params['num_fragments'] = 2
        self.assertEqual(self.bcc._fragment_seq('ACCGTTATATTT'),\
            ['ACCGTT','ATATTT'])
            
        self.bcc.Params['num_fragments'] = 4
        self.assertEqual(self.bcc._fragment_seq('ACCGTTATATTT'),\
            ['ACC','GTT','ATA','TTT'])
            
    def test_fragment_seq_uneven_len_frags(self):
        """BlastFragmentsChimeraChecker:_fragment_seq fns when frags aren't evenly divisible
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        
        # default: 3 frags
        self.assertEqual(self.bcc._fragment_seq('ACCGTTATATTTC'),\
            ['ACCGT','TATA','TTTC'])
        self.assertEqual(self.bcc._fragment_seq('ACCGTTATATTTCC'),\
            ['ACCGT','TATAT','TTCC'])
        
    def test_is_chimeric(self):
        """BlastFragmentsChimeraChecker: _is_chimeric functions as expected
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        # Basic cases, fewer than 4 assignments deep
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD','AB;CE','AB;CD']))
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD','AB;CD','AB;CD']))
        
        # Basic cases, exactly 4 assignments deep
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;GH']))
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
            
        # Leading/trialing spaces in taxa don't affect results
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF; GH ','AB ;CD;EF;GH','AB;CD; EF;GH']))
            
        # More complex cases -- only first four levels are considered
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH;IJ','AB;CD;EF;GH']))
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH;LM;XY','AB;CD;EF;GH;IJ','AB;CD;EF;GH;JK']))
            
        # unlikely case (but possible) where 5th level is identical,
        # but 4th is not
        self.bcc.Params['taxonomy_depth'] = 3
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH;XY','AB;CD;EF;JK;XY','AB;CD;EF;GH;XY']))
        self.bcc.Params['taxonomy_depth'] = 4
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH;XY','AB;CD;EF;JK;XY','AB;CD;EF;GH;XY']))
        self.bcc.Params['taxonomy_depth'] = 5
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH;XY','AB;CD;EF;JK;XY','AB;CD;EF;GH;XY']))
            
    def test_is_chimeric_no_blast_hit(self):
        """BlastFragmentsChimeraChecker: no blast hit is lack of evidence for chimeric
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        # 'No blast hit' values ignored (lack of evidence to claim sequence
        # is chimeric) 
        self.bcc.Params['taxonomy_depth'] = 3
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF','No blast hit','AB;CD;EF']))
        self.assertFalse(self.bcc._is_chimeric(\
            ['No blast hit','No blast hit','No blast hit']))
        self.assertFalse(self.bcc._is_chimeric(\
            ['No blast hit','AB;CD;EF','AB;CD;EF']))
            
        # Still called a chimera though if other assignments suggest that it is
        self.assertTrue(self.bcc._is_chimeric(\
            ['No blast hit','AB;CD;EF','AB;CD;XY']))
            
    def test_is_chimeric_alt_depth(self):
        """BlastFragmentsChimeraChecker: _is_chimeric functions as expected with alt depth
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        
        # chimeric at depth 5
        self.bcc.Params['taxonomy_depth'] = 5
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
        # chimeric at depth 4
        self.bcc.Params['taxonomy_depth'] = 4
        self.assertTrue(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
        # non-chimeric at depth 3
        self.bcc.Params['taxonomy_depth'] = 3
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
        # non-chimeric at depth 2
        self.bcc.Params['taxonomy_depth'] = 2
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
        # non-chimeric at depth 1
        self.bcc.Params['taxonomy_depth'] = 1
        self.assertFalse(self.bcc._is_chimeric(\
            ['AB;CD;EF;GH','AB;CD;EF;GH','AB;CD;EF;XY']))
        
        
    def test_get_taxonomy(self):
        """BlastFragmentsChimeraChecker: getTaxonomy functions with full-length sequence
        
            (Just testing the input/output here as the core functionality
             is tested in Qiime/tests/test_assign_taxonomy.py.)
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp}
        self.bcc = BlastFragmentsChimeraChecker(params)
        s1 = test_seq_coll.getSeq('s1')
        actual = self.bcc._get_taxonomy(str(s1))
        expected = "Archaea;Euryarchaeota;Halobacteriales;uncultured"
        self.assertEqual(actual,expected)
        
    def test_call(self):
        """BlastFragmentsChimeraChecker: correct assignments on known results
        """
        params = {'id_to_taxonomy_fp':self.id_to_taxonomy_fp,\
                  'reference_seqs_fp':self.reference_seqs_fp,\
                  'num_fragments':2}
        self.bcc = BlastFragmentsChimeraChecker(params)
        
        actual = list(self.bcc(self.input_seqs_fp))
        self.assertEqual(actual,self.expected1)
        

class ChimeraSlayerChimeraCheckerTests(TestCase):
    """ """
    
    def setUp(self):
        """ """   

        self.files_to_remove=[] 
        test_seqs_fp = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp)
        fh = open(test_seqs_fp,"w")
        fh.write(chimeras)
        fh.close()

        self.query_fp= test_seqs_fp
         
    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_chimeraSlayer_identify_chimeras_default_ref_db(self):
        """chimeraSlayer_identify_chimeras works with default reference DB"""

        # Real chimeras are identified as such
        observed = chimeraSlayer_identify_chimeras(self.query_fp)
        self.assertEqual(observed, [(chimera_id, parent_ids)])

        
    def test_chimeraSlayer_identify_chimeras_user_db(self):
        """chimeraSlayer_identify_chimeras works with user provided DB"""

        # set up DB
        ref_db_fp = get_tmp_filename(prefix="test_chimera_slayer_ref_db_")
        fh = open(ref_db_fp,"w")
        fh.write(ref_db)
        fh.close()

        ref_db_nast_fp = get_tmp_filename(prefix="test_chimera_slayer_nast_db_")
        fh_nast = open(ref_db_nast_fp,"w")
        fh_nast.write(ref_db_nast)
        fh_nast.close()

        self.files_to_remove.extend([ref_db_fp, ref_db_nast_fp])

        # Real chimeras are identified as such
        observed = chimeraSlayer_identify_chimeras(self.query_fp)
        self.assertEqual(observed, [(chimera_id, parent_ids)])

class ChimeraSlayer_app_tests(TestCase):
    
    def setUp(self):
        self.files_to_remove = []

    def tearDown(self):
        remove_files(self.files_to_remove,error_on_missing=False)

    def test_parse_CPS_files(self):
        """parse_CPS_files parse CPS format correctly"""
        
        lines=["\t".join(["ChimeraSlayer","chimera_X92624","7000004131495956",
                         "S000469847","1.0360","99.35","100","0.9354","89.70",
                         "0","YES","NAST:4595-4596","ECO:941-942"]),
                "\t".join(["ChimeraSlayer","Test","7000004131495956",
                         "S000469847","1.0360","99.35","100","0.9354","89.70",
                          "0","NO","NAST:4595-4596","ECO:941-942"])
               ]
        
        expected = [('chimera_X92624', ['7000004131495956', 'S000469847'])]
        self.assertEqual(parse_CPS_file(lines), expected)

        #Bad file raises error
        self.assertRaises(ValueError,parse_CPS_file, [""])
        self.assertRaises(ValueError,parse_CPS_file, ["Broken"])
        
    def test_get_chimeras_from_Nast_aligned(self):
        """get_chimeras_from_Nast_aligned identifies chimeras using default_db."""

        #empty input gives empty output
        seqs = ""
        test_seqs_fp = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp)
        fh = open(test_seqs_fp,"w")
        fh.write(seqs)
        fh.close()

        observed = get_chimeras_from_Nast_aligned(test_seqs_fp)
        self.assertEqual(observed, [])
        
        # no chimeras give empty output
        seqs = """>test1
GTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCGACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACCTCTTTCAGCAGGGACGAAGCGTAAGTGACGGTACCTGCAGAAGAAGCGCCGGCCAACTACGTGCCAGCAGCCGCGGTAAGAC
"""
        test_seqs_fp2 = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp2)
        fh = open(test_seqs_fp2, "w")
        fh.write(seqs)
        fh.close()
        observed = get_chimeras_from_Nast_aligned(test_seqs_fp2)
        self.assertEqual(observed, [])
        
        # Real chimeras are identified as such
        test_seqs_fp3 = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp3)
        fh = open(test_seqs_fp3,"w")
        fh.write(chimeras)
        fh.close()

        observed = get_chimeras_from_Nast_aligned(test_seqs_fp3)
        self.assertEqual(observed, [(chimera_id, parent_ids)])

    def test_get_chimeras_from_Nast_aligned_user_db(self):
        """get_chimeras_from_Nast_aligned identifies chimeras using user db.
        """

        # set up DB
        ref_db_fp = get_tmp_filename(prefix="test_chimera_slayer_ref_db_")
        fh = open(ref_db_fp,"w")
        fh.write(ref_db)
        fh.close()

        ref_db_nast_fp = get_tmp_filename(prefix="test_chimera_slayer_nast_db_")
        fh_nast = open(ref_db_nast_fp,"w")
        fh_nast.write(ref_db_nast)
        fh_nast.close()

        self.files_to_remove.extend([ref_db_fp, ref_db_nast_fp])

        #empty input gives empty output
        seqs = ""
        test_seqs_fp = get_tmp_filename(prefix="test_chimera_slayer")

        self.files_to_remove.append(test_seqs_fp)
        fh = open(test_seqs_fp,"w")
        fh.write(seqs)
        fh.close()

        observed = get_chimeras_from_Nast_aligned(test_seqs_fp, ref_db_nast_fp, ref_db_fp)
        self.assertEqual(observed, [])
        
        # no chimeras give empty output
        seqs = """>test1
GTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCGACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACCTCTTTCAGCAGGGACGAAGCGTAAGTGACGGTACCTGCAGAAGAAGCGCCGGCCAACTACGTGCCAGCAGCCGCGGTAAGAC
"""

        test_seqs_fp2 = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp2)
        fh = open(test_seqs_fp2,"w")
        fh.write(seqs)
        fh.close()
        observed = get_chimeras_from_Nast_aligned(test_seqs_fp2, ref_db_nast_fp, ref_db_fp)
        self.assertEqual(observed, [])
        
        # Real chimeras are identified as such
        test_seqs_fp3 = get_tmp_filename(prefix="test_chimera_slayer")
        self.files_to_remove.append(test_seqs_fp3)
        fh = open(test_seqs_fp3,"w")
        fh.write(chimeras)
        fh.close()

        observed = get_chimeras_from_Nast_aligned(test_seqs_fp3, ref_db_nast_fp, ref_db_fp)
        self.assertEqual(observed, [(chimera_id, parent_ids)])
            

id_to_taxonomy_string = \
"""AY800210\tArchaea;Euryarchaeota;Halobacteriales;uncultured
EU883771\tArchaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.
EF503699\tArchaea;Crenarchaeota;uncultured;uncultured
DQ260310\tArchaea;Euryarchaeota;Methanobacteriales;Methanobacterium
EF503697\tArchaea;Crenarchaeota;uncultured;uncultured"""

        
test_refseq_coll = LoadSeqs(data=[\
 ('AY800210','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('EU883771','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('EF503699','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('DQ260310','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('EF503697','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC')],aligned=False)
 
test_seq_coll = LoadSeqs(data=[\
 ('s1','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('s2','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('s3','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('s4','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('s5','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC'),\
 ('s6','ATAGTAGGTGATTGCGAAGACCGCGGAACCGGGACCTAGCACCCAGCCTGTACCGAGGGATGGGGAGCTGTGGCGGTCCACCGACGACCCTTTGTGACAGCCGATTCCTACAATCCCAGCAACTGCAATGATCCACTCTAGTCGGCATAACCGGGAATCGTTAACCTGGTAGGGTTCTCTACGTCTGAGTCTACAGCCCAGAGCAGTCAGGCTACTATACGGTTTGCTGCATTGCATAGGCATCGGTCGCGGGCACTCCTCGCGGTTTCAGCTAGGGTTTAAATGGAGGGTCGCTGCATGAGTATGCAAATAGTGCCACTGCTCTGATACAGAGAAGTGTTGATATGACACCTAAGACCTGGTCACAGTTTTAACCTGCCTACGCACACCAGTGTGCTATTGATTAACGATATCGGTAGACACGACCTTGGTAACCTGACTAACCTCATGGAAAGTGACTAGATAAATGGACCGGAGCCAACTTTCACCCGGAAAACGGACCGACGAATCGTCGTAGACTACCGATCTGACAAAATAAGCACGAGGGAGCATGTTTTGCGCAGGCTAGCCTATTCCCACCTCAAGCCTCGAGAACCAAGACGCCTGATCCGGTGCTGCACGAAGGGTCGCCTCTAGGTAAGGAGAGCTGGCATCTCCAGATCCGATATTTTACCCAACCTTTGCGCGCTCAGATTGTTATAGTGAAACGATTTAAGCCTGAACGGAGTTCCGCTCCATATGTGGGTTATATATGTGAGATGTATTAACTTCCGCAGTTGTCTCTTTCGGTGCAGTACGCTTGGTATGTGTCTCAAATAATCGGTATTATAGTGATCTGAGAGGTTTTAAG'),\
 ('c1','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG')],aligned=False)



#Test data taken from the ChimeraSlayer sample data

chimera_id = 'chimera_X92624|S000015368_d10.86_AF282889|S000390572_d11.18_nc4580_ec939'
parent_ids =['7000004131495956','S000469847']


ref_db_nast=""">7000004131495956
............................................................
......................................................T-GA--
T-CC-T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----T-T--A
ACACA-T-GC-A-AGT-CGA-GCGGAAA-------G-G---CC-----------------
------------------------------------------------------------
--------------CTTCGGGG--------------------------------------
------------------------------------------------------------
--------T--A-CTC--G--AG-C-GG-C-GA-A--C-------------GGG-TGAGT
-A--AC-AC-G-T-G-AG---CAA--C-CT-G--C-C-CTA--GG-C-------------
-----------------------------------------------------T-TT---
-GGG-AT-AA-CCC-------------------------C-G-G----------------
-------GAA-A---CCG-GGG-CTAA-TA---CC-G--A-AT-A---------------
------------------G-GACC-T--T--C-----------------G----GT-C--
------------------------------------------------------------
---------------------------------------------------------G-C
A-T---------------------------------------------------------
------------------------------------------------------------
-------------------G-A--C-T---------------G-TT-G-G-T-G------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-------GAAA--G-T--------------------------------------------
------------------------------------------------------------
--------------------------------TTT-------------------------
------------------------------------------------------------
-----------------------------------------------T---C-G------
--------G----C-C---T-GG-G---AT---G-G-----G-CTC-GCG--G-CCT--A
------TC--A--G-CT-T----G---TTGG-T-G-GG-G-T----GAT-GG-C-C-T-A
CCA--A-GG-C-G--A-CG-A------------CGG-G-T------AG-CC-G-G-CCT-
G-AG----A--GG-GC--G-AC-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G
-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TT
GCA-C-AA-T-GG--GC-GG-A----A-G-CC-T-GA-TG-CA-GCGA-CGCC-G-CG-T
---G-A-G--G--GA-T-G--A--C-G-G-CC-----TT-CG---------G-G-T-T-G
-T--A---AA-C-CTC--------TT-TC-A-G--C-AGG----GA-C--G--AAGCGTA
A--------------------G--T---------------------------GA----CG
-GT-A-C-CT-G-CA-G---------AA-----------GAAGC-GCC-GG-C-CAA---
C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AG--AC---GT-AG-GGC
-GCG-A-G-CG-TTGT-C-CGG-AT-TT-A--T-T--GGGC-GTA----AA-GAGC-TC-
-G-TA-G-G-C-G------------G--C-TT-G-T-C-GC----G-T-C-G---A-CTG
-TG-A-AA-AC--CC-GCA-G---------------------------------------
-----------------------------CT-C-AA------------------------
-------------------------------------------------CT-G-C-GG-G
C-C----T-G-C-A-G-T--------C--GA-T-A-C-G-G-GCA--G-G-C--------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-T-A-G-A-G-T-T-----C-GG--TA-G-G------------G-GA-G-AC-T----GG
--AATT-CCT-G-GT--GT-A-GCG-GTGAAA-TG-CGC-AGAT-A-TC-A-GGA--GG-
A-AC-A-CC-GG--T--G--GC-GAA-G--G-C---G----G--G-T-CTCTG------G
-GC-CG------------------------------------------------------
--------AT-A-C-T--GA--CG-----CT-GA-GG--A-G-CGA--AA-G-C------
--------G-TGGG-GAG-C-G-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-
G-T----C-CA--C-G-CTG-T-AAA--C-GTTG-GG--CG-CT---------A-GG--T
--G-T-GG-G-GG-G--C------------------------------------------
--------------------------------------------CTC-TC-C--------
------------------------------------------------------------
------------------------------------------------------------
----------------G-G-TTCT--C-T-G-T-GC-C------GC--A----GC-TAA-
-CG-C-A-T--T--AA-GC--G----C-CCC-GCC-T-G-GG-GAG-TA---CGG-----
C-C--G-C-A-A-GGC-T--AAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-
CCCG----C-A--C-A-A-GCG-GC-G--G--AG-CA-T--GC-GGA-TT-AATT-C-G-
ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-TGGGT-TT-G-AC-A-T-G-------
-----G-CCG-C---------------A-AA---A-C--CG--G-CA-G-A-G--A-TGT
C-G-G-G-T----------CC-------------------------------------T-
-TC-G------------------------------------------GG-----------
--GG-CGG---T--CA--------------------------------------------
-------C-A-G-G-T-GGTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G
--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC
-C-T-CG--TT--C-GATG--T-T-G-C-C---AG-C-G--C------------------
------------------------------------------------------------
------------------------------------------------------------
---------------------------------------------GTTATG---------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------G----C
-G------------G----G---G-A--CT---------------C-A-T-C-G-AA-G-
-AC-T-G-CCG--G-G------------------------------------G-T---CA
A----------------------------------C-T-C-G--G-A-GG-A--AGG-T-
-GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-C-C-C-CTT----AT-G--TC-C
-A-GG-GC-TT-CAC-GCATG-C--TA--CAATG---G-CCGG-T-A--C-AAT-GG-GC
------------------------------------------------------------
--------------------------------------T-G-C-G-A--T-ACCG-T--G
---------------------------------------A-GG-T-G-----------G-
-A-G-CG---A----------A--TCC-C------A-A-AAAGC-CG-G-T-C-T-CAG-
TTC--------GGA-T-CGGGG-TC--T-GCAA-CT-C----------------------
------------------------------------------------------------
---------------G-ACCCC-G-T-G-AA-G-TC-GGAGT-CG-C-TA--G-TA-AT-
C-G-C----AGA-TC-A-G-CA------AC--GCT-GC-G-GT-G-AAT-ACGT-T-CCC
GGGCCT-TGTA----CACACCG-CCC-GTC-----A---CG--TCA-CG-AA-A--G---
TCG-G-CA-AC-ACC--C-GAA------G--C-CGG-TG-G-C-C-C-AA-C-C-C----
-------------------------------------------------------T-TG-
T-----------------------------------------------------------
----------------------------------------G--GA----GG-G--A---G
C-CGT--CG--AAG-G----T-GGG-GC-TGG------------------------CG--
ATT-GGGA-CG-AAG-TCGTAACAA-GGTAG-CCGT-ACCGG..................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
..
>S000469847
............................................................
............................................................
............................................................
.....................................a---ag-g---------------
----------------c-c-ct--------------------------------------
--------------tcg-g-----------------------------------------
-------------------------------------------g-g--------------
-------ca--c--ac--g--aggc-gg-c-ga-a--c-------------ggg-tgagt
-a--ac-actg-t-g-gg---tga--t-ct-g--c-c-tcg--ca-c-------------
-----------------------------------------------------t-tc---
-ggg-at-aa-gcc-------------------------t-g-g----------------
-------gaa-a---ctg-ggt-ctaa-ta---cc-g--g-at-a---------------
------------------t-ga-c-c--t--g-----------------ct---gt-c--
------------------------------------------------------------
---------------------------------------------------------g-c
a-t---------------------------------------------------------
------------------------------------------------------------
-----------------g-g-c--g-g--------------tg--g-g---t-g------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-------gaaa--g-a--------------------------------------------
------------------------------------------------------------
--------------------------------ttt-------------------------
------------------------------------------------------------
-----------------------------------------------a--tc-g------
--------g----t-g---c-ga-g---at---g-g-----g-ccc-gcg--g-cct--a
------tc--a--g-ct-t----g---ttgg-t-g-gg---t----aat-gg-c-c-t-a
cca--a-gg-c-g--a-cg-a------------cgg-g-t------ag-cc-g-a-cct-
g-ag----a--ggggt--g-ac-c-gg-ccac-a-ctggg--a-c-tg-a-ga-c-ac-g
-g-cccaga-ctcc-tac-g--g-g-a-g-gc-a-gc-a-g-tg---gg-g-a-ata-tt
gca-c-aa-t-gg--gc-ga-a----a-g-cc-t-ga-tg-ca-gcga-cgcc-g-cg-t
---g-a-g--g--ga-t-g--a--c-g-g-cc-----tt-cg---------g-g-t-t-g
-t--a---aa-c-ctc--------tt-tc-g-a--c-agg----ga-c--g---aagcgc
a---a-------------g--t---------------------------------ga-cg
-gt-a-c-ct-g-ta-g---------aa-----------gaagc-acc-gg-c-caa---
c--t-acgt--gcca--g-c---a--gccg---c-gg--ta-at--ac---gt-ag-ggt
-gcg-a-g-cg-ttgt-c-cgg-aa-tt-a--c-t--gggc-gta----aa-gagc-tt-
-g-ta-g-g-c-g------------g--t-tt-g-t-c-gc----g-t-c-g---t-ctg
-tg-a-aa-ac--tc-aca-g---------------------------------------
-----------------------------ct-c-aa------------------------
-------------------------------------------------ct-g-t-ga-g
c-t----t-g-c-a-g-g--------c--ga-t-a-c-g-g-gcag---a-c--------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-t-t-g-a-g-t-a-----c-tt--ca-g-g------------g-ga-g-ac-t----gg
--aatt-cct-g-gt--gt-a-gcg-gtgaaa-tg-cgc-agat-a-tc-a-gga--gg-
a-ac-a-cc-gg--t--g--gc-gaa-g--g-c---g----g--g-t-ctctg------g
-ga-ag------------------------------------------------------
--------ta-a-c-t--ga--cg-----ct-ga-ga--a-g-cga--aa-g-c------
--------g-tggg-tag-c-g-aaca--gg-atta-g-ata-c-----cc-t-g-gta-
g-t----c-ca--c-g-ccg-t-aaa--c-ggtg-gg--ta-ct---------a-gg--t
--g-t-gg-g-tt-t-----c---------------------------------------
-------------------------------------------cttc-ca----------
------------------------------------------------------------
------------------------------------------------------------
------------c---g-g-g-at--c-c-g-t-gc-c------gt--a----gc-taa-
-cg-c-a-t--t--aa-gt--a----c-ccc-gcc-t-g-gg-gag-ta---cgg-----
c-c--g-c-a-a-ggc-t--aaa-actc-aaa---------ggaa-ttg-acggg-g-g-
cccg----c-a--c-a-a-gcg-gc-g--g--ag-ca-t--gt-gga-tt-aatt-c-g-
atg-caac-g-cg-a-ag-a-a-cc-tt-a-cc-tgggt-tt-g-ac-a-t-a-------
-----c-acc-g-g-------------a-aa-c-c-t--gc--a-ga-g-a-t--g-t-a
--g-g-c-c----------cc-------------------------------------c-
-tt-g-----------------------------------------tgg-----------
--tc-ggt---g--ta--------------------------------------------
-------c-a-g-g-t-ggtg-ca-tgg-ct--gtc-gtc-a-gc-tc---g-tg-tc-g
--tga-ga-tgt-t-gg-g-tt-aa-gt-cccgc-aa--------c-gag-cgc-a-acc
-c-t-ta--tc--t-tatg--t-t-g-c-c---ag-c-g--c--gt-aa-----------
------------------------------------------------------------
------------------------------------------------------------
----------------------------------------------t-------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
--------------------------------------------------g---g----c
-g------------g----g---g-a--ct---------------c-g-t-g-a-ga-g-
-ac-t-g-ccg--g-g------------------------------------g-t---ca
a----------------------------------c-t-c-g--g-a-gg-a--agg-t-
-gggg-a-cgac-gtc--aagt-c---atc-a-t-g-c-c-c-ctt----at-g--tc-c
-a-gg-gc-tt-cac-acatg-c--ta--caatg---g-ccgg-t-a--c-aga-gg-gc
------------------------------------------------------------
--------------------------------------t-g-c-g-a--t-accg-t--g
---------------------------------------a-gg-t-g-----------g-
-a-g-cg---a----------a--tcc-c------t-t-aaagc-cg-g-t-c-t-cag-
ttc--------gga-t-cgggg-tc--t-gcaa-ct-c----------------------
------------------------------------------------------------
---------------g-acccc-g-t-g-aa-g-tt-ggagt-cg-c-ta--g-ta-at-
c-g-c----aga-tc-a-g-ca------ac--gct-gc-g-gt-g-aat-acgt-t-ccc
gggcct-tgta----cacaccg-ccc-gtc-----a---cg--tca-tg-aa-a--g---
tcg-g-ta-ac-acc--c-gaa------g--c-cgg-tg-g-c-c-t-aa-c-c------
-----------------------------------------------------cc-ttg-
tg----------------------------------------------------------
-------------------------------------------gg-a--gg-g--a---g
c-cgt--cg--aag-g----t-ggg-at-cgg------------------------cg--
att-ggga-cg-aag-tcgtaacaa-ggtag-ccgt-accggaa-ggtg-cggc-tggat
cccct.......................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
............................................................
..
"""

ref_db = """>7000004131495956
TGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGAAAGGCCCTTCGGGGTACTCGAGCGGCGAACGGGTGAGTAACACGTGAGCAACCTGCCCTAGGCTTTGGGATAACCCCGGGAAACCGGGGCTAATACCGAATAGGACCTTCGGTCGCATGACTGTTGGTGGAAAGTTTTTCGGCCTGGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGGGGTGATGGCCTACCAAGGCGACGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCGACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACCTCTTTCAGCAGGGACGAAGCGTAAGTGACGGTACCTGCAGAAGAAGCGCCGGCCAACTACGTGCCAGCAGCCGCGGTAAGACGTAGGGCGCGAGCGTTGTCCGGATTTATTGGGCGTAAAGAGCTCGTAGGCGGCTTGTCGCGTCGACTGTGAAAACCCGCAGCTCAACTGCGGGCCTGCAGTCGATACGGGCAGGCTAGAGTTCGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGCCGATACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGTTGGGCGCTAGGTGTGGGGGGCCTCTCCGGTTCTCTGTGCCGCAGCTAACGCATTAAGCGCCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGTTTGACATGGCCGCAAAACCGGCAGAGATGTCGGGTCCTTCGGGGGCGGTCACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGTTCGATGTTGCCAGCGCGTTATGGCGGGGACTCATCGAAGACTGCCGGGGTCAACTCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTATGTCCAGGGCTTCACGCATGCTACAATGGCCGGTACAATGGGCTGCGATACCGTGAGGTGGAGCGAATCCCAAAAAGCCGGTCTCAGTTCGGATCGGGGTCTGCAACTCGACCCCGTGAAGTCGGAGTCGCTAGTAATCGCAGATCAGCAACGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACGTCACGAAAGTCGGCAACACCCGAAGCCGGTGGCCCAACCCTTGTGGAGGGAGCCGTCGAAGGTGGGGCTGGCGATTGGGACGAAGTCGTAACAAGGTAGCCGTACCGG
>S000469847
aaggcccttcggggcacacgaggcggcgaacgggtgagtaacactgtgggtgatctgcctcgcacttcgggataagcctgggaaactgggtctaataccggatatgacctgctgtcgcatggcggtgggtggaaagatttatcggtgcgagatgggcccgcggcctatcagcttgttggtgggtaatggcctaccaaggcgacgacgggtagccgacctgagaggggtgaccggccacactgggactgagacacggcccagactcctacgggaggcagcagtggggaatattgcacaatgggcgaaagcctgatgcagcgacgccgcgtgagggatgacggccttcgggttgtaaacctctttcgacagggacgaagcgcaagtgacggtacctgtagaagaagcaccggccaactacgtgccagcagccgcggtaatacgtagggtgcgagcgttgtccggaattactgggcgtaaagagcttgtaggcggtttgtcgcgtcgtctgtgaaaactcacagctcaactgtgagcttgcaggcgatacgggcagacttgagtacttcaggggagactggaattcctggtgtagcggtgaaatgcgcagatatcaggaggaacaccggtggcgaaggcgggtctctgggaagtaactgacgctgagaagcgaaagcgtgggtagcgaacaggattagataccctggtagtccacgccgtaaacggtgggtactaggtgtgggtttccttccacgggatccgtgccgtagctaacgcattaagtaccccgcctggggagtacggccgcaaggctaaaactcaaaggaattgacgggggcccgcacaagcggcggagcatgtggattaattcgatgcaacgcgaagaaccttacctgggtttgacatacaccggaaacctgcagagatgtaggcccccttgtggtcggtgtacaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaacccttatcttatgttgccagcgcgtaatggcggggactcgtgagagactgccggggtcaactcggaggaaggtggggacgacgtcaagtcatcatgccccttatgtccagggcttcacacatgctacaatggccggtacagagggctgcgataccgtgaggtggagcgaatcccttaaagccggtctcagttcggatcggggtctgcaactcgaccccgtgaagttggagtcgctagtaatcgcagatcagcaacgctgcggtgaatacgttcccgggccttgtacacaccgcccgtcacgtcatgaaagtcggtaacacccgaagccggtggcctaaccccttgtgggagggagccgtcgaaggtgggatcggcgattgggacgaagtcgtaacaaggtagccgtaccggaaggtgcggctggatcccct
"""

chimeras=""">chimera_X92624|S000015368_d10.86_AF282889|S000390572_d11.18_nc4580_ec939 
.................................................................................................................TT-GA--T-CC-T-G-GCGC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----T-T--AACACA-T-GC-A-AGT-CGA-GCG
GAAA-------G-G---CC-------------------------------------------------------------------------------------------CTTCGGGG--------------------------------------------------------------------------------------
--------------------T--A-CTC--G--AG-C-GG-C-GA-A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-AG---CAA--C-CT-C--C-C-CTA--GG-C------------------------------------------------------------------T-TT----GGG-AT-AA-C
CC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAA-TA---CC-G--A-AT-A---------------------------------G-GACC-G--T--C-----------------G----AT-C--------------------------------------
---------------------------------------------------------------------------------G-CA-T---------------------------------------------------------------------------------------------------------------------
-------------------G-A--T-C---------------G-TT-G-G-T-G------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------GAAA--G-T--------------------------------------------------------------------------------------------
--------------------------------------------TTT------------------------------------------------------------------------------------------------------------------------------------T---C-G--------------G---
-C-C---T-GG-G---AT---G-G-----G-CTC-GCG--G-CCT--A------TC--A--G-CT-T----G---TTGG-T-G-GG-G-T----GAT-GG-C-C-T-ACCA--A-GG-C-G--A-CG-A------------CGG-G-T------AG-CC-G-G-CCT-G-AG----A--GG-GC--G-AC-C-GC-CCAC-A-C
TGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GC-GG-A----A-G-CC-T-GA-TG-CA-GCGA-CGCC-G-CG-T---G-A-G--G--GA-T-G--A--C-G-G-CC-----TT-CG---------G-G-T-T-G
-T--A---AA-C-CTC--------TT-TC-A-G--C-AGG----GA-C--G--AAGCGTAA--------------------G--T---------------------------GA----CG-GT-A-C-CT-G-CA-G---------AA-----------GAAGC-GCC-GG-C-CAA---C--T-ACGT--GCCA--G-C---A
--GCCG---C-GG--TA-AG--AC---GT-AG-GGC-GCG-A-G-CG-TTGT-C-CGG-AT-TT-A--T-T--GGGC-GTA----AA-GAGC-TC--G-TA-G-G-C-G------------G--C-TT-G-T-C-GC----G-T-C-G---A-CTG-TG-A-AA-AC--CC-GCA-G---------------------------
-----------------------------------------CT-C-AA-------------------------------------------------------------------------CT-G-C-GG-GC-C----T-G-C-A-G-T--------C--GA-T-A-C-G-G-GCA--G-G-C--------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------T-A-G-A-G-T-T-----C-GG--TA-G-G------------G-GA-G-AC-T----GG--AATT-CCT-G-GT--GT-A-GCG-GTGAAA-TG-CGC-AGAT-A-TC-A-GGA--GG-A-AC-A-CC-GG--T--G--GC-GAA-G--G-C---
G----G--G-T-CTCTG------G-GC-CG--------------------------------------------------------------AT-A-C-T--GA--CG-----CT-GA-GG--A-G-CGA--AA-G-C--------------G-TGGG-GAG-C-G-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-
G-C----C-CA--C-G-CTG-C-AAA--C-GTTG-GG--CG-CT---------A-GG--T--G-T-GG-G-GG-G--C--------------------------------------------------------------------------------------CTC-TC-C--------------------------------
----------------------------------------------------------------------------------------------------------------G-G-TTCT--C-T-G-T-GC-C------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--G----C-CCC-GCC-T-G-GG-GAG-T
A---CGG-----C-C--G-C-A-A-GGC-T--AAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GC-G--G--AG-CA-T--GT-GGA-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-TGGGT-TT-G-AC-A-T-A------------C-ACC-G
-G-------------A-AA-C-C-T--GC--A-GA-G-A-T--G-T-A--G-GCC-C----------CC-------------------------------------C--TT-G-----------------------------------------TGG-------------TC-GGT---G--TA--------------------
-------------------------------C-A-G-G-T-GGTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TC--T-TATG--T-T-G-C-C---AG-C-G--C--GT-AA-----------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------T-------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------G---G----C-G------------G----G---G-A--CT---------------C-G-T-G-A-GA-G--AC-T-G-CCG-
-G-G------------------------------------G-T---CAA----------------------------------C-T-C-G--G-A-GG-A--AGG-T--GGGG-A-CGAC-GTC--AAGT-C---ATC-A-T-G-C-C-C-CTT----AT-G--TC-C-A-GG-GC-TT-CAC-ACATG-C--TA--CAATG--
-G-CCGG-T-A--C-AGA-GG-GC--------------------------------------------------------------------------------------------------T-G-C-G-A--T-ACCG-T--G---------------------------------------A-GG-T-G-----------G-
-A-G-CG---A----------A--TCC-C------T-T-AAAGC-CG-G-T-C-T-CAG-TTC--------GGA-T-CGGGG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACCCC-G
-T-G-AA-G-TT-GGAGT-CG-C-TA--G-TA-AT-C-G-C----AGA-TC-A-G-CA------AC--GCT-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CG--TCA-TG-AA-A--G---TCG-G-TA-AC-ACC--C-GAA------G--C-CGG-TG-G-C-C-T-
AA-C-C-----------------------------------------------------------CC-TTG-TG-----------------------------------------------------------------------------------------------------GG-A--GG-G--A---GC-CGT--CG--A
AG-G----T-GGG-AT-CGG------------------------CG--ATT-GGGA-CG-AAG-TCGTAACAA-GGTAG-CCGT-ACCGG..................................................................................................................
............................................................................................................................................................................................................
............................................................................................................................................................................................................
............................................................................................................................................................................................................
......................................................................................................................................
"""



if __name__ == "__main__":
    main()
