#!/usr/bin/env python
#file test_check_id_map.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","William Walters"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "William Walters"
__email__ = "rob@spot.colorado.edu","william.a.walters@colorado.edu"
__status__ = "Release"


from collections import defaultdict
from numpy import array
from string import strip
from cogent.util.unit_test import TestCase, main
from StringIO import StringIO
from qiime.check_id_map import (find_diff_length, CharFilter, lwu, 
    DupChecker, SameChecker, 
    run_checks, filename_has_space, run_description_missing, adapt_dupchecker,
    sampleid_missing, blank_header, bad_char_in_header, pad_rows, 
    barcode_missing, description_missing, wrap_arrays,
    check_vals_by_type, check_vals_by_contains, check_field_types,
    check_same_length, check_bad_chars, check_mixed_caps,
    check_missing_descriptions, check_duplicate_descriptions,
    check_description_chars, process_id_map, get_primers_barcodes,
    check_primers_barcodes, check_missing_sampleIDs,
    check_dup_var_barcodes_primers
    )



class TopLevelTests(TestCase):
    """Tests of top-level functions"""
    def test_find_diff_length(self):
        """find_diff_length should find wrong length item"""
        fdl = find_diff_length
        self.assertEqual(fdl(['a','b','c']), [])
        self.assertEqual(fdl(['a','b','ccc','d']), [[2, 'ccc',3,'a',1]])

    def test_lwu(self):
        """lwu should case-convert, strip whitespace/underscores"""
        x = '  aBc_D e Fg\t g'
        self.assertEqual(lwu(x), 'abcdefgg')
        self.assertEqual(lwu('abc'), 'abc')
        self.assertEqual(lwu(''), '')

    def test_run_checks(self):
        """run_checks should run a series of checks on data"""
        def bad_if_upper(x, raw_data=None):
            if x.upper() == x:
                return x.lower(), 'X is uppercase'
            else:
                return x, ''
        def bad_if_short(x, raw_data=None):
            if len(x) < 10:
                return x+'-----', 'X is short'
            else:
                return x, ''
        def bad_if_lower(x, raw_data=None):
            if x.lower() == x:
                return x, 'X is lowercase'
            else:
                return x, ''

        checks = [(bad_if_lower, 'warning'), 
            (bad_if_upper, 'error'),
            (bad_if_short, 'error'), 
            ]
        
        problems = defaultdict(list)
        result = run_checks('ABC', checks, problems)
        self.assertEqual(problems, {'error':['X is uppercase', 'X is short']})
        self.assertEqual(result, 'abc-----')

    def test_filename_has_space(self):
        """filename_has_space should complain if space in filename"""
        self.assertEqual(filename_has_space('x.txt'), ('x.txt',''))
        self.assertEqual(filename_has_space('x .txt'), ('x_.txt', 
            'Filename may not contain spaces. '+ \
            'Please re-upload without spaces, e.g. x .txt -> x_.txt.'))

    def test_run_description_missing(self):
        """run_description_missing should complain if no run description"""
        self.assertEqual(run_description_missing('x'), ('x', ''))
        self.assertEqual(run_description_missing(''), 
            ('No run description supplied.',  \
             'Run description was not supplied, using default value.'))

    def test_adapt_dupchecker(self):
        """adapt_dupchecker should adapt DupChecker to correct API"""
        raw_dup_checker = adapt_dupchecker(lambda x:x, 'DupCheck')
        self.assertEqual(raw_dup_checker(['x','y']), (['x','y'],''))
        self.assertEqual(raw_dup_checker(['x','x']), (['x','x'], \
            "DupChecker 'DupCheck' found the following possible duplicates. If these metadata should have the same name, please correct.:\nGroup\tOriginal names\nx\tx, x\n"))

    def test_sampleid_missing(self):
        """sampleid_missing should complain if sampleid missing"""
        self.assertEqual(sampleid_missing(['#SampleID', 'x','y']),\
            (['#SampleID', 'x','y'],''))
        res = sampleid_missing(['x','y','z'])
        self.assertEqual(res[1][:10], 'SampleID f')

    def test_blank_header(self):
        """blank_header should complain if whitespace fields in header"""
        self.assertEqual(blank_header(['x','y','z']), (['x','y','z'],''))
        self.assertEqual(blank_header(['x',' ', 'z']), (['x',' ','z'], \
            'Found an empty or all whitespace header. Please check the input '+\
            'file for missing headers or headers consisting only of '+\
            'forbidden characters.'))

    def test_bad_char_in_header(self):
        """bad_char_in_header should complain if bad char in header"""
        self.assertEqual(bad_char_in_header(['x','y','z']), (['x','y','z'], ''))
        header, msg = bad_char_in_header(['x','\\','z'])
        self.assertNotEqual(msg, '')

    def test_pad_rows(self):
        """pad_rows should produce correct table"""
        good_table = [['a','b'],['c','d']]
        bad_too_long = [['a','b'],['c','d','e']]
        bad_too_short = [['a','b'],['c']]
        self.assertEqual(pad_rows(good_table), good_table)
        self.assertEqual(pad_rows(bad_too_long),good_table)
        self.assertEqual(pad_rows(bad_too_short), [['a','b'],['c','']])


    def test_barcode_missing(self):
        """barcode_missing should complain if barcode missing"""
        fields = ['x','BarcodeSequence', 'y']
        self.assertEqual(barcode_missing(fields), (fields, ''))
        fields = ['x','y']
        self.assertEqual(barcode_missing(fields), (fields, 
        "Second field should be barcode field: "+
        "expected BarcodeSequence but got y.  Correct header errors before attempting to address warnings."))
        fields = ['x']
        self.assertEqual(barcode_missing(fields), (fields,
        "Second field should be barcode field but got < 2 fields.  Correct header errors before attempting to address warnings."))

    def test_description_missing(self):
        """description_missing should complain if description missing"""
        fields = ['x','Description']
        self.assertEqual(description_missing(fields), (fields, ''))
        fields = ['x','y']
        self.assertEqual(description_missing(fields),\
         (['x', 'y', 'Description'], \
         'Last field should be description field: expected Description but got y.  Correct header errors before attempting to fix warnings.'))

    def test_wrap_arrays(self):
        """wrap_arrays should return correct headers and dict"""
        good_data = array([
            ['SampleID','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','No','x'],
            ])
        sample_descriptions = array(['Description', 'xtest', 'ytest'])
        header, sample_desc, data_as_dict = wrap_arrays(
            sample_descriptions, good_data)
        self.assertEqual(header, ['bc', 'ph', 'ctl', 'x'])
        self.assertEqual(sample_desc, {'x':'xtest', 'y':'ytest'})
        self.assertEqual(data_as_dict, 
            {'x':{'bc':'x','ph':'3','ctl':'Yes','x':'x'},
             'y':{'bc':'y','ph':'4','ctl':'No','x':'x'},
             })

    def test_check_vals_by_type(self):
        """check_vals_by_type should return indices that can't convert"""
        good_vals = [1,'2',3.0]
        bad_vals = [1, 'x', [3,4], 4, None]
        self.assertEqual(check_vals_by_type(good_vals, int), [])
        self.assertEqual(check_vals_by_type(bad_vals, int), [1,2,4])

    def test_check_vals_by_contains(self):
        """check_vals_by_contains should return indices not in supplied object"""
        good_vals = ['a','b','a']
        bad_vals = [None, 'x', 'a', 3]
        self.assertEqual(check_vals_by_contains(good_vals, 'ab'), [])
        self.assertEqual(check_vals_by_contains(bad_vals, 'ab'), [0,1,3])

    def test_check_field_types(self):
        """check_field_types should return string of errors for invalid fields"""
        field_types = {'bc':'uid','ph':float,'sample':'uid','ctl':['Yes','No']}
        good_data = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','No','x'],
            ])
        self.assertEqual(check_field_types((good_data, field_types)), 
                ((good_data, field_types),''))
        bad_ctl = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','Nx','x'],
            ])
        self.assertEqual(check_field_types((bad_ctl, field_types)),
            ((bad_ctl, field_types),
            "Could not find Nx (sample id y, col ctl) in allowed vals "+
            "['Yes', 'No']"))

        bad_ph = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','x','Yes','x'],
            ['y','y','4','No','x'],
            ])
        self.assertEqual(check_field_types((bad_ph, field_types)),
            ((bad_ph, field_types),
            "Could not convert x (sample id x, col ph) to right type"))

        bad_bc = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','x','4','No','x'],
            ])
        self.assertEqual(check_field_types((bad_bc, field_types)),
            ((bad_bc, field_types),
                "DupChecker 'bc' found the following possible duplicates. If these metadata should have the same name, please correct.:\nGroup\tOriginal names\nx\tx, x\n"))

    def test_check_lengths(self):
        """check_lengths should return string of errors for invalid fields"""
        field_types = {'bc':'uid','ph':float,'sample':'uid','ctl':['Yes','No']}
        good_data = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','No','x'],
            ])
        self.assertEqual(check_same_length((good_data, field_types), 'bc'), 
                ((good_data, field_types),''))
        bad_bc = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','yx','4','No','x'],
            ])
        self.assertEqual(check_same_length((bad_bc, field_types),'bc'),
            ((bad_bc, field_types),
            "In field bc, item yx (sample id y) differs in length from "+
            "first item x (2 and 1).Location (row, column):\t1,1"))

    def test_check_bad_chars(self):
        """check_bad_chars should return string of errors for invalid fields"""
        field_types = {'bc':'uid','ph':float,'sample':'uid','ctl':['Yes','No']}
        good_data = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','No','x'],
            ])
        self.assertEqual(check_bad_chars((good_data, field_types)), 
                ((good_data, field_types),''))
        bad_vals = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x!','3','Yes','x'],
            ['y','y','>','No','x'],
            ])
        
        self.assertEqual(check_bad_chars((bad_vals, field_types)),((array([['sample', 'bc', 'ph', 'ctl', 'x'],['x', 'x_', '3', 'Yes', 'x'], ['y', 'y', '_', 'No', 'x']], dtype='|S6'), {'sample': 'uid', 'ctl': ['Yes', 'No'], 'ph': float, 'bc': 'uid'}), 'Removed bad chars from cell x! (now x_) in sample id x, col bc. Location (row, column):\t0,1\nRemoved bad chars from cell > (now _) in sample id y, col ph. Location (row, column):\t1,2'))

    def test_check_mixed_caps(self):
        """check_mixed_caps should return string of errors for invalid fields"""
        field_types = {'bc':'uid','ph':float,'sample':'uid','ctl':['Yes','No']}
        good_data = array([
            ['sample','bc','ph','ctl','x'],
            ['x','x','3','Yes','x'],
            ['y','y','4','No','x'],
            ])
        self.assertEqual(check_mixed_caps((good_data, field_types)), 
                ((good_data, field_types),''))
        bad_vals = array([
            ['sample','bc','ph','ctl','x'],
            ['x','Y','3','Yes','x'],
            ['y','y','>','  yes_ ','x'],
            ])
        self.assertEqual(check_mixed_caps((bad_vals, field_types)),
            ((bad_vals, field_types),
            "DupChecker 'Caps and Whitespace' found the following possible duplicates. If these metadata should have the same name, please correct. Found in field bc:\nGroup\tOriginal names\ny\tY, y\n\nDupChecker 'Caps and Whitespace' found the following possible duplicates. If these metadata should have the same name, please correct. Found in field ctl:\nGroup\tOriginal names\nyes\tYes,   yes_ \n"
        ))

    def test_check_missing_descriptions(self):
        """check_missing_descriptions should add run description and warn"""
        cmd = check_missing_descriptions
        good_sd = ['x','y','z']
        bad_sd = ['x', ' ', '']
        sample_ids = ['1','2','3']
        rd = 'test'
        self.assertEqual(cmd((good_sd, sample_ids, rd)), 
            ((good_sd, sample_ids, rd), ''))
        self.assertEqual(cmd((bad_sd, sample_ids, rd)),\
         ((['x', 'missing_description', 'missing_description'], \
         ['1', '2', '3'], 'test'), \
         "These sample ids lack descriptions (replaced with 'missing_description'): 2,3"))


    def test_check_duplicate_descriptions(self):
        """check_duplicate_descriptions should warn about duplicates"""
        cdd = check_duplicate_descriptions
        good_sd = ['Description','x','y','z']
        dup_sd = ['Description','x', 'y', 'x']
        sample_ids = ['#SampleID','1','2','3']
        raw_data_good = [['#SampleID','Description'],['1','x'],['2','y'],['3','z']]
        raw_data_dup = [['#SampleID','Description'],['1','x'],['2','y'],['3','x']]
        rd = 'test'
        self.assertEqual(cdd((good_sd, sample_ids, rd),raw_data=raw_data_good), 
            ((good_sd, sample_ids, rd), ''))
        self.assertEqual(cdd((dup_sd, sample_ids, rd),raw_data=raw_data_dup), \
         ((dup_sd, sample_ids, rd), 'These sample ids have duplicate descriptions:\n1,3: x\nRow, column for all duplicate descriptions:\nLocation (row, column):\t0,1\nLocation (row, column):\t2,1'))
        # ((bad_sd, sample_ids, rd) removed

    def test_check_description_chars(self):
        """check_description_chars should warn about bad chars"""
        cdc = check_description_chars
        good_sd = ['Description','x' , 'y' , 'z']
        bad_sd = ['Description','<' , 'y' , 'x>']
        sample_ids = ['#SampleID','1','2','3']
        raw_data_good = [['#SampleID','Description'],['1','x'],['2','y'],\
         ['3','z']]
        raw_data_bad = [['#SampleID','Description'],['1','<'],['2','y'],\
         ['3','x>']]
        rd = 'test'
        self.assertEqual(cdc((good_sd, sample_ids, rd), raw_data=raw_data_good), 
            ((good_sd, sample_ids, rd), ''))
        self.assertEqual(cdc((bad_sd, sample_ids, rd), raw_data=raw_data_bad),
            ((['Description', '_', 'y', 'x_'], ['#SampleID', '1', '2', '3'], 'test'), "These sample ids have bad characters in their descriptions:\n1: changed '<' to '_'\n3: changed 'x>' to 'x_'\nRow, column for all descriptions with bad characters:\nLocation (row, column):\t0,1\nLocation (row, column):\t2,1"))

    def test_process_id_map(self):
        """process_id_map should return correct results on small test map"""
        s = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tX\tDescription
#fake data
x\tAA\tACGT\t3\tsample_x
y\t"AC"\tACGT\t4\t"sample_y"
z\tGG\tACGT\t5\tsample_z"""
        f = StringIO(s)
        f.name='test.xls'
        headers, id_map, description_map, run_description, errors, warnings = \
            process_id_map(f)

        self.assertEqual(headers, ['BarcodeSequence', 'LinkerPrimerSequence', \
         'X'])
        self.assertEqual(id_map, {'y': {'X': '4', 'LinkerPrimerSequence': \
         'ACGT', 'BarcodeSequence': 'AC'}, 'x': {'X': '3', \
         'LinkerPrimerSequence': 'ACGT', 'BarcodeSequence': 'AA'}, 'z': \
        {'X': '5', 'LinkerPrimerSequence': 'ACGT', 'BarcodeSequence': 'GG'}})
        self.assertEqual(description_map, {
            'x':'sample_x',
            'y':'sample_y',
            'z':'sample_z',
        })
        self.assertEqual(run_description, ['fake data'])
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])
        
    def test_get_primers_barcodes(self):
        """ get_primers_barcodes should properly return primers, barcodes """
        
        expected_barcodes = ['ATCG','ATTA']
        expected_primers = ['CCTT','CCTAT']
        
        mapping_data = \
         [['#SampleID','BarcodeSequence','LinkerPrimerSequence','DOB'],
         ['PC.354','ATCG','CCTT','20061218'],
         ['PC.356','ATTA','CCTAT','20061216']]
         
        primers, barcodes = get_primers_barcodes(mapping_data, \
         is_barcoded=True, disable_primer_check=False)
         
        self.assertEqual(barcodes, expected_barcodes)
        self.assertEqual(primers, expected_primers)
        
        # Should return empty list if barcodes disabled
        expected_barcodes = []
        
        primers, barcodes = get_primers_barcodes(mapping_data, \
         is_barcoded=False, disable_primer_check=False)
         
        self.assertEqual(barcodes, expected_barcodes)
        self.assertEqual(primers, expected_primers)
        
        # Should return empty list if primers disabled
        expected_barcodes = ['ATCG','ATTA']
        expected_primers = []
        
        primers, barcodes = get_primers_barcodes(mapping_data, \
         is_barcoded=True, disable_primer_check=True)
         
        self.assertEqual(barcodes, expected_barcodes)
        self.assertEqual(primers, expected_primers)
        
        # Both primers and barcodes should be empty lists if both disabled
        expected_barcodes = []
        expected_primers = []
        
        primers, barcodes = get_primers_barcodes(mapping_data, \
         is_barcoded=False, disable_primer_check=True)
         
        self.assertEqual(barcodes, expected_barcodes)
        self.assertEqual(primers, expected_primers)
        
        
        


class CharFilterTests(TestCase):
    """Tests of CharFilter class"""
    def test_init(self):
        """CharFilter object should init and set properties right"""
        c = CharFilter('abc')
        self.assertEqual(c.Chars, 'abc')
        self.assertEqual(c.Name, None)
        self.assertEqual(c.Invert, False)
        self.assertEqual(c.stripF, strip)

    def test_call(self):
        """CharFilter object should omit or keep desired chars"""
        c = CharFilter('abc')
        d = CharFilter('abc',invert_charset=True)
        e = CharFilter('abc', default_char='#')
        f = CharFilter('abc', invert_charset=True, default_char='#')
        s = ('aceace')
        self.assertEqual(c(s), 'acac')
        self.assertEqual(d(s), 'ee')
        self.assertEqual(e(s), 'ac#ac#')
        self.assertEqual(f(s), '##e##e')
        #should automatically strip whitespace
        s2 = '   aceacef   '
        self.assertEqual(c(s2), 'acac')
        self.assertEqual(d(s2), 'eef')

    def test_badChars(self):
        """CharFilter badChars should report bad chars"""
        c = CharFilter('abc')
        d = CharFilter('abc',invert_charset=True)
        #should automatically strip whitespace
        s2 = '   aceacef   '
        self.assertEqual(c.badChars(s2), set('ef'))
        self.assertEqual(d.badChars(s2), set('ac'))

    def test_errMsg(self):
        """CharFilter errMsg should produce useful error message"""
        c = CharFilter('abc')
        d = CharFilter('abc', name='XYZ', invert_charset=True)
        #should automatically strip whitespace
        s2 = '  aceacef  '
        self.assertEqual(c.errMsg(s2), 
            "Filter 'None' found bad chars 'e,f' in input '  aceacef  '")
        self.assertEqual(d.errMsg(s2),
            "Filter 'XYZ' found bad chars 'a,c' in input '  aceacef  '")

    def test_resultAndError(self):
        """CharFilter resultAndError should return both result and error"""
        c = CharFilter('abc')
        s2 = '  aceacef  '
        self.assertEqual(c.resultAndError(s2), (c(s2), c.errMsg(s2)))

class DupCheckerTests(TestCase):
    """DupChecker object should check duplicates correctly."""
    def test_init(self):
        """DupChecker init should occur without errors"""
        d = DupChecker(lwu)
        self.assertEqual(d.Name, None)
        self.assertEqual(d.CanonicalF, lwu)

    def test_call(self):
        """DupChecker call should report dups"""
        d = DupChecker(lwu)
        self.assertEqual(d(['a','b','c']), {})
        self.assertEqual(d(['a','a','b']), {'a':['a','a']})
        self.assertEqual(d(['a','A','b']), {'a':['a','A']})
        self.assertEqual(d(['AA','A_a','b']), {'aa':['AA','A_a']})
        d2 = DupChecker()
        self.assertEqual(d2(['a','b']),{})
        self.assertEqual(d2(['a','A']),{})
        self.assertEqual(d2(['a','a']),{'a':['a','a']})

        d3 = DupChecker(lwu, allow_exact_dup=True)
        self.assertEqual(d3(['a','b','c']), {})
        self.assertEqual(d3(['a','a','b']), {})
        self.assertEqual(d3(['a','A','b']), {'a':['a','A']})
        self.assertEqual(d3(['AA','A_a','b']), {'aa':['AA','A_a']})


    def test_errMsg(self):
        """DupChecker errMsg should return useful error message"""
        d = DupChecker(lwu, 'Test')
        self.assertEqual(d.errMsg(['a','b','c']), '')
        self.assertEqual(d.errMsg(['a','A','b']), 
            "DupChecker 'Test' found the following possible duplicates. If these metadata should have the same name, please correct.:\nGroup\tOriginal names\na\ta, A\n")

    def test_dupIndices(self):
        """DupChecker dupIndices should report dup indices"""
        d = DupChecker(lwu)
        self.assertEqual(d.dupIndices(['a','b','c']), {})
        self.assertEqual(d.dupIndices(['a','a','b']), {'a':[0,1]})
        self.assertEqual(d.dupIndices(['a','A','b']), {'a':[0,1]})
        self.assertEqual(d.dupIndices(['AA','A_a','b']), {'aa':[0,1]})
        d2 = DupChecker()
        self.assertEqual(d2.dupIndices(['a','b']),{})
        self.assertEqual(d2.dupIndices(['a','A']),{})
        self.assertEqual(d2.dupIndices(['a','a']),{'a':[0,1]})

        d3 = DupChecker(lwu, allow_exact_dup=True)
        self.assertEqual(d3.dupIndices(['a','b','c']), {})
        self.assertEqual(d3.dupIndices(['a','a','b']), {})
        self.assertEqual(d3.dupIndices(['a','A','b']), {'a':[0,1]})
        self.assertEqual(d3.dupIndices(['AA','A_a','b']), {'aa':[0,1]})


class SameCheckerTests(TestCase):
    """SameChecker object should enforce similarity constraints correctly."""
    def test_init(self):
        """SameChecker init should occur without errors"""
        s = SameChecker(lwu)
        self.assertEqual(s.Name, None)
        self.assertEqual(s.CanonicalF, lwu)

    def test_call(self):
        """SameChecker call should report mismatches"""
        s = SameChecker(lwu)
        self.assertEqual(s(['a','A','A ']), [])
        self.assertEqual(s(['A','a','B']), [[2,'B','b','A','a']])
        s2 = SameChecker()
        self.assertEqual(s2(['a','a']), [])
        self.assertEqual(s2(['a','A']), [[1,'A','A','a','a']])

    def test_errMsg(self):
        """DupChecker errMsg should return useful error message"""
        s = SameChecker(lwu, 'Test')
        self.assertEqual(s.errMsg(['a','A','A ']), '')
        self.assertEqual(s.errMsg(['a','A','B']), 
            "SameChecker 'Test' found the following values different from the first:\nIndex\tVal\tf(Val)\tFirst\tf(First)\n2\tB\tb\ta\ta\n")

    def test_check_primers_barcodes(self):
        """ Should give warnings for invalid or missing primers/barcodes """
        
        
        problems = defaultdict(list)
        barcodes_good = ['CACGC','CCACG','GGTTA']
        # The linker sequence, usually two base pairs is considered to be
        # part of the primer.
        primers_good = ['GGATTCG','AATRCGG','CANGCRT']
        # Should append nothing to problems with valid barcodes, primers.
        self.assertEqual(check_primers_barcodes(primers_good, barcodes_good, \
         problems), defaultdict(list))
        barcodes_bad = ['CAC1C','','GGAAT']
        primers_bad = ['1GGATTCG','ATCCATCG','']
        # Should create warnings about invalid characters and missing barcode
        # and primer
        self.assertEqual(check_primers_barcodes(primers_bad, barcodes_bad, \
         problems),  defaultdict(list, {'warning': ['The primer 1GGATTCG has invalid characters.  Location (row, column):\t0,2', 'Missing primer.  Location (row, column):\t2,2', 'The barcode CAC1C has invalid characters.  Location (row, column):\t0,1', 'Missing barcode. Location (row, column):\t1,1']}))
         
        # Should not raise errors if barcodes missing and disabled
        problems = defaultdict(list)
        barcodes_absent = ['','','']
        primers_good = ['GGATTCG','AATRCGG','CANGCRT']
        self.assertEqual(check_primers_barcodes(primers_good, barcodes_absent, \
         problems, is_barcoded=False), defaultdict(list))
         
        # Should not raise errors if primers missing and disabled
        problems = defaultdict(list)
        barcodes_absent = ['CACGC','CCACG','GGTTA']
        primers_good = ['','','']
        self.assertEqual(check_primers_barcodes(primers_good, barcodes_absent, \
         problems, is_barcoded=True, disable_primer_check=True),\
         defaultdict(list))
         
    def test_check_missing_sampleIDs(self):
        """ Should give warnings if missing sample IDs from given list """
        
        problems = defaultdict(list)
        sample_IDs_good = ['#SampleID','Sample_1','Sample_2']
        # Should not create any warnings
        self.assertEqual(check_missing_sampleIDs(sample_IDs_good, problems), \
         defaultdict(list))
        # Should raise errors for empty/whitespace sample ID list items
        sample_IDs_bad = ['#SampleID','','Sample_2']
        self.assertEqual(check_missing_sampleIDs(sample_IDs_bad, problems), \
         defaultdict(list, {'warning': ['Missing Sample ID.  Location (row, column):\t0,0']}))
        sample_IDs_bad = ['#SampleID','Sample_1','\t']
        problems = defaultdict(list)
        self.assertEqual(check_missing_sampleIDs(sample_IDs_bad, problems), \
         defaultdict(list, {'warning': ['Missing Sample ID.  Location (row, column):\t1,0']}))
         
    def test_check_dup_var_barcodes_primers(self):
        """ Should give warnings/location of duplicate barcode+primer seqs """
       
        test_barcodes = ['AATCGA' , 'TACCGT' , 'ATCCGTAT']
        test_primers = ['CCGGAT' , 'CCGGAT' , 'CCGGAT']
        problems = defaultdict(list)
        
        # Since there are no duplicates when the barcodes and primers are
        # concatenated, there should be nothing added to problems
              
        self.assertEqual(check_dup_var_barcodes_primers(test_primers,\
         test_barcodes, problems), defaultdict(list))
         
        test_barcodes = ['AATCGA' , 'AATCGAC' , 'ATCCGTAT']
        test_primers = ['CCGGAT' , 'CGGAT' , 'CCGGAT']
        
        # The first and second barcode+primers should be equal, should
        # append a warning, give duplicate sequences, and location of problem
        
        self.assertEqual(check_dup_var_barcodes_primers(test_primers,\
         test_barcodes, problems), defaultdict(list,  {'warning': ['The barcode + primer sequence "AATCGACCGGAT" has duplicate results.  Location (row, column):\t0,1', 'The barcode + primer sequence "AATCGACCGGAT" has duplicate results.  Location (row, column):\t1,1']}))
         
         
            


if __name__ =='__main__':
    main()
