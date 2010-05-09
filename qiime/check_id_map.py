#!/usr/bin/env python

"""Parse mapping file, checking for a number of undesirable characteristics.

Specifically, we check that:
    - the filename does not contain spaces (warn + rewrite if it does)
    - there is an overall run description at the start (warn + add if missing)
    - there is a SampleID field and, if barcoded, a BarcodeSequence field
      (warn + correct reasonable variants of both)
    - there are not duplicate header fields (error)
    - there are not duplicate near-unique but not exactly unique values 
      within each column (warning)

Overall strategy:
    - maintain list of errors and warnings (initially empty).
    - for each check in the following classes:
      - check_name
      - check_run_description
      - check_col_headers
      - check_cols
    - run the check f(data) -> msg
    - collect the messages, assigning to error or warning
    - return errors and warnings

Returns both errors and warnings as lists of formatted strings: should not
raise exceptions itself under normal circumstances (e.g. if file is
mis-formatted).

SampleID column is required - must contain unique values. 
BarcodeSequence field required (if is barcoded) and must contain unique values

It is somewhat inefficient to read in the whole table, but on the other hand
if reading in the mapping file is the bottleneck the downstream analysis is
likely to prove somewhat challenging as well...
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","William Walters"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "William Walters"
__email__ = "rob@spot.colorado.edu","william.a.walters@colorado.edu"
__status__ = "Development"

from collections import defaultdict
from string import strip, letters, digits, strip, translate
from os.path import basename, isdir
from os import makedirs
from cogent.util.transform import (keep_chars, exclude_chars, trans_except,
    trans_all)
from numpy import array
from qiime.parse import parse_mapping_file
from qiime.format import format_map_file
from optparse import OptionParser
from cogent.seqsim.sequence_generators import IUPAC_DNA


DESC_KEY = "Description"
SAMPLE_ID_KEY = "SampleID"
BARCODE_KEY = "BarcodeSequence"
LINKER_PRIMER_KEY = "LinkerPrimerSequence"
NEGATIVE_CONTROL_KEY = "NegativeControl"



STANDARD_FIELD_TYPES = {SAMPLE_ID_KEY:'uid', BARCODE_KEY:'uid', \
      LINKER_PRIMER_KEY:'uid',NEGATIVE_CONTROL_KEY:['Yes','No']}

# Header category most restrictive
ALLOWED_CHARS_HEADER = '_' + digits + letters
# Metadata should be fairly restricted as well
EXTRA_DESC_CHARS = "+-%."
#EXTRA_COMMENT_CHARS =  ''
ALLOWED_DESC_CHARS = ALLOWED_CHARS_HEADER + EXTRA_DESC_CHARS
# Comment line(s) following header should not be used for parsing data
# Should allow more characters.
ALLOWED_SUBCAT_CHARS = ALLOWED_CHARS_HEADER + EXTRA_DESC_CHARS

def find_diff_length(items):
    """Return items that differ in length from the first, with indices.
    
    (also returns length of current item and of first item, and the items).
    Useful for identifying element in a list that is wrong length, e.g. in
    a table where all the rows are supposed to have the same number of fields.

    Returns [] if there was no problem (note: if None is the element with a
    diff length, you'll get a TypeError when you try to calculate its length
    rather than returning it as a value, so the apparent conflict here
    doesn't actually exist in practice).
    """
    result = []
    orig_len = None
    for i, item in enumerate(items):
        if orig_len is None:
            orig_len = len(item)
            orig_item = item
        else:
            if len(item) != orig_len:
                result.append([i, item, len(item), orig_item, orig_len])
    return result

class CharFilter(object):
    """Checks string for allowed chars.
    
    Returns cleaned result. Other methods return error message etc."""

    def __init__(self, chars, name=None, invert_charset=False,strip_f=strip,
        default_char=None):

        """Returns new CharFilter object."""
        self.Chars = chars
        self.Name = name
        self.Invert = invert_charset
        if invert_charset:
            if default_char:
                trans_table = trans_all(chars, default_char)
                self.Filter = lambda s: s.translate(trans_table)
            else:
                self.Filter = exclude_chars(chars)
        else:
            if default_char:
                trans_table = trans_except(chars, default_char)
                self.Filter = lambda s: s.translate(trans_table)
            else:
                self.Filter = keep_chars(chars)
        self.stripF = strip_f

    def __call__(self, input):
        """Returns cleaned version of input."""
        return self.Filter(self.stripF(input))

    def badChars(self, input):
        """Returns set of bad chars in input (empty if OK)."""
        orig = self.stripF(input)
        new = self(orig)
        return set(orig) - set(new)

    def errMsg(self, input):
        """Generates error message from input ('' if OK)."""
        bad = self.badChars(input)
        if bad:
            return "Filter '%s' found bad chars '%s' in input '%s'" % \
                    (str(self.Name), ','.join(sorted(bad)), input)
        else:
            return ''

    def resultAndError(self, input):
        """Convenience wrapper: returns result and error message."""
        return self(input), self.errMsg(input)
        
header_filter = CharFilter(ALLOWED_CHARS_HEADER, "Header Filter", 
    default_char='#')
descr_filter = CharFilter(ALLOWED_DESC_CHARS, "Description Filter", 
    default_char='_')
subcat_filter = CharFilter(ALLOWED_SUBCAT_CHARS, "Subcat Filter", 
    default_char='_')

class DupChecker(object):
    """Checks set of objects for duplicates in canonical representation.
    
    Returns dict of {collision:[originals]}, empty if input OK.
    """

    def __init__(self, canonical_f=None, name=None, allow_exact_dup=False,\
     raw_data=None):
        """Returns new DupChecker, will use canonical_f for test."""
        self.CanonicalF = canonical_f
        self.Name = name
        self.AllowExactDup = allow_exact_dup
        self.Raw_Data=raw_data
        # Need to allow exact duplicates for linker-primer sequences
        if self.Name=="LinkerPrimerSequence":
            self.AllowExactDup = True

    def __call__(self, input):
        """Checks input for non-uniqueness."""
        result = defaultdict(list)
        if self.CanonicalF:
            for i in input:
                result[self.CanonicalF(i)].append(i)
        else:
            for i in input:
                result[i].append(i)

        bad_pairs = []
        for k, v in result.items():
            if self.AllowExactDup:
                v_to_test = set(v)
            else:
                v_to_test = v
            if len(v_to_test) > 1:
                bad_pairs.append((k, v))
        
        return dict(bad_pairs)

    def errMsg(self, input, field_name=None):
        """Generates error message from input ('' if OK)."""
        bad = sorted(self(input).items())
        if bad:
            # Record column index for row, column locations
            if self.Raw_Data:
                column_index = self.Raw_Data[0].index(field_name)
            elif self.Name=="BarcodeSequence":
                column_index = 1
            else:
                column_index = None

            res = "DupChecker '%s' found the following possible duplicates" \
                % str(self.Name) +". If these metadata should have the same "+\
                "name, please correct."
            if field_name:
                res += " Found in field %s" % field_name
            res += ':\n'
            res += "Group\tOriginal names\n"
            for k, v in sorted(bad):
                res += '%s\t%s\n' % (k, ', '.join(map(str,v)))
            # Get the row indices of bad data
            row_indices = []
            index_counter = 0
            bad_data = []
            for bad_d in bad:
                bad_data += bad_d[1]
            for row_datum in input:
                if row_datum in bad_data:
                    row_indices.append(index_counter)
                index_counter += 1
            # Only append data if raw_data passed
            if column_index:
                # Append row, column indices messages
                res += "Row, column for all possible duplicate descriptions:\n"
                for row in row_indices:
                    res += "Location (row, column):\t%d,%d\n" % (row, column_index)
                

            return res
        else:
            return ''

    def dupIndices(self, input):
        """Checks input for non-uniqueness, returning indices of collisions."""
        result = defaultdict(list)
        if self.CanonicalF:
            for i, item in enumerate(input):
                result[self.CanonicalF(item)].append(i)
        else:
            for i, item in enumerate(input):
                result[item].append(i)

        bad_pairs = []
        for k, v in result.items():
            if self.AllowExactDup:
                v_to_test = set([input[i] for i in v])
            else:
                v_to_test = v
            if len(v_to_test) > 1:
                bad_pairs.append((k, v))
        return dict(bad_pairs)


def lowercase_whitespace_underscore(s):
    """Returns s lowercase and stripped of whitespace and underscores"""
    return s.lower().strip().replace('_','').replace(' ','',).replace('\t','')

lwu = lowercase_whitespace_underscore

class SameChecker(object):
    """Checks that set of objects is same in canonical representation.
    
    Returns list of index, input, f(input), first, f(first); empty if input OK.
    """

    def __init__(self, canonical_f=None, name=None):
        """Returns new SameChecker, will use canonical_f for test."""
        self.CanonicalF = canonical_f
        self.Name = name

    def __call__(self, input):
        """Checks input for uniqueness."""
        overall_result = []
        orig_val = None
        first_val = None
        for i, val in enumerate(input):
            if self.CanonicalF:
                res = self.CanonicalF(val)
            else:
                res = val
            if not i:
                first_res = res
                first_val = val
            elif res != first_res:
                overall_result.append([i, val, res, first_val, first_res])
        return overall_result

    def errMsg(self, input, field_name=None):
        """Generates error message from input ('' if OK)."""
        bad = self(input)
        if bad:
            res = "SameChecker '%s' " % str(self.Name) +\
                "found the following values different from the first" \
                
            if field_name:
                res += "in field %s" % field_name
            res += ':\n'
            res += "Index\tVal\tf(Val)\tFirst\tf(First)\n"
            for i in bad:
                res += '\t'.join(map(str,i))+'\n'
            return res
        else:
            return ''

def run_checks(data, checks, problems, all_mapping_data=None):
    """Runs checks on data, reports issues in problems and returns clean data.

    checks should be list of (check_f, type) tuples.

    WARNING: checks are performed only once, so if a later check _introduces_
    a problem you would have seen in an earlier check you won't see it.
    """
    
    for check, type_ in checks:
        data, problem = check(data, raw_data=all_mapping_data)
        if problem:
            problems[type_].append(problem)
    return data

def filename_has_space(fname, raw_data=None):
    """Returns message if filename contains space character"""
    if ' ' in fname:
        return fname.replace(' ','_'), \
            "Filename may not contain spaces. "+\
            "Please re-upload without spaces, e.g. %s -> %s." \
            % (fname, fname.replace(' ', '_'))
    return fname, ''

RUN_DESCRIPTION_DEFAULT = "No run description supplied."

def run_description_missing(desc, default_value=RUN_DESCRIPTION_DEFAULT):
    """Returns message if run description missing"""
    if not desc:
        return default_value, \
            "Run description was not supplied, using default value."
    return desc, ''

def adapt_dupchecker(f, name, field_name=None):
    """Returns function that adapts DupChecker to API for checkers."""
    dup_checker = DupChecker(f, name)
    def inner_f(data, field_name=field_name, raw_data=None):
        """returns f(data) -> (clean_data, error_msg)"""
        return data, dup_checker.errMsg(data, field_name=field_name)
    return inner_f

raw_dup_checker = adapt_dupchecker(lambda x:x, 'Duplicate checker')
space_dup_checker = adapt_dupchecker(lwu, 
    'Duplicate checker including whitespace and capitalization')
space_dup_checker_header = adapt_dupchecker(lwu, 
    'Duplicate checker including whitespace and capitalization', 'Header')

#checks for valid headers
def sampleid_missing(fields, field_name=SAMPLE_ID_KEY, raw_data=None):
    """Returns error message if sample id field doesn't start with #"""


    try:
        if fields[0].strip() == "#" + field_name:
            return fields, ''
    except (TypeError, IndexError):
        pass
    return fields, 'SampleID field must start with #SampleID: '+\
    'please ensure that this is not a binary (e.g. Excel) file.' +\
    ' and that the %s field is first.  Found %s' %\
     (SAMPLE_ID_KEY,fields[0].strip())

def blank_header(fields, raw_data=None):
    """Returns error message if any header is empty"""
    stripped_fields = map(strip, fields)
    if '' in stripped_fields:
        return fields, 'Found an empty or all whitespace header. ' +\
            'Please check the input file for missing headers or ' +\
            'headers consisting only of forbidden characters.'
    else:
        return fields, ''

def bad_char_in_header(fields, raw_data=None):
    """Returns error message if bad char in header"""
    bad_chars = []
    filtered_fields=[]
    for f in fields:
        filtered, bad = header_filter.resultAndError(f)
        if bad:
            bad_chars.append([f, bad])
        filtered_fields.append(filtered)
    if bad_chars:
        return filtered_fields, "Found bad characters in these headers:\n" +\
            '\n'.join(["%s\t%s" % (bad, f) for bad, f in bad_chars])
    return fields, ''

def barcode_missing(fields, raw_data=None):
    """Returns error message if second field is not barcode field"""
    if len(fields) < 2:
        return fields, \
            'Second field should be barcode field but got < 2 fields.  '+\
            'Correct header errors before attempting to address warnings.'
    elif fields[1] == BARCODE_KEY:
        return fields, ''
    else:
        return fields, "Second field should be barcode field:" + \
            " expected %s but got %s." % (BARCODE_KEY, fields[1]) +\
            "  Correct header errors before attempting to address warnings."
            
def linker_primer_missing(fields, raw_data=None):
    """Returns error message if third field is not linker_primer field"""
    if len(fields) < 3:
        return fields, \
            'Third field should be linker_primer field but got < 3 fields.  '+\
            'Correct header errors before attempting to address warnings.'
    elif fields[2] == LINKER_PRIMER_KEY or fields[1]== LINKER_PRIMER_KEY:
        return fields, ''
    else:
        return fields, "Third field should be linker_primer field:" + \
            " expected %s but got %s." % (LINKER_PRIMER_KEY, fields[2]) +\
            " Correct header errors before attempting to address warnings."

def description_missing(fields, raw_data=None):
    """Returns error message if last field is not description field"""
    if fields[-1] == DESC_KEY:
        return fields, ''
    else:
        return fields + [DESC_KEY], "Last field should be description field:"+\
            " expected %s but got %s." % (DESC_KEY, fields[-1]) +\
            "  Correct header errors before attempting to fix warnings."
            

def pad_rows(table):
    """Ensures that table has missing fields padded with empty string."""
    num_cols = len(table[0])
    result = []
    for row in table:
        if len(row) == num_cols:
            result.append(row)
        elif len(row) > num_cols:
            result.append(row[:num_cols])
        else: #must be too short
            result.append(row + ['']*(num_cols-len(row)))
    return result
    

def wrap_arrays(sample_descriptions, data):
    """Wraps sample descriptions and data into appropriate dicts.

    Assumptions:
    - first field in sample_descriptions is 'Description'
    - sample_descriptions exist for each row
    - first row in data is the header
    - first col in data is the sample id
    - descriptions already removed from the data matrix by this point
    
    Results:
    - header (= first row as list of strings)
    - sample_desc (= dict of {sample_id: description} for each sample)
    - data (= dict of {sample_id:{field_name:val}} for each sample/field)
    """
    sample_descriptions = list(sample_descriptions)[1:]
    header = list(data[0, 1:])  #remove first field = SampleID
    body = data[1:]
    sample_ids = body[:,0]
    body = body[:,1:]
    if len(sample_ids) != len(sample_descriptions):
        raise ValueError, "Didn't get same number of sample ids and "+\
         "descriptions!"
    sample_desc = dict(zip(sample_ids, sample_descriptions))
    data_as_dict = {}
    for sample_id, fields in zip(sample_ids, body):
        data_as_dict[sample_id] = dict(zip(header, fields))
    return header, sample_desc, data_as_dict

def check_vals_by_type(vals, type_):
    """Checks each val can be converted to type, returns index if fails.
    """
    result = []
    for i, v in enumerate(vals):
        try:
            type_(v)
        except (TypeError, ValueError):
            result.append(i)
    return result

def check_vals_by_contains(vals, contains):
    """Checks each val is in contains, returns index if fails."""
    result = []
    for i, v in enumerate(vals):
        try:
            if v not in contains:
                result.append(i)
        except (KeyError, ValueError, TypeError):
            result.append(i)
    return result

def check_field_types((data, field_types), raw_data=None):
    """Checks that field types match data"""
    errors = []
    col_headers = list(data[0])
    body = data[1:]
    for col, type_ in field_types.items():
        bad_indices = []
        if col in col_headers:
            index = col_headers.index(col)
            vals = body[:,index]
            if isinstance(type_, type):
                bad_indices = check_vals_by_type(vals, type_)
                for i in bad_indices:
                    errors.append(
                    "Could not convert %s (sample id %s, col %s) to right type"
                        % (vals[i], body[i,0], col))
            elif type_ == 'uid':
                dup_checker = DupChecker(name=col)
                err_msg = dup_checker.errMsg(vals)
                if err_msg:
                    errors.append(err_msg)
            else:
                bad_indices = check_vals_by_contains(vals, type_)
                for i in bad_indices:
                    errors.append(
                "Could not find %s (sample id %s, col %s) in allowed vals %s"%
                    (vals[i], body[i,0], col, type_))
    return (data, field_types), '\n'.join(errors)

def check_same_length((data, field_types),col_name=BARCODE_KEY, raw_data=None):
    """Checks field lengths, reporting mismatch: assumes column present."""
    col_headers = list(data[0])
    errors = []
    
    
    if col_name not in col_headers:
        errors.append(("The required field %s is missing from the mapping file"
            % col_name))
    else:
        col_index = col_headers.index(col_name)
        result = find_diff_length(data[1:,col_index])
        for index, item, len_item, orig_item, len_orig_item in result:
            errors.append(("In field %s, item %s (sample id %s) "+
            "differs in length from first item %s (%s and %s).") %
            (col_name,item,data[index+1,0],orig_item,len_item,len_orig_item)+\
            "Location (row, column):\t%d,1" % index)
    return (data, field_types), '\n'.join(errors)



def check_bad_chars((data, field_types), filter_f=subcat_filter, raw_data=None):
    """Checks all fields for bad chars, removing and warning."""
    problems = []
    headers, body = data[0], data[1:]

    for i, row in enumerate(body):
        for j, val in enumerate(row):
            new_val, e = filter_f.resultAndError(val)
            if e:
                row[j] = new_val
                problems.append(
            "Removed bad chars from cell %s (now %s) in sample id %s, col %s." %
             (val, new_val, row[0], headers[j]) + " Location (row, column):"+\
             "\t%d,%d" % (i,j))
    return (data, field_types), '\n'.join(problems)

def check_mixed_caps((data, field_types), dup_f=lwu, 
    dup_name='Caps and Whitespace', allow_exact_dup=True, raw_data=None):
    """Checks all fields for mixed caps, warning."""
    
    dup_checker = DupChecker(dup_f, dup_name, allow_exact_dup,\
     raw_data)
    problems = []
    headers, body = data[0], data[1:]
    for i, col in enumerate(body.T):
        errmsg = dup_checker.errMsg(col, field_name=headers[i])
        if errmsg:
            problems.append(errmsg)
    return (data, field_types), '\n'.join(problems)

def check_missing_descriptions((sample_descriptions, sample_ids, 
    run_description), column_index=None, raw_data=None):
    """Returns warnings for sample ids with missing descriptions."""
    missing = []
    for i, desc in enumerate(sample_descriptions):
        if not desc.strip():
            missing.append(i)
    if missing:
        err =  "These sample ids lack descriptions (replaced with "+\
         "'missing_description'): %s" %\
         ','.join(sorted([sample_ids[i] for i in missing]))
        for i in missing:
            sample_descriptions[i] = 'missing_description'
    else:
        err = ''
    return (sample_descriptions, sample_ids, run_description), err

def check_missing_sampleIDs(sample_ids, problems):
    """Returns warnings for missing sample IDs"""

    # sample IDs will always be in the first column (0)
    column = 0
    
    row = 0
    
    for sample_ID in sample_ids:
        # skip header
        if sample_ID == "#SampleID":
            continue
        if not sample_ID.strip():
            problems['warning'].append('Missing Sample ID.  ' +\
             'Location (row, column):\t%d,%d' % (row, column))
        row += 1
            

    return problems

def check_duplicate_descriptions((sample_descriptions, sample_ids,
    run_description), raw_data=None):
    """Returns warnings for duplicate descriptions"""
    
    d = DupChecker()
    dup_indices = d.dupIndices(sample_descriptions)
    problems = []
    for k, v in dup_indices.items():
        problems.append("%s: %s" % (','.join([sample_ids[i] for i in v]), k))
    #if there were any problems, insert a useful header
    if problems:
        problems.insert(0, 
            "These sample ids have duplicate descriptions:")
    #Insert a list of (rows,columns) locations for later error correction
    #Get column index for Description column
    if raw_data:
        column_index = len(raw_data[0]) - 1
    #Get list of the row indices of duplications
    row_indices = []
    for dups in dup_indices.values():
        row_indices += dups

    #If there are problems, append list of rows,columns in the standard way
    if problems:
        problems.append("Row, column for all duplicate descriptions:")
        for row in row_indices:
            # Row requires a correction to start at zero
            problems.append("Location (row, column):\t%d,%d" % \
             (row-1, column_index))


    return (sample_descriptions, sample_ids, run_description), \
        '\n'.join(problems)
        
def check_duplicate_sample_ids((sample_descriptions, sample_ids,
    run_description), raw_data=None):
    """Returns warnings for duplicate sample_ids"""
    

    d = DupChecker()
    dup_indices = d.dupIndices(sample_ids)
    problems = []
    for k, v in dup_indices.items():
        problems.append("%s: %s" % (','.join([sample_ids[i] for i in v]), k))
    #if there were any problems, insert a useful header
    if problems:
        problems.insert(0, 
            "These sample ids have duplicate ids:")
    # Column_index for sample_ids is always 0
    column_index = 0
    #Get list of the row indices of duplications
    row_indices = []
    for dups in dup_indices.values():
        row_indices += dups

    #If there are problems, append list of rows,columns in the standard way
    if problems:
        problems.append("Row, column for all duplicate ids:")
        for row in row_indices:
            # Row requires a correction to start at zero
            problems.append("Location (row, column):\t%d,%d" % \
             (row-1, column_index))


    return (sample_descriptions, sample_ids, run_description), \
        '\n'.join(problems)
        
def check_primers_barcodes(primers, barcodes, problems, is_barcoded=True,
 disable_primer_check=False):
    """Returns warnings for primers/barcodes that have invalid characters 
    
    The check_primers_barcodes function only tests for valid IUPAC DNA
    characters and for the presence of a primer or barcode.  No testing
    for valid Golay/Hamming barcodes or duplicates are performed in this
    function."""

    
    for row in range(len(primers)):

        for base in primers[row]:
            try:
                IUPAC_DNA[base]
            except KeyError:
                # The primers are always located in the third column
                problems['warning'].append('The primer %s ' % primers[row] +\
                'has invalid characters.  Location (row, column):\t' +\
                '%d,2' % row)
        if len(primers[row])==0 and not disable_primer_check:
            problems['warning'].append('Missing primer.  ' +\
             'Location (row, column):\t%d,2' % row)
    
    if is_barcoded:
        for row in range(len(barcodes)):
            for base in barcodes[row]:
                try:
                    IUPAC_DNA[base]
                except KeyError:
                    # The barcodes are always located in the second column
                    problems['warning'].append('The barcode %s ' % barcodes[row] +\
                     'has invalid characters.  Location (row, column):\t' +\
                     '%d,1' % row)
            if len(barcodes[row])==0:
                problems['warning'].append('Missing barcode. '+\
                'Location (row, column):\t%d,1' % row)

    return problems

def check_description_chars((sample_descriptions, sample_ids, run_description),\
 filter_f=descr_filter, raw_data=None):
    """Returns warnings for descriptions with bad chars, replacing them."""
    new_descr = map(filter_f, sample_descriptions)
    errors = []

    # Store descriptions that have changed to find indices
    changed_desc_indices = []
    # Use a counter to record row indices of bad chars
    row_index = 0
    for id_, old, new in zip(sample_ids, sample_descriptions, new_descr):
        if old != new:
            errors.append("%s: changed '%s' to '%s'" % (id_, old, new))
            changed_desc_indices.append(row_index)
        row_index += 1
    if errors:
        errors.insert(0, \
         "These sample ids have bad characters in their descriptions:")
            
    #Insert a list of (rows,columns) locations for later error correction
    #Get column index for Description column
    column_index = len(raw_data[0]) - 1


    #If there are errors, append list of rows,columns in the standard way
    if errors:
        errors.append("Row, column for all descriptions with bad characters:")
        for row in changed_desc_indices:
            # Row requires a correction to start at zero
            errors.append("Location (row, column):\t%d,%d" % \
             (row-1, column_index))

    
    
    return (new_descr, sample_ids, run_description), '\n'.join(errors)

def get_sample_description_column(data):
    """ Returns column of sample_description, used for indexing errors """
    
    # Test to ensure that column referenced is the sample_description
    '''if not(data[0][-1]=='Description'):
        raise ValueError,('Incorrect mapping data passed.  Final column '+\
        'should be "Description"')
    else:
        # Return len of column corrected by -1 for proper indexing
        return (len(data[0])-1) '''
    
    return (len(data[0])-1)

STANDARD_FILENAME_CHECKS = [(filename_has_space, 'error')]
# Removed run description checks, no longer using to fill in missing desc.
''' STANDARD_RUN_DESCRIPTION_CHECKS = [(run_description_missing, 'warning'), 
    (descr_filter.resultAndError, 'warning')] '''
STANDARD_SAMPLE_DESCRIPTION_CHECKS = [
    (check_missing_descriptions, 'warning'),
    (check_duplicate_descriptions, 'warning'),
    (check_duplicate_sample_ids, 'warning'),
    (check_description_chars, 'warning'),
    ]
STANDARD_COL_HEADER_CHECKS = [(sampleid_missing, 'error'),
    (bad_char_in_header, 'error'),
    (space_dup_checker_header, 'error'), 
    (blank_header, 'error'),
    (description_missing, 'error'),
    ]
BARCODE_COL_HEADER_CHECKS = [(barcode_missing, 'error')]
STANDARD_COL_CHECKS = [
        (check_field_types, 'error'),
        (check_bad_chars, 'warning'),
        (check_mixed_caps, 'warning'),
        ]
        

BARCODE_COL_CHECKS = [(check_same_length, 'warning')]
PRIMER_COL_CHECKS = [(linker_primer_missing, 'error')]

def get_primers_barcodes(data, is_barcoded, disable_primer_check):
    """ Returns list of primers, barcodes from mapping file """
    
    primers=[]
    barcodes=[]
    
    if not is_barcoded and disable_primer_check:
        return primers, barcodes
    # Convert DNA characters to uppercase, as IUPAC_DNA dictionary only
    # contains uppercase characters
    for sample in data:
        if sample[1]=="BarcodeSequence":
            continue
        if not is_barcoded:
            if sample[1]== "LinkerPrimerSequence":
                continue
        if is_barcoded and not disable_primer_check:
            barcodes.append(sample[1].upper())
            primers.append(sample[2].upper())
        elif not is_barcoded and not disable_primer_check:
            primers.append(sample[2].upper())
        elif is_barcoded and disable_primer_check:
            barcodes.append(sample[1].upper())

    
    return primers, barcodes
    
def check_dup_var_barcodes_primers(primers, barcodes, problems):
    """ Checks that no duplicate seqs occur when barcodes/primers appended """
    
    # Get list of concatenated barcodes + primers
    concat_barcodes_primers = []
    
    for primer, barcode in map(None, primers, barcodes):
        concat_barcodes_primers.append(barcode+primer)
    
    
    
    for seq_index in range(len(concat_barcodes_primers)):
        # Check for any duplicates, append warning to problems
        if concat_barcodes_primers.count(concat_barcodes_primers[seq_index])>1:
            problems['warning'].append('The barcode + primer sequence '+\
            '"%s"' % concat_barcodes_primers[seq_index] + ' has '+\
             'duplicate results.  Location (row, column):\t' +\
             '%d,1' % seq_index)
    
    return problems
    

def process_id_map(infile, disable_primer_check=False, is_barcoded=True, \
    char_replace="_", var_len_barcodes = False, 
    filename_checks=STANDARD_FILENAME_CHECKS, 
    #run_description_checks=STANDARD_RUN_DESCRIPTION_CHECKS,
    sample_description_checks=STANDARD_SAMPLE_DESCRIPTION_CHECKS,
    col_header_checks=STANDARD_COL_HEADER_CHECKS,
    col_checks=STANDARD_COL_CHECKS,
    field_types=STANDARD_FIELD_TYPES):
    """ Parse ID mapping file.
   
    Returns the following:

    headers: list of header values
    id_map: {sample_id:{header:value}}
    description_map: {sample_id:description}
    run_description: run description as a string
    errors: list of error messages generated
    warnings: list of warning messages generated
    """

    problems = defaultdict(list)
    errors, warnings = [], []
    col_headers, id_map, description_map, run_description = \
        None, None, None, None

    #check infile name
    infile_name = basename(infile.name)
    run_checks(infile_name, filename_checks, problems)
    
    #read data
    try:
        data, headers, run_description = parse_mapping_file(infile, \
        suppress_stripping=True)
        headers[0] = "#" + headers[0]
        col_headers = headers
        data.insert(0, headers)
        if run_description:
        	for n in range(len(run_description)):
        		run_description[n] = run_description[n].replace('\n','')
        # Need to replace newline characters in data
        for row in range(len(data)):
            for column in range(len(data[row])):
                data[row][column] = data[row][column].replace("\n","")
    except (TypeError, ValueError), e:
        problems['error'].append(
            "Couldn't read map file '%s': failed with error message %s" 
            % (infile_name, e))
        #Note: this error is fatal so we have to bail out if we get it
        return col_headers, id_map, description_map, run_description, errors, \
            warnings
            
    # Save raw data for referencing source 'cells' in log file.
    raw_data = data
    
    #check col values
    data = array(pad_rows(data))
    


    #add barcode checks if needed
    if is_barcoded:
        col_header_checks.extend(BARCODE_COL_HEADER_CHECKS)
        if not var_len_barcodes:
            col_checks.extend(BARCODE_COL_CHECKS)

        
    
    #check col headers
    col_headers = run_checks(col_headers, col_header_checks, problems, raw_data)

    #check col values
    data = array(pad_rows(data))
    

    sample_description_column = get_sample_description_column(data)
    
    #Note: last field should be description if default checks are applied, so
    #need to remove. However, we are not making any assumptions here, so if the
    #last col isn't a description we will add one.
    #here, data and description both include the column headers.
    if col_headers[-1] == DESC_KEY:
        data, sample_descriptions = data[:,:-1], data[:,-1]
    else:
        data, sample_descriptions = data, array([DESC_KEY] + ['']*(len(data)-1))

    data, field_types = run_checks((data, field_types), col_checks, problems, \
     raw_data)
    sample_ids = data[:,0]
    




    #check sample descriptions
    sample_descriptions, sample_ids, run_description = \
    run_checks((sample_descriptions, sample_ids,run_description, \
     ), sample_description_checks, problems, raw_data)
     
    #check primers,barcodes for valid IUPAC DNA characters
    primers, barcodes = get_primers_barcodes(data, is_barcoded, \
     disable_primer_check)
    problems = check_primers_barcodes(primers, barcodes, problems, \
     is_barcoded, disable_primer_check)
     
    if var_len_barcodes:
        problems = check_dup_var_barcodes_primers(primers, barcodes, problems)
        
    
    #check for missing sample_IDs
    problems = check_missing_sampleIDs(sample_ids, problems)

    #return formatted output
    headers, description_map, id_map = wrap_arrays(sample_descriptions, data)
    errors = problems['error']
    warnings = problems['warning']
    
    return headers, id_map, description_map, run_description, errors, warnings

def write_corrected_file(headers, id_map, description_map, run_description, \
output_filepath, chars_replaced=False):
    """ Writes corrected mapping file with illegal characters replaced """
    
    outf = open(output_filepath, 'w')
    
    if chars_replaced:
        outfile_data=format_map_file(headers, id_map, DESC_KEY, SAMPLE_ID_KEY,\
         description_map, run_description)
    else:
        outfile_data="# No invalid characters were found and replaced.\n"+\
        "# Note that non-IUPAC DNA characters found in primer or barcode\n"+\
        "# sequences are not replaced but will be listed in the .log file." 
     
    for data in outfile_data:
        outf.write(data)
        
    return
    
def test_for_replacement_chars(warnings):
    """ Checks for replacement character warnings, returns true if so """
    
    
    for warning in warnings:
        if warning.startswith("Removed ") or \
         warning.startswith("These sample ids have bad characters") or \
         warning.startswith("These sample ids lack descriptions (replaced "):
            return True
    
    return False
    
def write_logfile(errors, warnings, log_filepath, mapping_filepath):
    """ Writes errors/warnings or lack thereof to log filepath """
    
    try:
        log_f = open(log_filepath,"w")
    except IOError:
        # Attempt to create log file name based on mapping file name
        log_filepath += mapping_filepath.replace(".txt",".log")
        log_f = open(log_filepath,"w")
    
    if not (errors or warnings):
        log_f.write("No errors or warnings for mapping file %s" % \
         mapping_filepath)
    else:
        log_f.write("#Listed locations of errors/warnings in row/column "+\
        "format have an index that is \n#in reference to the beginning of "+\
        "the sample IDs and metadata.\n#Location (0,0) is the first SampleID "+\
        "in a given data set.")
         
    if errors:
        log_f.write('\nERRORS-------------------\n')
        for ix, f in enumerate(errors):
            log_f.write("%d: %s\n" % (ix, f))
    if warnings:
        log_f.write('\nWARNINGS -------------------\n')
        for ix, f in enumerate(warnings):
            log_f.write("%d: %s\n" % (ix, f))
            
    return log_filepath


def check_mapping_file(infile_name, output_dir, has_barcodes, char_replace, \
 verbose, var_len_barcodes, disable_primer_check):
    """ Central program function for checking mapping file """
	
    headers, id_map, description_map, run_description, errors, warnings = \
     process_id_map(open(infile_name, 'U'), has_barcodes, char_replace,\
     var_len_barcodes, disable_primer_check)

    chars_replaced = test_for_replacement_chars(warnings)

    
    mapping_root_name = infile_name.split("/")[-1].replace(".txt","")
    
    if not output_dir.endswith("/"):
        output_dir += "/"
    try:
        if not isdir(output_dir):
            makedirs(output_dir)
    except IOError:
        raise IOError,('Unable to create output directory %s ' % output_dir)
    try:
        corrected_output_filepath = output_dir + mapping_root_name + \
         '_corrected.txt'
        outf = open(corrected_output_filepath, 'w')
        outf.close()
    except IOError:
        raise IOError, ('Unable to create corrected output file %s ' %\
         corrected_output_filepath)
    try:
        log_filepath = output_dir + mapping_root_name + '.log'
        outf = open(log_filepath, 'w')
        outf.close()
    except IOError:
        raise IOError, ('Unable to create log file %s ' % log_filepath)
            
    write_corrected_file(headers, id_map, description_map, run_description,\
     corrected_output_filepath, chars_replaced)
         
    log_filepath = write_logfile(errors, warnings, log_filepath, infile_name)
    
    if verbose and (errors or warnings):
        print('Errors and/or warnings occurred, see log file %s' % log_filepath)
    if verbose and not(errors or warnings):
        print('No errors or warnings for mapfile %s' % infile_name)
	
	


    
    
    
    
        

