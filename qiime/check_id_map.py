#!/usr/bin/env python
#check_id_map.py: checks that sample mapping file is ok
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

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
from collections import defaultdict
from string import strip, letters, digits, strip, translate
from os.path import basename
from cogent.util.transform import (keep_chars, exclude_chars, trans_except,
    trans_all)
from numpy import array
from pipe454.parse import parse_map
from optparse import OptionParser

DESC_KEY = "Description"
SAMPLE_ID_KEY = "SampleID"
BARCODE_KEY = "BarcodeSequence"
NEGATIVE_CONTROL_KEY = "NegativeControl"

STANDARD_FIELD_TYPES = {SAMPLE_ID_KEY:'uid', BARCODE_KEY:'uid', \
      NEGATIVE_CONTROL_KEY:['Yes','No']}

ALLOWED_CHARS = '+#.' + digits + letters
ALLOWED_SUBCAT_CHARS = ALLOWED_CHARS + '_'
EXTRA_DESC_CHARS = ", -():/'"
ALLOWED_DESC_CHARS = ALLOWED_SUBCAT_CHARS + EXTRA_DESC_CHARS

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

class DupChecker(object):
    """Checks set of objects for duplicates in canonical representation.
    
    Returns dict of {collision:[originals]}, empty if input OK.
    """

    def __init__(self, canonical_f=None, name=None, allow_exact_dup=False):
        """Returns new DupChecker, will use canonical_f for test."""
        self.CanonicalF = canonical_f
        self.Name = name
        self.AllowExactDup = allow_exact_dup

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
            res = "DupChecker '%s' found the following possible duplicates" \
                % str(self.Name)
            if field_name:
                res += " in field %s" % field_name
            res += ':\n'
            res += "Group\tOriginal names\n"
            for k, v in sorted(bad):
                res += '%s\t%s\n' % (k, ', '.join(map(str,v)))
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

def run_checks(data, checks, problems):
    """Runs checks on data, reports issues in problems and returns clean data.

    checks should be list of (check_f, type) tuples.

    WARNING: checks are performed only once, so if a later check _introduces_
    a problem you would have seen in an earlier check you won't see it.
    """
    for check, type_ in checks:
        data, problem = check(data)
        if problem:
            problems[type_].append(problem)
    return data

def filename_has_space(fname):
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
    def inner_f(data, field_name=field_name):
        """returns f(data) -> (clean_data, error_msg)"""
        return data, dup_checker.errMsg(data, field_name=field_name)
    return inner_f

raw_dup_checker = adapt_dupchecker(lambda x:x, 'Duplicate checker')
space_dup_checker = adapt_dupchecker(lwu, 
    'Duplicate checker including whitespace and capitalization')
space_dup_checker_header = adapt_dupchecker(lwu, 
    'Duplicate checker including whitespace and capitalization', 'Header')

#checks for valid headers
def sampleid_missing(fields, field_name=SAMPLE_ID_KEY):
    """Returns error message if sample id field doesn't start with #"""
    try:
        if fields[0] == '#' + field_name:
            return fields, ''
    except (TypeError, IndexError):
        pass
    return fields, 'SampleID field must start with #: '+\
    'please ensure that this is not a binary (e.g. Excel) file.' +\
    ' and that the %s field is first.' % SAMPLE_ID_KEY

def blank_header(fields):
    """Returns error message if any header is empty"""
    stripped_fields = map(strip, fields)
    if '' in stripped_fields:
        return fields, 'Found an empty or all whitespace header. ' +\
            'Please check the input file for missing headers or ' +\
            'headers consisting only of forbidden characters.'
    else:
        return fields, ''

def bad_char_in_header(fields):
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

def barcode_missing(fields):
    """Returns error message if second field is not barcode field"""
    if len(fields) < 2:
        return fields, \
            'Second field should be barcode field but got < 2 fields.'
    elif fields[1] == BARCODE_KEY:
        return fields, ''
    else:
        return fields, "Second field should be barcode field:" + \
            " expected %s but got %s." % (BARCODE_KEY, fields[1])

def description_missing(fields):
    """Returns error message if last field is not description field"""
    if fields[-1] == DESC_KEY:
        return fields, ''
    else:
        return fields + [DESC_KEY], "Last field should be description field:"+\
            " expected %s but got %s." % (DESC_KEY, fields[-1])

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
    
def format_map_file(headers, id_map, description_map=None, \
    run_description=None):
    """Generates string for formatted map file.
    
    Input:
        headers: list of strings corresponding to col headers
        id_map: dict of {id:{header:val}}
        description_map: dict of {id:description}
        run_description: either string, or list of strings
    """
    result = []
    if DESC_KEY in headers:
        headers.remove(DESC_KEY)
    if SAMPLE_ID_KEY in headers:
        headers.remove(SAMPLE_ID_KEY)
    header_line = '\t'.join([SAMPLE_ID_KEY] + headers + [DESC_KEY])
    if not header_line.startswith('#'):
        header_line = '#' + header_line
    result.append(header_line)
    if run_description:
        if not isinstance(run_description, str):
            run_description = '\n#'.join(run_description)
        if not run_description.startswith('#'):
            run_description = '#'+run_description
        result.append(run_description)
    for id_, fields in sorted(id_map.items()):
        curr_line = [id_]
        curr_line.extend([fields.get(h,'') for h in headers])
        curr_line.append(description_map.get(id_,''))
        result.append('\t'.join(map(str, curr_line)))
    return '\n'.join(result)

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
        raise ValueError, "Didn't get same # of sample ids and descriptions!"
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

def check_field_types((data, field_types)):
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

def check_same_length((data, field_types),col_name=BARCODE_KEY):
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
            (col_name,item,data[index+1,0],orig_item,len_item,len_orig_item))
    return (data, field_types), '\n'.join(errors)

descr_filter = CharFilter(ALLOWED_DESC_CHARS, "Description Filter", 
    default_char='#')
header_filter = CharFilter(ALLOWED_CHARS, "Header Filter", 
    default_char='#')
subcat_filter = CharFilter(ALLOWED_SUBCAT_CHARS, "Subcat Filter", 
    default_char='#')

def check_bad_chars((data, field_types), filter_f=subcat_filter):
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
                (val, new_val, row[i], headers[j]))
    return (data, field_types), '\n'.join(problems)

def check_mixed_caps((data, field_types), dup_f=lwu, 
    dup_name='Caps and Whitespace', allow_exact_dup=True):
    """Checks all fields for mixed caps, warning."""
    dup_checker = DupChecker(dup_f, dup_name, allow_exact_dup=allow_exact_dup)
    problems = []
    headers, body = data[0], data[1:]
    for i, col in enumerate(body.T):
        errmsg = dup_checker.errMsg(col, field_name=headers[i])
        if errmsg:
            problems.append(errmsg)
    return (data, field_types), '\n'.join(problems)

def check_missing_descriptions((sample_descriptions, sample_ids, 
    run_description)):
    """Returns warnings for sample ids with missing descriptions."""
    missing = []
    for i, desc in enumerate(sample_descriptions):
        if not desc.strip():
            missing.append(i)
    if missing:
        err =  \
        "These sample ids lack descriptions (used run description): %s" % \
        ','.join(sorted([sample_ids[i] for i in missing]))
        for i in missing:
            sample_descriptions[i] = run_description
    else:
        err = ''
    return (sample_descriptions, sample_ids, run_description), err

def check_duplicate_descriptions((sample_descriptions, sample_ids,
    run_description)):
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
    return (sample_descriptions, sample_ids, run_description), \
        '\n'.join(problems)

def check_description_chars((sample_descriptions, sample_ids, run_description),
    filter_f=descr_filter):
    """Returns warnings for descriptions with bad chars, replacing them."""
    new_descr = map(filter_f, sample_descriptions)
    errors = []
    for id_, old, new in zip(sample_ids, sample_descriptions, new_descr):
        if old != new:
            errors.append("%s: changed '%s' to '%s'" % (id_, old, new))
    if errors:
        errors.insert(0, 
            "These sample ids have bad characters in their descriptions:")
    return (new_descr, sample_ids, run_description), '\n'.join(errors)


STANDARD_FILENAME_CHECKS = [(filename_has_space, 'error')]
STANDARD_RUN_DESCRIPTION_CHECKS = [(run_description_missing, 'warning'), 
    (descr_filter.resultAndError, 'warning')]
STANDARD_SAMPLE_DESCRIPTION_CHECKS = [
    (check_missing_descriptions, 'warning'),
    (check_duplicate_descriptions, 'warning'),
    (check_description_chars, 'warning'),
    ]
STANDARD_COL_HEADER_CHECKS = [(sampleid_missing, 'error'),
    (bad_char_in_header, 'warning'),
    (space_dup_checker_header, 'error'), 
    (blank_header, 'error'),
    (description_missing, 'warning'),
    ]
BARCODE_COL_HEADER_CHECKS = [(barcode_missing, 'error')]
STANDARD_COL_CHECKS = [
        (check_field_types, 'error'),
        (check_bad_chars, 'warning'),
        (check_mixed_caps, 'warning'),
        ]
BARCODE_COL_CHECKS = [(check_same_length, 'warning')]

def parse_id_map(infile, is_barcoded=True,
    filename_checks=STANDARD_FILENAME_CHECKS, 
    run_description_checks=STANDARD_RUN_DESCRIPTION_CHECKS,
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
        data, run_description = parse_map(infile, return_header=True)
        col_headers = data[0]
    except (TypeError, ValueError), e:
        problems['error'].append(
            "Couldn't read map file '%s': failed with error message %s" 
            % (infile_name, e))
        #Note: this error is fatal so we have to bail out if we get it
        return col_headers, id_map, description_map, run_description, errors, \
            warnings

    #note: parse_map always returns run_description as list of lines
    run_description = ' '.join(map(strip, run_description))

    #check run description
    run_description = run_checks(run_description, run_description_checks, \
        problems) 

    #add barcode checks if needed
    if is_barcoded:
        col_header_checks.extend(BARCODE_COL_HEADER_CHECKS)
        col_checks.extend(BARCODE_COL_CHECKS)

    #check col headers
    col_headers = run_checks(col_headers, col_header_checks, problems)

    #check col values
    data = array(pad_rows(data))

    #Note: last field should be description if default checks are applied, so
    #need to remove. However, we are not making any assumptions here, so if the
    #last col isn't a description we will add one.
    #here, data and description both include the column headers.
    if col_headers[-1] == DESC_KEY:
        data, sample_descriptions = data[:,:-1], data[:,-1]
    else:
        data, sample_descriptions = data, array([DESC_KEY] + ['']*(len(data)-1))

    data, field_types = run_checks((data, field_types), col_checks, problems)
    sample_ids = data[:,0]

    #check sample descriptions
    sample_descriptions, sample_ids, run_description = \
    run_checks((sample_descriptions, sample_ids, 
        run_description), sample_description_checks, problems)

    #return formatted output
    headers, description_map, id_map = wrap_arrays(sample_descriptions, data)
    errors = problems['error']
    warnings = problems['warning']
    return headers, id_map, description_map, run_description, errors, warnings

def make_cmd_parser():
    """Returns command-line options"""
    parser = OptionParser()
    parser.add_option('-m', '--map', dest='map_fname',
        help='name of mapping file')
    parser.add_option('-b', '--is-barcoded', dest='has_barcodes',
        type=int, default=1,
        help='1 if barcoded (default), 0 otherwise: ' + \
            'must include BarcodeSequence field if barcoded.')
    options, args = parser.parse_args()
    return options

if __name__ == "__main__":
    from sys import argv, exit, stderr
    options = make_cmd_parser()
    infile_name = options.map_fname
    has_barcodes = options.has_barcodes
    stderr.write("\nReading id map:" + infile_name + '\n')
    headers, id_map, description_map, run_description, errors, warnings = \
        parse_id_map(open(infile_name, 'U'), has_barcodes)

    if errors:
        stderr.write('\nERRORS-------------------\n')
        for ix, f in enumerate(errors):
            stderr.write("%d: %s\n" % (ix, f))
    if warnings:
        stderr.write('\nWARNINGS -------------------\n')
        for ix, f in enumerate(warnings):
            stderr.write("%d: %s\n" % (ix, f))

    print format_map_file(headers, id_map, description_map, run_description)
