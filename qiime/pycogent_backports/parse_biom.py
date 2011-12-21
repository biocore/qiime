#!/usr/bin/env python
# File created on 19 dec 2011
from __future__ import division
import qiime.pycogent_backports.rich_otu_table as rt
import json
import numpy
from qiime.parse import process_otu_table_sample_ids, parse_otu_table
from string import strip

MATRIX_ELEMENT_TYPE = {'int':int,'float':float,'str':str,
                       u'int':int,u'float':float,u'str':str}
def parse_biom_table(json_fh):
    """parses a biom format otu table into a rich otu table object

    input is an open filehandle or compatable object (e.g. list of lines)

    sparse/dense will be determined by "matrix_type" in biom file, and 
    either a SparseOTUTable or DenseOTUTable object will be returned
    note that sparse here refers to the compressed format of [row,col,count]
    dense refers to the full / standard matrix representations
    """

    json_table = json.load(json_fh)
    if json_table['type'].lower() != 'otu table':
        raise ValueError('type not OTU table')

    sample_ids = [col['id'] for col in json_table['columns']]
    # null metadata -> None object in metadata list 
    sample_metadata = [col['metadata'] for col in json_table['columns']]
    ObservationIds = [row['id'] for row in json_table['rows']]
    obs_metadata = [row['metadata'] for row in json_table['rows']]

    dtype = MATRIX_ELEMENT_TYPE[json_table['matrix_element_type']]
    if json_table['matrix_type'].lower() == 'sparse': # lowercase is expected
        import pysparse
        dims = json_table['shape']
        data = pysparse.spmatrix.ll_mat(*dims)
        for entry in json_table['data']:
            data[entry[0],entry[1]] = entry[2]
        table_obj = rt.SparseOTUTable(Data=data, 
            SampleIds=sample_ids, ObservationIds=ObservationIds,
            SampleMetadata=sample_metadata, 
            ObservationMetadata=obs_metadata)

    elif json_table['matrix_type'].lower() == 'dense':
        data = numpy.asarray(json_table['data'],
            dtype=dtype)
        table_obj = rt.DenseOTUTable(Data=data, 
            SampleIds=sample_ids, ObservationIds=ObservationIds,
            SampleMetadata=sample_metadata, 
            ObservationMetadata=obs_metadata)
    else:
        raise ValueError( 'invalid matrix_type in biom table')

    return table_obj

def parse_otu_table_to_rich_otu_table(lines,count_map_f=int,dense=False):
    """parses an otu table (tab delimited) (sample ID x OTU ID map)

    Returns a rich otu table object (a subclass of Table), 
    sparse by default (or see parameter 'dense')
    """
    if dense:
        sample_ids, otu_ids, otu_table, metadata = parse_otu_table(lines,count_map_f=count_map_f)
        if len(metadata) > 0:
            metadata = [{'taxonomy':elem} for elem in metadata]
        else:
            metadata = None
        table_obj = rt.DenseOTUTable(Data=otu_table,
        SampleIds=sample_ids, ObservationIds=otu_ids,
        SampleMetadata=None, ObservationMetadata=metadata)
        return table_obj

    otu_ids = []
    metadata = []
    sample_ids = []
    # iterate over lines in the OTU table -- keep track of line number 
    # to support legacy (Qiime 1.2.0 and earlier) OTU tables
    two_d_dict = {} # {(row,col):value}, row is otu/observaiton
    otu_idx = 0 # keep track of observation/otu lines, (skip comments and headers, etc.)
    for i, line in enumerate(lines):
        line = line.strip()
        if line:
            if i == 1 and line.startswith('#OTU ID') and not sample_ids:
                # we've got a legacy OTU table
                try:
                    sample_ids, has_metadata = process_otu_table_sample_ids(
                     line.strip().split('\t')[1:])
                except ValueError:
                    raise ValueError, \
                     "Error parsing sample IDs in OTU table. Appears to be a"+\
                     " legacy OTU table. Sample ID line:\n %s" % line
            elif not line.startswith('#'):
                if not sample_ids:
                    # current line is the first non-space, non-comment line 
                    # in OTU table, so contains the sample IDs
                    try:
                        sample_ids, has_metadata = process_otu_table_sample_ids(
                         line.strip().split('\t')[1:])
                    except ValueError:
                        raise ValueError,\
                         "Error parsing sample IDs in OTU table."+\
                         " Sample ID line:\n %s" % line
                else:
                    # current line is OTU line in OTU table
                    fields = line.split('\t')

                    # if there is OTU metadata the last column gets appended
                    # to the metadata list
                    # otherwise all columns are appended to otu_table
                    if has_metadata:
                        abund_fields = fields[1:-1]
                        metadata.append({'taxonomy':map(strip, fields[-1].split(';'))})
                    else:
                        abund_fields = fields[1:]
                        metadata.append(None)

                    # added in a try/except to handle OTU tables containing
                    # floating numbers
                    try:
                        abunds = map(count_map_f,abund_fields)
                    except ValueError:
                        abunds = map(float,abund_fields)

                    if all([abund==0 for abund in abunds]):
                        continue #don't increment otu_idx counter
                    else:
                        # grab the OTU ID
                        otu_id = fields[0].strip()
                        otu_ids.append(otu_id)

                    for j, abund in enumerate(abunds):
                        if abund != 0:
                            two_d_dict[(otu_idx,j)] = abund
                    otu_idx += 1 # this is needed for indexing into two_d_dict
                    # this sets dimensions of matrix, so it must be accurate

    data = rt.to_ll_mat(two_d_dict)
    table_obj = rt.SparseOTUTable(Data=data, 
        SampleIds=sample_ids, ObservationIds=otu_ids,
        SampleMetadata=None, ObservationMetadata=metadata)
    return(table_obj)

