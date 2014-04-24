#!/usr/bin/env python

from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh", "Jose Carlos Clemente Litran",
               "Jai Ram Rideout"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from numpy import array, concatenate, asarray, transpose, log, invert, asarray,\
    float32, float64, minimum, inf
from optparse import OptionParser
from qiime.util import MissingFileError
import os
from qiime.filter import filter_otus_from_otu_table
from biom.parse import parse_biom_table


def make_html_doc(js_filename):
    """Create the basic framework for the OTU table heatmap"""
    html_script = \
        r'''
    <html>
    <head>
    \t<script type="text/javascript" src="js/overlib.js"></script>
        <script type="text/javascript" src="%s"></script>
    \t<script type="text/javascript" src="js/otu_count_display.js"></script>
    \t<script type="text/javascript" src="./js/jquery.js"></script>
    \t<script type="text/javascript" src="./js/jquery.tablednd_0_5.js"></script>
        <script type="text/javascript">


        $(document).ready(function() {

        \t$('#otu_table_body').tableDnD({
        \t\tonDragStart: function(table, new_row) {
        \t\t\tif (row==new_row.parentNode.rowIndex && is_selected==1){
        \t\t\t\tchange_sel_row=1;
        \t\t\t}else{
        \t\t\t\told_row=new_row.parentNode.rowIndex;
        \t\t\t\tchange_sel_row=0;
        \t\t\t}
        \t\t},
        \t\tonDrop: function(table, new_row) {
        \t\t\tif (change_sel_row==1){
        \t\t\t\trow=new_row.rowIndex;
        \t\t\t}else if(old_row<row && new_row.rowIndex>row){
        \t\t\t\trow=row-1;
        \t\t\t}else if(old_row>row && new_row.rowIndex<row){
        \t\t\t\trow=row+1;
        \t\t\t}
        \t\t},
            \tdragHandle: "dragHandle"
        \t});
            var otu_cutoff=document.getElementById("otu_count_cutoff");
            otu_cutoff.value=otu_num_cutoff;
        });
        </script>
    \t<style type="text/css">
    \t    th.rotate{
    \t\t\twhite-space : nowrap;
    \t\t\t-webkit-transform: rotate(-90deg) translate(20px, 0px);
    \t\t\t-moz-transform: rotate(-90deg) translate(20px, 0px);
    \t\t\tfont-family:arial;
    \t\t\tfont-size:9px;
    \t\t}
    \t\tth.lineage{
        \t    white-space : nowrap;
        \t    text-align:left;
        \t    font-family:arial;
        \t    font-size:10px;
        \t    font-weight: bolder;
        \t}
        \ttd.dragHandle{
            \twhite-space : nowrap;
            \ttext-align:left;
            \tfont-family:arial;
            \tfont-size:10px;
            \tfont-weight: bolder;
        \t}
        \ttd{
            \twhite-space : nowrap;
            \tfont-family:arial;
            \tfont-size:10px;
            \ttext-align:center;
            \tfont-weight: bolder;
        \t}
        \ttable{
            \tborder-spacing: 0;
            \ttext-align:center;
        \t}
        \tp{
            \t\ttext-align:left;
            \t\tfont-weight: normal;
        \t}
    \t</style>
    </head>
    <body>
    \t<p>
    \t\tFilter by Counts per OTU: <input type="text" id="otu_count_cutoff" value="">
    \t\t<input type="button" onclick="javascript:create_OTU_intervals();" value="Sample ID">
    \t\t<input type="button" onclick="javascript:write_taxon_heatmap();" value="Taxonomy">
    \t</p>
    \t<br><br><br><br><br><br>
    \t<table id='otu_table_html'>
    \t\t<thead id='otu_table_head'>
    \t\t</thead>
    \t\t<tbody id='otu_table_body'>
    \t\t<tr><td class="dragHandle"></td>
    \t\t</tr>
    \t\t<tr><td class="dragHandle"></td>
    \t\t</tr>
    \t\t</tbody>
    \t</table>

    </body>
    </html>''' % (js_filename)
    return html_script


def create_javascript_array(otu_table, use_floats=False):
    """Convert the OTU table counts into a javascript array"""
    # Build up list of strings and concatenate at end, as this string can be
    # very large and have many concatenations.
    js_array = ['var OTU_table=new Array();\n'
                'var i=0;\n'
                'for (i==0;i<%i;i++) {\n'
                'OTU_table[i]=new Array();}\n' %
                (len(otu_table.sample_ids) + 2)]

    # 0 ['#OTU ID', 'OTU2', 'OTU3']
    #1 ['Sample1', 1, 2]
    #2 ['Sample2', 5, 4]
    #3 ['Consensus Lineage', 'Archaea', 'Bacteria']

    # OTU ids first
    js_array.append("OTU_table[0][0]='#OTU ID';\n")
    for (idx, otu_id) in enumerate(otu_table.observation_ids):
        js_array.append("OTU_table[0][%i]='%s';\n" % (idx + 1, otu_id))

    # Sample ids and values in the table
    i = 1
    for (sam_val, sam_id, meta) in otu_table.iter_samples():
        js_array.append("OTU_table[%i][0]='%s';\n" % (i, sam_id))
        for (idx, v) in enumerate(sam_val):
            if use_floats:
                js_array.append("OTU_table[%i][%i]=%.4f;\n" %
                                (i, idx + 1, float(v)))
            else:
                # don't quite understand why int(float(v)), rather than int(v)
                js_array.append("OTU_table[%i][%i]=%d;\n" %
                                (i, idx + 1, int(float(v))))
        i += 1

    # Consensus lineages for each OTU
    last_idx = len(otu_table.sample_ids) + 1
    js_array.append("OTU_table[%i][0]='Consensus Lineage';\n" % last_idx)
    i = 1
    for (otu_val, otu_id, meta) in otu_table.iter_observations():
        js_array.append("OTU_table[%i][%i]='%s';\n" %
                        (last_idx, i, ";".join(meta['taxonomy']).strip('"')))
        i += 1

    return ''.join(js_array)


def filter_by_otu_hits(num_otu_hits, otu_table):
    """Filter the OTU table by the number of otus per sample"""
    # Filter out rows with sum > num_otu_hits
    new_otu_table = filter_otus_from_otu_table(
        otu_table, otu_table.observation_ids,
        num_otu_hits, inf, 0, inf)

    return new_otu_table


def get_log_transform(otu_table, eps=None):
    """ This function and the one in make_otu_heatmap.py are essentially the same except
    the non-negative transform at the end of this function. Dan Knights suggests this might
    be due to this script not being able to handle negative values, hence the transform.
    """
    # explicit conversion to float: transform
    def f(s_v, s_id, s_md):
        return float64(s_v)
    float_otu_table = otu_table.transform_samples(f)

    if eps is None:
        # get the minimum among nonzero entries and divide by two
        eps = inf
        for (obs, sam) in float_otu_table.nonzero():
            eps = minimum(eps, float_otu_table.get_value_by_ids(obs, sam))
        if eps == inf:
            raise ValueError('All values in the OTU table are zero!')

    # set zero entries to eps/2 using a transform

    def g2(x):
        return [i if i != 0 else eps / 2 for i in x]

    # do we have map in OTU object?
    g = lambda x: x if (x != 0) else eps / 2

    def g_m(s_v, s_id, s_md):
        return asarray(map(g, s_v))

    eps_otu_table = float_otu_table.transform_samples(g_m)

    # take log of all values with transform
    def h(s_v, s_id, s_md):
        return log(s_v)
    log_otu_table = eps_otu_table.transform_samples(h)

    # one more transform
    min_val = inf
    for val in log_otu_table.iterSampleData():
        min_val = minimum(min_val, val.min())

    def i(s_v, s_id, s_md):
        return s_v - min_val

    res_otu_table = log_otu_table.transform_samples(i)

    return res_otu_table


def get_otu_counts(fpath):
    """Reads the OTU table file into memory"""

    try:
        otu_table = parse_biom_table(open(fpath, 'U'))
    except (TypeError, IOError):
        raise MissingFileError('OTU table file required for this analysis')

    if (otu_table.observation_metadata is None or
            otu_table.observation_metadata[0]['taxonomy'] is None):
        raise ValueError(
            '\n\nThe lineages are missing from the OTU table. Make sure you included the lineages for the OTUs in your OTU table. \n')

    return otu_table


def generate_heatmap_plots(
        num_otu_hits, otu_table, otu_sort, sample_sort, dir_path,
        js_dir_path, filename, fractional_values=False):
    """Generate HTML heatmap and javascript array for OTU counts"""

    # Filter by number of OTU hits
    # rows come transposed in the original code
    filtered_otu_table = filter_by_otu_hits(num_otu_hits, otu_table)

    if otu_sort:
        # Since the BIOM object comes back with fewer Observation_ids, we need to
        # remove those from the original sort_order
        actual_observations = filtered_otu_table.observation_ids
        new_otu_sort_order = []
        for i in otu_sort:
            if i in actual_observations:
                new_otu_sort_order.append(i)

        filtered_otu_table = filtered_otu_table.sortObservationOrder(
            new_otu_sort_order)

    # This sorts the samples by the order supplied
    if sample_sort:
        # Since the BIOM object may come back with fewer Sampleids, we need to
        # remove those from the original sample_sort
        actual_samples = filtered_otu_table.sample_ids
        new_sample_sort_order = []
        for i in sample_sort:
            if i in actual_samples:
                new_sample_sort_order.append(i)

        filtered_otu_table = filtered_otu_table.sortSampleOrder(
            new_sample_sort_order)

    # Convert OTU counts into a javascript array
    js_array = create_javascript_array(filtered_otu_table, fractional_values)

    # Write otu filter number
    js_otu_cutoff = 'var otu_num_cutoff=%d;\n' % num_otu_hits

    # Write js array to file
    js_filename = os.path.join(js_dir_path, filename) + '.js'
    jsfile = open(js_filename, 'w')
    jsfile.write(js_otu_cutoff)
    jsfile.write(js_array)
    jsfile.close()

    # Write html file
    html_filename = os.path.join(dir_path, filename) + '.html'
    js_file_location = 'js/' + filename + '.js'
    table_html = make_html_doc(js_file_location)
    ofile = open(html_filename, 'w')
    ofile.write(table_html)
    ofile.close()
