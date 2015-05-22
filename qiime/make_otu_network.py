#!/usr/bin/env python
# file make_otu_network.py

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Julia Goodrich"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"


"""
This script generates the otu networks and statistics

Author: Julia Goodrich (julia.goodrich@colorado.edu)
Status: Prototype

Requirements:
Python 2.5
"""

from collections import defaultdict
from numpy import nonzero, arange
from string import strip
import os
from time import strftime
from random import choice, randrange
from qiime.stats import G_2_by_2
from qiime.colors import iter_color_groups, Color, data_colors
from qiime.parse import parse_mapping_file
from biom import load_table


def get_sample_info(lines):
    """Collects information from mapping file for easy use in following steps

        Input: lines - the mapping file as a list of lines

        Output: cat_by_sample - a dictionary where the key is a id and the value
                    is a list of the metadata from the mapping file that
                    corresponds to the sample id
                sample_by_cat - a dictionary where the key is a tuple
                    (label,value) label is the heading for one of the columns in
                    the mapping file, and value is one of the items in the
                    column label
                len(category_labels) - the number of columns in the mapping file
                meta_dict - keyed by sample id, values are the lines from
                    mapping file without the sample id
                labels - headings of columns for edge file output
                node_labels - headings of columns for node file output
                label_list - a list containing 1 list for each header in the
                    mapping file, each list contains all of the possible values
                    for that header.
    """
    mapping_data, header, comments = parse_mapping_file(lines)
    labels = ["from", "to", "eweight", "consensus_lin"]
    node_labels = ["node_name", "node_disp_name", "ntype", "degree",
                   "weighted_degree", "consensus_lin"]
    cat_by_sample = {}
    sample_by_cat = defaultdict(list)
    meta_dict = {}
    category_labels = header[1:-1]
    labels.extend(category_labels)
    node_labels.extend(category_labels)
    label_list = [[] for c in category_labels]
    for r in mapping_data:
        categories = r[0:len(category_labels) + 1]
        sample = categories[0]
        meta_dict[sample] = ['\t'.join(categories[1:]), 0]

        cat_by_sample[sample] = [(l.strip(), c.strip())
                                 for l, c in zip(category_labels, categories[1:])]

        cat_list = []
        for i, (l, c) in enumerate(zip(category_labels, categories[1:])):
            if c not in label_list[i]:
                label_list[i].append(c)
            l = l.strip()
            c = c.strip()
            cat_list.append((l, c))
            sample_by_cat[(l, c)].append(sample)

        cat_by_sample[sample] = cat_list

    return cat_by_sample, sample_by_cat, len(category_labels), meta_dict,\
        labels, node_labels, label_list


def get_connection_info(otu_table_fp, num_meta, meta_dict):
    """Collects information from OTU table for easy use in following steps

        Input: otu_table_fh - the file path of a biom-formatted OTU table
               num_meta - the number of columns of metadata in the mapping file
               meta_dict - meta_dict returned by get_sample_info

        Output:
               con_by_sample - dictionary where the key is a sample id and the
                    value for that key is a set containing all of the sample ids
                    that the key sample shares an OTU with
               node_file - a list containing the lines for the file that defines
                    the network nodes (sample nodes and OTU nodes)
               edge_file - a list containing the lines for the file that defines
                    the network edges (connection between sample and OTU)
               red_node_file - a list containing the lines for the file that
                    defines the network nodes, where a reduced node is an OTU
                    node that is many OTU's collapsed into one node this happens
                    when more than one OTU is only connected to the same one
                    sample, then they are collapsed into one node.
               red_edge_file - a list containing the lines for the edge file
                    using the collapsed nodes.
               otu_dc - dictionary where key is a node degree and the value is
                    the number of otu nodes that have that degree
               sample_dc - dictionary where key is a node degree and the value is
                    the number of sample nodes that have that degree
               degree_counts - dictionary with all node degree counts
    """
    con_by_sample = defaultdict(set)
    node_file = []
    edge_file = []
    red_nodes = defaultdict(int)
    red_node_file = []
    red_edge_file = []
    multi = defaultdict(list)
    edge_from = []
    to = []
    otu_dc = defaultdict(int)
    degree_counts = defaultdict(int)
    sample_dc = defaultdict(int)
    sample_num_seq = defaultdict(int)
    con_list = []

    otu_table = load_table(otu_table_fp)

    is_con = False
    # This could be moved to OTU table sub-class
    if (otu_table.metadata(axis='observation') is not None and
            'taxonomy' in otu_table.metadata(axis='observation')[0]):
        is_con = True

    for otu_values, otu_id, otu_md in otu_table.iter(axis='observation'):
        con = ''
        if is_con:
            con = ':'.join(otu_md['taxonomy'][:6])
            con = con.replace(" ", "_")
            con = con.replace("\t", "_")
        if con not in con_list:
            con_list.append(con)
        non_zero_counts = otu_values.nonzero()[0]
        degree = len(non_zero_counts)
        weighted_degree = sum(otu_values)
        node_file_line = [otu_id, '', 'otu_node', str(degree),
                          str(weighted_degree), con]
        node_file_line.extend(['otu'] * num_meta)
        node_file.append('\t'.join(node_file_line))

        if len(non_zero_counts) != 1:
            red_node_file.append('\t'.join(node_file_line))

        otu_dc[degree] += 1
        degree_counts[degree] += 1
        sample_ids = otu_table.ids()
        samples = [sample_ids[i] for i in non_zero_counts]
        for i, s in enumerate(samples):
            if s not in meta_dict.keys():
                continue
            con_by_sample[s].update(samples[0:i])
            con_by_sample[s].update(samples[i + 1:])
            sample_num_seq[s] += float(otu_values[non_zero_counts[i]])

            edge_from.append(s)
            to.append(otu_id)
            meta = meta_dict[s]
            meta[1] += 1
            data_num = str(otu_values[non_zero_counts[i]])
            edge_file.append('\t'.join([s, otu_id,
                                        data_num, con, meta[0]]))
            multi[otu_id].append(
                (s, float(otu_values[non_zero_counts[i]]), meta[0]))
            if len(non_zero_counts) == 1:
                red_nodes[(sample_ids[non_zero_counts[0]], meta[0])] += degree
            else:
                red_edge_file.append('\t'.join([s, otu_id,
                                                data_num, con, meta[0]]))

    num_otu_nodes = len(node_file)
    for s in meta_dict:
        meta = meta_dict[s]
        degree = meta[1]
        sample_dc[degree] += 1
        degree_counts[degree] += 1
        weighted_degree = sample_num_seq[s]
        node_file_line = '\t'.join([s, s, 'user_node', str(meta[1]),
                                    str(weighted_degree), 'other', meta[0]])
        node_file.append(node_file_line)
        red_node_file.append(node_file_line)

    for n, d in red_nodes.items():
        red_node_file_line = ['@' + n[0], '',
                              'otu_collapsed', str(d), str(float(d)), 'other']
        red_node_file_line.extend(['otu'] * num_meta)
        red_node_file.append('\t'.join(red_node_file_line))
        red_edge_file.append(
            '\t'.join([n[0], '@' + n[0], "1.0", "missed", n[1]]))

    return con_by_sample, node_file, edge_file, red_node_file,\
        red_edge_file, otu_dc, degree_counts, sample_dc


def get_num_con_cat(con_by_sample, cat_by_sample):
    """finds the number of samples connected to the same OTU split by metadata
    """
    num_con_cat = defaultdict(float)
    num_con = 0
    for sample, connects in con_by_sample.items():
        sample_categories = cat_by_sample[sample]
        for s_con in connects:
            if s_con not in cat_by_sample.keys():
                continue
            for s_cat, con_cat in zip(sample_categories, cat_by_sample[s_con]):
                if s_cat == con_cat:
                    num_con_cat[s_cat[0]] += 0.5
            num_con += 0.5

    return num_con_cat, num_con


def get_num_cat(sample_by_cat, samples_in_otus):
    """Builds a dictionary of numbers of samples keyed by metadata value
    """
    num_cat = defaultdict(int)
    for cat, samples in sample_by_cat.items():
        num_samples = len(set(samples_in_otus) & set(samples))
        num_cat[cat[0]] += (num_samples * (num_samples - 1)) / 2
    return num_cat


def make_table_file(lines, labels, dir_path, filename):
    """Generic table writer, given lines, headers, and filename it writes to file
    """
    lines.sort()
    lines.insert(0, '\t'.join(labels))

    output = open(os.path.join(dir_path, filename), 'w')
    output.write('\n'.join(lines))
    output.close()


def make_stats_files(sample_dc, otu_dc, degree_counts, num_con_cat, num_con,
                     num_cat, cat_by_sample, dir_path):
    """Creates and writes the network statistic files including degree info

        Input: sample_dc, otu_dc, degree_counts - output from get_connection_info
               num_con_cat, num_con - output from get_num_con_cat
               num_cat - output from get_num_cat
               cat_by_sample - output from get_sample_info
               dir_path - directory path that the files should be written to
    """
    output = open(os.path.join(dir_path,
                               "stats/real_dc_sample_degree.txt"), 'w')
    sample_dc_out = sorted(sample_dc.items())
    sample_dc_str = '\n'.join(['\t'.join(map(str, t)) for t in sample_dc_out])
    output.write(''.join(["# Just Sample degree counts\n",
                          "Degree	Sample Count\n", sample_dc_str]))
    output.close()

    output = open(os.path.join(dir_path,
                               "stats/real_dc_otu_degree.txt"), 'w')
    otu_dc_out = sorted(otu_dc.items())
    otu_dc_str = '\n'.join(['\t'.join(map(str, t)) for t in otu_dc_out])
    output.write(''.join(["# Just OTU degree counts\n",
                          "Degree	OTU Count\n", otu_dc_str]))
    output.close()

    output = open(os.path.join(dir_path,
                               "stats/real_dc_sample_otu_degree.txt"), 'w')
    dc_out = sorted(degree_counts.items())
    dc_str = '\n'.join(['\t'.join(map(str, t)) for t in dc_out])
    output.write(''.join(["# Sample and OTU degree counts\n",
                          "Degree	Both Count \n", dc_str]))
    output.close()

    num_pairs = len(cat_by_sample) * (len(cat_by_sample) - 1) / 2

    num_pairs_line = "NUM PAIRS: %s" % str(num_pairs)
    num_cat_pairs_line = "NUM SAME CAT PAIRS: %s"
    num_con_pairs_line = "NUM CONNECTED PAIRS: %s" % int(num_con)

    for cat, num in num_con_cat.items():
        filename = "stats/real_cat_stats_%s.txt" % cat
        output = open(os.path.join(dir_path, filename), 'w')
        num_neither = int((num_pairs - num_con) - (num_cat[cat] - num))
        stats_line = ''.join(['(', str(int(num)), ', ', str(int(num_con - num)),
                              ', ', str(int(num_cat[cat] - num)), ', ',
                              str(num_neither), ')'])
        G_stat = G_2_by_2(int(num), int(num_con - num),
                          int(num_cat[cat] - num), num_neither)
        output.write(
            '\n'.join([num_pairs_line,
                       num_cat_pairs_line % num_cat[
                           cat],
                       num_con_pairs_line,
                       stats_line,
                       str(G_stat)]))
        output.close()


def make_props_files(labels, label_list, dir_path, data,
                     background_color, label_color, prefs):
    """Creates and writes the props file for cytoscape to automate coloring

        Input: labels - the metadata header labels
               label_list - a list containing 1 list for each header in the
                   mapping file, each list contains all of the possible values
                   for that header, this is output from get_sample_info.
               dir_path - directory path that the files should be written to
               prefs - dictionary of user color choices
               data - returned by sample_color_prefs_and_map_data_from_options
                      from qiime.colors
               background_color - color to be used for background of network
               label_color - color to be used for the labels on the network

    """
    cat_connected_num = 0
    mapping = data['map']
    groups_and_colors = iter_color_groups(mapping, prefs)
    for params in groups_and_colors:
        l = params[0]
        if l == "SampleID" or l == "Description":
            continue
        m = params[2]
        c = params[3]
        output = open(os.path.join(dir_path, "props/custom.%s.props" % l), 'w')
        props_str_list = [l] * 5
        props_str_list.append(','.join(map(str, label_color.toRGB())))
        props_str_list.extend([l] * 22)
        props_str_list.append(','.join(map(str, label_color.toRGB())))
        props_str_list.extend([l] * 16)
        props_str_list.append(props_edge % (l, l))
        props_str_list.append(l)
        props_str_list.append(
            '\n'.join([props_edge_meta % (l, s, ','.join(map(str, c[n].toRGB())))
                       for s, n in m.items()]))
        props_str_list.extend([l] * 109)
        props_str_list.append(props_node % (l, l))
        props_str_list.append(l)
        props_str_list.append(
            '\n'.join([props_node_meta % (l, s, ','.join(map(str, c[n].toRGB())))
                       for s, n in m.items()]))
        props_str_list.extend([l] * 48)
        props_str_list[98] = ','.join(map(str, background_color.toRGB()))
        props_str_list[109] = ','.join(map(str, label_color.toRGB()))
        props_str_list[132] = ','.join(map(str, label_color.toRGB()))
        output.write(props_file_str % tuple(props_str_list))
        output.close()


def create_network_and_stats(
        dir_path, map_lines, otu_table_fp, prefs, data, background_color, label_color):
    """Creates and writes the edge, node, props, and stats files for the network

        Input: dir_path - directory path that the files should be written to
               map_lines - list of lines from mapping file
               otu_table - file path of a biom-formatted OTU table file
               prefs - dictionary of user color choices
               data - returned by sample_color_prefs_and_map_data_from_options
                      from qiime.colors
               background_color - color to be used for background of network
               label_color - color to be used for the labels on the network

    """
    cat_by_sample, sample_by_cat, num_meta, meta_dict, labels, node_labels,\
        label_list = get_sample_info(map_lines)
    con_by_sample, node_file, edge_file, red_node_file,\
        red_edge_file, otu_dc, degree_counts, sample_dc, \
        = get_connection_info(otu_table_fp, num_meta, meta_dict)
    num_con_cat, num_con = get_num_con_cat(con_by_sample, cat_by_sample)
    num_cat = get_num_cat(sample_by_cat, con_by_sample.keys())
    dir_path = os.path.join(dir_path, "otu_network")
    make_table_file(edge_file, labels, dir_path, "real_edge_table.txt")
    make_table_file(node_file, node_labels, dir_path, "real_node_table.txt")
    make_table_file(red_edge_file, labels, dir_path,
                    "real_reduced_edge_table.txt")
    make_table_file(red_node_file, node_labels, dir_path,
                    "real_reduced_node_table.txt")
    make_stats_files(
        sample_dc,
        otu_dc,
        degree_counts,
        num_con_cat,
        num_con,
        num_cat,
        cat_by_sample,
        dir_path)
    if background_color == 'white':
        background_color = Color('white', (255, 255, 255))
    elif background_color == 'black':
        background_color = Color('black', (0, 0, 0))
    else:
        try:
            background_color = data_colors[background_color]
        except KeyError:
            raise KeyError("background_color unknown")

    if label_color == 'white':
        label_color = Color('white', (255, 255, 255))
    elif label_color == 'black':
        label_color = Color('black', (0, 0, 0))
    else:
        try:
            label_color = data_colors[label_color]
        except KeyError:
            raise KeyError("label_color unknown")

    make_props_files(
        labels,
        label_list,
        dir_path,
        data,
        background_color,
        label_color,
        prefs)


props_node_meta = "nodeFillColorCalculator.%s-Node Color-Discrete\ Mapper.mapping.map.%s=%s"
props_edge_meta = "edgeColorCalculator.%s-Edge\ Color-Discrete\ Mapper.mapping.map.%s=%s"
props_node = "nodeFillColorCalculator.%s-Node\ Color-Discrete\ Mapper.mapping.controller=%s"
props_edge = "edgeColorCalculator.%s-Edge\ Color-Discrete\ Mapper.mapping.controller=%s"


props_file_str = """# Properties generated by Qiime
edgeAppearanceCalculator.%s.defaultEdgeColor=255,255,0
edgeAppearanceCalculator.%s.defaultEdgeFont=SanSerif,plain,10
edgeAppearanceCalculator.%s.defaultEdgeFontSize=10.0
edgeAppearanceCalculator.%s.defaultEdgeLabel=
edgeAppearanceCalculator.%s.defaultEdgeLabelColor=%s
edgeAppearanceCalculator.%s.defaultEdgeLabelOpacity=255
edgeAppearanceCalculator.%s.defaultEdgeLabelPosition=C,C,c,0,0
edgeAppearanceCalculator.%s.defaultEdgeLineStyle=SOLID
edgeAppearanceCalculator.%s.defaultEdgeLineWidth=0.1
edgeAppearanceCalculator.%s.defaultEdgeOpacity=200
edgeAppearanceCalculator.%s.defaultEdgeSourceArrow=NONE
edgeAppearanceCalculator.%s.defaultEdgeSourceArrowColor=0,0,0
edgeAppearanceCalculator.%s.defaultEdgeSourceArrowOpacity=255
edgeAppearanceCalculator.%s.defaultEdgeSourceArrowShape=NONE
edgeAppearanceCalculator.%s.defaultEdgeTargetArrow=NONE
edgeAppearanceCalculator.%s.defaultEdgeTargetArrowColor=0,0,0
edgeAppearanceCalculator.%s.defaultEdgeTargetArrowOpacity=255
edgeAppearanceCalculator.%s.defaultEdgeTargetArrowShape=NONE
edgeAppearanceCalculator.%s.defaultEdgeToolTip=
edgeAppearanceCalculator.%s.defaultNodeBorderColor=0,0,0
edgeAppearanceCalculator.%s.defaultNodeBorderOpacity=255
edgeAppearanceCalculator.%s.defaultNodeFillColor=255,255,255
edgeAppearanceCalculator.%s.defaultNodeFont=Default,plain,12
edgeAppearanceCalculator.%s.defaultNodeFontSize=12.0
edgeAppearanceCalculator.%s.defaultNodeHight=30.0
edgeAppearanceCalculator.%s.defaultNodeLabel=
edgeAppearanceCalculator.%s.defaultNodeLabelColor=%s
edgeAppearanceCalculator.%s.defaultNodeLabelOpacity=255
edgeAppearanceCalculator.%s.defaultNodeLabelPosition=C,C,c,0,0
edgeAppearanceCalculator.%s.defaultNodeLineStyle=SOLID
edgeAppearanceCalculator.%s.defaultNodeLineWidth=1.0
edgeAppearanceCalculator.%s.defaultNodeOpacity=255
edgeAppearanceCalculator.%s.defaultNodeShape=rect
edgeAppearanceCalculator.%s.defaultNodeSize=35.0
edgeAppearanceCalculator.%s.defaultNodeToolTip=
edgeAppearanceCalculator.%s.defaultNodeWidth=70.0
edgeAppearanceCalculator.%s.edgeColorCalculator=%s-Edge Color-Discrete Mapper
edgeAppearanceCalculator.%s.edgeLineWidthCalculator=%s-Edge Line Width-Continuous Mapper
edgeAppearanceCalculator.%s.edgeOpacityCalculator=%s-Edge Opacity-Continuous Mapper
edgeAppearanceCalculator.%s.nodeSizeLocked=true
edgeColorCalculator.BasicDiscrete.mapping.controller=interaction
edgeColorCalculator.BasicDiscrete.mapping.controllerType=4
edgeColorCalculator.BasicDiscrete.mapping.map.pd=255,0,51
edgeColorCalculator.BasicDiscrete.mapping.map.pp=0,204,0
edgeColorCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeColorCalculator.BasicDiscrete.visualPropertyType=EDGE_COLOR
%s
edgeColorCalculator.%s-Edge\ Color-Discrete\ Mapper.mapping.controllerType=3
%s
edgeColorCalculator.%s-Edge\ Color-Discrete\ Mapper.mapping.type=DiscreteMapping
edgeColorCalculator.%s-Edge\ Color-Discrete\ Mapper.visualPropertyType=EDGE_COLOR
edgeFontFaceCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeFontFaceCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeFontFaceCalculator.BasicContinuous.mapping.bv0.equal=Default,bold,12
edgeFontFaceCalculator.BasicContinuous.mapping.bv0.greater=Default,plain,11
edgeFontFaceCalculator.BasicContinuous.mapping.bv0.lesser=Default,bold,12
edgeFontFaceCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeFontFaceCalculator.BasicContinuous.mapping.bv1.equal=Default,plain,11
edgeFontFaceCalculator.BasicContinuous.mapping.bv1.greater=Default,italic,10
edgeFontFaceCalculator.BasicContinuous.mapping.bv1.lesser=Default,plain,11
edgeFontFaceCalculator.BasicContinuous.mapping.controller=pvalues
edgeFontFaceCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeFontFaceCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeFontFaceCalculator.BasicContinuous.visualPropertyType=EDGE_FONT_FACE
edgeFontFaceCalculator.BasicDiscrete.mapping.controller=interaction
edgeFontFaceCalculator.BasicDiscrete.mapping.controllerType=4
edgeFontFaceCalculator.BasicDiscrete.mapping.map.pd=Serif,italic,11
edgeFontFaceCalculator.BasicDiscrete.mapping.map.pp=Serif,bold,12
edgeFontFaceCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeFontFaceCalculator.BasicDiscrete.visualPropertyType=EDGE_FONT_FACE
edgeFontSizeCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeFontSizeCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeFontSizeCalculator.BasicContinuous.mapping.bv0.equal=12.0
edgeFontSizeCalculator.BasicContinuous.mapping.bv0.greater=11.0
edgeFontSizeCalculator.BasicContinuous.mapping.bv0.lesser=12.0
edgeFontSizeCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeFontSizeCalculator.BasicContinuous.mapping.bv1.equal=11.0
edgeFontSizeCalculator.BasicContinuous.mapping.bv1.greater=10.0
edgeFontSizeCalculator.BasicContinuous.mapping.bv1.lesser=11.0
edgeFontSizeCalculator.BasicContinuous.mapping.controller=pvalues
edgeFontSizeCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeFontSizeCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeFontSizeCalculator.BasicContinuous.visualPropertyType=EDGE_FONT_SIZE
edgeFontSizeCalculator.BasicDiscrete.mapping.controller=interaction
edgeFontSizeCalculator.BasicDiscrete.mapping.controllerType=4
edgeFontSizeCalculator.BasicDiscrete.mapping.map.pd=11.0
edgeFontSizeCalculator.BasicDiscrete.mapping.map.pp=12.0
edgeFontSizeCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeFontSizeCalculator.BasicDiscrete.visualPropertyType=EDGE_FONT_SIZE
edgeLabelCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeLabelCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeLabelCalculator.BasicContinuous.mapping.bv0.equal=good
edgeLabelCalculator.BasicContinuous.mapping.bv0.greater=ok
edgeLabelCalculator.BasicContinuous.mapping.bv0.lesser=good
edgeLabelCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeLabelCalculator.BasicContinuous.mapping.bv1.equal=ok
edgeLabelCalculator.BasicContinuous.mapping.bv1.greater=poor
edgeLabelCalculator.BasicContinuous.mapping.bv1.lesser=ok
edgeLabelCalculator.BasicContinuous.mapping.controller=pvalues
edgeLabelCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeLabelCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeLabelCalculator.BasicContinuous.visualPropertyType=EDGE_LABEL
edgeLabelCalculator.BasicDiscrete.mapping.controller=interaction
edgeLabelCalculator.BasicDiscrete.mapping.controllerType=4
edgeLabelCalculator.BasicDiscrete.mapping.map.pd=pd
edgeLabelCalculator.BasicDiscrete.mapping.map.pp=pp
edgeLabelCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeLabelCalculator.BasicDiscrete.visualPropertyType=EDGE_LABEL
edgeLabelCalculator.EdgeLabel.mapping.controller=interaction
edgeLabelCalculator.EdgeLabel.mapping.type=PassThroughMapping
edgeLabelCalculator.EdgeLabel.visualPropertyType=EDGE_LABEL
edgeLabelCalculator.testPassThrough.mapping.controller=label
edgeLabelCalculator.testPassThrough.mapping.type=PassThroughMapping
edgeLabelCalculator.testPassThrough.visualPropertyType=EDGE_LABEL
edgeLineWidthCalculator.Edge\ Line\ Width-Continuous\ Mapper.mapping.boundaryvalues=0
edgeLineWidthCalculator.Edge\ Line\ Width-Continuous\ Mapper.mapping.controller=Order
edgeLineWidthCalculator.Edge\ Line\ Width-Continuous\ Mapper.mapping.interpolator=LinearNumberToNumberInterpolator
edgeLineWidthCalculator.Edge\ Line\ Width-Continuous\ Mapper.mapping.type=ContinuousMapping
edgeLineWidthCalculator.Edge\ Line\ Width-Continuous\ Mapper.visualPropertyType=EDGE_LINE_WIDTH
edgeLineWidthCalculator.Edge\ Line\ Width-Discrete\ Mapper.mapping.controller=eweight
edgeLineWidthCalculator.Edge\ Line\ Width-Discrete\ Mapper.mapping.controllerType=3
edgeLineWidthCalculator.Edge\ Line\ Width-Discrete\ Mapper.mapping.type=DiscreteMapping
edgeLineWidthCalculator.Edge\ Line\ Width-Discrete\ Mapper.visualPropertyType=EDGE_LINE_WIDTH
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.boundaryvalues=3
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.domainvalue=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.equal=2000.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.greater=2000.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.lesser=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.domainvalue=116.45779901742935
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.equal=0.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.greater=0.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.lesser=0.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv2.domainvalue=595.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv2.equal=5.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv2.greater=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv2.lesser=5.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.controller=eweight
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.interpolator=LinearNumberToNumberInterpolator
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.type=ContinuousMapping
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.visualPropertyType=EDGE_LINE_WIDTH
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.boundaryvalues=2
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.domainvalue=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.equal=2.753411
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.greater=2.753411
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv0.lesser=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.domainvalue=1135.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.equal=4.1048427
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.greater=1.0
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.bv1.lesser=4.1048427
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.controller=eweight
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.interpolator=LinearNumberToNumberInterpolator
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.mapping.type=ContinuousMapping
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Continuous\ Mapper.visualPropertyType=EDGE_LINE_WIDTH
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Discrete\ Mapper.mapping.controller=eweight
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Discrete\ Mapper.mapping.controllerType=3
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Discrete\ Mapper.mapping.type=DiscreteMapping
edgeLineWidthCalculator.%s-Edge\ Line\ Width-Discrete\ Mapper.visualPropertyType=EDGE_LINE_WIDTH
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.boundaryvalues=2
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv0.domainvalue=3.900255730841309
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv0.equal=45.23877
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv0.greater=45.23877
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv0.lesser=30.0
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv1.domainvalue=1135.0
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv1.equal=45.02375
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv1.greater=30.0
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.bv1.lesser=45.02375
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.controller=eweight
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.interpolator=LinearNumberToNumberInterpolator
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.mapping.type=ContinuousMapping
edgeOpacityCalculator.%s-Edge\ Opacity-Continuous\ Mapper.visualPropertyType=EDGE_OPACITY
edgeSourceArrowCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeSourceArrowCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeSourceArrowCalculator.BasicContinuous.mapping.bv0.equal=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.bv0.greater=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.bv0.lesser=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeSourceArrowCalculator.BasicContinuous.mapping.bv1.equal=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.bv1.greater=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.bv1.lesser=NONE
edgeSourceArrowCalculator.BasicContinuous.mapping.controller=pvalues
edgeSourceArrowCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeSourceArrowCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeSourceArrowCalculator.BasicContinuous.visualPropertyType=EDGE_SRCARROW
edgeSourceArrowCalculator.BasicDiscrete.mapping.controller=interaction
edgeSourceArrowCalculator.BasicDiscrete.mapping.controllerType=4
edgeSourceArrowCalculator.BasicDiscrete.mapping.map.pd=NONE
edgeSourceArrowCalculator.BasicDiscrete.mapping.map.pp=NONE
edgeSourceArrowCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeSourceArrowCalculator.BasicDiscrete.visualPropertyType=EDGE_SRCARROW
edgeTargetArrowCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeTargetArrowCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeTargetArrowCalculator.BasicContinuous.mapping.bv0.equal=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.bv0.greater=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.bv0.lesser=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeTargetArrowCalculator.BasicContinuous.mapping.bv1.equal=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.bv1.greater=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.bv1.lesser=NONE
edgeTargetArrowCalculator.BasicContinuous.mapping.controller=pvalues
edgeTargetArrowCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeTargetArrowCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeTargetArrowCalculator.BasicContinuous.visualPropertyType=EDGE_TGTARROW
edgeTargetArrowCalculator.BasicDiscrete.mapping.controller=interaction
edgeTargetArrowCalculator.BasicDiscrete.mapping.controllerType=4
edgeTargetArrowCalculator.BasicDiscrete.mapping.map.pd=NONE
edgeTargetArrowCalculator.BasicDiscrete.mapping.map.pp=NONE
edgeTargetArrowCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeTargetArrowCalculator.BasicDiscrete.visualPropertyType=EDGE_TGTARROW
edgeTooltipCalculator.BasicContinuous.mapping.boundaryvalues=2
edgeTooltipCalculator.BasicContinuous.mapping.bv0.domainvalue=0.0010
edgeTooltipCalculator.BasicContinuous.mapping.bv0.equal=smallest
edgeTooltipCalculator.BasicContinuous.mapping.bv0.greater=small
edgeTooltipCalculator.BasicContinuous.mapping.bv0.lesser=smallest
edgeTooltipCalculator.BasicContinuous.mapping.bv1.domainvalue=0.01
edgeTooltipCalculator.BasicContinuous.mapping.bv1.equal=small
edgeTooltipCalculator.BasicContinuous.mapping.bv1.greater=large
edgeTooltipCalculator.BasicContinuous.mapping.bv1.lesser=small
edgeTooltipCalculator.BasicContinuous.mapping.controller=pvalues
edgeTooltipCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
edgeTooltipCalculator.BasicContinuous.mapping.type=ContinuousMapping
edgeTooltipCalculator.BasicContinuous.visualPropertyType=EDGE_TOOLTIP
edgeTooltipCalculator.BasicDiscrete.mapping.controller=interaction
edgeTooltipCalculator.BasicDiscrete.mapping.controllerType=4
edgeTooltipCalculator.BasicDiscrete.mapping.map.pd=pd Tip
edgeTooltipCalculator.BasicDiscrete.mapping.map.pp=pp Tip
edgeTooltipCalculator.BasicDiscrete.mapping.type=DiscreteMapping
edgeTooltipCalculator.BasicDiscrete.visualPropertyType=EDGE_TOOLTIP
globalAppearanceCalculator.%s.defaultBackgroundColor=%s
globalAppearanceCalculator.%s.defaultEdgeReverseSelectionColor=255,0,0
globalAppearanceCalculator.%s.defaultEdgeSelectionColor=255,0,0
globalAppearanceCalculator.%s.defaultNodeReverseSelectionColor=0,255,0
globalAppearanceCalculator.%s.defaultNodeSelectionColor=255,255,0
globalAppearanceCalculator.%s.defaultSloppySelectionColor=128,128,128
nodeAppearanceCalculator.%s.defaultEdgeColor=0,0,0
nodeAppearanceCalculator.%s.defaultEdgeFont=SanSerif,plain,10
nodeAppearanceCalculator.%s.defaultEdgeFontSize=10.0
nodeAppearanceCalculator.%s.defaultEdgeLabel=
nodeAppearanceCalculator.%s.defaultEdgeLabelColor=%s
nodeAppearanceCalculator.%s.defaultEdgeLabelOpacity=255
nodeAppearanceCalculator.%s.defaultEdgeLabelPosition=C,C,c,0,0
nodeAppearanceCalculator.%s.defaultEdgeLineStyle=SOLID
nodeAppearanceCalculator.%s.defaultEdgeLineWidth=1.0
nodeAppearanceCalculator.%s.defaultEdgeOpacity=255
nodeAppearanceCalculator.%s.defaultEdgeSourceArrow=NONE
nodeAppearanceCalculator.%s.defaultEdgeSourceArrowColor=0,0,0
nodeAppearanceCalculator.%s.defaultEdgeSourceArrowOpacity=255
nodeAppearanceCalculator.%s.defaultEdgeSourceArrowShape=NONE
nodeAppearanceCalculator.%s.defaultEdgeTargetArrow=NONE
nodeAppearanceCalculator.%s.defaultEdgeTargetArrowColor=0,0,0
nodeAppearanceCalculator.%s.defaultEdgeTargetArrowOpacity=255
nodeAppearanceCalculator.%s.defaultEdgeTargetArrowShape=NONE
nodeAppearanceCalculator.%s.defaultEdgeToolTip=
nodeAppearanceCalculator.%s.defaultNodeBorderColor=0,0,0
nodeAppearanceCalculator.%s.defaultNodeBorderOpacity=100
nodeAppearanceCalculator.%s.defaultNodeFillColor=255,0,51
nodeAppearanceCalculator.%s.defaultNodeFont=Arial Narrow Bold,plain,12
nodeAppearanceCalculator.%s.defaultNodeFontSize=50.0
nodeAppearanceCalculator.%s.defaultNodeHight=30.0
nodeAppearanceCalculator.%s.defaultNodeLabel=label
nodeAppearanceCalculator.%s.defaultNodeLabelColor=%s
nodeAppearanceCalculator.%s.defaultNodeLabelOpacity=255
nodeAppearanceCalculator.%s.defaultNodeLabelPosition=C,C,c,0,0
nodeAppearanceCalculator.%s.defaultNodeLineStyle=SOLID
nodeAppearanceCalculator.%s.defaultNodeLineWidth=0.1
nodeAppearanceCalculator.%s.defaultNodeOpacity=155
nodeAppearanceCalculator.%s.defaultNodeShape=ellipse
nodeAppearanceCalculator.%s.defaultNodeSize=35.0
nodeAppearanceCalculator.%s.defaultNodeToolTip=
nodeAppearanceCalculator.%s.defaultNodeWidth=70.0
nodeAppearanceCalculator.%s.nodeFillColorCalculator=%s-Node Color-Discrete Mapper
nodeAppearanceCalculator.%s.nodeFontFaceCalculator=%s-Node Font Face-Discrete Mapper
nodeAppearanceCalculator.%s.nodeFontSizeCalculator=%s-Node Font Size-Discrete Mapper
nodeAppearanceCalculator.%s.nodeLabelCalculator=%s-Node Label-Passthrough Mapper
nodeAppearanceCalculator.%s.nodeOpacityCalculator=%s-Node Opacity-Discrete Mapper
nodeAppearanceCalculator.%s.nodeShapeCalculator=%s-Node Shape-Discrete Mapper
nodeAppearanceCalculator.%s.nodeSizeLocked=true
nodeAppearanceCalculator.%s.nodeUniformSizeCalculator=%s-Node Size-Continuous Mapper
nodeBorderColorCalculator.RedGreen.mapping.boundaryvalues=3
nodeBorderColorCalculator.RedGreen.mapping.bv0.domainvalue=-2.5
nodeBorderColorCalculator.RedGreen.mapping.bv0.equal=255,0,0
nodeBorderColorCalculator.RedGreen.mapping.bv0.greater=255,0,0
nodeBorderColorCalculator.RedGreen.mapping.bv0.lesser=0,0,255
nodeBorderColorCalculator.RedGreen.mapping.bv1.domainvalue=0.0
nodeBorderColorCalculator.RedGreen.mapping.bv1.equal=255,255,255
nodeBorderColorCalculator.RedGreen.mapping.bv1.greater=255,255,255
nodeBorderColorCalculator.RedGreen.mapping.bv1.lesser=255,255,255
nodeBorderColorCalculator.RedGreen.mapping.bv2.domainvalue=2.1
nodeBorderColorCalculator.RedGreen.mapping.bv2.equal=0,255,102
nodeBorderColorCalculator.RedGreen.mapping.bv2.greater=0,0,0
nodeBorderColorCalculator.RedGreen.mapping.bv2.lesser=0,255,102
nodeBorderColorCalculator.RedGreen.mapping.controller=gal1RGexp
nodeBorderColorCalculator.RedGreen.mapping.interpolator=LinearNumberToColorInterpolator
nodeBorderColorCalculator.RedGreen.mapping.type=ContinuousMapping
nodeBorderColorCalculator.RedGreen.visualPropertyType=NODE_BORDER_COLOR
nodeFillColorCalculator.Node\ Color-Discrete\ Mapper.mapping.controller=Order
nodeFillColorCalculator.Node\ Color-Discrete\ Mapper.mapping.controllerType=-1
nodeFillColorCalculator.Node\ Color-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeFillColorCalculator.Node\ Color-Discrete\ Mapper.visualPropertyType=NODE_FILL_COLOR
nodeFillColorCalculator.RedGreen.mapping.boundaryvalues=3
nodeFillColorCalculator.RedGreen.mapping.bv0.domainvalue=-2.5
nodeFillColorCalculator.RedGreen.mapping.bv0.equal=255,0,0
nodeFillColorCalculator.RedGreen.mapping.bv0.greater=255,0,0
nodeFillColorCalculator.RedGreen.mapping.bv0.lesser=0,0,255
nodeFillColorCalculator.RedGreen.mapping.bv1.domainvalue=0.0
nodeFillColorCalculator.RedGreen.mapping.bv1.equal=255,255,255
nodeFillColorCalculator.RedGreen.mapping.bv1.greater=255,255,255
nodeFillColorCalculator.RedGreen.mapping.bv1.lesser=255,255,255
nodeFillColorCalculator.RedGreen.mapping.bv2.domainvalue=2.1
nodeFillColorCalculator.RedGreen.mapping.bv2.equal=0,255,102
nodeFillColorCalculator.RedGreen.mapping.bv2.greater=0,0,0
nodeFillColorCalculator.RedGreen.mapping.bv2.lesser=0,255,102
nodeFillColorCalculator.RedGreen.mapping.controller=gal1RGexp
nodeFillColorCalculator.RedGreen.mapping.interpolator=LinearNumberToColorInterpolator
nodeFillColorCalculator.RedGreen.mapping.type=ContinuousMapping
nodeFillColorCalculator.RedGreen.visualPropertyType=NODE_FILL_COLOR
%s
nodeFillColorCalculator.%s-Node\ Color-Discrete\ Mapper.mapping.controllerType=3
%s
nodeFillColorCalculator.%s-Node\ Color-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeFillColorCalculator.%s-Node\ Color-Discrete\ Mapper.visualPropertyType=NODE_FILL_COLOR
nodeFontFaceCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeFontFaceCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeFontFaceCalculator.BasicContinuous.mapping.bv0.equal=Serif,italic,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv0.greater=Serif,italic,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv0.lesser=Serif,italic,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeFontFaceCalculator.BasicContinuous.mapping.bv1.equal=null,plain,10
nodeFontFaceCalculator.BasicContinuous.mapping.bv1.greater=Serif,bold,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv1.lesser=Serif,italic,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeFontFaceCalculator.BasicContinuous.mapping.bv2.equal=Serif,bold,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv2.greater=Serif,bold,12
nodeFontFaceCalculator.BasicContinuous.mapping.bv2.lesser=Serif,bold,12
nodeFontFaceCalculator.BasicContinuous.mapping.controller=expression
nodeFontFaceCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
nodeFontFaceCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeFontFaceCalculator.BasicContinuous.visualPropertyType=NODE_FONT_FACE
nodeFontFaceCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeFontFaceCalculator.BasicDiscrete.mapping.controllerType=-1
nodeFontFaceCalculator.BasicDiscrete.mapping.map.autophagy=Serif,italic,11
nodeFontFaceCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=Serif,plain,10
nodeFontFaceCalculator.BasicDiscrete.mapping.map.cell\ cycle=Serif,bold,12
nodeFontFaceCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeFontFaceCalculator.BasicDiscrete.visualPropertyType=NODE_FONT_FACE
nodeFontFaceCalculator.%s-Node\ Font\ Face-Discrete\ Mapper.mapping.controller=ntype
nodeFontFaceCalculator.%s-Node\ Font\ Face-Discrete\ Mapper.mapping.controllerType=4
nodeFontFaceCalculator.%s-Node\ Font\ Face-Discrete\ Mapper.mapping.map.user_node=Arial Narrow Bold,plain,12
nodeFontFaceCalculator.%s-Node\ Font\ Face-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeFontFaceCalculator.%s-Node\ Font\ Face-Discrete\ Mapper.visualPropertyType=NODE_FONT_FACE
nodeFontSizeCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeFontSizeCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv0.equal=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv0.greater=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv0.lesser=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv1.equal=10.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv1.greater=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv1.lesser=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv2.equal=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv2.greater=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.bv2.lesser=12.0
nodeFontSizeCalculator.BasicContinuous.mapping.controller=expression
nodeFontSizeCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
nodeFontSizeCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeFontSizeCalculator.BasicContinuous.visualPropertyType=NODE_FONT_SIZE
nodeFontSizeCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeFontSizeCalculator.BasicDiscrete.mapping.controllerType=-1
nodeFontSizeCalculator.BasicDiscrete.mapping.map.autophagy=11.0
nodeFontSizeCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=10.0
nodeFontSizeCalculator.BasicDiscrete.mapping.map.cell\ cycle=12.0
nodeFontSizeCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeFontSizeCalculator.BasicDiscrete.visualPropertyType=NODE_FONT_SIZE
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.controller=ntype
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.controllerType=4
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.map.otu_collapsed=10.0
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.map.otu_node=10.0
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.map.user_node=28.0
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeFontSizeCalculator.%s-Node\ Font\ Size-Discrete\ Mapper.visualPropertyType=NODE_FONT_SIZE
nodeHeightCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeHeightCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeHeightCalculator.BasicContinuous.mapping.bv0.equal=100.0
nodeHeightCalculator.BasicContinuous.mapping.bv0.greater=100.0
nodeHeightCalculator.BasicContinuous.mapping.bv0.lesser=100.0
nodeHeightCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeHeightCalculator.BasicContinuous.mapping.bv1.equal=50.0
nodeHeightCalculator.BasicContinuous.mapping.bv1.greater=50.0
nodeHeightCalculator.BasicContinuous.mapping.bv1.lesser=50.0
nodeHeightCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeHeightCalculator.BasicContinuous.mapping.bv2.equal=100.0
nodeHeightCalculator.BasicContinuous.mapping.bv2.greater=100.0
nodeHeightCalculator.BasicContinuous.mapping.bv2.lesser=100.0
nodeHeightCalculator.BasicContinuous.mapping.controller=expression
nodeHeightCalculator.BasicContinuous.mapping.interpolator=LinearNumberToNumberInterpolator
nodeHeightCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeHeightCalculator.BasicContinuous.visualPropertyType=NODE_HEIGHT
nodeHeightCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeHeightCalculator.BasicDiscrete.mapping.controllerType=-1
nodeHeightCalculator.BasicDiscrete.mapping.map.autophagy=50.0
nodeHeightCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=10.0
nodeHeightCalculator.BasicDiscrete.mapping.map.cell\ cycle=60.0
nodeHeightCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeHeightCalculator.BasicDiscrete.visualPropertyType=NODE_HEIGHT
nodeLabelCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeLabelCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeLabelCalculator.BasicContinuous.mapping.bv0.equal=under
nodeLabelCalculator.BasicContinuous.mapping.bv0.greater=under
nodeLabelCalculator.BasicContinuous.mapping.bv0.lesser=under
nodeLabelCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeLabelCalculator.BasicContinuous.mapping.bv1.equal=zero
nodeLabelCalculator.BasicContinuous.mapping.bv1.greater=over
nodeLabelCalculator.BasicContinuous.mapping.bv1.lesser=under
nodeLabelCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeLabelCalculator.BasicContinuous.mapping.bv2.equal=over
nodeLabelCalculator.BasicContinuous.mapping.bv2.greater=over
nodeLabelCalculator.BasicContinuous.mapping.bv2.lesser=over
nodeLabelCalculator.BasicContinuous.mapping.controller=expression
nodeLabelCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
nodeLabelCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeLabelCalculator.BasicContinuous.visualPropertyType=NODE_LABEL
nodeLabelCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeLabelCalculator.BasicDiscrete.mapping.controllerType=-1
nodeLabelCalculator.BasicDiscrete.mapping.map.autophagy=autophagy
nodeLabelCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=unknown
nodeLabelCalculator.BasicDiscrete.mapping.map.cell\ cycle=cell cycle
nodeLabelCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeLabelCalculator.BasicDiscrete.visualPropertyType=NODE_LABEL
nodeLabelCalculator.BasicPassThrough.mapping.controller=label
nodeLabelCalculator.BasicPassThrough.mapping.type=PassThroughMapping
nodeLabelCalculator.BasicPassThrough.visualPropertyType=NODE_LABEL
nodeLabelCalculator.NodeLabel.mapping.controller=ID
nodeLabelCalculator.NodeLabel.mapping.type=PassThroughMapping
nodeLabelCalculator.NodeLabel.visualPropertyType=NODE_LABEL
nodeLabelCalculator.%s-Node\ Label-Passthrough\ Mapper.mapping.controller=node_disp_name
nodeLabelCalculator.%s-Node\ Label-Passthrough\ Mapper.mapping.type=PassThroughMapping
nodeLabelCalculator.%s-Node\ Label-Passthrough\ Mapper.visualPropertyType=NODE_LABEL
nodeLabelCalculator.id\ as\ node\ labels.mapping.controller=ID
nodeLabelCalculator.id\ as\ node\ labels.mapping.type=PassThroughMapping
nodeLabelCalculator.id\ as\ node\ labels.visualPropertyType=NODE_LABEL
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.controller=ntype
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.controllerType=4
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.map.otu_collapsed=210.0
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.map.otu_node=220.0
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.map.user_node=200.0
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeOpacityCalculator.%s-Node\ Opacity-Discrete\ Mapper.visualPropertyType=NODE_OPACITY
nodeShapeCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeShapeCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeShapeCalculator.BasicContinuous.mapping.bv0.equal=ellipse
nodeShapeCalculator.BasicContinuous.mapping.bv0.greater=ellipse
nodeShapeCalculator.BasicContinuous.mapping.bv0.lesser=ellipse
nodeShapeCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeShapeCalculator.BasicContinuous.mapping.bv1.equal=diamond
nodeShapeCalculator.BasicContinuous.mapping.bv1.greater=rect
nodeShapeCalculator.BasicContinuous.mapping.bv1.lesser=ellipse
nodeShapeCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeShapeCalculator.BasicContinuous.mapping.bv2.equal=rect
nodeShapeCalculator.BasicContinuous.mapping.bv2.greater=rect
nodeShapeCalculator.BasicContinuous.mapping.bv2.lesser=rect
nodeShapeCalculator.BasicContinuous.mapping.controller=expression
nodeShapeCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
nodeShapeCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeShapeCalculator.BasicContinuous.visualPropertyType=NODE_SHAPE
nodeShapeCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeShapeCalculator.BasicDiscrete.mapping.controllerType=-1
nodeShapeCalculator.BasicDiscrete.mapping.map.autophagy=rect
nodeShapeCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=ellipse
nodeShapeCalculator.BasicDiscrete.mapping.map.cell\ cycle=diamond
nodeShapeCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeShapeCalculator.BasicDiscrete.visualPropertyType=NODE_SHAPE
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.controller=ntype
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.controllerType=4
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.map.otu_collapsed=diamond
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.map.otu_node=roundrect
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.map.user_node=ellipse
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.mapping.type=DiscreteMapping
nodeShapeCalculator.%s-Node\ Shape-Discrete\ Mapper.visualPropertyType=NODE_SHAPE
nodeTooltipCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeTooltipCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeTooltipCalculator.BasicContinuous.mapping.bv0.equal=less
nodeTooltipCalculator.BasicContinuous.mapping.bv0.greater=less
nodeTooltipCalculator.BasicContinuous.mapping.bv0.lesser=less
nodeTooltipCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeTooltipCalculator.BasicContinuous.mapping.bv1.equal=same
nodeTooltipCalculator.BasicContinuous.mapping.bv1.greater=more
nodeTooltipCalculator.BasicContinuous.mapping.bv1.lesser=less
nodeTooltipCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeTooltipCalculator.BasicContinuous.mapping.bv2.equal=more
nodeTooltipCalculator.BasicContinuous.mapping.bv2.greater=more
nodeTooltipCalculator.BasicContinuous.mapping.bv2.lesser=more
nodeTooltipCalculator.BasicContinuous.mapping.controller=expression
nodeTooltipCalculator.BasicContinuous.mapping.interpolator=FlatInterpolator
nodeTooltipCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeTooltipCalculator.BasicContinuous.visualPropertyType=NODE_TOOLTIP
nodeTooltipCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeTooltipCalculator.BasicDiscrete.mapping.controllerType=-1
nodeTooltipCalculator.BasicDiscrete.mapping.map.autophagy=autophagy Tip
nodeTooltipCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=unknown Tip
nodeTooltipCalculator.BasicDiscrete.mapping.map.cell\ cycle=cell cycle Tip
nodeTooltipCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeTooltipCalculator.BasicDiscrete.visualPropertyType=NODE_TOOLTIP
nodeUniformSizeCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeUniformSizeCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv0.equal=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv0.greater=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv0.lesser=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv1.equal=50.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv1.greater=50.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv1.lesser=50.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv2.equal=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv2.greater=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.bv2.lesser=100.0
nodeUniformSizeCalculator.BasicContinuous.mapping.controller=expression
nodeUniformSizeCalculator.BasicContinuous.mapping.interpolator=LinearNumberToNumberInterpolator
nodeUniformSizeCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeUniformSizeCalculator.BasicContinuous.visualPropertyType=NODE_SIZE
nodeUniformSizeCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeUniformSizeCalculator.BasicDiscrete.mapping.controllerType=-1
nodeUniformSizeCalculator.BasicDiscrete.mapping.map.autophagy=50.0
nodeUniformSizeCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=10.0
nodeUniformSizeCalculator.BasicDiscrete.mapping.map.cell\ cycle=60.0
nodeUniformSizeCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeUniformSizeCalculator.BasicDiscrete.visualPropertyType=NODE_SIZE
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.boundaryvalues=3
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv0.domainvalue=1.0
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv0.equal=9.968475
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv0.greater=9.968475
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv0.lesser=1.0
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv1.domainvalue=292.13042813539505
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv1.equal=24.050201
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv1.greater=24.050201
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv1.lesser=24.050201
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv2.domainvalue=460.0
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv2.equal=30.955704
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv2.greater=1.0
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.bv2.lesser=30.955704
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.controller=degree
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.interpolator=LinearNumberToNumberInterpolator
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.mapping.type=ContinuousMapping
nodeUniformSizeCalculator.%s-Node\ Size-Continuous\ Mapper.visualPropertyType=NODE_SIZE
nodeWidthCalculator.BasicContinuous.mapping.boundaryvalues=3
nodeWidthCalculator.BasicContinuous.mapping.bv0.domainvalue=-1.0
nodeWidthCalculator.BasicContinuous.mapping.bv0.equal=100.0
nodeWidthCalculator.BasicContinuous.mapping.bv0.greater=100.0
nodeWidthCalculator.BasicContinuous.mapping.bv0.lesser=100.0
nodeWidthCalculator.BasicContinuous.mapping.bv1.domainvalue=0.0
nodeWidthCalculator.BasicContinuous.mapping.bv1.equal=50.0
nodeWidthCalculator.BasicContinuous.mapping.bv1.greater=50.0
nodeWidthCalculator.BasicContinuous.mapping.bv1.lesser=50.0
nodeWidthCalculator.BasicContinuous.mapping.bv2.domainvalue=1.0
nodeWidthCalculator.BasicContinuous.mapping.bv2.equal=100.0
nodeWidthCalculator.BasicContinuous.mapping.bv2.greater=100.0
nodeWidthCalculator.BasicContinuous.mapping.bv2.lesser=100.0
nodeWidthCalculator.BasicContinuous.mapping.controller=expression
nodeWidthCalculator.BasicContinuous.mapping.interpolator=LinearNumberToNumberInterpolator
nodeWidthCalculator.BasicContinuous.mapping.type=ContinuousMapping
nodeWidthCalculator.BasicContinuous.visualPropertyType=NODE_WIDTH
nodeWidthCalculator.BasicDiscrete.mapping.controller=GO Biological Process (level 4)
nodeWidthCalculator.BasicDiscrete.mapping.controllerType=-1
nodeWidthCalculator.BasicDiscrete.mapping.map.autophagy=50.0
nodeWidthCalculator.BasicDiscrete.mapping.map.biological_process\ unknown=10.0
nodeWidthCalculator.BasicDiscrete.mapping.map.cell\ cycle=60.0
nodeWidthCalculator.BasicDiscrete.mapping.type=DiscreteMapping
nodeWidthCalculator.BasicDiscrete.visualPropertyType=NODE_WIDTH"""
