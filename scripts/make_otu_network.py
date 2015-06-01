#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Julia Goodrich", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"

"""
This script generates the otu networks and statistics

"""

from qiime.util import parse_command_line_parameters, create_dir
from qiime.util import make_option
from qiime.make_otu_network import create_network_and_stats
from tempfile import mkdtemp
from qiime.colors import sample_color_prefs_and_map_data_from_options
import os
import shutil

script_info = {}
script_info[
    'brief_description'] = """Make an OTU network and calculate statistics"""
script_info['script_description'] = """This script generates the otu network files to be passed into cytoscape and statistics for those networks. It uses the OTU fileand the user metadata mapping file.

Network-based analysis is used to display and analyze how OTUs are partitioned between samples. This is a powerful way to display visually large and highly complex datasets in such a way that similarities and differences between samples are emphasized. The visual output of this analysis is a clustering of samples according to their shared OTUs - samples that share more OTUs cluster closer together. The degree to which samples cluster is based on the number of OTUs shared between samples (when OTUs are found in more than one sample) and this is weighted according to the number of sequences within an OTU. In the network diagram, there are two kinds of "nodes" represented, OTU-nodes and sample-nodes. These are shown with symbols such as filled circles and filled squares. If an OTU is found within a sample, the two nodes are connected with a line (an "edge"). (OTUs found only in one sample are given a second, distinct OTU-node shape.) The nodes and edges can then be colored to emphasize certain aspects of the data. For instance, in the initial application of this analysis in a microbial ecology study, the gut bacteria of a variety of mammals was surveyed, and the network diagrams were colored according to the diets of the animals, which highlighted the clustering of hosts by diet category (herbivores, carnivores, omnivores). In a meta-analysis of bacterial surveys across habitat types, the networks were colored in such a way that the phylogenetic classification of the OTUs was highlighted: this revealed the dominance of shared Firmicutes in vertebrate gut samples versus a much higher diversity of phyla represented amongst the OTUs shared by environmental samples.

Not just pretty pictures: the connections within the network are analyzed statistically to provide support for the clustering patterns displayed in the network. A G-test for independence is used to test whether sample-nodes within categories (such as diet group for the animal example used above) are more connected within than a group than expected by chance. Each pair of samples is classified according to whether its members shared at least one OTU, and whether they share a category. Pairs are then tested for independence in these categories (this asks whether pairs that share a category also are equally likely to share an OTU). This statistical test can also provide support for an apparent lack of clustering when it appears that a parameter is not contributing to the clustering.

This OTU-based approach to comparisons between samples provides a counterpoint to the tree-based PCoA graphs derived from the UniFrac analyses. In most studies, the two approaches reveal the same patterns. They can reveal different aspects of the data, however. The network analysis can provide phylogenetic information in a visual manner, whereas PCoA-UniFrac clustering can reveal subclusters that may be obscured in the network. The PCs can be pulled out individually and regressed against other metadata; the network analysis can provide a visual display of shared versus unique OTUs. Thus, together these tools can be used to draw attention to disparate aspects of a dataset, as desired by the author.

In more technical language: OTUs and samples are designated as two types of nodes in a bipartite network in which OTU-nodes are connected via edges to sample-nodes in which their sequences are found. Edge weights are defined as the number of sequences in an OTU. To cluster the OTUs and samples in the network, a stochastic spring-embedded algorithm is used, where nodes act like physical objects that repel each other, and connections act a springs with a spring constant and a resting length: the nodes are organized in a way that minimized forces in the network. These algorithms are implemented in Cytoscape (Shannon et al., 2003)."""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Example:""",
     """Create network cytoscape and statistic files in a user-specified output directory. This example uses an OTU table (-i) and the metadata mapping file (-m), and the results are written to the "otu_network/" folder.""",
     """%prog -i otu_table.biom -m Fasting_Map.txt -o otu_network"""))

script_info[
    'output_description'] = """The result of make_otu_network.py consists of a folder which contains edge and node files to be loaded into cytoscape along with props files labeled by category, which can used for coloring."""


script_info['required_options'] = [
    make_option('-i', '--input_fp', type='existing_filepath',
                help='name of otu table file in biom format [REQUIRED]'),
    # note that the options list gets passed around, so it is required that
    # the option be called --map_fname for this value - an annoying name
    # but not refactoring now...
    make_option('-m', '--map_fname', type='existing_filepath',
                help='name of input map file [REQUIRED]'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='output directory for all analyses [REQUIRED]')
]

script_info['optional_options'] = [
    make_option('-b', '--colorby', dest='colorby', type='string',
                help='This is the categories to color by in the plots from the \
user-generated mapping file. The categories must match the name of a column \
header in the mapping file exactly and multiple categories can be list by comma \
separating them without spaces. The user can also combine columns in the \
mapping file by separating the categories by "&&" without spaces \
[default=%default]'),
    make_option('-p', '--prefs_path', help='This is the user-generated preferences \
file. NOTE: This is a file with a dictionary containing preferences for the \
analysis [default: %default]', type='existing_filepath'),
    make_option('-k', '--background_color', type='string', help='This is the background color to \
use in the plots. [default: %default]'),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    prefs, data, background_color, label_color, ball_scale, arrow_colors = \
        sample_color_prefs_and_map_data_from_options(opts)

    dir_path = opts.output_dir

    create_dir(dir_path)
    create_dir(os.path.join(dir_path, "otu_network"))
    create_dir(os.path.join(dir_path, "otu_network/props"))
    create_dir(os.path.join(dir_path, "otu_network/stats"))

    map_lines = open(opts.map_fname, 'U').readlines()
    create_network_and_stats(
        dir_path,
        map_lines,
        opts.input_fp,
        prefs,
        data,
        background_color,
        label_color)

if __name__ == "__main__":
    main()
