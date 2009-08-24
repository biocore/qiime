#!/usr/bin/env python

"""Collapses tree nodes based on similarity branchlength and taxonomy"""

from numpy import arange, average
from collections import defaultdict
from operator import itemgetter
from optparse import make_option, OptionParser
from cogent.parse.tree import DndParser
from cogent.seqsim.tree import RangeNode

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

options = [make_option('--taxonomy-file',dest='taxonomy_file',default=None),
           make_option('--tree-file',dest='tree_file',default=None),
           make_option('--start-range',dest='start_range',type='float',
               default=0.01),
           make_option('--end-range',dest='end_range',type='float', 
               default=0.15),
           make_option('--step-size',dest='step_size',type='float', 
               default=0.01),
           make_option('--output-name',dest='output_name', default=None),
           make_option('--branch-length-correction',dest='bl_correction', \
                   default=1.0, type='float')]

class CollapseNode(RangeNode):
    def __init__(self, *args, **kwargs):
        self.AtThreshold = 0.0
        self.TaxonId = None
        self.TaxonomyInformation = {}
        super(CollapseNode, self).__init__(*args, **kwargs)
    
    def setRepresentativeTaxonId(self):
        """Sets representative id"""
        # collect descending taxonomy information for a given type if assigned
        all_descending_info = []
        for t in self.tips():
            taxon_id = t.TaxonId
            if taxon_id is None:
                continue

            all_descending_info.append(taxon_id)

        # no assigned taxonomic information descends
        if not all_descending_info:
            return

        # pick the most representative taxonomic information
        counts = defaultdict(int)
        for i in all_descending_info:
            counts[i] += 1
        best_taxon, best_count = sorted(counts.items(), key=itemgetter(1))[-1]
        self.TaxonId = best_taxon

    def setTaxonInfoAtTips(self, tax_lookup, tax_filter_f):
        """Sets taxonomy id at tips"""
        for n in self.tips():
            n.TaxonomyInformation = tax_lookup.get(n.Name, None)
            if n.TaxonomyInformation is None:
                #print "no taxonomy information:", n.Name
                n.TaxonId = None
            else:
                n.TaxonId = tax_filter_f(n.TaxonomyInformation)

ignore_strings = ['environmental samples','unclassified',
                  'unknown','environmental sequence']

def ncbi_taxonid_filter(taxonomy_information):
    """Returns taxon id, or None if unassigned"""
    info = taxonomy_information['ncbi_tax_string_format_2'].lower()
    found = False
    for i in ignore_strings:
        if i in info:
            found = True
            break
    if found:
        return None
    else:
        return taxonomy_information['ncbi_tax_id']

def pick_tree_OTUs(tree, threshold):
    """groups sequences in tree based on BL threshold

    part of cathy's node collapsing code
    """
    OTUs = []
    assigned = []
    #do a tip to root tree traversal
    for tip in tree.tips():
        if tip in assigned:
            continue

        #get OTU for node
        OTU_node = get_tip_relatives(tip, threshold)
        #if OTU_node.isTip() and not OTU_node.TaxonId:
        #    continue
        OTU_node.AtThreshold = threshold
        OTU_node.setRepresentativeTaxonId()
        assigned.extend(OTU_node.tips())
        OTUs.append(OTU_node)
   
    # lop off children
    for n in OTUs:
        n.Children = []

def get_tip_relatives(tip, threshold):
    """gets deepest node that is not past the branch length threshold

    part of cathy's node collapsing code
    """
    ancestors = tip.ancestors()
    ancestors.insert(0, tip)
    for index, node in enumerate(ancestors):
        if node == tip:
            continue

        node_tip_dists = [node.distance(x) for x in node.tips()]

        if not len(node_tip_dists):
            continue

        avg_bl = average(node_tip_dists)

        if avg_bl > threshold and index == 0:
            return tip
        elif avg_bl > threshold and index > 0:
            return ancestors[index-1]
    # if we make it through the forloop...
    return tip

def parse_taxonomy_file(ids_at_tips, lines, id_type='prokMSA_id'):
    """ """
    # first line, toss #
    header = lines[0][1:].strip()
    header_fields = header.split('\t')
    ids = set(ids_at_tips)
    result = {}

    for line in lines[1:]:
        curr = {}
        fields = line.strip().split('\t')
        for h,f in zip(header_fields, fields):
            curr[h] = f
        curr_id = curr[id_type]
        if curr_id not in ids:
            continue
        else:
            result[curr_id] = curr

    return result

def main():
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args()

    tree = DndParser(open(opts.tree_file), constructor=CollapseNode)
    for n in tree.traverse():
        if n.Length:
            n.Length *= opts.bl_correction
    
    print "loaded tree"
    ids_at_tips = [n.Name for n in tree.tips()]
    taxonfile = open(opts.taxonomy_file)
    taxon_lookup = parse_taxonomy_file(ids_at_tips, taxonfile.readlines())
    print "parsed taxonomy file"

    tree.assignIds()
    print "assigned ids"

    tree.setTaxonInfoAtTips(taxon_lookup, ncbi_taxonid_filter)
    print "set taxon information"

    for threshold in arange(opts.start_range, opts.end_range, opts.step_size):
        pick_tree_OTUs(tree, threshold)
        print "set for threshold", threshold
        lines = []
        for t in tree.tips():
            #if not t.TaxonId:
            #    continue
            
            line = [t.Id, t.TaxonId, t.AtThreshold]
            print '\t'.join(map(str, line))
        ### need to spit out info from tips... blah blah...
        ### node id, taxon id, threshold from which it came
if __name__ == '__main__':
    main()

