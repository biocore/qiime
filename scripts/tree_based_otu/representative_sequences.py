#!/usr/bin/env python

from cogent.core.alignment import Alignment
from cogent.parse.tree import DndParser
from cogent.parse.fasta import MinimalFastaParser
from cogent.seqsim.tree import RangeNode
from optparse import OptionParser, make_option
from collections import defaultdict
from operator import add
from random import choice

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

options = [make_option('--fasta-file',dest='fasta_file',default=None),
           make_option('--tree-file',dest='tree_file',default=None),
           make_option('--otu-file',dest='otu_file',default=None),
           make_option('--degap',dest='degap',default=False,\
                   action='store_true')]

class NoRepresentativeSequence(Exception):
    pass

def random_tiebreaker_f(aln, seqids, random_f=choice):
    """Returns a random seqid"""
    if len(seqids) == 1:
        return seqids[0]
    return random_f(seqids)

class RepSeqNode(RangeNode):
    def setTipNames(self):
        """Sets and caches descending tip names"""
        for n in self.postorder():
            if n.isTip():
                n.TipNames = [n.Name]
            else:
                n.TipNames = reduce(add, [c.TipNames for c in n.Children])

class RepresentativeSequence(object):
    """Abstract base class for obtaining a representative sequence"""
    def __init__(self, aln, otu_map, tiebreaker_f=None, params=None):
        """Set input data

        aln : alignment contiaining all sequences
        otu_map : otu -> sequence ids
        tiebreaker_f : f(aln, seqids), returns "best" seqid in case of tie
        params : additional subclass params
        """
        self.Alignment = aln
        self.OTUmap = otu_map

        if tiebreaker_f is None:
            tiebreaker_f = random_tiebreaker_f
        self.TieBreaker = tiebreaker_f

        if params is None:
            params = {}
        self.Params = params

    def __call__(self):
        """Returns an alignment of representative sequences and otu mapping"""
        otu_to_seq_id = {}

        for otu in self.OTUmap:
            otu_to_seq_id[otu] = self[otu]
        rep_aln = self.Alignment.takeSeqs(otu_to_seq_id.values())

        return rep_aln, otu_to_seq_id

    def __getitem__(self, arg):
        """Wrapper for self._get_representative_sequence_id"""
        return self._get_representative_seq_id(arg)

    def _get_representative_seq_id(self, otu):
        """Returns the representative sequence id

        NoRepresentativeSequence is raise if a representative sequence could not
        be determined.

        NOTE: method must be defined in subclass
        """
        raise NotImplementedError

    def _get_otu_cluster(self, otu):
        """Returns a new SequenceCollection object for the OTU cluster"""
        res = self.Alignment.takeSeqs(self.OTUmap[otu])

        if not res:
            raise NoRepresentativeSequence, "No sequence ids map to %s" % otu
        else:
            return res

class LongestSequence(RepresentativeSequence):
    """Returns longest sequence for an OTU cluster as representative"""
    def _get_representative_seq_id(self, otu):
        """Returns longest sequence for a given cluster"""
        try:
            seq_cluster = self._get_otu_cluster(otu).degap()
        except:
            print otu
            print type(otu)
            print self.OTUmap[otu]
            print self.Alignment.Names[:100]
            raise SystemExit
        len_hist = defaultdict(list)
        for seqid, seq in seq_cluster.items():
            len_hist[len(seq)].append(seqid)

        longest = max(len_hist.keys())

        return self.TieBreaker(seq_cluster, len_hist[longest])

# not positive on the name here...
class BestSharedSubstring(RepresentativeSequence):
    """Returns the sequence id with the best shared substring

    - find length such that MinSeqPercent of seqs are >= that length
    - truncate all reads to that length
    - determine most identical truncated string
    - representative sequence is longest that contains the truncated string
    """
    MinSeqPercent = 0.95

    def _get_representative_seq_id(self, otu):
        """Returns most highly represented sequence by substring"""
        seq_cluster = self._get_otu_cluster(otu)

        length_idx = itemgetter(0)
        seqid_idx = itemgetter(1)

        min_len_percent = self.Params.get('MinSeqPercent',self.MinSeqPercent)

        seqs_lengths = []
        for seqid, seq in seq_cluster.items():
            seqs_lengths.append((len(seq), seqid))

        sorted_lengths = sorted(seqs_lengths, key=length_idx)[::-1]

        n_to_keep = len(sorted_lengths) * min_len_percent
        to_keep = sorted_lengths[:n_to_keep]
        min_len = length_idx(to_keep[-1])
        reduced_cluster = seq_cluster.takeSeqs(map(seqid_idx, to_keep))

        ### do we assume alignment here? gaps? do we need to realign?
        raise NotImplementedError
        
class RandomSequence(RepresentativeSequence):
    """Returns a random sequence for an OTU cluster"""
    RandomF = choice
    def _get_representative_seq_id(self, otu):
        """Returns a random sequence id from an otu cluster"""
        rand_f = self.Params.get('RandomF',self.RandomF)
        return rand_f(self.OTUmap[otu])

def parse_otu_lines(lines):
    """Parses otu file. Expected format: OTU\tTaxonId\tThreshold

    Returns a dict, Thresold -> OTUs
    """
    d = defaultdict(list)

    for line in lines:
#   
        otu_id, taxon_id, threshold = line.strip().split()
        d[threshold].append(otu_id)

    return d

def load_tree(tree_str):
    """Loads tree, assigns ids, returns tree and dict indexed by id"""
    tree = DndParser(tree_str, constructor=RepSeqNode)
    tree.assignIds()
    for n in tree.preorder(): # RangeNode traverse is just tips
        n.Id = str(n.Id)
    id_dict = tree.indexByAttr('Id')
    return tree, id_dict

get_first_field = lambda x: x.split()[0]
replace_bad_gap = lambda x: x.replace('.','-')

def load_alignment(lines, seqid_f=get_first_field, seq_f=replace_bad_gap):
    """Cleans seq_ids and sequences, returns Alignment"""
    if seqid_f is None:
        seqid_f = lambda x: x
    if seq_f is None:
        seq_f = lambda x: x

    seq_dict = {}
    for seqid, seq in MinimalFastaParser(lines):
        seq_dict[seqid_f(seqid)] = seq_f(seq)

    return Alignment(seq_dict)

def otus_to_seqids(OTUs, node_dict):
    """Returns mapping otu->[seqids]"""
    ##### no need to call tips. start high threshold -> low, cache tip names
    ##### @ nodes. only need a single call per high threshold otu to 
    ##### subsets?? and can then just query the cache attr

    res = {}
    for otu in OTUs:
        node = node_dict[otu]
        if not hasattr(node, 'TipNames'):
            node.setTipNames()
        res[otu] = node.TipNames
    return res

def main():
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args()

    tree, node_dict = load_tree(open(opts.tree_file))
    aln = load_alignment(open(opts.fasta_file))
    threshold_to_otus = parse_otu_lines(open(opts.otu_file))

    for threshold in sorted(threshold_to_otus)[::-1]:
        OTUs = threshold_to_otus[threshold]
        otu_map = otus_to_seqids(OTUs, node_dict)
        ls = LongestSequence(aln, otu_map)
        result_aln, result_map = ls()

        f = open(opts.fasta_file + '-OTUs_at_%s.fasta' % threshold, 'w')
        if opts.degap:
            f.write(result_aln.degap().toFasta())
        else:
            f.write(result_aln.toFasta())
        f.close()

        # same as seqid -> otu
        res_map = ['\t'.join([v,k]) for k,v in result_map.items()]
        f = open(opts.fasta_file + '-OTUs_at_%s.map' % threshold, 'w')
        f.write('\n'.join(res_map))
        f.close()

if __name__ == '__main__':
    main()
