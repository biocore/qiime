#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger","Jens Reeder","William Walters"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

"""Contains code for OTU picking, using several techniques.

This module has the responsibility for taking a set of sequences and
grouping those sequences by similarity.
"""

from copy import copy
from itertools import ifilter
from os.path import splitext, split
from os import makedirs
from itertools import imap
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.mothur import parse_otu_list as mothur_parse
from cogent.app.util import get_tmp_filename
from cogent.app.cd_hit import cdhit_clusters_from_seqs
from cogent.app.dotur import dotur_from_alignment
from cogent.app.mothur import Mothur
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.app.blast import blast_seqs, Blastall, BlastResult
from qiime.pycogent_backports.uclust import get_clusters_from_fasta_filepath
from cogent.core.sequence import DnaSequence
from cogent.util.misc import remove_files
from cogent import LoadSeqs, DNA, Alignment
from cogent.util.trie import build_prefix_map
from cogent.util.misc import flatten
from qiime.util import FunctionWithParams, sort_fasta_by_abundance
from qiime.parse import fields_to_dict

class OtuPicker(FunctionWithParams):
    """An OtuPicker dereplicates a set of sequences at a given similarity.

    This is an abstract class: subclasses should implement the __call__
    method.
    """
    
    Name = 'OtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        cdhit) params, so leaving it as a dict rather than setting
        attributes. Some standard entries in params are:

        Algorithm: algorithm used (e.g. nearest-neighbor, furthest-neighbor)
        Similarity: similarity threshold, e.g. 0.97
        Application: 3rd-party application used, if any, e.g. cdhit
        """
        self.Params = params

    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "OtuPicker is an abstract class"
        

    def _prefilter_exact_prefixes(self,seqs,prefix_length=100):
        """
        """
        unique_prefixes = {}
        for seq_id,seq in seqs:
            seq_len = len(seq)
            seq_id = seq_id.split()[0]
            current_prefix = seq[:prefix_length]
            try:
                prefix_data = unique_prefixes[current_prefix]
                if seq_len > prefix_data[2]:
                    # if this is the longest seq with this prefix so far,
                    # update the list of seq_ids, the best seq_len, and the
                    # best hit seq_id
                    prefix_data[0].append(seq_id)
                    prefix_data[1] = seq_id
                    prefix_data[2] = seq_len
                    prefix_data[3] = seq
                else:
                    # if longer have been seen, only update the list of seq_ids
                    prefix_data[0].append(seq_id)
            except KeyError:
                # list of seq_ids mapped to this prefix, best hit seq_id, best hit seq_len
                unique_prefixes[current_prefix] = [[seq_id],seq_id,seq_len,seq]

        # construct the result objects
        filtered_seqs = []
        seq_id_map = {}
        for data in unique_prefixes.values():
            filtered_seqs.append((data[1],data[3]))
            seq_id_map[data[1]] = data[0]

        return filtered_seqs, seq_id_map


    def _prefilter_with_trie(self, seq_path):

        trunc_id = lambda (a,b): (a.split()[0],b)
        # get the prefix map
        mapping=build_prefix_map(imap(trunc_id, MinimalFastaParser(open(seq_path))))
        for key in mapping.keys():
                mapping[key].append(key)

        # collect the representative seqs
        filtered_seqs=[]
        for (label,seq) in MinimalFastaParser(open(seq_path)):
            label=label.split()[0]
            if label in mapping:
                filtered_seqs.append((label,seq))
        return filtered_seqs, mapping
        
    def _map_filtered_clusters_to_full_clusters(self,clusters,filter_map):
        """
        """
        results = []
        for cluster in clusters:
            full_cluster = []
            for seq_id in cluster:
                full_cluster += filter_map[seq_id]
            results.append(full_cluster)
        return results
        
class BlastOtuPicker(OtuPicker):
    """Blast-based OTU picker: clusters sequence by their 'best' blast hit.
    
        The 'best blast hit' for a sequence is defined as the database 
         sequence which achieves the longest alignment with percent sequence
         identity greater than or equal to the OTU similarity threshold
         (default in Params['Similarity'] = 0.97). Database hits must have an
         e-value threshold less than or equal to the max_e_value threshold 
         (default in Params['max_e_value'] as 1e-10).
    """
    
    def __init__(self, params):
        """Return new BlastOtuPicker object with specified params.
        
        """
        _params = {'max_e_value':1e-10,\
                   'seqs_per_blast_run':1000,\
                   'Similarity':0.97}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__(self,seq_path,result_path=None,log_path=None,
        blast_db=None,refseqs_fp=None):
        
        self.log_lines = []
        
        if not blast_db:
            self.blast_db, self.db_files_to_remove = \
                build_blast_db_from_fasta_path(refseqs_fp)
            self.log_lines.append('Reference seqs fp (to build blast db): %s'%\
             refseqs_fp)
        else:
            self.blast_db = blast_db
            self.db_files_to_remove = []
             
        self.log_lines.append('Blast database: %s' % self.blast_db)
        
        clusters, failures = self._cluster_seqs(\
         MinimalFastaParser(open(seq_path)))
        self.log_lines.append('Num OTUs: %d' % len(clusters))
        
        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for cluster_id,cluster in clusters.items():
                of.write('%s\t%s\n' % (cluster_id,'\t'.join(cluster)))
            of.close()
            result = None
            self.log_lines.append('Result path: %s\n' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = clusters
            self.log_lines.append('Result path: None, returned as dict.')
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            self.log_lines = [str(self)] + self.log_lines
            log_file.write('\n'.join(self.log_lines))
            failures.sort()
            log_file.write('Num failures: %d\n' % len(failures))
            log_file.write('Failures: %s\n' % '\t'.join(failures))
    
        remove_files(self.db_files_to_remove,error_on_missing=False)
        # return the result (note this is None if the data was
        # written to file)
        return result
        
    def _cluster_seqs(self,seqs):
        """
        """
        # blast seqs seq_per_blast_run at a time
        # Build object to keep track of the current set of sequences to be
        # blasted, and the results (i.e., seq_id -> (taxonomy,quaility score) 
        # mapping)
        seqs_per_blast_run = self.Params['seqs_per_blast_run']
        current_seqs = []
        result = {}
        failures = []
        
        # Iterate over the (seq_id, seq) pairs
        for seq_id, seq in seqs:
            # append the current seq_id,seq to list of seqs to be blasted
            current_seqs.append((seq_id,seq))
            # When there are self.SeqsPerBlastRun in the list, blast them
            if len(current_seqs) == seqs_per_blast_run:
                # update the result object
                current_clusters, current_failures =\
                 self._blast_seqs(current_seqs)
                result = self._update_cluster_map(result,current_clusters)
                failures += current_failures
                # reset the list of seqs to be blasted
                current_seqs = []
        # Cluster the remaining sequences
        current_clusters, current_failures = self._blast_seqs(current_seqs)
        result = self._update_cluster_map(result,current_clusters)
        failures += current_failures
        return result, failures
         
    def _update_cluster_map(self,cluster_map,new_clusters):
        for cluster_id, seq_ids in new_clusters.items():
            try:
                cluster_map[cluster_id] += seq_ids
            except KeyError:
                cluster_map[cluster_id] = seq_ids
        return cluster_map
        
    def _blast_seqs(self,seqs):
        """
        """
        result = {}
        failures = []
        if not seqs: 
            return result, failures
        # Get the blast hits with e-values less than self.Params['max_e_value']
        # and percent identity greater than self.Params['Similarity']
        blast_hits = get_blast_hits(seqs,self.blast_db,\
         max_e_value=self.Params['max_e_value'],\
         min_pct_identity=self.Params['Similarity'])
        # Choose the longest alignment out of the acceptable blast hits -- 
        # the result will therefore be the blast hit with at least
        # self.Params['Similarity'] percent identity to the input sequence
        seq_id_to_best_blast_hit = \
         self._choose_longest_blast_hit(blast_hits)
        for seq_id, blast_hit in seq_id_to_best_blast_hit.items():
            if blast_hit == None:
                failures.append(seq_id)
            else:
                cluster_id = blast_hit['SUBJECT ID']
                try:
                    result[cluster_id].append(seq_id)
                except KeyError:
                    result[cluster_id] = [seq_id]
        return result, failures
    
    def _choose_longest_blast_hit(self,blast_hits):
        """ choose the longest blast match 
            
            This function assumes that the blast_hits below 
             self.Params['Similarity'] have already been filtered out, 
             and therefore the longest alignment is the best blast pick.
        """
        result = {}
        # iterate over the queries and their acceptable blast hits
        for query,blast_hits in blast_hits.items():
            choice = None
            len_longest = 0
            # iterate over the acceptable blast hits
            for blast_hit in blast_hits:
                # if the alignment is the longest we've seen so far (or 
                # the first), hold on to it as a possible best hit
                len_current = blast_hit['ALIGNMENT LENGTH']
                if len_current >= len_longest:
                    choice = blast_hit
                    len_longest = len_current 
            query = query.split()[0]    #get rid of spaces
            result[query] = choice
        return result
  
## START MOVE TO BLAST APP CONTROLLER
## The following two functions should be move to the blast application
## controller. When that's done, qiime.assign_taxonomy needs to be updated
## to use these functions rather that the member functions which these 
## are replicas of. Note that when moving to the blast app controller,
## tests should be extractable from test_assign_taxonomy.py.

# THIS FUNCTION SHOULD DO THE SeqsPerBlastRun splitting, would be _much_
# cleaner that way. 
def get_blast_hits(seqs,blast_db,max_e_value=1e-10,min_pct_identity=0.75):
    """ blast each seq in seqs against blast_db and retain good hits
    """
    max_evalue = max_e_value
    min_percent_identity = min_pct_identity
    seq_ids = [s[0] for s in seqs]
    result = {}
    
    blast_result = blast_seqs(\
     seqs,Blastall,blast_db=blast_db,\
     params={'-p':'blastn','-n':'F'},\
     add_seq_names=False)
     
    if blast_result['StdOut']:
        lines = [x for x in blast_result['StdOut']]
        blast_result = BlastResult(lines)
    else:
        return {}.fromkeys(seq_ids,[])
        
    for seq_id in seq_ids:
        blast_result_id = seq_id.split()[0]
        result[seq_id] = []
        if blast_result_id in blast_result:
            for e in blast_result[blast_result_id][0]:
                if (float(e['E-VALUE']) <= max_evalue and \
                    float(e['% IDENTITY']) / 100. >= min_percent_identity):
                    result[seq_id].append(e)

    return result
## END MOVE TO BLAST APP CONTROLLER


class PrefixSuffixOtuPicker(OtuPicker):
    
    Name = 'PrefixSuffixOtuPicker'
    
    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        params contains both generic and per-method (e.g. for
        cdhit application controller) params.
        
        Some generic entries in params are:
    
        Algorithm: algorithm used
        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '-c' parameter
         to the cd-hit application controllers)
        Application: 3rd-party application used
        """
        _params = {'Similarity':0.97,\
         'Algorithm':'Prefix/suffix exact matching'}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__ (self, seq_path, result_path=None, log_path=None, 
        prefix_length=50,suffix_length=50):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.
        prefix_prefilter_length: prefilters the sequence collection so 
         sequences whose first prefix_prefilter_length characters are 
         identical will automatically be grouped into the same OTU [off by 
         default, 100 is typically a good value if this filtering is 
         desired] -- useful for large sequence collections, when cdhit doesn't
         scale well

        """
        log_lines = []
        log_lines.append('Prefix length: %d' % suffix_length)
        log_lines.append('Suffix length: %d' % suffix_length)
        
        assert prefix_length >= 0, 'Prefix length (%d) must be >= 0' % prefix_length
        assert suffix_length >= 0, 'Suffix length (%d) must be >= 0' % suffix_length

        clusters = self._collapse_exact_matches(\
         MinimalFastaParser(open(seq_path)),prefix_length,suffix_length)
        log_lines.append('Num OTUs: %d' % len(clusters))
        
        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_lines.append('Result path: %s' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_lines.append('Result path: None, returned as dict.')
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_lines = [str(self)] + log_lines
            log_file.write('\n'.join(log_lines))
    
        # return the result (note this is None if the data was
        # written to file)
        return result
    
    def _build_seq_hash(self,seq,prefix_length,suffix_length): 
        """ Merge the prefix and suffix into a hash for the OTU
        """
        len_seq = len(seq)
        
        if len_seq <= prefix_length + suffix_length:
            return seq
            
        prefix = seq[:prefix_length]
        suffix = seq[len_seq-suffix_length:]
            
        return prefix + suffix
    
    def _collapse_exact_matches(self,seqs,prefix_length,suffix_length):
        """ Cluster sequences into sets with identical prefix/suffix
        """
        cluster_map = {}
        for seq_id, seq in seqs:
            seq_id = seq_id.split()[0]
            seq_hash = self._build_seq_hash(seq,prefix_length,suffix_length)
            try:
                cluster_map[seq_hash].append(seq_id)
            except KeyError:
                cluster_map[seq_hash] = [seq_id]
        
        return cluster_map.values()


class TrieOtuPicker(OtuPicker):
    
    Name = 'TrieOtuPicker'
    
    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        params contains both generic and per-method (e.g. for
        cdhit application controller) params.
        
        Some generic entries in params are:
    
        Algorithm: algorithm used
        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '-c' parameter
         to the cd-hit application controllers)
        Application: 3rd-party application used
        """
        _params = {'Similarity':0.97,\
         'Algorithm':'Trie prefix or suffix matching',\
         'Reverse':False}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.

        """
        log_lines = []
        
        # Get the appropriate sequence iterator
        if self.Params['Reverse']:
            # Reverse the sequences prior to building the prefix map. 
            # This effectively creates a suffix map.
            # Also removes descriptions from seq identifier lines
            seqs = imap(lambda s: (s[0].split()[0], s[1][::-1]),\
                        MinimalFastaParser(open(seq_path)))
            log_lines.append(\
             'Seqs reversed for suffix mapping (rather than prefix mapping).')
        else:
            # remove descriptions from seq identifier lines
            seqs = imap(lambda s: (s[0].split()[0], s[1]),\
                        MinimalFastaParser(open(seq_path)))
        
        # Build the mapping
        mapping=build_prefix_map(seqs)
        log_lines.append('Num OTUs: %d' % len(mapping))
        
        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,(otu_id,members) in enumerate(mapping.iteritems()):
                of.write('%s\t%s\n' % (i,'\t'.join([otu_id] + members)))
            of.close()
            result = None
            log_lines.append('Result path: %s' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            #add key to cluster_members
            for key in mapping.keys():
                mapping[key].append(key)
            result = dict(enumerate(mapping.values()))
            log_lines.append('Result path: None, returned as dict.')
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_lines = [str(self)] + log_lines
            log_file.write('\n'.join(log_lines))
    
        # return the result (note this is None if the data was
        # written to file)
        return result


class CdHitOtuPicker(OtuPicker):
    
    Name = 'CdHitOtuPicker'
    
    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        params contains both generic and per-method (e.g. for
        cdhit application controller) params.
        
        Some generic entries in params are:
    
        Algorithm: algorithm used
        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '-c' parameter
         to the cd-hit application controllers)
        Application: 3rd-party application used
        """
        _params = {'Similarity':0.97,\
         'Application':'cdhit',\
         'Algorithm':'cdhit: "longest-sequence-first list removal algorithm"'}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__ (self, seq_path, result_path=None, log_path=None, 
        id_len=0, prefix_prefilter_length=None, trie_prefilter=False):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.
        id_len: if set, truncates ids to n chars (you don't want this!)
        prefix_prefilter_length: prefilters the sequence collection so 
         sequences whose first prefix_prefilter_length characters are 
         identical will automatically be grouped into the same OTU [off by 
         default, 100 is typically a good value if this filtering is 
         desired] -- useful for large sequence collections, when cdhit doesn't
         scale well
        trie_prefilter: prefilter the sequence collection such that all sequences
         which are a prefix of another sequence are clustered with the other sequence.
         Togther with cd-hit this is a non-heuristic filter reduces run time a lot.
         Still a bit slower than the prefix_prefilter toggled with prefix_prefilter_length.
        """
        moltype = DNA
        log_lines = []
        
        # create the params dict to pass to cd-hit-est -- IS THERE A
        # BETTER WAY TO MAKE self.Params INTO THE params DICT TO PASS
        # TO THE APPLICATION CONTROLLERS?
        cd_hit_params = copy(self.Params)
        del cd_hit_params['Application']
        del cd_hit_params['Algorithm']
        cd_hit_params['-d'] = id_len  #turn off id truncation
        cd_hit_params['-g'] = "1"
        if (prefix_prefilter_length!=None and trie_prefilter):
            log_lines.append("Both prefilters selected. Deactivate \
            trie_prefilter")
            trie_prefilter=False

        if prefix_prefilter_length != None:
            log_lines.append(\
             'Prefix-based prefiltering, prefix length: %d' \
             % prefix_prefilter_length )
            seqs, filter_map = self._prefilter_exact_prefixes(\
              MinimalFastaParser(open(seq_path)),prefix_prefilter_length)
            log_lines.append(\
             'Prefix-based prefiltering, post-filter num seqs: %d' \
             % len(seqs))
            
        elif trie_prefilter:
            log_lines.append(\
                         'Trie-based prefiltering')
            seqs, filter_map = self._prefilter_with_trie(seq_path)

            log_lines.append(\
                         'Trie-based prefiltering, post-filter num seqs: %d' \
                         % len(seqs))
            
        else:
            log_lines.append('No prefix-based prefiltering.')
            # Load the seq path. Right now, cdhit_clusters_from_seqs
            # doesn't support being passed a file path even though the 
            # seqs do get written to a fasta file before being passed
            # to cd-hit-est. We may want to change that in the future 
            # to avoid the overhead of loading large sequence collections
            # during this step. 
            seqs = LoadSeqs(seq_path,moltype=moltype,aligned=False)
        
        
        # Get the clusters by running cd-hit-est against the
        # sequence collection
        clusters = cdhit_clusters_from_seqs(\
         seqs=seqs,moltype=moltype,params=cd_hit_params)
        
        if prefix_prefilter_length != None or trie_prefilter:
            clusters = self._map_filtered_clusters_to_full_clusters(\
             clusters,filter_map)
        
        
        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_lines.append('Result path: %s' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_lines.append('Result path: None, returned as dict.')
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_lines = [str(self)] + log_lines
            log_file.write('\n'.join(log_lines))
    
        # return the result (note this is None if the data was
        # written to file)
        return result


class UclustOtuPickerBase(OtuPicker):
    
    def _presort_by_abundance(self,seq_path):
        """ Preform pre-sorting of input by abundance """
        
        # Turn off uclust's sorting
        self.Params['suppress_sort'] = True
        
        # Get a temp file name for the sorted fasta file
        sorted_input_seqs_filepath = \
         get_tmp_filename(prefix=self.Name,suffix='.fasta')
        
        # Sort input seqs by abundance, and write to the temp
        # file
        sort_fasta_by_abundance(open(seq_path,'U'),
         open(sorted_input_seqs_filepath,'w'))
        
        # Return the sorted sequences filepath
        return sorted_input_seqs_filepath
    
    def _write_log(self,log_path,log_lines):
        # if the user provided a log file path, log the run
        log_file = open(log_path,'w')
        log_file.write('\n'.join([str(self)] + log_lines))
        log_file.close()
    
    def _prepare_results(self,result_path,clusters,log_lines):
        """
        """
        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for cluster_id,cluster in clusters:
                of.write('%s\t%s\n' % (cluster_id,'\t'.join(cluster)))
            of.close()
            result = None
            log_lines.append('Result path: %s' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(clusters)
            log_lines.append('Result path: None, returned as dict.')
            
        return result

class UclustOtuPicker(UclustOtuPickerBase):
    """ Uclust based OTU picker

    Important note - the default behaviour of uclust is to ignore
    sequences of 32 nucleotides or less.  These will be omitted
    in the clusters generated. """

    Name = 'UclustOtuPicker'
    
    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        params contains both generic and per-method (e.g. for
        uclust application controller) params.
        
        Some generic entries in params are:
    
        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '--id' parameter
         to the uclust application controllers)
        Application: 3rd-party application used
        """
        _params = {'Similarity':0.97,
         'Application':'uclust',
         'max_accepts':8,
         'max_rejects':32,
         'enable_rev_strand_matching':False,
         'optimal':False,
         'exact':False,
         'suppress_sort':True,
         'presort_by_abundance':True}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__(self,
                 seq_path,
                 result_path=None,
                 log_path=None,
                 HALT_EXEC=False):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.

        """
        if self.Params['presort_by_abundance']:
            # seq path will become the temporary sorted sequences
            # filepath, to be cleaned up after the run
            seq_path = self._presort_by_abundance(seq_path)
            files_to_remove = [seq_path]
        else:
            # create a dummy list of files to clean up
            files_to_remove = []
        
        # perform the clustering
        clusters, failures, seeds = get_clusters_from_fasta_filepath(
         seq_path,
         percent_ID = self.Params['Similarity'],
         optimal = self.Params['optimal'],
         exact = self.Params['exact'],
         suppress_sort = self.Params['suppress_sort'],
         enable_rev_strand_matching =
          self.Params['enable_rev_strand_matching'],
         max_accepts=self.Params['max_accepts'],
         max_rejects=self.Params['max_rejects'],
         HALT_EXEC=HALT_EXEC)
        
        # clean up any temp files that were created
        remove_files(files_to_remove)
        
        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))
        
        clusters = enumerate(clusters)
        result = self._prepare_results(result_path,clusters,log_lines)
        
        if log_path:
            self._write_log(log_path,log_lines)
    
        # return the result (note this is None if the data was
        # written to file)
        return result


class UclustReferenceOtuPicker(UclustOtuPickerBase):
    """Uclust reference OTU picker: clusters seqs by match to ref collection
    
    """
    
    def __init__(self, params):
        """Return new UclustReferenceOtuPicker object with specified params.
        
        """
        _params = {'Similarity':0.97,
                   'Application':'uclust',
                   'enable_rev_strand_matching':True,
                   'max_accepts':8,
                   'max_rejects':32,
                   'suppress_new_clusters':False,
                   'optimal':False,
                   'exact':False,
                   'suppress_sort':False,
                   'new_cluster_identifier':'qiime_otu_',
                   'next_new_cluster_number':1,
                   'presort_by_abundance':True}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__(self,
                 seq_fp,
                 refseqs_fp,
                 next_new_cluster_number=None,
                 new_cluster_identifier=None,
                 result_path=None,
                 log_path=None,
                 HALT_EXEC=False):
        
        if new_cluster_identifier:
            self.Params['new_cluster_identifier'] = new_cluster_identifier
        if next_new_cluster_number != None:
            self.Params['next_new_cluster_number'] = next_new_cluster_number
            
        if self.Params['presort_by_abundance']:
            # seq path will become the temporary sorted sequences
            # filepath, to be cleaned up after the run
            seq_fp = self._presort_by_abundance(seq_fp)
            files_to_remove = [seq_fp]
        else:
            # create a dummy list of files to clean up
            files_to_remove = []
        
        # perform the clustering
        cluster_map, failures, new_seeds = get_clusters_from_fasta_filepath(
            seq_fp,
            subject_fasta_filepath=refseqs_fp,
            percent_ID=self.Params['Similarity'],
            enable_rev_strand_matching=self.Params['enable_rev_strand_matching'],
            max_accepts=self.Params['max_accepts'],
            max_rejects=self.Params['max_rejects'],
            suppress_new_clusters=self.Params['suppress_new_clusters'],
            optimal=self.Params['optimal'],
            exact=self.Params['exact'],
            suppress_sort=self.Params['suppress_sort'],
            return_cluster_maps=True,
            HALT_EXEC=HALT_EXEC)
        
        self._rename_clusters(cluster_map,new_seeds)
        
        # clean up any temp files that were created
        remove_files(files_to_remove)
        
        log_lines = []
        log_lines.append('Reference seqs:%s' % refseqs_fp)
        log_lines.append('Num OTUs:%d' % len(cluster_map))
        log_lines.append('Num new OTUs:%d' % len(new_seeds))
        log_lines.append('Num failures:%d' % len(failures))
        log_lines.append('Failures:%s' % '\t'.join(failures))
        
        cluster_map = cluster_map.items()
        result = self._prepare_results(result_path,cluster_map,log_lines)
 
        if log_path:
            self._write_log(log_path,log_lines)
    
        # return the result (note this is None if the data was
        # written to file)
        return result
    
    def _rename_clusters(self,cluster_map,new_seeds):
        """ """
        next_new_cluster_number = self.Params['next_new_cluster_number']
        new_cluster_identifier = self.Params['new_cluster_identifier']
        new_seed_lookup = {}.fromkeys(new_seeds)
        
        for seed,cluster in cluster_map.items():
            del cluster_map[seed]
            if seed in new_seed_lookup:
                new_cluster_id = '%s%d' % (new_cluster_identifier, 
                                           next_new_cluster_number)
                next_new_cluster_number += 1
            else:
                new_cluster_id = seed.split()[0]
                
            cluster_map[new_cluster_id] = cluster
        
        self.Params['next_new_cluster_number'] = next_new_cluster_number

class DoturOtuPicker(OtuPicker):
    Name = 'DoturOtuPicker'

    def __init__(self, params):
        """Return new DoturOtuPicker object with specified params.

        The Behavior of a DoturOtuPicker is controlled by the
        following parameters:

        Similarity: similarity threshold corresponding to genus-level
        OTUs (default 0.97)

        Distance Function: a function which takes two sequences and
        returns a distance ranging between 0 and 1 (default
        DnaSequence.fracDiff())

        Moltype: type of molecule (default: DNA)
        """
        params['Application'] = 'dotur'
        params['Algorithm'] = 'Furthest-neighbor clustering algorithm'

        if 'Similarity' not in params:
            params['Similarity'] = 0.97
        if 'Distance Function' not in params:
            def default_distance_function(first, second):
                first = DnaSequence(first)
                return first.fracDiff(second)
            params['Distance Function'] = default_distance_function 
        if 'Moltype' not in params:
            params['Moltype'] = DNA

        OtuPicker.__init__(self, params)

    def __call__(self, seq_path, result_path=None, log_path=None,\
        prefix_prefilter_length=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.
        """
        if prefix_prefilter_length != None:
            raise NotImplementedError,\
             "DOTUR otu-picking does not currently support prefix-based pre-filtering."
        moltype = self.Params['Moltype']
        distance_function = self.Params['Distance Function']
        similarity_threshold = self.Params['Similarity']

        seqs = LoadSeqs(seq_path, moltype=moltype, aligned=True)
        parsed_results = dotur_from_alignment(seqs, moltype, distance_function)
        clusters = self.__pick_clusters(parsed_results, similarity_threshold)

        # From here down, this is all copied straight from
        # CdHitOtuPicker, and prime for refactoring into a private
        # method of OtuPicker

        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_str = 'Result path: %s' % result_path
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_str = 'Result path: None, returned as dict.'
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_file.write(str(self))
            log_file.write('\n')
            log_file.write('%s\n' % log_str)
    
        # return the result (note this is None if the data was
        # written to file)
        return result

    @staticmethod
    def __pick_clusters(dotur_results, similarity_threshold):
        """Returns OTU's that satisfy the given similarity threshold.

        A higher similarity threshold results in a larger number of
        distinct OTU's.

        Since DOTUR returns all sets of OTU's at various degrees of
        similarity, this process is sufficiently involved to warrant
        its own static method.  This method should probably live with
        the DOTUR application controller in PyCogent, and has been set
        up for easy refactoring.
        """
        # In dotur, a lower similarity score means more otu's.  To
        # find otu's that satisfy a similarity threshold of 0.9, we
        # must find the largest dotur score less than (1 - 0.9 =) 0.1.
        dotur_score_threshold = 1 - similarity_threshold

        # info stored in various positions of the results list
        i_score = 0
        i_otu_list = 2

        # 'unique' in first result set indicates score of 0.0
        # change to numerical value
        dotur_results[0][i_score] = 0.0

        prev_res = dotur_results[0]
        for res in dotur_results:
            score = res[i_score]
            if score <= dotur_score_threshold:
                prev_res = res
            else:
                break
        
        # return only the list of otu's
        retval = prev_res[i_otu_list]

        return retval


class MothurOtuPicker(OtuPicker):
    Name = 'MothurOtuPicker'
    ClusteringAlgorithms = ['furthest', 'nearest', 'average']

    def __init__(self, params):
        """Return new MothurOtuPicker object with specified params.
        
        Valid params are:

        Algorithm
            Algorithm used for clustering (valid choices are nearest,
            furthest, average)
        Similarity
            Similarity threshold for OTUs (default 0.97)
        """
        params['Application'] = 'mothur'
        if 'Algorithm' not in params:
            params['Algorithm'] = 'furthest'
        if 'Similarity' not in params:
            params['Similarity'] = 0.97
        if params['Algorithm'] not in self.ClusteringAlgorithms:
            raise ValueError('Unsupported algorithm %s.  Choices are %s' % \
                             (params['Algorithm'], self.ClusteringAlgorithms))
        super(MothurOtuPicker, self).__init__(params)

    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        app = Mothur(InputHandler='_input_as_path')
        app.Parameters['method'].on(self.Params['Algorithm'])
        results = app(seq_path)
        parsed_otus = mothur_parse(results['otu list'])
        clusters = self.__pick_clusters(parsed_otus)
        results.cleanUp()

        # From here down, this is all copied straight from
        # CdHitOtuPicker, and prime for refactoring into a private
        # method of OtuPicker

        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_str = 'Result path: %s' % result_path
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_str = 'Result path: None, returned as dict.'
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_file.write(str(self))
            log_file.write('\n')
            log_file.write('%s\n' % log_str)

        # return the result (note this is None if the data was
        # written to file)
        return result

    def __pick_clusters(self, mothur_results):
        """Returns OTU's that satisfy the given similarity threshold.
        """
        # Sanity check
        if not 0 <= self.Params['Similarity'] <= 1:
            raise ValueError(
                'Similarity threshold must be number between 0 and 1 '
                '(received %)' % similarity_threshold)

        # A lower mothur score means more otu's.  To find otu's that
        # satisfy a similarity threshold of 0.9, we must find the
        # largest score less than or equal to (1 - 0.9 =) 0.1.
        score_threshold = 1 - self.Params['Similarity']

        my_score, my_otus = mothur_results.next()
        for score, otus in mothur_results:

            # Sanity check
            if score < my_score:
                raise ValueError(
                    'Mothur results not in ascending order.  This is an error '
                    'in the Mothur application controller, and it should be '
                    'reported to the PyCogent developers.')

            if score <= score_threshold:
                my_score, my_otus = score, otus
            else:
                # Scores are only getting larger, so bail out now
                break
        return my_otus

# Some functions to support merging OTU tables
# generated one after another. This functionality is currently available
# via Qiime/scripts/merge_otu_maps.py and will be incorporated into the
# MetaPickOtus or ChainedPickOtus class when that comes into existence.

def expand_otu_map_seq_ids(otu_map,seq_id_map):
    for otu_id, seq_ids in otu_map.items():
        mapped_seq_ids = flatten(\
         [seq_id_map[seq_id] for seq_id in seq_ids])
        otu_map[otu_id] = mapped_seq_ids
    return otu_map
                
def map_otu_map_files(otu_files):
    # passing delim=None splits on any whitespace, so can handle mixed tabs
    # and spaces
    result = fields_to_dict(otu_files[0],delim=None)
    for otu_file in otu_files[1:]:
        current_otu_map = fields_to_dict(otu_file,delim=None)
        result = expand_otu_map_seq_ids(current_otu_map,result)
    return result

def write_otu_map(otu_map,output_fp):
    """
    """
    of = open(output_fp,'w')
    for otu_id in sorted(otu_map.keys()):
        of.write('\t'.join([otu_id] + otu_map[otu_id]))
        of.write('\n')
    of.close()
# End functions to support merging OTU tables


otu_picking_method_constructors = {
    'cdhit': CdHitOtuPicker,
    'prefix_suffix': PrefixSuffixOtuPicker,
    'mothur': MothurOtuPicker,
    'trie':TrieOtuPicker,
    'blast':BlastOtuPicker,
    'uclust': UclustOtuPicker,
    'uclust_ref':UclustReferenceOtuPicker
    }
    
otu_picking_method_choices = otu_picking_method_constructors.keys()


        

