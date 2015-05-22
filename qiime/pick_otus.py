#!/usr/bin/env python
__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Jens Reeder",
               "William Walters", "Jose Carlos Clemente Litran",
               "Adam Robbins-Pianka", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

"""Contains code for OTU picking, using several techniques.

This module has the responsibility for taking a set of sequences and
grouping those sequences by similarity.
"""

from copy import copy
from itertools import ifilter
from os.path import splitext, split, abspath, join, dirname
from os import makedirs, close, rename
from itertools import imap
from tempfile import mkstemp

from bfillings.mothur import parse_otu_list as mothur_parse

from skbio.util import remove_files, flatten
from skbio.tree import CompressedTrie, fasta_to_pairlist
from skbio.parse.sequences import parse_fasta
from skbio.alignment import SequenceCollection
from skbio.sequence import DNA

from qiime.util import FunctionWithParams, get_qiime_temp_dir
from qiime.sort import sort_fasta_by_abundance
from qiime.parse import fields_to_dict

from bfillings.blast import blast_seqs, Blastall, BlastResult
from bfillings.formatdb import build_blast_db_from_fasta_path
from bfillings.mothur import Mothur
from bfillings.cd_hit import cdhit_clusters_from_seqs
from bfillings.uclust import get_clusters_from_fasta_filepath
from bfillings.sortmerna_v2 import (build_database_sortmerna,
                                    sortmerna_ref_cluster)
from bfillings.usearch import (usearch_qf,
                            usearch61_denovo_cluster,
                            usearch61_ref_cluster)
from bfillings.sumaclust_v1 import sumaclust_denovo_cluster
from bfillings.swarm_v127 import swarm_denovo_cluster


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

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError("OtuPicker is an abstract class")

    def _prefilter_exact_prefixes(self, seqs, prefix_length=100):
        """
        """
        unique_prefixes = {}
        for seq_id, seq in seqs:
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
                # list of seq_ids mapped to this prefix, best hit seq_id, best
                # hit seq_len
                unique_prefixes[current_prefix] = [[seq_id],
                                                   seq_id,
                                                   seq_len,
                                                   seq]

        # construct the result objects
        filtered_seqs = []
        seq_id_map = {}
        for data in unique_prefixes.values():
            filtered_seqs.append((data[1], data[3]))
            seq_id_map[data[1]] = data[0]

        return filtered_seqs, seq_id_map

    def _prefilter_exact_matches(self, seqs):
        """
        """
        unique_sequences = {}
        seq_id_map = {}
        filtered_seqs = []
        for seq_id, seq in seqs:
            seq_id = seq_id.split()[0]
            try:
                temp_seq_id = unique_sequences[seq]
            except KeyError:
                # unseen sequence so create a new temp_seq_id,
                # a new unique_sequence entry, and new seq_id_map
                # entry, and add the sequence to the list of
                # filtered seqs -- this will retain the order
                # of the input sequences too
                temp_seq_id = 'QiimeExactMatch.%s' % seq_id
                unique_sequences[seq] = temp_seq_id
                seq_id_map[temp_seq_id] = []
                filtered_seqs.append((temp_seq_id, seq))
            seq_id_map[temp_seq_id].append(seq_id)
        return filtered_seqs, seq_id_map

    def _prefilter_with_trie(self, seq_path):

        trunc_id = lambda a_b: (a_b[0].split()[0], a_b[1])
        # get the prefix map
        with open(seq_path, 'U') as seq_lines:
            t = CompressedTrie(fasta_to_pairlist(imap(trunc_id,
                                                      parse_fasta(seq_lines))))
        mapping = t.prefix_map
        for key in mapping.keys():
                mapping[key].append(key)

        # collect the representative seqs
        filtered_seqs = []
        for (label, seq) in parse_fasta(open(seq_path)):
            label = label.split()[0]
            if label in mapping:
                filtered_seqs.append((label, seq))
        return filtered_seqs, mapping

    def _map_filtered_clusters_to_full_clusters(self, clusters, filter_map):
        """
            Input:  clusters, a list of cluster lists
                    filter_map, the seq_id in each clusters
                                is the key to the filter_map
                                containing all seq_ids with
                                duplicate FASTA sequences
            Output: an extended list of cluster lists
        """
        results = []
        for cluster in clusters:
            full_cluster = []
            for seq_id in cluster:
                full_cluster += filter_map[seq_id]
            results.append(full_cluster)
        return results


class SortmernaV2OtuPicker(OtuPicker):

    """ SortMeRNA-based version 2 OTU picker: clusters queries by their 'best'
        alignment to a reference seed.

        The 'best' alignment for a query is the one with:
            1. the lowest E-value score (at most 1)
            2. percent sequence identity greater than or equal to the OTU
                similarity threshold (default in Params['similarity'] = 0.97)
            3. percent query coverage greater than or equal to the OTU
                coverage threshold (default in Params['coverage'] = 0.97)
    """

    def __init__(self, params):
        """ Return a new SortmernaV2OtuPicker object with specified params.
        """

        OtuPicker.__init__(self, params)

    def __call__(self, seq_path, result_path=None, log_path=None,
                 sortmerna_db=None, refseqs_fp=None, failure_path=None):

        """ Purpose   : Call to construct the reference database (if not provided)
                        and to launch sortmerna.
            Parameters: seq_path, path to reads file;
                        result_path, path to OTU mapping file;
                        log_path, path to QIIME log file;
                        sortmerna_db, path to sortmerna indexed database;
                        refseqs_fp, path to reference sequences;
                        failure_path, path to text file of reads failing
                          to align with similarity & coverage thresholds;
            Return    : None (output is always written to file)
        """

        self.log_lines = []
        prefilter_identical_sequences =\
            self.Params['prefilter_identical_sequences']

        # Indexed database not provided, build it
        if not sortmerna_db:
            # write index to output directory
            self.sortmerna_db, self.files_to_remove = \
                build_database_sortmerna(abspath(refseqs_fp),
                                         max_pos=self.Params['max_pos'],
                                         output_dir=dirname(abspath(result_path)))

            self.log_lines.append('Reference seqs fp (to build '
                                  'sortmerna database): %s' %
                                  abspath(refseqs_fp))
        # Indexed database provided
        else:
            self.sortmerna_db = sortmerna_db
            self.files_to_remove = []

        self.log_lines.append('SortMeRNA database: %s' % self.sortmerna_db)

        original_fasta_path = seq_path

        # Collapse identical sequences to a new file
        if prefilter_identical_sequences:
            exact_match_id_map, seq_path =\
                self._apply_identical_sequences_prefilter(seq_path)

        # Call sortmerna for reference clustering
        cluster_map, failures, smr_files_to_remove =\
            sortmerna_ref_cluster(seq_path=seq_path,
                                  sortmerna_db=self.sortmerna_db,
                                  refseqs_fp=refseqs_fp,
                                  result_path=result_path,
                                  tabular=self.Params['blast'],
                                  max_e_value=self.Params['max_e_value'],
                                  similarity=self.Params['similarity'],
                                  coverage=self.Params['coverage'],
                                  threads=self.Params['threads'],
                                  best=self.Params['best'],
                                  HALT_EXEC=False)

        # Remove temporary files
        self.files_to_remove.extend(smr_files_to_remove)
        remove_files(self.files_to_remove, error_on_missing=False)

        # Expand identical sequences to create full OTU map
        if prefilter_identical_sequences:
            cluster_names = cluster_map.keys()
            clusters = [cluster_map[c] for c in cluster_names]
            clusters =\
                self._map_filtered_clusters_to_full_clusters(
                    clusters, exact_match_id_map)
            cluster_map = dict(zip(cluster_names, clusters))

            # Expand failures
            temp_failures = []
            for fa in failures:
                temp_failures.extend(exact_match_id_map[fa])
            failures = temp_failures

        self.log_lines.append('Num OTUs: %d' % len(cluster_map))
        self.log_lines.append('Num failures: %d' % len(failures))

        # Write failures to file
        if failure_path is not None:
            failure_file = open(failure_path, 'w')
            failure_file.write('\n'.join(failures))
            failure_file.write('\n')
            failure_file.close()

        # Write OTU map
        if result_path:
            # If the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster (this will over-write the default SortMeRNA
            # OTU map with the extended OTU map)
            of = open(result_path, 'w')
            for cluster_id, cluster in cluster_map.items():
                of.write('%s\t%s\n' % (cluster_id, '\t'.join(cluster)))
            of.close()
            result = None
            self.log_lines.append('Result path: %s\n' % result_path)
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = cluster_map
            self.log_lines.append('Result path: None, returned as dict.')

        # Log the run
        if log_path:
            log_file = open(log_path, 'w')
            self.log_lines = [str(self)] + self.log_lines
            log_file.write('\n'.join(self.log_lines))
            log_file.write('\n')

        return result

    def _apply_identical_sequences_prefilter(self, seq_path):
        """
            Input : a filepath to input FASTA reads
            Method: prepares and writes de-replicated reads
                    to a temporary FASTA file, calls
                    parent method to do the actual
                    de-replication
            Return: exact_match_id_map, a dictionary storing
                    de-replicated amplicon ID as key and
                    all original FASTA IDs with identical
                    sequences as values;
                    unique_seqs_fp, filepath to FASTA file
                    holding only de-replicated sequences
        """
        # Creating mapping for de-replicated reads
        with open(seq_path, 'U') as s_path:
            seqs_to_cluster, exact_match_id_map =\
                self._prefilter_exact_matches(parse_fasta(s_path))

        # Create temporary file for storing the de-replicated reads
        fd, unique_seqs_fp = mkstemp(
            prefix='SortMeRNAExactMatchFilter', suffix='.fasta')
        close(fd)

        self.files_to_remove.append(unique_seqs_fp)

        # Write de-replicated reads to file
        unique_seqs_f = open(unique_seqs_fp, 'w')
        for seq_id, seq in seqs_to_cluster:
            unique_seqs_f.write('>%s count=%d;\n%s\n' %
                                (seq_id,
                                 len(exact_match_id_map[seq_id]),
                                 seq))
        unique_seqs_f.close()

        # Clean up the seqs_to_cluster as we don't need
        # it again
        del(seqs_to_cluster)

        return exact_match_id_map, unique_seqs_fp


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
        _params = {'max_e_value': 1e-10,
                   'seqs_per_blast_run': 1000,
                   'Similarity': 0.97,
                   'min_aligned_percent': 0.50,
                   'blast_program': 'blastn',
                   'is_protein': False}
        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None,
                 blast_db=None, refseqs_fp=None):

        self.log_lines = []

        if not blast_db:
            self.blast_db, self.db_files_to_remove = \
                build_blast_db_from_fasta_path(abspath(refseqs_fp),
                                               is_protein=self.Params[
                                                   'is_protein'],
                                               output_dir=get_qiime_temp_dir())
            self.log_lines.append('Reference seqs fp (to build blast db): %s' %
                                  abspath(refseqs_fp))
        else:
            self.blast_db = blast_db
            self.db_files_to_remove = []

        self.log_lines.append('Blast database: %s' % self.blast_db)

        clusters, failures = self._cluster_seqs(parse_fasta(open(seq_path)))
        self.log_lines.append('Num OTUs: %d' % len(clusters))

        if result_path:
            # if the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for cluster_id, cluster in clusters.items():
                of.write('%s\t%s\n' % (cluster_id, '\t'.join(cluster)))
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
            log_file = open(log_path, 'w')
            self.log_lines = [str(self)] + self.log_lines
            log_file.write('\n'.join(self.log_lines))
            failures.sort()
            log_file.write('Num failures: %d\n' % len(failures))
            log_file.write('Failures: %s\n' % '\t'.join(failures))

        remove_files(self.db_files_to_remove, error_on_missing=False)
        # return the result (note this is None if the data was
        # written to file)
        return result

    def _cluster_seqs(self, seqs):
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
            current_seqs.append((seq_id, seq))
            # When there are self.SeqsPerBlastRun in the list, blast them
            if len(current_seqs) == seqs_per_blast_run:
                # update the result object
                current_clusters, current_failures =\
                    self._blast_seqs(current_seqs)
                result = self._update_cluster_map(result, current_clusters)
                failures += current_failures
                # reset the list of seqs to be blasted
                current_seqs = []
        # Cluster the remaining sequences
        current_clusters, current_failures = self._blast_seqs(current_seqs)
        result = self._update_cluster_map(result, current_clusters)
        failures += current_failures
        return result, failures

    def _update_cluster_map(self, cluster_map, new_clusters):
        for cluster_id, seq_ids in new_clusters.items():
            try:
                cluster_map[cluster_id] += seq_ids
            except KeyError:
                cluster_map[cluster_id] = seq_ids
        return cluster_map

    def _blast_seqs(self, seqs):
        """
        """
        result = {}
        failures = []
        if not seqs:
            return result, failures
        # Get the blast hits with e-values less than self.Params['max_e_value']
        # and percent identity greater than self.Params['Similarity']
        blast_hits = get_blast_hits(seqs, self.blast_db,
                                    max_e_value=self.Params['max_e_value'],
                                    min_pct_identity=self.Params['Similarity'],
                                    min_aligned_percent=self.Params[
                                        'min_aligned_percent'],
                                    blast_program=self.Params['blast_program'])
        # Choose the longest alignment out of the acceptable blast hits --
        # the result will therefore be the blast hit with at least
        # self.Params['Similarity'] percent identity to the input sequence
        seq_id_to_best_blast_hit = \
            self._choose_longest_blast_hit(blast_hits)
        for seq_id, blast_hit in seq_id_to_best_blast_hit.items():
            if blast_hit is None:
                failures.append(seq_id)
            else:
                cluster_id = blast_hit['SUBJECT ID']
                try:
                    result[cluster_id].append(seq_id)
                except KeyError:
                    result[cluster_id] = [seq_id]
        return result, failures

    def _choose_longest_blast_hit(self, blast_hits):
        """ choose the longest blast match

            This function assumes that the blast_hits below
             self.Params['Similarity'] have already been filtered out,
             and therefore the longest alignment is the best blast pick.
        """
        result = {}
        # iterate over the queries and their acceptable blast hits
        for query, blast_hits in blast_hits.items():
            choice = None
            len_longest = 0
            # iterate over the acceptable blast hits
            for blast_hit in blast_hits:
                # if the alignment is the longest we've seen so far (or
                # the first), hold on to it as a possible best hit
                len_current = blast_hit['ALIGNMENT LENGTH']
                if len_current > len_longest:
                    choice = blast_hit
                    len_longest = len_current
            query = query.split()[0]  # get rid of spaces
            result[query] = choice
        return result


class BlastxOtuPicker(BlastOtuPicker):

    """Blastx-based OTU picker: clusters sequence by their 'best' blast hit.

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
        _params = {'max_e_value': 1e-3,
                   'seqs_per_blast_run': 1000,
                   'Similarity': 0.75,
                   'min_aligned_percent': 0.50,
                   'blast_program': 'blastx',
                   'is_protein': True}
        _params.update(params)
        OtuPicker.__init__(self, _params)

# START MOVE TO BLAST APP CONTROLLER
# The following two functions should be move to the blast application
# controller. When that's done, qiime.assign_taxonomy needs to be updated
# to use these functions rather that the member functions which these
# are replicas of. Note that when moving to the blast app controller,
# tests should be extractable from test_assign_taxonomy.py.

# THIS FUNCTION SHOULD DO THE SeqsPerBlastRun splitting, would be _much_
# cleaner that way.


def get_blast_hits(seqs,
                   blast_db,
                   max_e_value=1e-10,
                   min_pct_identity=0.75,
                   min_aligned_percent=0.50,
                   blast_program='blastn'):
    """ blast each seq in seqs against blast_db and retain good hits
    """
    max_evalue = max_e_value
    min_percent_identity = min_pct_identity
    seq_ids = [s[0] for s in seqs]
    result = {}

    blast_result = blast_seqs(
        seqs, Blastall, blast_db=blast_db,
        params={'-p': blast_program, '-n': 'F'},
        add_seq_names=False)

    if blast_result['StdOut']:
        lines = [x for x in blast_result['StdOut']]
        blast_result = BlastResult(lines)
    else:
        return {}.fromkeys(seq_ids, [])

    for seq_id, seq in seqs:
        blast_result_id = seq_id.split()[0]
        max_alignment_length = len(seq)
        if blast_program == 'blastx':
            # if this is a translated blast search, the max alignment
            # length is the number of 3mers in seq
            max_alignment_length /= 3
        min_alignment_length = max_alignment_length * min_aligned_percent
        result[seq_id] = []
        if blast_result_id in blast_result:
            for e in blast_result[blast_result_id][0]:
                if (float(e['E-VALUE']) <= max_evalue and
                        float(e['% IDENTITY']) / 100. >= min_percent_identity and
                        int(e['ALIGNMENT LENGTH']) >= min_alignment_length):
                    result[seq_id].append(e)

    return result
# END MOVE TO BLAST APP CONTROLLER


class SumaClustOtuPicker(OtuPicker):

    """ SumaClust is a de novo OTU picker, following the same clustering
        algorithm as Uclust. It is open source and supports multithreading,
        both SIMD and OpenMP.

        Clusters are created by their similarity threshold (default 0.97).
        If a query does not match any seed with this similarity threshold,
        it is used to create a new seed.

        Exact clustering (with parameter -e) assigns queries to the
        best-matching seed, rather than to the first seed with similarity
        threshold.
    """

    def __init__(self, params):
        """ Return a new SumaClustOtuPicker object with specified params.
            The defaults are set in the SumaClust API (see bfillings)
        """

        OtuPicker.__init__(self, params)

    def _apply_identical_sequences_prefilter(self, seq_path):
        """
            Input : a filepath to input FASTA reads
            Method: prepares and writes de-replicated reads
                    to a temporary FASTA file, calls
                    parent method to do the actual
                    de-replication
            Return: exact_match_id_map, a dictionary storing
                    de-replicated amplicon ID as key and
                    all original FASTA IDs with identical
                    sequences as values;
                    unique_seqs_fp, filepath to FASTA file
                    holding only de-replicated sequences
        """
        # creating mapping for de-replicated reads
        seqs_to_cluster, exact_match_id_map =\
            self._prefilter_exact_matches(parse_fasta(open(seq_path, 'U')))

        # create temporary file for storing the de-replicated reads
        fd, unique_seqs_fp = mkstemp(
            prefix='SumaClustExactMatchFilter', suffix='.fasta')
        close(fd)

        self.files_to_remove.append(unique_seqs_fp)

        # write de-replicated reads to file
        unique_seqs_f = open(unique_seqs_fp, 'w')
        for seq_id, seq in seqs_to_cluster:
            unique_seqs_f.write('>%s count=%d;\n%s\n'
                                % (seq_id,
                                   len(exact_match_id_map[seq_id]),
                                   seq))
        unique_seqs_f.close()
        # clean up the seqs_to_cluster list as it can be big and we
        # don't need it again
        del(seqs_to_cluster)

        return exact_match_id_map, unique_seqs_fp

    def __call__(self, seq_path=None, result_path=None, log_path=None):

        self.log_lines = []
        self.files_to_remove = []

        prefilter_identical_sequences =\
            self.Params['prefilter_identical_sequences']

        original_fasta_path = seq_path

        # Collapse idetical sequences to a new file
        if prefilter_identical_sequences:
            exact_match_id_map, seq_path =\
                self._apply_identical_sequences_prefilter(seq_path)

        # Run SumaClust, return a dict of output files
        clusters = sumaclust_denovo_cluster(
            seq_path=abspath(seq_path),
            result_path=abspath(result_path),
            shortest_len=self.Params['l'],
            similarity=self.Params['similarity'],
            threads=self.Params['threads'],
            exact=self.Params['exact'],
            HALT_EXEC=False)

        # Clean up any temp files that were created
        remove_files(self.files_to_remove)

        # Create file for expanded OTU map
        if prefilter_identical_sequences:
            clusters = self._map_filtered_clusters_to_full_clusters(
                clusters, exact_match_id_map)

        self.log_lines.append('Num OTUs: %d' % len(clusters))

        # Add prefix ID to de novo OTUs
        otu_id_prefix = self.Params['denovo_otu_id_prefix']
        if otu_id_prefix is None:
            clusters = dict(enumerate(clusters))
        else:
            clusters = dict(('%s%d' % (otu_id_prefix, i), c)
                            for i, c in enumerate(clusters))

        if result_path:
            # If the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster (this will over-write the default SumaClust
            # OTU map with the extended OTU map)
            of = open(result_path, 'w')
            for cluster_id, cluster in clusters.items():
                of.write('%s\t%s\n' % (cluster_id, '\t'.join(cluster)))
            of.close()
            result = None
            self.log_lines.append('Result path: %s\n' % result_path)
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = clusters
            self.log_lines.append('Result path: None, returned as dict.')

        # Log the run
        if log_path:
            log_file = open(log_path, 'w')
            self.log_lines.insert(0, str(self))
            log_file.write('\n'.join(self.log_lines))
            log_file.close()

        return result


class SwarmOtuPicker(OtuPicker):

    """ Swarm is a de novo OTU picker, an exact clustering method based
        on a single-linkage algorithm. It is open source and supports
        SSE2 multithreading.

        Clusters are created by their local clustering threshold 'd',
        which is computed as the number of nucleotide differences
        (substitution, insertion or deletion) between two amplicons
        in the optimal pairwise global alignment.

        This class is compatible with Swarm v.1.2.7
    """

    def __init__(self, params):
        """ Return a new SwarmOtuPicker object with specified params.
            The defaults are set in the Swarm API
        """

        OtuPicker.__init__(self, params)

    def __call__(self, seq_path=None, result_path=None, log_path=None):

        self.log_lines = []

        # Run Swarm, return a list of lists (clusters)
        clusters = swarm_denovo_cluster(
            seq_path=seq_path,
            d=self.Params['resolution'],
            threads=self.Params['threads'],
            HALT_EXEC=False)

        self.log_lines.append('Num OTUs: %d' % len(clusters))

        # Add prefix ID to de novo OTUs
        otu_id_prefix = self.Params['denovo_otu_id_prefix']
        if otu_id_prefix is None:
            clusters = dict(enumerate(clusters))
        else:
            clusters = dict(('%s%d' % (otu_id_prefix, i), c)
                            for i, c in enumerate(clusters))

        if result_path:
            # If the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for cluster_id, cluster in clusters.items():
                of.write('%s\t%s\n' % (cluster_id, '\t'.join(cluster)))
            of.close()
            result = None
            self.log_lines.append('Result path: %s\n' % result_path)
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = clusters
            self.log_lines.append('Result path: None, returned as dict.')

        # Log the run
        if log_path:
            log_file = open(log_path, 'w')
            self.log_lines.insert(0, str(self))
            log_file.write('\n'.join(self.log_lines))
            log_file.close()

        return result


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
        _params = {'Similarity': 0.97,
                   'Algorithm': 'Prefix/suffix exact matching'}
        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None,
                 prefix_length=50, suffix_length=50):
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
        log_lines.append('Prefix length: %d' % prefix_length)
        log_lines.append('Suffix length: %d' % suffix_length)

        assert prefix_length >= 0, 'Prefix length (%d) must be >= 0' % prefix_length
        assert suffix_length >= 0, 'Suffix length (%d) must be >= 0' % suffix_length

        clusters = self._collapse_exact_matches(parse_fasta(open(seq_path)),
                                                prefix_length, suffix_length)
        log_lines.append('Num OTUs: %d' % len(clusters))

        if result_path:
            # if the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for i, cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i, '\t'.join(cluster)))
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
            log_file = open(log_path, 'w')
            log_lines = [str(self)] + log_lines
            log_file.write('\n'.join(log_lines))

        # return the result (note this is None if the data was
        # written to file)
        return result

    def _build_seq_hash(self, seq, prefix_length, suffix_length):
        """ Merge the prefix and suffix into a hash for the OTU
        """
        len_seq = len(seq)

        if len_seq <= prefix_length + suffix_length:
            return seq

        prefix = seq[:prefix_length]
        suffix = seq[len_seq - suffix_length:]

        return prefix + suffix

    def _collapse_exact_matches(self, seqs, prefix_length, suffix_length):
        """ Cluster sequences into sets with identical prefix/suffix
        """
        cluster_map = {}
        for seq_id, seq in seqs:
            seq_id = seq_id.split()[0]
            seq_hash = self._build_seq_hash(seq, prefix_length, suffix_length)
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
        _params = {'Similarity': 0.97,
                   'Algorithm': 'Trie prefix or suffix matching',
                   'Reverse': False}
        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None):
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
            seqs = imap(lambda s: (s[0].split()[0], s[1][::-1]),
                        parse_fasta(open(seq_path)))
            log_lines.append(
                'Seqs reversed for suffix mapping (rather than prefix mapping).')
        else:
            # remove descriptions from seq identifier lines
            seqs = imap(lambda s: (s[0].split()[0], s[1]),
                        parse_fasta(open(seq_path)))

        # Build the mapping
        t = CompressedTrie(fasta_to_pairlist(seqs))
        mapping = t.prefix_map
        log_lines.append('Num OTUs: %d' % len(mapping))

        if result_path:
            # if the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for i, (otu_id, members) in enumerate(mapping.iteritems()):
                of.write('%s\t%s\n' % (i, '\t'.join([otu_id] + members)))
            of.close()
            result = None
            log_lines.append('Result path: %s' % result_path)
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            # add key to cluster_members
            for key in mapping.keys():
                mapping[key].append(key)
            result = dict(enumerate(mapping.values()))
            log_lines.append('Result path: None, returned as dict.')

        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path, 'w')
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
        _params = {'Similarity': 0.97,
                   'Application': 'cdhit',
                   'Algorithm': 'cdhit: "longest-sequence-first list removal algorithm"'}
        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None,
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
        log_lines = []

        # create the params dict to pass to cd-hit-est -- IS THERE A
        # BETTER WAY TO MAKE self.Params INTO THE params DICT TO PASS
        # TO THE APPLICATION CONTROLLERS?
        cd_hit_params = copy(self.Params)
        del cd_hit_params['Application']
        del cd_hit_params['Algorithm']
        cd_hit_params['-d'] = id_len  # turn off id truncation
        cd_hit_params['-g'] = "1"
        if (prefix_prefilter_length is not None and trie_prefilter):
            log_lines.append("Both prefilters selected. Deactivate trie_prefilter")
            trie_prefilter = False

        if prefix_prefilter_length is not None:
            log_lines.append(
                'Prefix-based prefiltering, prefix length: %d'
                % prefix_prefilter_length)
            with open(seq_path) as seq_f:
                seqs, filter_map = self._prefilter_exact_prefixes(
                    parse_fasta(seq_f, label_to_name=lambda x: x.split()[0]),
                    prefix_prefilter_length)
            log_lines.append(
                'Prefix-based prefiltering, post-filter num seqs: %d' % len(seqs))
        elif trie_prefilter:
            log_lines.append(
                'Trie-based prefiltering')
            seqs, filter_map = self._prefilter_with_trie(seq_path)
            log_lines.append(
                'Trie-based prefiltering, post-filter num seqs: %d' % len(seqs))
        else:
            log_lines.append('No prefix-based prefiltering.')
            # Load the seq path. Right now, cdhit_clusters_from_seqs
            # doesn't support being passed a file path even though the
            # seqs do get written to a fasta file before being passed
            # to cd-hit-est. We may want to change that in the future
            # to avoid the overhead of loading large sequence collections
            # during this step.
            with open(seq_path) as seq_f:
                seqs = SequenceCollection.from_fasta_records(
                    parse_fasta(seq_f, label_to_name=lambda x: x.split()[0]),
                    DNA)
            seqs = dict(seqs.iteritems())

        # Get the clusters by running cd-hit-est against the
        # sequence collection
        clusters = cdhit_clusters_from_seqs(
            seqs=seqs, params=cd_hit_params)
        if prefix_prefilter_length is not None or trie_prefilter:
            clusters = self._map_filtered_clusters_to_full_clusters(
                clusters, filter_map)

        if result_path:
            # if the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for i, cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i, '\t'.join(cluster)))
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
            log_file = open(log_path, 'w')
            log_lines = [str(self)] + log_lines
            log_file.write('\n'.join(log_lines))

        # return the result (note this is None if the data was
        # written to file)
        return result


class UclustOtuPickerBase(OtuPicker):

    def _presort_by_abundance(self, seq_path):
        """ Preform pre-sorting of input by abundance """

        # Turn off uclust's sorting - if doing our presort by
        # abundance we _always_ need to disable uclust's sorting.
        self.Params['suppress_sort'] = True

        # Get a temp file name for the sorted fasta file
        fd, sorted_input_seqs_filepath = \
            mkstemp(prefix=self.Name, suffix='.fasta')
        close(fd)
        # Sort input seqs by abundance, and write to the temp
        # file
        sort_fasta_by_abundance(open(seq_path, 'U'),
                                open(sorted_input_seqs_filepath, 'w'))

        # Return the sorted sequences filepath
        return sorted_input_seqs_filepath

    def _write_log(self, log_path, log_lines):
        # if the user provided a log file path, log the run
        log_file = open(log_path, 'w')
        log_file.write('\n'.join([str(self)] + log_lines))
        log_file.close()

    def _prepare_results(self, result_path, clusters, log_lines):
        """
        """
        if result_path:
            # if the user provided a result_path, write the
            # results to file with one tab-separated line per
            # cluster
            of = open(result_path, 'w')
            for cluster_id, cluster in clusters:
                of.write('%s\t%s\n' % (cluster_id, '\t'.join(cluster)))
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

    def _apply_identical_sequences_prefilter(self, seq_path):
        """ """
        fd, unique_seqs_fp = mkstemp(
            prefix='UclustExactMatchFilter', suffix='.fasta')
        close(fd)
        seqs_to_cluster, exact_match_id_map =\
            self._prefilter_exact_matches(parse_fasta(open(seq_path, 'U')))
        self.files_to_remove.append(unique_seqs_fp)
        unique_seqs_f = open(unique_seqs_fp, 'w')
        for seq_id, seq in seqs_to_cluster:
            unique_seqs_f.write('>%s\n%s\n' % (seq_id, seq))
        unique_seqs_f.close()
        # clean up the seqs_to_cluster list as it can be big and we
        # don't need it again
        # del(seqs_to_cluster)
        return exact_match_id_map, unique_seqs_fp


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
        _params = {'Similarity': 0.97,
                   'Application': 'uclust',
                   'max_accepts': 1,
                   'max_rejects': 8,
                   'stepwords': 8,
                   'word_length': 8,
                   'enable_rev_strand_matching': False,
                   'optimal': False,
                   'exact': False,
                   'suppress_sort': True,
                   'presort_by_abundance': True,
                   'new_cluster_identifier': None,
                   'stable_sort': True,
                   'save_uc_files': True,
                   'output_dir': '.',
                   'prefilter_identical_sequences': True}
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
        prefilter_identical_sequences =\
            self.Params['prefilter_identical_sequences']
        original_fasta_path = seq_path
        self.files_to_remove = []

        if self.Params['presort_by_abundance']:
            # seq path will become the temporary sorted sequences
            # filepath, to be cleaned up after the run
            seq_path = self._presort_by_abundance(seq_path)
            self.files_to_remove.append(seq_path)

        # Collapse idetical sequences to a new file
        if prefilter_identical_sequences:
            exact_match_id_map, seq_path =\
                self._apply_identical_sequences_prefilter(seq_path)

        # perform the clustering
        clusters, failures, seeds = get_clusters_from_fasta_filepath(
            seq_path,
            original_fasta_path,
            percent_ID=self.Params['Similarity'],
            optimal=self.Params['optimal'],
            exact=self.Params['exact'],
            suppress_sort=self.Params['suppress_sort'],
            enable_rev_strand_matching=
            self.Params['enable_rev_strand_matching'],
            max_accepts=self.Params['max_accepts'],
            max_rejects=self.Params['max_rejects'],
            stepwords=self.Params['stepwords'],
            word_length=self.Params['word_length'],
            stable_sort=self.Params['stable_sort'],
            save_uc_files=self.Params['save_uc_files'],
            output_dir=self.Params['output_dir'],
            HALT_EXEC=HALT_EXEC)

        # clean up any temp files that were created
        remove_files(self.files_to_remove)

        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))

        # expand identical sequences to create full OTU map
        if prefilter_identical_sequences:
            clusters = self._map_filtered_clusters_to_full_clusters(
                clusters, exact_match_id_map)

        otu_id_prefix = self.Params['new_cluster_identifier']
        if otu_id_prefix is None:
            clusters = enumerate(clusters)
        else:
            clusters = [('%s%d' % (otu_id_prefix, i), c)
                        for i, c in enumerate(clusters)]
        result = self._prepare_results(result_path, clusters, log_lines)

        if log_path:
            self._write_log(log_path, log_lines)

        # return the result (note this is None if the data was
        # written to file)
        return result


class UsearchOtuPicker(UclustOtuPickerBase):

    """ Usearch based OTU picker

    """

    Name = 'UsearchOtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.

        params contains both generic and per-method (e.g. for
        usearch application controller) params.

        Some generic entries in params are:

        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '--id' parameter
         to the uclust application controllers)
        Application: 3rd-party application used
        """

        _params = {
            'percent_id': 0.97,
            'percent_id_err': 0.97,
            'Application': 'usearch',
            'minsize': 4,
            'abundance_skew': 2,
            'db_filepath': None,
            'rev': False,
            'label_prefix': "",
            'label_suffix': "",
            'retain_label_as_comment': False,
            'count_start': 0,
            'perc_id_blast': 0.97,
            'save_intermediate_files': False,
            'global_alignment': True,
            'sizein': True,
            'sizeout': True,
            'w': 64,
            'slots': 16769023,
            'maxrejects': 64,
            'minlen': 64,
            'de_novo_chimera_detection': True,
            'reference_chimera_detection': True,
            'cluster_size_filtering': True,
            'output_dir': '.',
            'remove_usearch_logs': False,
            'derep_fullseq': False,
            'chimeras_retention': 'union',
            'verbose': False}

        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self,
                 seq_path,
                 output_dir='.',
                 log_path=None,
                 HALT_EXEC=False,
                 failure_path=None,
                 result_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu, and a list
         of seq ids that failed the filters.

        Parameters:
        seq_path: path to file of sequences
        output_dir: directory to output results, including log files and
         intermediate files if flagged for.
        log_path: path to log, which includes dump of params.

        """

        original_fasta_path = seq_path
        self.files_to_remove = []

        if self.Params['db_filepath'] is None:
            db_fp = None
        else:
            db_fp = abspath(self.Params['db_filepath'])

        # perform the filtering/clustering
        clusters, failures = usearch_qf(
            seq_path,
            output_dir=self.Params['output_dir'],
            percent_id=self.Params['percent_id'],
            percent_id_err=self.Params['percent_id_err'],
            minsize=self.Params['minsize'],
            abundance_skew=self.Params['abundance_skew'],
            db_filepath=db_fp,
            rev=self.Params['rev'],
            label_prefix=self.Params['label_prefix'],
            label_suffix=self.Params['label_suffix'],
            retain_label_as_comment=self.Params['retain_label_as_comment'],
            count_start=self.Params['count_start'],
            perc_id_blast=self.Params['perc_id_blast'],
            save_intermediate_files=self.Params['save_intermediate_files'],
            global_alignment=self.Params['global_alignment'],
            sizein=self.Params['sizein'],
            sizeout=self.Params['sizeout'],
            w=self.Params['w'],
            slots=self.Params['slots'],
            maxrejects=self.Params['maxrejects'],
            minlen=self.Params['minlen'],
            de_novo_chimera_detection=self.Params[
                'de_novo_chimera_detection'],
            reference_chimera_detection=self.Params[
                'reference_chimera_detection'],
            cluster_size_filtering=self.Params['cluster_size_filtering'],
            remove_usearch_logs=self.Params['remove_usearch_logs'],
            derep_fullseq=self.Params['derep_fullseq'],
            chimeras_retention=self.Params['chimeras_retention'],
            verbose=self.Params['verbose'],
            HALT_EXEC=HALT_EXEC)

        # clean up any temp files that were created
        remove_files(self.files_to_remove)

        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))
        log_lines.append('Num failures:%d' % len(failures))

        if failure_path:
            failure_file = open(failure_path, 'w')
            failure_file.write('\n'.join(failures))
            failure_file.close()

        if log_path:
            self._write_log(log_path, log_lines)

        if result_path:

            result_out = open(result_path, "w")
            for cluster_id in clusters:
                result_out.write(cluster_id + "\t" +
                                 "\t".join(clusters[cluster_id]) + '\n')

            result = None

        else:

            result = clusters

        return result


class UsearchReferenceOtuPicker(UclustOtuPickerBase):

    """ Usearch reference based OTU picker

    """

    Name = 'UsearchReferenceOtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.

        params contains both generic and per-method (e.g. for
        usearch application controller) params.

        Some generic entries in params are:

        Similarity: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '--id' parameter
         to the uclust application controllers)
        Application: 3rd-party application used
        """

        _params = {
            'percent_id': 0.97,
            'percent_id_err': 0.97,
            'Application': 'usearch',
            'minsize': 4,
            'abundance_skew': 2,
            'db_filepath': None,
            'rev': False,
            'label_prefix': "",
            'label_suffix': "",
            'retain_label_as_comment': False,
            'count_start': 0,
            'perc_id_blast': 0.97,
            'save_intermediate_files': False,
            'global_alignment': True,
            'sizein': True,
            'sizeout': True,
            'w': 64,
            'slots': 16769023,
            'maxrejects': 64,
            'minlen': 64,
            'de_novo_chimera_detection': True,
            'reference_chimera_detection': True,
            'cluster_size_filtering': True,
            'output_dir': '.',
            'remove_usearch_logs': False,
            'suppress_new_clusters': False,
            'derep_fullseq': False,
            'chimeras_retention': 'union',
            'verbose': False}

        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self,
                 seq_path,
                 refseqs_fp,
                 output_dir='.',
                 log_path=None,
                 HALT_EXEC=False,
                 failure_path=None,
                 result_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu, and a list
         of seq ids that failed the filters.

        Parameters:
        seq_path: path to file of sequences
        output_dir: directory to output results, including log files and
         intermediate files if flagged for.
        log_path: path to log, which includes dump of params.

        """

        original_fasta_path = seq_path
        self.files_to_remove = []

        # perform the filtering/clustering
        clusters, failures = usearch_qf(
            seq_path,
            refseqs_fp,
            output_dir=self.Params['output_dir'],
            percent_id=self.Params['percent_id'],
            percent_id_err=self.Params['percent_id_err'],
            minsize=self.Params['minsize'],
            abundance_skew=self.Params['abundance_skew'],
            db_filepath=self.Params['db_filepath'],
            rev=self.Params['rev'],
            label_prefix=self.Params['label_prefix'],
            label_suffix=self.Params['label_suffix'],
            retain_label_as_comment=self.Params['retain_label_as_comment'],
            count_start=self.Params['count_start'],
            perc_id_blast=self.Params['perc_id_blast'],
            save_intermediate_files=self.Params['save_intermediate_files'],
            global_alignment=self.Params['global_alignment'],
            sizein=self.Params['sizein'],
            sizeout=self.Params['sizeout'],
            w=self.Params['w'],
            slots=self.Params['slots'],
            maxrejects=self.Params['maxrejects'],
            minlen=self.Params['minlen'],
            de_novo_chimera_detection=self.Params[
                'de_novo_chimera_detection'],
            reference_chimera_detection=self.Params[
                'reference_chimera_detection'],
            cluster_size_filtering=self.Params['cluster_size_filtering'],
            remove_usearch_logs=self.Params['remove_usearch_logs'],
            suppress_new_clusters=self.Params['suppress_new_clusters'],
            derep_fullseq=self.Params['derep_fullseq'],
            chimeras_retention=self.Params['chimeras_retention'],
            verbose=self.Params['verbose'],
            HALT_EXEC=HALT_EXEC)

        # clean up any temp files that were created
        remove_files(self.files_to_remove)

        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))
        log_lines.append('Num failures:%d' % len(failures))
        log_lines.append('Reference database for OTU picking: %s' %
                         abspath(refseqs_fp))

        if failure_path is not None:
            failure_file = open(failure_path, 'w')
            failure_file.write('\n'.join(failures))
            failure_file.close()

        if log_path:
            self._write_log(log_path, log_lines)

        if result_path:

            result_out = open(result_path, "w")
            for cluster_id in clusters:
                result_out.write(cluster_id + "\t" +
                                 "\t".join(clusters[cluster_id]) + '\n')

            result = None

        else:

            result = clusters

        return result


class Usearch610DeNovoOtuPicker(UclustOtuPickerBase):

    """ Usearch based OTU picker, de novo clustering only

    """

    Name = 'Usearch610DeNovoOtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.

        params contains both generic and per-method (e.g. for
        usearch61 application controller) params.

        Some generic entries in params are:

        percent_id: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '--id' parameter
         to the uclust application controllers)
        Application: 3rd-party application used
        """

        _params = {
            'percent_id': 0.97,
            'Application': 'usearch61',
            'rev': False,
            'save_intermediate_files': False,
            'minlen': 64,
            'output_dir': '.',
            'remove_usearch_logs': False,
            'verbose': False,
            'wordlength': 8,
            'usearch_fast_cluster': False,
            'usearch61_sort_method': 'abundance',
            'usearch61_maxrejects': 32,
            'usearch61_maxaccepts': 1,
            'sizeorder': False,
            'threads': 1.0
        }

        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self,
                 seq_path,
                 output_dir='.',
                 log_path=None,
                 HALT_EXEC=False,
                 result_path=None,
                 otu_prefix="denovo"):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu

        Parameters:
        seq_path: path to file of sequences
        output_dir: directory to output results, including log files and
         intermediate files if flagged for.
        log_path: path to log, which includes dump of params.
        HALT_EXEC: Setting for halting execution of application controller,
         should only be True when debugging.
        result_path: If supplied will write out results (e.g. OTU mapping file),
         otherwise a dict is returned with data.

        """

        # perform de novo clustering
        clusters = usearch61_denovo_cluster(
            seq_path,
            percent_id=self.Params['percent_id'],
            rev=self.Params['rev'],
            save_intermediate_files=self.Params['save_intermediate_files'],
            minlen=self.Params['minlen'],
            output_dir=self.Params['output_dir'],
            remove_usearch_logs=self.Params['remove_usearch_logs'],
            verbose=self.Params['verbose'],
            wordlength=self.Params['wordlength'],
            usearch_fast_cluster=self.Params['usearch_fast_cluster'],
            usearch61_sort_method=self.Params['usearch61_sort_method'],
            otu_prefix=otu_prefix,
            usearch61_maxrejects=self.Params['usearch61_maxrejects'],
            usearch61_maxaccepts=self.Params['usearch61_maxaccepts'],
            sizeorder=self.Params['sizeorder'],
            threads=self.Params['threads'],
            HALT_EXEC=HALT_EXEC
        )

        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))

        if log_path:
            self._write_log(log_path, log_lines)

        if result_path:
            result_out = open(result_path, "w")
            for cluster_id in clusters:
                result_out.write(cluster_id + "\t" +
                                 "\t".join(clusters[cluster_id]) + '\n')
            result_out.close()
            result = None
        else:
            result = clusters

        return result


class Usearch61ReferenceOtuPicker(UclustOtuPickerBase):

    """ Usearch based OTU picker, supports closed or open reference OTU picking

    """

    Name = 'Usearch61ReferenceOtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.

        params contains both generic and per-method (e.g. for
        usearch61 application controller) params.

        Some generic entries in params are:

        percent_id: similarity threshold, default 0.97, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '--id' parameter
         to the uclust application controllers)
        Application: 3rd-party application used
        """

        _params = {
            'percent_id': 0.97,
            'Application': 'usearch61',
            'rev': False,
            'save_intermediate_files': False,
            'minlen': 64,
            'output_dir': '.',
            'remove_usearch_logs': False,
            'verbose': False,
            'wordlength': 8,
            'usearch_fast_cluster': False,
            'usearch61_sort_method': 'abundance',
            'usearch61_maxrejects': 32,
            'usearch61_maxaccepts': 1,
            'sizeorder': False,
            'suppress_new_clusters': False,
            'threads': 1.0
        }

        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self,
                 seq_path,
                 refseqs_fp,
                 output_dir='.',
                 log_path=None,
                 HALT_EXEC=False,
                 result_path=None,
                 failure_path=None,
                 otu_prefix="denovo"):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu

        Parameters:
        seq_path: path to file of sequences
        refseqs_fp: Reference database to pick OTUs against
        output_dir: directory to output results, including log files and
         intermediate files if flagged for.
        log_path: path to log, which includes dump of params.
        HALT_EXEC: Setting for halting execution of application controller,
         should only be True when debugging.
        result_path: If supplied will write out results (e.g. OTU mapping file),
         otherwise a dict is returned with data.

        """

        # perform reference clustering
        clusters, failures = usearch61_ref_cluster(
            seq_path,
            refseqs_fp,
            percent_id=self.Params['percent_id'],
            rev=self.Params['rev'],
            save_intermediate_files=self.Params['save_intermediate_files'],
            minlen=self.Params['minlen'],
            output_dir=self.Params['output_dir'],
            remove_usearch_logs=self.Params['remove_usearch_logs'],
            verbose=self.Params['verbose'],
            wordlength=self.Params['wordlength'],
            usearch_fast_cluster=self.Params['usearch_fast_cluster'],
            usearch61_sort_method=self.Params['usearch61_sort_method'],
            otu_prefix=otu_prefix,
            usearch61_maxrejects=self.Params['usearch61_maxrejects'],
            usearch61_maxaccepts=self.Params['usearch61_maxaccepts'],
            sizeorder=self.Params['sizeorder'],
            suppress_new_clusters=self.Params['suppress_new_clusters'],
            threads=self.Params['threads'],
            HALT_EXEC=HALT_EXEC
        )

        log_lines = []
        log_lines.append('Num OTUs:%d' % len(clusters))

        if log_path:
            self._write_log(log_path, log_lines)

        if result_path:
            result_out = open(result_path, "w")
            for cluster_id in clusters:
                result_out.write(cluster_id + "\t" +
                                 "\t".join(clusters[cluster_id]) + '\n')
            result_out.close()
            result = None
        else:
            result = clusters

        if failure_path:
            self._write_failures(failure_path, failures)

        return result, failures

    def _write_failures(self, failure_path, failures):
        failure_file = open(failure_path, 'w')
        failure_file.write('\n'.join(failures))
        failure_file.close()


class UclustReferenceOtuPicker(UclustOtuPickerBase):

    """Uclust reference OTU picker: clusters seqs by match to ref collection

    """

    def __init__(self, params):
        """Return new UclustReferenceOtuPicker object with specified params.

        """
        _params = {'Similarity': 0.97,
                   'Application': 'uclust',
                   'enable_rev_strand_matching': False,
                   'max_accepts': 1,
                   'max_rejects': 8,
                   'stepwords': 8,
                   'word_length': 8,
                   'suppress_new_clusters': False,
                   'optimal': False,
                   'exact': False,
                   'suppress_sort': False,
                   'new_cluster_identifier': 'QiimeOTU',
                   'next_new_cluster_number': 1,
                   'presort_by_abundance': True,
                   'stable_sort': True,
                   'save_uc_files': True,
                   'output_dir': '.',
                   'prefilter_identical_sequences': True}
        _params.update(params)
        OtuPicker.__init__(self, _params)

    def __call__(self,
                 seq_fp,
                 refseqs_fp,
                 next_new_cluster_number=None,
                 new_cluster_identifier=None,
                 result_path=None,
                 log_path=None,
                 failure_path=None,
                 HALT_EXEC=False):

        original_fasta_path = seq_fp
        prefilter_identical_sequences =\
            self.Params['prefilter_identical_sequences']

        if new_cluster_identifier:
            self.Params['new_cluster_identifier'] = new_cluster_identifier
        if next_new_cluster_number is not None:
            self.Params['next_new_cluster_number'] = next_new_cluster_number
        self.files_to_remove = []

        if self.Params['presort_by_abundance']:
            # seq path will become the temporary sorted sequences
            # filepath, to be cleaned up after the run
            seq_fp = self._presort_by_abundance(seq_fp)
            self.files_to_remove.append(seq_fp)

        # Collapse idetical sequences to a new file
        if prefilter_identical_sequences:
            exact_match_id_map, seq_fp =\
                self._apply_identical_sequences_prefilter(seq_fp)

        # perform the clustering
        cluster_map, failures, new_seeds = get_clusters_from_fasta_filepath(
            seq_fp,
            original_fasta_path,
            subject_fasta_filepath=refseqs_fp,
            percent_ID=self.Params['Similarity'],
            enable_rev_strand_matching=self.Params[
                'enable_rev_strand_matching'],
            max_accepts=self.Params['max_accepts'],
            max_rejects=self.Params['max_rejects'],
            stepwords=self.Params['stepwords'],
            word_length=self.Params['word_length'],
            suppress_new_clusters=self.Params['suppress_new_clusters'],
            optimal=self.Params['optimal'],
            exact=self.Params['exact'],
            suppress_sort=self.Params['suppress_sort'],
            return_cluster_maps=True,
            stable_sort=self.Params['stable_sort'],
            save_uc_files=self.Params['save_uc_files'],
            output_dir=self.Params['output_dir'],
            HALT_EXEC=HALT_EXEC)

        # expand identical sequences to create full OTU map
        if prefilter_identical_sequences:
            # expand the clusters (while retaining the names of
            # the clusters so we know which are new OTUs and
            # which are reference OTUs)
            cluster_names = cluster_map.keys()
            clusters = [cluster_map[c] for c in cluster_names]
            clusters = self._map_filtered_clusters_to_full_clusters(
                clusters, exact_match_id_map)
            cluster_map = dict(zip(cluster_names, clusters))

            # expand failures
            temp_failures = []
            for fa in failures:
                temp_failures.extend(exact_match_id_map[fa])
            failures = temp_failures

        self._rename_clusters(cluster_map, new_seeds)

        # clean up any temp files that were created
        remove_files(self.files_to_remove)

        log_lines = []
        log_lines.append('Reference seqs:%s' % abspath(refseqs_fp))
        log_lines.append('Num OTUs:%d' % len(cluster_map))
        log_lines.append('Num new OTUs:%d' % len(new_seeds))
        log_lines.append('Num failures:%d' % len(failures))

        cluster_map = cluster_map.items()
        result = self._prepare_results(result_path, cluster_map, log_lines)

        if log_path:
            self._write_log(log_path, log_lines)

        if failure_path:
            self._write_failures(failure_path, failures)

        # return the result (note this is None if the data was
        # written to file)
        return result

    def _rename_clusters(self, cluster_map, new_seeds):
        """ """
        next_new_cluster_number = self.Params['next_new_cluster_number']
        new_cluster_identifier = self.Params['new_cluster_identifier']
        new_seed_lookup = {}.fromkeys(new_seeds)

        for seed, cluster in cluster_map.items():
            del cluster_map[seed]
            if seed in new_seed_lookup:
                new_cluster_id = '%s%d' % (new_cluster_identifier,
                                           next_new_cluster_number)
                next_new_cluster_number += 1
            else:
                new_cluster_id = seed.split()[0]

            cluster_map[new_cluster_id] = cluster

        self.Params['next_new_cluster_number'] = next_new_cluster_number

    def _write_failures(self, failure_path, failures):
        # if the user provided a log file path, log the run
        failure_file = open(failure_path, 'w')
        failure_file.write('\n'.join(failures))
        failure_file.close()


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
            raise ValueError('Unsupported algorithm %s.  Choices are %s' %
                             (params['Algorithm'], self.ClusteringAlgorithms))
        super(MothurOtuPicker, self).__init__(params)

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        app = Mothur(
            InputHandler='_input_as_path',
            TmpDir=get_qiime_temp_dir())
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
            of = open(result_path, 'w')
            for i, cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i, '\t'.join(cluster)))
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
            log_file = open(log_path, 'w')
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


def expand_otu_map_seq_ids(otu_map, seq_id_map):
    for otu_id, seq_ids in otu_map.items():
        mapped_seq_ids = flatten(
            [seq_id_map[seq_id] for seq_id in seq_ids])
        otu_map[otu_id] = mapped_seq_ids
    return otu_map


def expand_failures(failures, seq_id_map):
    result = []
    for failure in failures:
        failure = failure.strip()
        result += seq_id_map[failure]
    return result


def map_otu_map_files(otu_files, failures_file=None):
    # passing delim=None splits on any whitespace, so can handle mixed tabs
    # and spaces
    result = fields_to_dict(otu_files[0], delim=None)
    for otu_file in otu_files[1:]:
        current_otu_map = fields_to_dict(otu_file, delim=None)
        result = expand_otu_map_seq_ids(current_otu_map, result)
    if failures_file:
        result = expand_failures(failures_file, result)
    return result

# End functions to support merging OTU tables


otu_picking_method_constructors = {
    'cdhit': CdHitOtuPicker,
    'prefix_suffix': PrefixSuffixOtuPicker,
    'mothur': MothurOtuPicker,
    'trie': TrieOtuPicker,
    'blast': BlastOtuPicker,
    'uclust': UclustOtuPicker,
    'uclust_ref': UclustReferenceOtuPicker,
    'usearch': UsearchOtuPicker,
    'usearch_ref': UsearchReferenceOtuPicker,
    'usearch61': Usearch610DeNovoOtuPicker,
    'usearch61_ref': Usearch61ReferenceOtuPicker,
    'sumaclust': SumaClustOtuPicker,
    'sortmerna': SortmernaV2OtuPicker,
    'swarm': SwarmOtuPicker
}

otu_picking_method_choices = otu_picking_method_constructors.keys()
