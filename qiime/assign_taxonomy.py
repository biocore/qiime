#!/usr/bin/env python

from __future__ import division

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger",
               "Antonio Gonzalez Pena", "David Soergel", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


import logging
import os
import re
from os import remove
from os.path import abspath, dirname
from itertools import count
from string import strip
from tempfile import NamedTemporaryFile, mkdtemp
from cStringIO import StringIO
from collections import Counter, defaultdict
from shutil import rmtree

from skbio.parse.sequences import parse_fasta
from skbio.util import remove_files

from bfillings.blast import blast_seqs, Blastall, BlastResult
from bfillings.formatdb import build_blast_db_from_fasta_path
from bfillings.uclust import Uclust
from bfillings.sortmerna_v2 import (build_database_sortmerna,
                                 sortmerna_map)
from bfillings import rdp_classifier
from bfillings import mothur
from bfillings import rtax

from qiime.util import FunctionWithParams, get_rdp_jarpath, get_qiime_temp_dir

"""Contains code for assigning taxonomy, using several techniques.

This module has the responsibility for taking a set of sequences and
providing a taxon assignment for each sequence."""


def validate_rdp_version(rdp_jarpath=None):
    if rdp_jarpath is None:
        rdp_jarpath = get_rdp_jarpath()
    if rdp_jarpath is None:
        raise RuntimeError(
            "RDP classifier is not installed or not accessible to QIIME. "
            "See install instructions here: "
            "http://qiime.org/install/install.html#rdp-install"
        )

    rdp_jarname = os.path.basename(rdp_jarpath)
    version_match = re.search("\d\.\d", rdp_jarname)
    if version_match is None:
        raise RuntimeError(
            "Unable to detect RDP Classifier version in file %s" % rdp_jarname
        )

    version = float(version_match.group())
    if version < 2.1:
        raise RuntimeError(
            "RDP Classifier does not look like version 2.2 or greater."
            "Versions of the software prior to 2.2 have different "
            "formatting conventions and are no longer supported by QIIME. "
            "Detected version %s from file %s" % (version, rdp_jarpath)
        )
    return version


class TaxonAssigner(FunctionWithParams):

    """A TaxonAssigner assigns a taxon to each of a set of sequences.

    This is an abstract class: subclasses should implement the __call__
    method.
    """
    Name = 'TaxonAssigner'

    def __init__(self, params):
        """Return new TaxonAssigner object with specified params.

        Note: expect params to contain both generic and per-method (e.g. for
        RDP classifier w/ Hugenholtz taxonomy) params, so leaving it as a dict
        rather than setting attributes. Some standard entries in params are:

        Taxonomy: taxonomy used (e.g. RDP, Hugenholtz)
        Similarity: similarity threshold for assignment, e.g. 0.97
        Bootstrap: bootstrap support for assignment, e.g. 0.80
        Application: 3rd-party application used, if any, e.g. RDP classifier
        """
        self.Params = params

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError("TaxonAssigner is an abstract class")

    @staticmethod
    def _parse_id_to_taxonomy_file(f):
        """ parse the id_to_taxonomy file into a dict mapping id -> taxonomy
        """
        result = {}
        for line in f:
            line = line.strip()
            if line:
                identifier, taxonomy = map(strip, line.split('\t'))
                result[identifier] = taxonomy
        return result

    def _tax_assignments_to_consensus_assignments(self,
                                                  query_to_assignments):
        """ For each query id and list of assignments,
            call _get_consensus_assigment to compute the
            consensus assignment.

            Parameters
            ----------
            query_to_assignments : dict of list of lists
                The keys in the dict correspond to query IDs and
                the values are a list of lists holding associated
                taxonomies.

            Returns
            -------
            query_to_assignments: dict
                The keys in the dict correspond to query IDs and
                the values carry a single consensus taxonomy
                assignment.
        """
        for query_id, assignments in query_to_assignments.iteritems():
            consensus_assignment = self._get_consensus_assignment(assignments)
            query_to_assignments[query_id] = consensus_assignment

        return query_to_assignments

    def _get_consensus_assignment(self, assignments):
        """ compute the consensus assignment from a list of assignments
            (method applied to SortMeRNATaxonAssigner and UclustConsensusTaxonAssigner)
        """
        num_input_assignments = len(assignments)
        consensus_assignment = []

        # if the assignments don't all have the same number
        # of levels, the resulting assignment will have a max number
        # of levels equal to the number of levels in the assignment
        # with the fewest number of levels. this is to avoid
        # a case where, for example, there are n assignments, one of
        # which has 7 levels, and the other n-1 assignments have 6 levels.
        # A 7th level in the result would be misleading because it
        # would appear to the user as though it was the consensus
        # across all n assignments.
        num_levels = min([len(a) for a in assignments])

        # iterate over the assignment levels
        for level in range(num_levels):
            # count the different taxonomic assignments at the current level.
            # the counts are computed based on the current level and all higher
            # levels to reflect that, for example, 'p__A; c__B; o__C' and
            # 'p__X; c__Y; o__C' represent different taxa at the o__ level (since
            # they are different at the p__ and c__ levels).
            current_level_assignments = \
                Counter([tuple(e[:level + 1]) for e in assignments])
            # identify the most common taxonomic assignment, and compute the
            # fraction of assignments that contained it. it's safe to compute the
            # fraction using num_assignments because the deepest level we'll
            # ever look at here is num_levels (see above comment on how that
            # is decided).
            tax, max_count = current_level_assignments.most_common(1)[0]
            max_consensus_fraction = max_count / num_input_assignments
            # check whether the most common taxonomic assignment is observed
            # in at least min_consensus_fraction of the sequences
            if max_consensus_fraction >= self.Params['min_consensus_fraction']:
                # if so, append the current level only (e.g., 'o__C' if tax is
                # 'p__A; c__B; o__C', and continue on to the next level
                consensus_assignment.append((tax[-1], max_consensus_fraction))
            else:
                # if not, there is no assignment at this level, and we're
                # done iterating over levels
                break

        # construct the results
        # determine the number of levels in the consensus assignment
        consensus_assignment_depth = len(consensus_assignment)
        if consensus_assignment_depth > 0:
            # if it's greater than 0, generate a list of the
            # taxa assignments at each level
            assignment_result = [a[0] for a in consensus_assignment]
            # and assign the consensus_fraction_result as the
            # consensus fraction at the deepest level
            consensus_fraction_result = \
                consensus_assignment[consensus_assignment_depth - 1][1]
        else:
            # if there are zero assignments, indicate that the taxa is
            # unknown
            assignment_result = [self.Params['unassignable_label']]
            # and assign the consensus_fraction_result to 1.0 (this is
            # somewhat arbitrary, but could be interpreted as all of the
            # assignments suggest an unknown taxonomy)
            consensus_fraction_result = 1.0

        return (
            assignment_result, consensus_fraction_result, num_input_assignments
        )


class SortMeRNATaxonAssigner(TaxonAssigner):
    """ Assign taxonomy using SortMeRNA
    """

    Name = 'SortMeRNATaxonAssigner'
    Application = "SortMeRNA"
    Citation = ("SortMeRNA is hosted at:\n"
                "http://bioinfo.lifl.fr/RNA/sortmerna\n"
                "https://github.com/biocore/sortmerna\n\n"
                "The following paper should be cited if this resource is "
                "used:\n\n"
                "Kopylova, E., Noe L. and Touzet, H.,\n"
                "SortMeRNA: fast and accurate filtering of ribosomal RNAs "
                "in\n"
                "metatranscriptomic data, Bioinformatics (2012) 28(24)\n"
                )
    _tracked_properties = ['Application', 'Citation']

    def __init__(self, params):
        _params = {
            # id to taxonomy filepath
            'id_to_taxonomy_fp': None,
            # reference sequences filepath
            'reference_sequences_fp': None,
            # reference sequences indexed database
            'sortmerna_db': None,
            # Fraction of sequence hits that a taxonomy assignment
            # must show up in to be considered the consensus assignment
            'min_consensus_fraction': 0.51,
            # minimum identity to consider a hit
            'min_percent_id': 90.0,
            # minimum query coverage to consider a hit
            'min_percent_cov': 90.0,
            # output 10 best alignments
            'best_N_alignments': 10,
            # E-value
            'e_value': 1,
            # threads
            'threads': 1,
            # label to apply for queries that cannot be assigned
            'unassignable_label': 'Unassigned'
        }
        _params.update(params)
        super(SortMeRNATaxonAssigner, self).__init__(_params)

    def __call__(self,
                 seq_path,
                 result_path=None,
                 log_path=None,
                 HALT_EXEC=False):
        """Returns mapping of each seq to (taxonomy, consensus fraction, n).

        Parameters
        ----------
        seq_path : str, mandatory
            The filepath to input sequences.
        result_path : str, optional
            The filepath to store resulting alignments.
        log_path : str, optional
            The filepath to store logging information.
        HALT_EXEC : bool, debugging parameter
            If passed, will exit just before the sortmerna command in issued
            and will print out the command that would have been called
            to stdout.

        Returns
        -------
        dict if result_path=None
            The results will be stored in a dict:
                dict{query_id:[tax, consensus fraction, n]}

        None if result_path
            The results will be written to result_path as tab-separated
            lines of:
                query_id <tab> tax <tab> consensus fraction <tab> n

            The values represent:
            tax: the consensus taxonomy assignment
            consensus fraction: the fraction of the assignments for the
            query that contained the lowest level tax assignment that is
            included in tax (e.g., if the assignment goes to genus level,
            this will be the fraction of assignments that had the consensus
            genus assignment)
            n: the number of assignments that were considered when
            constructing the consensus
        """
        # Check input reference sequence and taxonomy are provided
        if self.Params['reference_sequences_fp'] is None:
            raise ValueError("Filepath for reference sequences is mandatory.")
        if self.Params['id_to_taxonomy_fp'] is None:
            raise ValueError("Filepath for id to taxonomy map is mandatory.")

        # initialize the logger
        logger = self._get_logger(log_path)
        logger.info(str(self))
        self.dirs_to_remove = []

        # Indexed database not provided, build it
        if not self.Params['sortmerna_db']:
            output_dir = mkdtemp()
            self.sortmerna_db, files_to_remove = \
                build_database_sortmerna(abspath(self.Params[
                    'reference_sequences_fp']),
                    output_dir=output_dir)
            self.dirs_to_remove.append(output_dir)

        # Indexed database provided
        else:
            self.sortmerna_db = self.Params['sortmerna_db']

        # Set SortMeRNA's output directory
        if result_path is None:
            output_dir = mkdtemp()
            self.dirs_to_remove.append(output_dir)
        else:
            output_dir = dirname(abspath(result_path))

        # Call sortmerna mapper
        app_result =\
            sortmerna_map(seq_path=seq_path,
                          output_dir=output_dir,
                          sortmerna_db=self.sortmerna_db,
                          refseqs_fp=self.Params['reference_sequences_fp'],
                          e_value=self.Params['e_value'],
                          threads=self.Params['threads'],
                          best=self.Params['best_N_alignments'],
                          HALT_EXEC=False)

        with open(self.Params['id_to_taxonomy_fp'], "U") as id_to_taxonomy_f:
            self.id_to_taxonomy_map =\
                self._parse_id_to_taxonomy_file(id_to_taxonomy_f)

        blast_tabular_fp = app_result['BlastAlignments'].name
        query_to_assignments = self._blast_to_tax_assignments(blast_tabular_fp)
        result = self._tax_assignments_to_consensus_assignments(
            query_to_assignments)

        # Write results to file
        if result_path is not None:
            with open(result_path, 'w') as of:
                of.write('#OTU ID\ttaxonomy\tconfidence\tnum hits\n')
                for seq_id, (assignment, consensus_fraction, n) in result.items():
                    assignment_str = ';'.join(assignment)
                    of.write('%s\t%s\t%1.2f\t%d\n' % (
                        seq_id, assignment_str, consensus_fraction, n))
            result = None
            logger.info('Result path: %s' % result_path)
        else:
            # If no result_path was provided, the result dict is
            # returned as-is.
            logger.info('Result path: None, returned as dict.')

        # clean up
        map(rmtree, self.dirs_to_remove)

        return result

    def _get_logger(self, log_path=None):
        if log_path is not None:
            handler = logging.FileHandler(log_path, mode='w')
        else:
            class NullHandler(logging.Handler):

                def emit(self, record):
                    pass
            handler = NullHandler()
        logger = logging.getLogger("SortMeRNATaxonAssigner logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _blast_to_tax_assignments(self,
                                  blast_output_fp):
        """ Parse SortMeRNA's Blast-like tabular format for query
            IDs and the references they map to, use the reference IDs
            to find the associated taxonomies in the id_to_taxonomy_map.

            Three types of alignments are possible,

            1. The Null alignment (E-value threshold failed):
            not16S.1_130\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t*\t0

            2. All alignments for a query pass the E-value threshold
            but fail the %id threshold (3rd column is %id):
            f1_4866\t426848\t85.4\t121\t15\t3\t1\t121\t520\t641\t4.79e-32\t131\t72M1D7M1I13M1D28M31S\t79.6
            f1_4866\t342684\t84\t91\t9\t6\t1\t91\t522\t612\t2.8e-19\t89\t55M1D4M1I12M1D3M1D4M1I1M1I9M61S\t59.9

            3. Some/all alignments for a query pass both E-value and %id
            thresholds:
            f2_1271\t295053\t100\t128\t0\t0\t1\t128\t520\t647\t1.15e-59\t223\t128M\t100
            f2_1271\t42684\t84.8\t124\t17\t2\t1\t124\t527\t650\t2.63e-32\t132\t101M1D6M1I16M4S\t96.9

            Parameters
            ----------
            blast_output_fp : str
                Filepath to Blast-like tabular alignments.

            Returns
            -------
            result : dict of list of lists
                The keys in the dict correspond to query IDs and
                the values are a list of lists holding associated
                taxonomies.
        """
        min_percent_id = self.Params['min_percent_id']
        result = defaultdict(list)
        with open(blast_output_fp, "U") as blast_output:
            for line in blast_output:
                fields = line.split('\t')
                query_id = fields[0]
                subject_id = fields[1]
                percent_id = float(fields[2])
                # sequence was not aligned
                if subject_id == "*":
                    result[query_id].append([])
                # sequence was aligned, passing %id threshold
                elif percent_id >= min_percent_id:
                    # if exists, remove the empty alignment (failing %id
                    # threshold) for this sequence (Blast tabular output
                    # will list all alignments passing E-value threshold,
                    # not necessarily the %id threshold). It should happen
                    # rarely that an alignment passing the %id threshold
                    # comes after an alignment that failed the threshold,
                    # but it can happen (Blast alignments are often ordered
                    # from highest %id to lowest), though as sortmerna uses
                    # a heuristic, this isn't always guaranteed.
                    if [] in result[query_id]:
                        result[query_id].remove([])
                    # add alignment passing %id threshold
                    subject_tax = self.id_to_taxonomy_map[
                        subject_id].strip().split(';')
                    result[query_id].append(subject_tax)
                # sequence was aligned, however failing %id threshold
                # if no alignment results have been recorded for this
                # sequence up to now, add an empty list
                elif not result[query_id]:
                    result[query_id].append([])

        return result


class BlastTaxonAssigner(TaxonAssigner):

    """ Assign taxon best on best blast hit above a threshold
    """
    Name = 'BlastTaxonAssigner'
    SeqsPerBlastRun = 1000

    def __init__(self, params):
        """ Initialize the object
        """
        _params = {
            'Min percent identity': 90.0,
            'Max E value': 1e-30,
            'Application': 'blastn/megablast'
        }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)

    def __call__(self, seq_path=None, seqs=None,
                 result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.
        """
        assert seq_path or seqs, \
            "Must provide either seqs or seq_path when calling a BlastTaxonAssigner."

        # initialize the logger
        logger = self._get_logger(log_path)
        logger.info(str(self))

        # assign the blast database, either as a pre-exisiting database
        # specified as self.Params['blast_db'] or by creating a
        # temporary database from the sequence file specified
        # as self.Params['reference_seqs_filepath']
        try:
            blast_db = self.Params['blast_db']
        except KeyError:
            # build a temporary blast_db
            reference_seqs_path = self.Params['reference_seqs_filepath']
            refseqs_dir, refseqs_name = os.path.split(reference_seqs_path)

            blast_db_dir = mkdtemp(prefix='bltax-', dir=get_qiime_temp_dir())
            blast_db, db_files_to_remove = build_blast_db_from_fasta_path(
                abspath(reference_seqs_path), output_dir=blast_db_dir)

        # build the mapping of sequence identifier
        # (wrt to the blast db seqs) to taxonomy
        id_to_taxonomy_map = self._parse_id_to_taxonomy_file(
            open(self.Params['id_to_taxonomy_filepath'], 'U'))

        # Iterate over the input self.SeqsPerBlastRun seqs at a time.
        # There are two competing issues here when dealing with very large
        # inputs. If all sequences are read in at once, the containing object
        # can be very large, causing the system to page. On the other hand,
        # in such cases it would be very slow to treat each sequence
        # individually, since blast requires a filepath. Each call would
        # therefore involve writing a single sequence to file, opening/closing
        # and removing the file. To balance this, sequences are read in and
        # blasted in chunks of self.SeqsPerBlastRun (defualt: 1000) at a time.
        # This appears to solve the problem with the largest sets I've worked
        # with so far.

        if seq_path:
            # Get a seq iterator
            seqs = parse_fasta(open(seq_path))
        # Build object to keep track of the current set of sequence to be
        # blasted, and the results (i.e., seq_id -> (taxonomy,quaility score)
        # mapping)
        current_seqs = []
        result = {}

        # Iterate over the (seq_id, seq) pairs
        for seq_id, seq in seqs:
            # append the current seq_id,seq to list of seqs to be blasted
            current_seqs.append((seq_id, seq))

            # When there are 1000 in the list, blast them
            if len(current_seqs) == self.SeqsPerBlastRun:
                # update the result object
                result.update(self._seqs_to_taxonomy(
                    current_seqs, blast_db, id_to_taxonomy_map))
                # reset the list of seqs to be blasted
                current_seqs = []
        # Assign taxonomy to the remaining sequences
        result.update(self._seqs_to_taxonomy(
            current_seqs, blast_db, id_to_taxonomy_map))
        # End iteration over the input self.SeqsPerBlastRun seqs at a time.

        # Write log data if we have a path (while the logger can handle
        # being called if we are not logging, some of these steps are slow).
        if log_path is not None:
            num_inspected = len(result)
            logger.info('Number of sequences inspected: %s' % num_inspected)
            num_null_hits = [r[1] for r in result.values()].count(None)
            logger.info('Number with no blast hits: %s' % num_null_hits)

        if result_path:
            # if the user provided a result_path, write the
            # results to file
            of = open(result_path, 'w')
            for seq_id, (lineage, confidence, blast_hit_id) in result.items():
                of.write('%s\t%s\t%s\t%s\n' %
                         (seq_id, lineage, confidence, blast_hit_id))
            of.close()
            result = None
            logger.info('Result path: %s' % result_path)
        else:
            # Returning the data as a dict, so no modification to result
            # is necessary.
            pass

            # if no result_path was provided, return the data as a dict
            logger.info('Result path: None, returned as dict.')

        # clean-up temp blastdb files, if a temp blastdb was created
        if 'reference_seqs_filepath' in self.Params:
            remove_files(db_files_to_remove)
            rmtree(blast_db_dir)

        # return the result
        return result

    def _seqs_to_taxonomy(self, seqs, blast_db, id_to_taxonomy_map):
        """ Assign taxonomy to (seq_id,seq) pairs
        """
        # Handle the case of no seqs passed in
        if not seqs:
            return {}
        # blast the seqs
        blast_hits = self._get_blast_hits(blast_db, seqs)

        # select the best blast hit for each query sequence
        best_blast_hit_ids = self._get_first_blast_hit_per_seq(blast_hits)

        # map the identifier of the best blast hit to (taxonomy, e-value)
        return self._map_ids_to_taxonomy(
            best_blast_hit_ids, id_to_taxonomy_map)

    def _get_logger(self, log_path=None):
        if log_path is not None:
            handler = logging.FileHandler(log_path, mode='w')
        else:
            class NullHandler(logging.Handler):

                def emit(self, record):
                    pass
            handler = NullHandler()
        logger = logging.getLogger("BlastTaxonAssigner logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _map_ids_to_taxonomy(self, hits, id_to_taxonomy_map):
        """ map {query_id:(best_blast_seq_id,e-val)} to {query_id:(tax,e-val,best_blast_seq_id)}
        """
        for query_id, hit in hits.items():
            query_id = query_id.split()[0]
            try:
                hit_id, e_value = hit
                hits[query_id] = \
                    (id_to_taxonomy_map.get(hit_id, None), e_value, hit_id)
            except TypeError:
                hits[query_id] = ('No blast hit', None, None)

        return hits

    def _get_blast_hits(self, blast_db, seqs):
        """ blast each seq in seqs against blast_db and retain good hits
        """
        max_evalue = self.Params['Max E value']
        min_percent_identity = self.Params['Min percent identity']
        if min_percent_identity < 1.0:
            min_percent_identity *= 100.0
        seq_ids = [s[0] for s in seqs]
        result = {}

        blast_result = blast_seqs(
            seqs, Blastall, blast_db=blast_db,
            params={'-p': 'blastn', '-n': 'T'},
            add_seq_names=False)

        if blast_result['StdOut']:
            lines = [x for x in blast_result['StdOut']]
            blast_result = BlastResult(lines)
        else:
            return {}.fromkeys(seq_ids, [])

        for seq_id in seq_ids:
            blast_result_id = seq_id.split()[0]
            try:
                result[seq_id] = [(e['SUBJECT ID'], float(e['E-VALUE']))
                                  for e in blast_result[blast_result_id][0]
                                  if (float(e['E-VALUE']) <= max_evalue and
                                      float(e['% IDENTITY']) >= min_percent_identity)]
            except KeyError:
                result[seq_id] = []

        return result

    def _get_first_blast_hit_per_seq(self, blast_hits):
        """ discard all blast hits except the best for each query sequence
        """
        result = {}
        for k, v in blast_hits.items():
            k = k.split()[0]  # get rid of spaces
            try:
                result[k] = v[0]
            except IndexError:
                # If there is no good blast hit, do we want to
                # leave the key out, or have it point to None?
                result[k] = None

        return result


class MothurTaxonAssigner(TaxonAssigner):

    """Assign taxonomy using Mothur's naive Bayes implementation
    """
    Name = 'MothurTaxonAssigner'
    Application = "Mothur"
    Citation = (
        "Schloss, P.D., et al., Introducing mothur: Open-source, platform-"
        "independent, community-supported software for describing and "
        "comparing microbial communities. Appl Environ Microbiol, 2009. "
        "75(23):7537-41."
    )
    _tracked_properties = ['Application', 'Citation']

    def __init__(self, params):
        _params = {
            'Confidence': 0.50,
            'Iterations': None,
            'KmerSize': None,
            'id_to_taxonomy_fp': None,
            'reference_sequences_fp': None,
        }
        _params.update(params)
        super(MothurTaxonAssigner, self).__init__(_params)

    def _format_id_to_taxonomy(self, id_to_taxonomy_file):
        """Reformat taxa to comply with Mothur formatting requirements.

        Mothur requires lineages to be semicolon-separated with no space
        following the semicolon.  (QIIME convention is to include a
        space.)  Taxa may have no internal spaces.  Furthermore, each
        lineage must end with a semi-colon.

        Returns the re-formatted id-to-taxonomy file as an open file
        object.
        """
        mothur_tax_file = NamedTemporaryFile(
            prefix='MothurTaxonAssigner_',
            suffix='.txt',
            dir=get_qiime_temp_dir())
        original_taxonomy = self._parse_id_to_taxonomy_file(id_to_taxonomy_file)
        for seq_id, lineage in original_taxonomy.iteritems():
            mothur_tax_file.write(seq_id)
            mothur_tax_file.write('\t')
            taxa = [t.strip() for t in lineage.split(';')]
            for taxon in taxa:
                mothur_tax_file.write(self._format_taxon(taxon))
                mothur_tax_file.write(';')
            mothur_tax_file.write('\n')
        mothur_tax_file.seek(0)
        return mothur_tax_file

    def _unformat_result(self, result):
        """Transform results to remove any changes introduced by formatting.
        """
        unformatted_result = {}
        for seq_id, (taxa, conf) in result.iteritems():
            unformatted_taxa = [self._unformat_taxon(t) for t in taxa]
            unformatted_result[seq_id] = (unformatted_taxa, conf)
        return unformatted_result

    def _format_taxon(self, taxon):
        """Format taxon for MOTHUR, removing internal spaces.

        Original taxon names are saved to self._original_taxa for later lookup.
        """
        # Create private attribute to store unformatted taxon names.
        # If _unformat_taxon() is called without first calling
        # _format_taxon(), this attribute will be missing, and an
        # AttributeError will be raised.
        if not hasattr(self, "_original_taxa"):
            self._original_taxa = {}
        # Escape backslashes
        mothur_taxon = taxon.replace("\\", "\\\\")
        # Escape underscores
        mothur_taxon = mothur_taxon.replace("_", "\\_")
        # Now we can safely replace spaces with underscores
        mothur_taxon = mothur_taxon.replace(' ', '_')
        if mothur_taxon != taxon:
            previously_registered_taxon = self._original_taxa.get(mothur_taxon)
            # If we have not yet registered the escaped taxon name, add it now.
            if previously_registered_taxon is None:
                self._original_taxa[mothur_taxon] = taxon
            # Otherwise, check that the previously registered taxon is
            # consistent with the current taxon.  If we have not
            # escaped the taxon names properly, two distinct taxa may
            # be registered under the same name.  This should probably
            # never happen, but I can't prove it, so we check and
            # raise an error if the taxa are inconsistent.
            elif taxon != previously_registered_taxon:
                raise ValueError(
                    "Taxon %s conflicts with another taxon, %s. "
                    "Please change one of the names." % (
                        taxon, previously_registered_taxon))
        return mothur_taxon

    def _unformat_taxon(self, taxon):
        """Recover original taxon names that were altered due to formatting.

        Looks up taxon names in the attribute self._original_taxa.  If
        self._format_taxon() was never called, this attribute will be
        missing, and an AttributeError will be raised.
        """
        return self._original_taxa.get(taxon, taxon)

    def __call__(self, seq_path, result_path=None, log_path=None):
        seq_file = open(seq_path)
        percent_confidence = int(self.Params['Confidence'] * 100)
        with open(self.Params['id_to_taxonomy_fp'], "U") as tax_file:
            mothur_tax_file = self._format_id_to_taxonomy(tax_file)
        try:
            result = mothur.mothur_classify_file(
                query_file=seq_file,
                ref_fp=self.Params['reference_sequences_fp'],
                tax_fp=mothur_tax_file.name,
                cutoff=percent_confidence,
                iters=self.Params['Iterations'],
                ksize=self.Params['KmerSize'],
                output_fp=None,
                tmp_dir=get_qiime_temp_dir()
            )
        finally:
            mothur_tax_file.close()
        result = self._unformat_result(result)
        if result_path is not None:
            with open(result_path, "w") as f:
                for seq_id, (taxa, conf) in result.iteritems():
                    lineage = ';'.join(taxa)
                    f.write("%s\t%s\t%.2f\n" % (seq_id, lineage, conf))
            return None
        if log_path:
            self.writeLog(log_path)
        return result


class RdpTaxonAssigner(TaxonAssigner):

    """Assign taxon using RDP's naive Bayesian classifier
    """
    Name = "RdpTaxonAssigner"
    Application = "RDP classfier"
    Citation = "Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7."
    Taxonomy = "RDP"
    _tracked_properties = ['Application', 'Citation', 'Taxonomy']

    def __init__(self, params):
        """Return new RdpTaxonAssigner object with specified params.

        Standard entries in params are:

        Taxonomy: taxonomy used (e.g. RDP, Hugenholtz)
        """
        _params = {
            'Confidence': 0.50,
            'id_to_taxonomy_fp': None,
            'reference_sequences_fp': None,
            'training_data_properties_fp': None,
            'max_memory': None
        }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for
        each seq.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, dumps the
            result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        tmp_dir = get_qiime_temp_dir()
        min_conf = self.Params['Confidence']
        training_data_properties_fp = self.Params[
            'training_data_properties_fp']
        reference_sequences_fp = self.Params['reference_sequences_fp']
        id_to_taxonomy_fp = self.Params['id_to_taxonomy_fp']
        max_memory = self.Params['max_memory']

        seq_file = open(seq_path, 'U')
        if reference_sequences_fp and id_to_taxonomy_fp:
            # Train and assign taxonomy
            taxonomy_file, training_seqs_file = self._generate_training_files()
            results = rdp_classifier.train_rdp_classifier_and_assign_taxonomy(
                training_seqs_file, taxonomy_file, seq_file,
                min_confidence=min_conf,
                classification_output_fp=result_path,
                max_memory=max_memory, tmp_dir=tmp_dir)

            if result_path is None:
                results = self._training_set.fix_results(results)
            else:
                self._training_set.fix_output_file(result_path)
        else:
            # Just assign taxonomy, using properties file if passed
            if training_data_properties_fp:
                fix_ranks = False
            else:
                fix_ranks = True
            results = rdp_classifier.assign_taxonomy(
                seq_file, min_confidence=min_conf, output_fp=result_path,
                training_data_fp=training_data_properties_fp,
                max_memory=max_memory, fixrank=fix_ranks, tmp_dir=tmp_dir)

        if log_path:
            self.writeLog(log_path)

        return results

    def _generate_training_files(self):
        """Returns a tuple of file objects suitable for passing to the
        RdpTrainer application controller.
        """
        tmp_dir = get_qiime_temp_dir()
        training_set = RdpTrainingSet()
        reference_seqs_file = open(self.Params['reference_sequences_fp'], 'U')
        id_to_taxonomy_file = open(self.Params['id_to_taxonomy_fp'], 'U')

        for seq_id, seq in parse_fasta(reference_seqs_file):
            training_set.add_sequence(seq_id, seq)

        for line in id_to_taxonomy_file:
            seq_id, lineage_str = map(strip, line.split('\t'))
            training_set.add_lineage(seq_id, lineage_str)

        training_set.dereplicate_taxa()

        rdp_taxonomy_file = NamedTemporaryFile(
            prefix='RdpTaxonAssigner_taxonomy_', suffix='.txt', dir=tmp_dir)
        rdp_taxonomy_file.write(training_set.get_rdp_taxonomy())
        rdp_taxonomy_file.seek(0)

        rdp_training_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssigner_training_seqs_', suffix='.fasta',
            dir=tmp_dir)
        for rdp_id, seq in training_set.get_training_seqs():
            rdp_training_seqs_file.write('>%s\n%s\n' % (rdp_id, seq))
        rdp_training_seqs_file.seek(0)

        self._training_set = training_set

        return rdp_taxonomy_file, rdp_training_seqs_file


class RdpTrainingSet(object):

    def __init__(self):
        self._tree = RdpTree()
        self.sequences = {}
        self.sequence_nodes = {}
        self.lineage_depth = None

    def add_sequence(self, seq_id, seq):
        self.sequences[seq_id] = seq

    def add_lineage(self, seq_id, lineage_str):
        for char, escape_str in _QIIME_RDP_ESCAPES:
            lineage_str = re.sub(char, escape_str, lineage_str)
        lineage = self._parse_lineage(lineage_str)
        seq_node = self._tree.insert_lineage(lineage)
        self.sequence_nodes[seq_id] = seq_node

    def dereplicate_taxa(self):
        return self._tree.dereplicate_taxa()

    def _parse_lineage(self, lineage_str):
        """Returns a list of taxa from the semi-colon-separated
        lineage string of an id_to_taxonomy file.
        """
        lineage = lineage_str.strip().split(';')
        if self.lineage_depth is None:
            self.lineage_depth = len(lineage)
        if len(lineage) != self.lineage_depth:
            raise ValueError(
                'Because the RDP Classifier operates in a bottom-up manner, '
                'each taxonomy assignment in the id-to-taxonomy file must have '
                'the same number of ranks.  Detected %s ranks in the first '
                'item of the file, but detected %s ranks later in the file. '
                'Offending taxonomy string: %s' %
                (self.lineage_depth, len(lineage), lineage_str))
        return lineage

    def get_training_seqs(self):
        """Returns an iterator of valid training sequences in
        RDP-compatible format

        Each training sequence is represented by a tuple (rdp_id,
        seq).  The rdp_id consists of two items: the original sequence
        ID with whitespace replaced by underscores, and the lineage
        with taxa separated by semicolons.
        """
        # Rdp requires unique sequence IDs without whitespace.  Can't
        # trust user IDs to not have whitespace, so we replace all
        # whitespace with an underscore.  Classification may fail if
        # the replacement method generates a name collision.
        for seq_id, node in self.sequence_nodes.iteritems():
            seq = self.sequences.get(seq_id)
            if seq is not None:
                lineage = node.get_lineage()
                rdp_id = '%s %s' % (
                    re.sub('\s',
                           '_',
                           seq_id),
                    ';'.join(lineage))
                yield rdp_id, seq

    def get_rdp_taxonomy(self):
        return self._tree.get_rdp_taxonomy()

    def fix_output_file(self, result_path):
        # Ultimate hack to replace mangled taxa names
        temp_results = StringIO()
        for line in open(result_path):
            line = re.sub(
                _QIIME_RDP_TAXON_TAG + "[^;\n\t]*", '', line)
            for char, escape_str in _QIIME_RDP_ESCAPES:
                line = re.sub(escape_str, char, line)
            temp_results.write(line)
        open(result_path, 'w').write(temp_results.getvalue())

    def fix_results(self, results_dict):
        for seq_id, assignment in results_dict.iteritems():
            lineage, confidence = assignment
            lineage = re.sub(
                _QIIME_RDP_TAXON_TAG + "[^;\n\t]*", '', lineage)
            for char, escape_str in _QIIME_RDP_ESCAPES:
                lineage = re.sub(escape_str, char, lineage)
            results_dict[seq_id] = (lineage, confidence)
        return results_dict


class RdpTree(object):

    """Simple, specialized tree class used to generate a taxonomy
    file for the Rdp Classifier.
    """
    taxonomic_ranks = ' abcdefghijklmnopqrstuvwxyz'

    def __init__(self, name='Root', parent=None, counter=None):
        if counter is None:
            self.counter = count(0)
        else:
            self.counter = counter
        self.id = self.counter.next()
        self.name = name
        self.parent = parent
        self.seq_ids = []
        if parent is None:
            self.depth = 0
        else:
            self.depth = parent.depth + 1
        self.children = dict()  # name => subtree

    def insert_lineage(self, lineage):
        """Inserts an assignment into the taxonomic tree.

        Lineage must support the iterator interface, or provide an
        __iter__() method that returns an iterator.
        """
        lineage = iter(lineage)
        try:
            taxon = lineage.next()
            if taxon not in self.children:
                self.children[taxon] = self.__class__(
                    name=taxon, parent=self, counter=self.counter)
            retval = self.children[taxon].insert_lineage(lineage)
        except StopIteration:
            retval = self
        return retval

    def get_lineage(self):
        if self.parent is not None:
            return self.parent.get_lineage() + [self.name]
        else:
            return [self.name]

    def get_nodes(self):
        yield self
        for child in self.children.values():
            child_nodes = child.get_nodes()
            for node in child_nodes:
                yield node

    def dereplicate_taxa(self):
        # We check that there are no duplicate taxon names (case insensitive)
        # at a given depth. We must do a case insensitive check because the RDP
        # classifier converts taxon names to lowercase when it checks for
        # duplicates, and will throw an error otherwise.
        taxa_by_depth = {}
        for node in self.get_nodes():
            name = node.name
            depth = node.depth
            current_names = taxa_by_depth.get(depth, set())
            if name.lower() in current_names:
                node.name = name + _QIIME_RDP_TAXON_TAG + str(node.id)
            else:
                current_names.add(name.lower())
                taxa_by_depth[depth] = current_names

    def get_rdp_taxonomy(self):
        """Returns a string, in Rdp-compatible format.
        """
        # RDP uses 0 for the parent ID of the root node
        if self.parent is None:
            parent_id = 0
        else:
            parent_id = self.parent.id

        # top rank name must be norank, and bottom rank must be genus
        if self.depth == 0:
            rank_name = "norank"
        elif self.children:
            rank_name = self.taxonomic_ranks[self.depth]
        else:
            rank_name = "genus"

        fields = [
            self.id, self.name, parent_id, self.depth, rank_name]
        taxonomy_str = '*'.join(map(str, fields)) + "\n"

        # Recursively append lines from sorted list of subtrees
        child_names = sorted(self.children.keys())
        subtrees = [self.children[name] for name in child_names]
        for subtree in subtrees:
            taxonomy_str += subtree.get_rdp_taxonomy()
        return taxonomy_str


_QIIME_RDP_TAXON_TAG = "_qiime_unique_taxon_tag_"
_QIIME_RDP_ESCAPES = [
    ("&", "_qiime_ampersand_escape_"),
    (">", "_qiime_greaterthan_escape_"),
    ("<", "_qiime_lessthan_escape_"),
]


class RtaxTaxonAssigner(TaxonAssigner):

    """Assign taxon using RTAX
    """
    Name = "RtaxTaxonAssigner"
    # ", version 0.98"  # don't hardcode the version number, as it may change, and then the log output test would fail
    Application = "RTAX classifier"
    Citation = "Soergel D.A.W., Dey N., Knight R., and Brenner S.E.  2012.  Selection of primers for optimal taxonomic classification of environmental 16S rRNA gene sequences.  ISME J (6), 1440-1444"
    _tracked_properties = ['Application', 'Citation']

    def __init__(self, params):
        """Return new RtaxTaxonAssigner object with specified params.
        """
        _params = {
            'id_to_taxonomy_fp': None,
            'reference_sequences_fp': None,
            # 'delimiter': ","
            # use the amplicon ID, not including /1 or /3, as the primary key
            # for the query sequences
            'header_id_regex': "\\S+\\s+(\\S+?)\/",
            # OTU clustering produces ">clusterID read_1_id"
            'read_id_regex': "\\S+\\s+(\\S+)",
            # split_libraries produces >read_1_id ampliconID/1 .   This makes a
            # map between read_1_id and ampliconID.
            'amplicon_id_regex': "(\\S+)\\s+(\\S+?)\/",
            'read_1_seqs_fp': None,
            'read_2_seqs_fp': None,
            'single_ok': False,
            'no_single_ok_generic': False
        }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for
        each seq.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, dumps the
            result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """

        if log_path:
            self.writeLog(log_path)

        reference_sequences_fp = self.Params['reference_sequences_fp']
        assert reference_sequences_fp, \
            "Must provide reference_sequences_fp when calling an RtaxTaxonAssigner."

        id_to_taxonomy_fp = self.Params['id_to_taxonomy_fp']
        assert id_to_taxonomy_fp, \
            "Must provide id_to_taxonomy_fp when calling an RtaxTaxonAssigner."

        # delimiter = self.Params['delimiter']
        read_1_seqs_fp = self.Params['read_1_seqs_fp']
        assert read_1_seqs_fp, \
            "Must provide read_1_seqs_fp when calling an RtaxTaxonAssigner."

        # following params may all be null

        read_2_seqs_fp = self.Params['read_2_seqs_fp']
        single_ok = self.Params['single_ok']
        no_single_ok_generic = self.Params['no_single_ok_generic']
        header_id_regex = self.Params['header_id_regex']
        assert header_id_regex, \
            "Must not provide empty header_id_regex when calling an RtaxTaxonAssigner; leave unset"\
            "to use default if in doubt."

        read_id_regex = self.Params['read_id_regex']
        amplicon_id_regex = self.Params['amplicon_id_regex']

        # seq_file = open(seq_path, 'r')

        results = rtax.assign_taxonomy(
            seq_path, reference_sequences_fp, id_to_taxonomy_fp,
            read_1_seqs_fp, read_2_seqs_fp, single_ok=single_ok, no_single_ok_generic=no_single_ok_generic,
            header_id_regex=header_id_regex, read_id_regex=read_id_regex,
            amplicon_id_regex=amplicon_id_regex, output_fp=result_path,
            log_path=log_path, base_tmp_dir=get_qiime_temp_dir())

        return results


class UclustConsensusTaxonAssigner(TaxonAssigner):

    """Assign taxonomy using uclust
    """
    Name = "UclustConsensusTaxonAssigner"
    Application = "uclust"
    Citation = """uclust citation: Search and clustering orders of magnitude faster than BLAST. Edgar RC. Bioinformatics. 2010 Oct 1;26(19):2460-1.

uclust-based consensus taxonomy assigner by Greg Caporaso, citation: QIIME allows analysis of high-throughput community sequencing data. Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Pena AG, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. Nat Methods. 2010 May;7(5):335-6.
"""

    def __init__(self, params):
        """Returns a new UclustConsensusTaxonAssigner object with specified params
        """
        _params = {
            # Required, mapping of reference sequence to taxonomy
            'id_to_taxonomy_fp': None,
            # Required, reference sequence fasta file
            'reference_sequences_fp': None,
            # max-accepts parameter, as passed to uclust
            'max_accepts': 3,
            # Fraction of sequence hits that a taxonomy assignment
            # must show up in to be considered the consensus assignment
            'min_consensus_fraction': 0.51,
            # minimum identity to consider a hit (passed to uclust as --id)
            'similarity': 0.90,
            # label to apply for queries that cannot be assigned
            'unassignable_label': 'Unassigned'
        }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)

        if self.Params['id_to_taxonomy_fp'] is None:
            raise ValueError(
                "id_to_taxonomy_fp must be provided when instantiating a UclustConsensusTaxonAssigner")
        if self.Params['reference_sequences_fp'] is None:
            raise ValueError(
                "reference_sequences_fp must be provided when instantiating a UclustConsensusTaxonAssigner")

        id_to_taxonomy_f = open(self.Params['id_to_taxonomy_fp'], 'U')
        self.id_to_taxonomy = self._parse_id_to_taxonomy_file(id_to_taxonomy_f)

    def __call__(self,
                 seq_path,
                 result_path=None,
                 uc_path=None,
                 log_path=None,
                 HALT_EXEC=False):
        """Returns mapping of each seq to (tax, consensus fraction, n)

        Results:
        If result_path is specified, the results will be written to file
         as tab-separated lines of:
          query_id <tab> tax <tab> consensus fraction <tab> n
        If result_path is None (default), the results will be returned
         as a dict of:
          {'query_id': (tax, consensus fraction, n)}
        In both cases, the values are:
         tax: the consensus taxonomy assignment
         consensus fraction: the fraction of the assignments for the
          query that contained the lowest level tax assignment that is
          included in tax (e.g., if the assignment goes to genus level,
          this will be the fraction of assignments that had the consensus
          genus assignment)
         n: the number of assignments that were considered when constructing
          the consensus

        Parameters:
        seq_path: path to file of query sequences
        result_path: path where results should be written. If None (default),
         returns results as a dict
        uc_path: path where .uc file should be saved. If None (default), and
         log_path is specified, the .uc contents will be written to appended to
         the log file.
        log_path: path where run log should be written. If None (default), no
         log file is written.
        HALT_EXEC: debugging paramter. If pass, will exit just before the
         uclust command is issued, and will print the command that would have
         been called to stdout.
        """

        # initialize the logger
        logger = self._get_logger(log_path)
        logger.info(str(self))

        # set the user-defined parameters
        params = {'--id': self.Params['similarity'],
                  '--maxaccepts': self.Params['max_accepts']}

        # initialize the application controller object
        app = Uclust(params,
                     TmpDir=get_qiime_temp_dir(),
                     HALT_EXEC=HALT_EXEC)

        # Configure for consensus taxonomy assignment
        app.Parameters['--rev'].on()
        app.Parameters['--lib'].on(self.Params['reference_sequences_fp'])
        app.Parameters['--libonly'].on()
        app.Parameters['--allhits'].on()

        if uc_path is None:
            uc = NamedTemporaryFile(prefix='UclustConsensusTaxonAssigner_',
                                    suffix='.uc',
                                    dir=get_qiime_temp_dir())
            uc_path = uc.name
            store_uc_in_log = True
        else:
            store_uc_in_log = False

        app_result = app({'--input': seq_path,
                          '--uc': uc_path})
        # get map of query id to all assignments
        result = self._uc_to_assignments(app_result['ClusterFile'])
        # get consensus assignment
        query_to_assignments = self._tax_assignments_to_consensus_assignments(result)
        if result_path is not None:
            # if the user provided a result_path, write the
            # results to file
            of = open(result_path, 'w')
            for seq_id, (assignment, consensus_fraction, n) in query_to_assignments.iteritems():
                assignment_str = ';'.join(assignment)
                of.write('%s\t%s\t%1.2f\t%d\n' %
                         (seq_id, assignment_str, consensus_fraction, n))
            of.close()
            result = None
            logger.info('Result path: %s' % result_path)
        else:
            # If no result_path was provided, the result dict is
            # returned as-is.
            logger.info('Result path: None, returned as dict.')

        if store_uc_in_log:
            # This is a little hackish, but we don't have a good way
            # to pass the uc_path value right now through the
            # assign_taxonomy.py script, so writing the contents to the
            # user-specified log file (since this is being stored for logging
            # purposes).
            app_result['ClusterFile'].seek(0)
            logger.info('\n.uc file contents:\n')
            for line in app_result['ClusterFile']:
                logger.info(line.strip())

        return result

    def _get_logger(self, log_path=None):
        if log_path is not None:
            handler = logging.FileHandler(log_path, mode='w')
        else:
            class NullHandler(logging.Handler):

                def emit(self, record):
                    pass
            handler = NullHandler()
        logger = logging.getLogger("UclustConsensusTaxonAssigner logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _uc_to_assignments(self, uc):
        """ return dict mapping query id to all taxonomy assignments
        """
        results = defaultdict(list)
        for line in uc:
            line = line.strip()
            if line.startswith('#') or line == "":
                continue
            elif line.startswith('H'):
                fields = line.split('\t')
                query_id = fields[8].split()[0]
                subject_id = fields[9].split()[0]
                tax = self.id_to_taxonomy[subject_id].split(';')
                results[query_id].append(tax)
            elif line.startswith('N'):
                fields = line.split('\t')
                query_id = fields[8].split()[0]
                results[query_id].append([])
        return results
