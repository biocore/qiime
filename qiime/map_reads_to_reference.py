#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import join, splitext, exists
from bfillings.blat import MinimalBlatParser9

from bfillings.blat import (assign_dna_reads_to_protein_database
                         as blat_assign_dna_reads_to_protein_database,
                         assign_dna_reads_to_dna_database
                         as blat_assign_dna_reads_to_dna_database)
from bfillings.usearch import (clusters_from_blast_uc_file,
                            assign_dna_reads_to_database
                            as usearch_assign_dna_reads_to_database)
from bfillings.bwa import (assign_dna_reads_to_dna_database
                        as bwa_assign_dna_reads_to_dna_database)

from qiime.format import format_observation_map
from qiime.parse import parse_taxonomy, MinimalSamParser
from qiime.make_otu_table import make_otu_table
from qiime.util import get_qiime_temp_dir, create_dir, write_biom_table


class DatabaseMapper(object):

    def __call__(self,
                 query_fasta_fp,
                 database_fasta_fp,
                 output_dir,
                 observation_metadata_fp=None,
                 params=None,
                 HALT_EXEC=False):

        if params is None:
            params = {}

        """ Call the DatabaseMapper """
        create_dir(output_dir)
        raw_output_fp = self._get_raw_output_fp(output_dir,
                                                params)
        output_observation_map_fp = '%s/observation_map.txt' % output_dir
        output_biom_fp = '%s/observation_table.biom' % output_dir
        log_fp = '%s/observation_table.log' % output_dir

        self._assign_dna_reads_to_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            raw_output_fp=raw_output_fp,
            temp_dir=get_qiime_temp_dir(),
            params=params,
            HALT_EXEC=HALT_EXEC)

        self._process_raw_output(raw_output_fp,
                                 log_fp,
                                 output_observation_map_fp)

        self._generate_biom_output(output_observation_map_fp,
                                   output_biom_fp,
                                   observation_metadata_fp)

    def _generate_biom_output(self,
                              observation_map_fp,
                              output_biom_fp,
                              observation_metadata_fp):
        if observation_metadata_fp is not None:
            observation_metadata = \
                parse_taxonomy(open(observation_metadata_fp, 'U'))
        else:
            observation_metadata = None

        biom_table = make_otu_table(open(observation_map_fp, 'U'),
                                    observation_metadata)
        write_biom_table(biom_table, output_biom_fp)

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      observation_metadata_fp,
                                      params,
                                      HALT_EXEC):
        raise NotImplementedError(
            "DatabaseMapper subclasses must override _assign_dna_reads_to_database")

    def _get_raw_output_fp(self, output_dir, params):
        """ Generate filepath for raw output

            subclasses will generally want to override this method

        """
        return join(output_dir, 'raw_output.txt')

    def _process_raw_output(self,
                            raw_output_fp,
                            log_fp,
                            output_observation_map_fp):
        raise NotImplementedError(
            "DatabaseMapper subclasses must override _process_raw_output")


class UsearchDatabaseMapper(DatabaseMapper):

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      temp_dir,
                                      params,
                                      HALT_EXEC):
        usearch_assign_dna_reads_to_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            output_fp=raw_output_fp,
            temp_dir=temp_dir,
            params=params,
            HALT_EXEC=HALT_EXEC)

    def _get_raw_output_fp(self,
                           output_dir,
                           params):
        """ Generate filepath for .uc file """
        return join(output_dir, 'out.uc')

    def _process_raw_output(self,
                            raw_output_fp,
                            log_fp,
                            output_observation_map_fp):
        """ Generate observation map and biom table from .uc file
        """
        hits, failures = clusters_from_blast_uc_file(
            open(raw_output_fp, 'U'), 9)
        observation_map_f = open(output_observation_map_fp, 'w')
        for line in format_observation_map(hits.items(), ''):
            observation_map_f.write(line)
        observation_map_f.close()


class BlatDatabaseMapper(DatabaseMapper):

    MaxEvalue = 1e-10
    MinId = 0.97

    def _get_raw_output_fp(self,
                           output_dir,
                           params):
        """ Generate filepath for .bl9 (blast9) file """
        return join(output_dir, 'out.bl9')

    def _process_raw_output(self,
                            raw_output_fp,
                            log_fp,
                            output_observation_map_fp):
        """ Generate observation map and biom table from .bl9 file
        """
        result = {}
        pct_id_field = 2
        evalue_field = 10
        output_observation_map_f = open(output_observation_map_fp, 'w')
        log_f = open(log_fp, 'w')
        for summary, blat_results in MinimalBlatParser9(
                open(raw_output_fp, 'U'),
                include_column_names=False):
            for e in blat_results:
                if (float(e[evalue_field]) <= self.MaxEvalue and
                        float(e[pct_id_field]) / 100. >= self.MinId):
                    query_id = e[0]
                    subject_id = e[1]
                    try:
                        result[subject_id].append(query_id)
                    except KeyError:
                        result[subject_id] = [query_id]
                    log_f.write('\t'.join(e))
                    log_f.write('\n')
                    break
        log_f.close()
        for e in result.items():
            output_observation_map_f.write(
                '%s\t%s\n' %
                (e[0], '\t'.join(e[1])))
        output_observation_map_f.close()
        return result

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      temp_dir,
                                      params,
                                      HALT_EXEC):
        blat_assign_dna_reads_to_protein_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            output_fp=raw_output_fp,
            temp_dir=temp_dir,
            params=params)


class BlatNtDatabaseMapper(BlatDatabaseMapper):

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      temp_dir,
                                      params,
                                      HALT_EXEC):
        blat_assign_dna_reads_to_dna_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            output_fp=raw_output_fp,
            params=params)


class BwaSwDatabaseMapper(DatabaseMapper):

    def _get_raw_output_fp(self,
                           output_dir,
                           params):
        """ Generate filepath for .bl9 (blast9) file """
        return join(output_dir, 'bwa_raw_out.sam')

    def _process_raw_output(self,
                            raw_output_fp,
                            log_fp,
                            output_observation_map_fp):
        """ Generate observation map and biom table from .bl9 file
        """
        result = {}
        query_id_field = 0
        flag_field = 1
        subject_id_field = 2
        output_observation_map_f = open(output_observation_map_fp, 'w')
        log_f = open(log_fp, 'w')
        for e in MinimalSamParser(open(raw_output_fp, 'U')):
            query_id = e[query_id_field]
            subject_id = e[subject_id_field]
            flag = int(e[flag_field])
            if (flag != 4):
                try:
                    result[subject_id].append(query_id)
                except KeyError:
                    result[subject_id] = [query_id]
                log_f.write('\t'.join(e))
                log_f.write('\n')

        log_f.close()
        for e in result.items():
            output_observation_map_f.write(
                '%s\t%s\n' %
                (e[0], '\t'.join(e[1])))
        output_observation_map_f.close()
        return result

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      temp_dir,
                                      params,
                                      HALT_EXEC):
        _params = {}
        _params.update(params)
        bwa_assign_dna_reads_to_dna_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            out_fp=raw_output_fp,
            params=_params)


class BwaShortDatabaseMapper(BwaSwDatabaseMapper):

    def _assign_dna_reads_to_database(self,
                                      query_fasta_fp,
                                      database_fasta_fp,
                                      raw_output_fp,
                                      temp_dir,
                                      params,
                                      HALT_EXEC):
        _aln_params = {'-f': splitext(raw_output_fp)[0] + '.sai'}
        if 'aln_params' in params:
            _aln_params.update(params['aln_params'])
        params['algorithm'] = 'bwa-short'
        params['aln_params'] = _aln_params
        bwa_assign_dna_reads_to_dna_database(
            query_fasta_fp=query_fasta_fp,
            database_fasta_fp=database_fasta_fp,
            out_fp=raw_output_fp,
            params=params)


def usearch_database_mapper(query_fp,
                            refseqs_fp,
                            output_dir,
                            evalue,
                            min_id,
                            queryalnfract,
                            targetalnfract,
                            maxaccepts,
                            maxrejects,
                            observation_metadata_fp=None,
                            HALT_EXEC=False):

        params = {}
        params['--evalue'] = evalue
        params['--id'] = min_id
        params['--queryalnfract'] = queryalnfract
        params['--targetalnfract'] = targetalnfract
        params['--maxaccepts'] = maxaccepts
        params['--maxrejects'] = maxrejects

        usearch_db_mapper = UsearchDatabaseMapper()
        usearch_db_mapper(query_fp,
                          refseqs_fp,
                          output_dir,
                          params=params,
                          observation_metadata_fp=observation_metadata_fp,
                          HALT_EXEC=HALT_EXEC)


def blat_database_mapper(query_fp,
                         refseqs_fp,
                         output_dir,
                         evalue,
                         min_id,
                         genetic_code,
                         observation_metadata_fp=None,
                         HALT_EXEC=False):

    params = {'-minIdentity': min_id,
              'genetic_code': genetic_code}

    blat_db_mapper = BlatDatabaseMapper()
    blat_db_mapper.MinId = min_id
    blat_db_mapper.MaxEvalue = evalue
    blat_db_mapper(query_fp,
                   refseqs_fp,
                   output_dir,
                   params=params,
                   observation_metadata_fp=observation_metadata_fp,
                   HALT_EXEC=HALT_EXEC)


def blat_nt_database_mapper(query_fp,
                            refseqs_fp,
                            output_dir,
                            evalue,
                            min_id,
                            observation_metadata_fp=None,
                            HALT_EXEC=False):

    params = {'-minIdentity': min_id}

    blat_db_mapper = BlatNtDatabaseMapper()
    blat_db_mapper.MinId = min_id
    blat_db_mapper.MaxEvalue = evalue
    blat_db_mapper(query_fp,
                   refseqs_fp,
                   output_dir,
                   params=params,
                   observation_metadata_fp=observation_metadata_fp,
                   HALT_EXEC=HALT_EXEC)


def bwa_sw_database_mapper(query_fp,
                           refseqs_fp,
                           output_dir,
                           observation_metadata_fp=None,
                           HALT_EXEC=False):

    bwa_db_mapper = BwaSwDatabaseMapper()
    params = {}
    bwa_db_mapper(query_fp,
                  refseqs_fp,
                  output_dir,
                  params=params,
                  observation_metadata_fp=observation_metadata_fp,
                  HALT_EXEC=HALT_EXEC)


def bwa_short_database_mapper(query_fp,
                              refseqs_fp,
                              output_dir,
                              max_diff,
                              observation_metadata_fp=None,
                              HALT_EXEC=False):

    bwa_db_mapper = BwaShortDatabaseMapper()
    if max_diff is not None:
        params = {'aln_params': {'-n': max_diff}}
    else:
        params = {}
    bwa_db_mapper(query_fp,
                  refseqs_fp,
                  output_dir,
                  params=params,
                  observation_metadata_fp=observation_metadata_fp,
                  HALT_EXEC=HALT_EXEC)
