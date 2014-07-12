#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import abspath, basename, splitext, join
from shutil import copyfile
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path
from skbio.parse.sequences import parse_fasta

from qiime.identify_chimeric_seqs import make_cidx_file
from qiime.util import write_degapped_fasta_to_file
from qiime.parallel.util import (ParallelWrapper, input_fasta_splitter,
                                 concatenate_files)
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp


class ParallelChimericSequenceIdentifier(ParallelWrapper):

    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        self._job_graph = nx.DiGraph()

        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        if jobs_to_start is None:
            # default to the number of workers
            jobs_to_start = context.get_number_of_workers()
        method = params['chimera_detection_method']
        self._log_file = generate_log_fp(output_dir)

        # Get a temp folder
        working_dir = mkdtemp(prefix='CSI_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        if method == 'blast_fragments':
            ref_seqs_fp = abspath(params['reference_seqs_fp'])
            # We first need to create the blast database
            handler = context.submit_async(build_blast_db_from_fasta_path,
                                           ref_seqs_fp, output_dir=working_dir)
            blast_db, db_files_to_remove = handler.get()
            self._filepaths_to_remove.extend(db_files_to_remove)
            params['blast_db'] = blast_db
            # Split the input_fp in multiple files
            handler = context.submit_async(input_fasta_splitter, input_fp,
                                           working_dir, jobs_to_start)
            fasta_fps = handler.get()
            self._filepaths_to_remove.extend(fasta_fps)
            # Now we create all the commands
            temp_out_fps = []
            ics_nodes = []
            for i, fasta_fp in enumerate(fasta_fps):
                # Create the output path
                temp_out_fp = join(
                    working_dir,
                    "%s_%d.txt" % (splitext(basename(fasta_fp))[0], i))
                temp_out_fps.append(temp_out_fp)
                cmd = (
                    "identify_chimeric_seqs.py -i %s -t %s -m blast_fragments "
                    "-o %s -n %s -d %s -e %s -b %s"
                    % (fasta_fp, params['id_to_taxonomy_fp'], temp_out_fp,
                       params['num_fragments'], params['taxonomy_depth'],
                       params['max_e_value'], params['blast_db']))
                node_name = "ICS_%d" % i
                ics_nodes.append(node_name)
                self._job_graph.add_node(node_name, job=(cmd,))

            # Adding temp_out_fps to the clean-up variable
            self._filepaths_to_remove.extend(temp_out_fps)
            # Merge the results by concatenating the output files
            output_fp = join(output_dir, "chimeric_seqs.txt")
            self._job_graph.add_node("CONCAT", job=(concatenate_files,
                                                    output_fp, temp_out_fps))
            # Make sure that the "CONCAT" node is the latest executed node
            for node in ics_nodes:
                self._job_graph.add_edge(node, "CONCAT")

        elif method == 'ChimeraSlayer':
            # We first need to copy the reference files to the working
            # directory because ChimeraSlayer creates an index file of the
            # reference and it will crash without write permission in the
            # reference sequences directory
            aln_ref_seqs_fp = abspath(params['aligned_reference_seqs_fp'])
            dest_aln_ref_seqs_fp = join(working_dir, basename(aln_ref_seqs_fp))
            copyfile(aln_ref_seqs_fp, dest_aln_ref_seqs_fp)
            aln_ref_seqs_fp = dest_aln_ref_seqs_fp
            # Adding the copy to the clean up variable
            self._filepaths_to_remove.append(aln_ref_seqs_fp)

            ref_seqs_fp = params['reference_seqs_fp']
            if ref_seqs_fp:
                ref_seqs_fp = abspath(ref_seqs_fp)
                # Copy the unaligned reference database
                dest_ref_seqs_fp = join(working_dir, basename(ref_seqs_fp))
                copyfile(ref_seqs_fp, dest_ref_seqs_fp)
                ref_seqs_fp = dest_ref_seqs_fp
            else:
                # The unaligned reference database does not exists, create it
                with open(aln_ref_seqs_fp, 'U') as f:
                    ref_seqs_fp = write_degapped_fasta_to_file(
                        parse_fasta(f), tmp_dir=working_dir)

            # Add the copied/new ref_seqs_fp to the clean up variable
            self._filepaths_to_remove.append(ref_seqs_fp)

            # Build the blast db of reference, otherwise ChimeraSlayer will do
            # it and parallel jobs clash
            handler = context.submit_async(build_blast_db_from_fasta_path,
                                           ref_seqs_fp, output_dir=working_dir)
            blast_db, db_files_to_remove = handler.get()
            self._filepaths_to_remove.extend(db_files_to_remove)

            # Make the index file. ChimeraSlayer first checks to see if the
            # index file exists. If not, it tries to create it. This can lead
            # to race conditions if several parallel jobs try to create it at
            # the same time
            self._job_graph.add_node("MAKE_CIDX", job=(make_cidx_file,
                                                       aln_ref_seqs_fp))
            self._filepaths_to_remove.append("%s.cidx" % aln_ref_seqs_fp)

            # Split the input_fp in multiple files
            handler = context.submit_async(input_fasta_splitter, input_fp,
                                           working_dir, jobs_to_start)
            fasta_fps = handler.get()
            self._filepaths_to_remove.extend(fasta_fps)

            min_div_ratio_str = ("--min_div_ratio %s" % params['min_div_ratio']
                                 if params['min_div_ratio'] else "")

            # Create the commands
            temp_out_fps = []
            ics_nodes = []
            for i, fasta_fp in enumerate(fasta_fps):
                # Create the output path
                temp_out_fp = join(
                    working_dir,
                    "%s_%d.txt" % (splitext(basename(fasta_fp))[0], i))
                temp_out_fps.append(temp_out_fp)
                cmd = (
                    "identify_chimeric_seqs.py -i %s -a %s -m ChimeraSlayer "
                    "-o %s -r %s %s"
                    % (fasta_fp, aln_ref_seqs_fp, temp_out_fp,
                       ref_seqs_fp, min_div_ratio_str))
                node_name = "ICS_%d" % i
                ics_nodes.append(node_name)
                self._job_graph.add_node(node_name, job=(cmd,))
                # This should be executed after the "MAKE_CIDX" node
                self._job_graph.add_edge("MAKE_CIDX", node_name)

            # Adding temp_out_fps to the clean_up variable
            self._filepaths_to_remove.extend(temp_out_fps)
            # Merge the results by concatenating the output files
            output_fp = join(output_dir, "chimeric_seqs.txt")
            self._job_graph.add_node("CONCAT", job=(concatenate_files,
                                                    output_fp, temp_out_fps))
            # Make sure that the "CONCAT" node is the latest executed node
            for node in ics_nodes:
                self._job_graph.add_edge(node, "CONCAT")

        else:
            raise ValueError(
                "Unrecognized chimera detection method '%s'." % method)
