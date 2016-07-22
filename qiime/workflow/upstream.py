#!/usr/bin/env python
# File created on 20 Feb 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Justin Kuczynski",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import split, splitext, abspath
from qiime.util import create_dir
from qiime.workflow.util import (print_to_stdout,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 log_input_md5s,
                                 get_params_str)


def run_pick_de_novo_otus(input_fp,
                          output_dir,
                          command_handler,
                          params,
                          qiime_config,
                          parallel=False,
                          logger=None,
                          suppress_md5=False,
                          status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime

        The steps performed by this function are:
          1) Pick OTUs;
          2) Pick a representative set;
          3) Align the representative set;
          4) Assign taxonomy;
          5) Filter the alignment prior to tree building - remove positions
             which are all gaps, and specified as 0 in the lanemask
          6) Build a phylogenetic tree;
          7) Build an OTU table.

    """

    # Prepare some variables for the later steps
    input_dir, input_filename = split(input_fp)
    input_basename, input_ext = splitext(input_filename)
    create_dir(output_dir)
    commands = []
    cluster_failures = False
    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [input_fp])

    # Prep the OTU picking command
    try:
        otu_picking_method = params['pick_otus']['otu_picking_method']
    except KeyError:
        otu_picking_method = 'uclust'
    pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir, input_basename)
    if parallel and (otu_picking_method == 'blast' or
                     otu_picking_method == 'uclust_ref'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''

        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --otu_picking_method
            # option. This works for now though.
            d = params['pick_otus'].copy()
            del d['otu_picking_method']
        except KeyError:
            pass

        if otu_picking_method == 'uclust_ref':
            try:
                suppress_new_clusters = d['suppress_new_clusters']
                del d['suppress_new_clusters']
                cluster_failures = False
            except KeyError:
                cluster_failures = True
                failure_otu_picking_method = 'uclust'

        params_str += ' %s' % get_params_str(d)
        otu_picking_script = 'parallel_pick_otus_%s.py' % otu_picking_method
        # Build the OTU picking command
        pick_otus_cmd = '%s -i %s -o %s -T %s' % (otu_picking_script,
                                                  input_fp,
                                                  pick_otu_dir,
                                                  params_str)
    else:
        try:
            params_str = get_params_str(params['pick_otus'])
        except KeyError:
            params_str = ''
        # Build the OTU picking command
        pick_otus_cmd = 'pick_otus.py -i %s -o %s %s' %\
            (input_fp, pick_otu_dir, params_str)

    commands.append([('Pick OTUs', pick_otus_cmd)])

    if cluster_failures:
        reference_otu_fp = otu_fp
        clustered_failures_dir = '%s/failure_otus/' % pick_otu_dir

        try:
            d = params['pick_otus'].copy()
            del d['otu_picking_method']
        except KeyError:
            pass

        if 'denovo_otu_id_prefix' not in d:
            d['denovo_otu_id_prefix'] = 'DeNovoOTU'
        params_str = ' %s' % get_params_str(d)

        failures_list_fp = '%s/%s_failures.txt' % \
            (pick_otu_dir, input_basename)
        failures_fasta_fp = '%s/%s_failures.fasta' % \
            (pick_otu_dir, input_basename)

        filter_fasta_cmd = 'filter_fasta.py -f %s -s %s -o %s' %\
            (input_fp, failures_list_fp, failures_fasta_fp)

        commands.append([('Generate failures fasta file',
                          filter_fasta_cmd)])

        # Prep the OTU picking command for
        failure_otu_fp = '%s/%s_failures_otus.txt' % (clustered_failures_dir,
                                                      input_basename)
        # Build the OTU picking command
        pick_otus_cmd = 'pick_otus.py -i %s -o %s -m %s %s' %\
            (failures_fasta_fp, clustered_failures_dir,
             failure_otu_picking_method, params_str)

        commands.append(
            [('Pick de novo OTUs for new clusters', pick_otus_cmd)])

        merged_otu_map_fp = '%s/merged_otu_map.txt' % clustered_failures_dir
        cat_otu_tables_cmd = 'cat %s %s >> %s' %\
            (reference_otu_fp, failure_otu_fp, merged_otu_map_fp)
        commands.append([('Merge OTU maps', cat_otu_tables_cmd)])
        otu_fp = merged_otu_map_fp

    # Prep the representative set picking command
    rep_set_dir = '%s/rep_set/' % output_dir
    create_dir(rep_set_dir)
    rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir, input_basename)
    rep_set_log_fp = '%s/%s_rep_set.log' % (rep_set_dir, input_basename)

    try:
        params_str = get_params_str(params['pick_rep_set'])
    except KeyError:
        params_str = ''
    # Build the representative set picking command
    pick_rep_set_cmd = 'pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
        (otu_fp, input_fp, rep_set_log_fp, rep_set_fp, params_str)
    commands.append([('Pick representative set', pick_rep_set_cmd)])

    # Prep the taxonomy assignment command
    try:
        assignment_method = params['assign_taxonomy']['assignment_method']
    except KeyError:
        assignment_method = 'uclust'
    assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
        (output_dir, assignment_method)
    taxonomy_fp = '%s/%s_rep_set_tax_assignments.txt' % \
        (assign_taxonomy_dir, input_basename)
    if parallel and (assignment_method == 'rdp' or
                     assignment_method == 'blast' or
                     assignment_method == 'uclust'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''

        # Grab the taxonomy assignment parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --assignment_method
            # option. This works for now though.
            d = params['assign_taxonomy'].copy()
            if 'assignment_method' in d:
                del d['assignment_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass

        # Build the parallel taxonomy assignment command
        assign_taxonomy_cmd = \
            'parallel_assign_taxonomy_%s.py -i %s -o %s -T %s' %\
            (assignment_method, rep_set_fp, assign_taxonomy_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['assign_taxonomy'])
        except KeyError:
            params_str = ''
        # Build the taxonomy assignment command
        assign_taxonomy_cmd = 'assign_taxonomy.py -o %s -i %s %s' %\
            (assign_taxonomy_dir, rep_set_fp, params_str)

    commands.append([('Assign taxonomy', assign_taxonomy_cmd)])

    # Prep the OTU table building command
    otu_table_fp = '%s/otu_table.biom' % output_dir
    try:
        params_str = get_params_str(params['make_otu_table'])
    except KeyError:
        params_str = ''
    # Build the OTU table building command
    make_otu_table_cmd = 'make_otu_table.py -i %s -t %s -o %s %s' %\
        (otu_fp, taxonomy_fp, otu_table_fp, params_str)

    commands.append([('Make OTU table', make_otu_table_cmd)])

    if cluster_failures:
        reference_otu_table_fp = '%s/reference_only_otu_table.biom' % output_dir
        # Build the OTU table building command
        make_otu_table_cmd = 'make_otu_table.py -i %s -t %s -o %s %s' %\
            (reference_otu_fp, taxonomy_fp, reference_otu_table_fp, params_str)

        commands.append(
            [('Make reference-only OTU table', make_otu_table_cmd)])

    # Prep the pynast alignment command
    try:
        alignment_method = params['align_seqs']['alignment_method']
    except KeyError:
        alignment_method = 'pynast'
    pynast_dir = '%s/%s_aligned_seqs' % (output_dir, alignment_method)
    aln_fp = '%s/%s_rep_set_aligned.fasta' % (pynast_dir, input_basename)
    if parallel and alignment_method == 'pynast':
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''

        # Grab the alignment parameters
        # Want to find a cleaner strategy for this: the parallel script
        # is method-specific, so doesn't take a --alignment_method
        # option. This works for now though.
        try:
            d = params['align_seqs'].copy()
        except KeyError:
            d = {}
        try:
            del d['alignment_method']
        except KeyError:
            pass
        params_str += ' %s' % get_params_str(d)

        # Build the parallel pynast alignment command
        align_seqs_cmd = 'parallel_align_seqs_pynast.py -i %s -o %s -T %s' %\
            (rep_set_fp, pynast_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['align_seqs'])
        except KeyError:
            params_str = ''
        # Build the pynast alignment command
        align_seqs_cmd = 'align_seqs.py -i %s -o %s %s' %\
            (rep_set_fp, pynast_dir, params_str)
    commands.append([('Align sequences', align_seqs_cmd)])

    # Prep the alignment filtering command
    filtered_aln_fp = '%s/%s_rep_set_aligned_pfiltered.fasta' %\
        (pynast_dir, input_basename)
    try:
        params_str = get_params_str(params['filter_alignment'])
    except KeyError:
        params_str = ''
    # Build the alignment filtering command
    filter_alignment_cmd = 'filter_alignment.py -o %s -i %s %s' %\
        (pynast_dir, aln_fp, params_str)
    commands.append([('Filter alignment', filter_alignment_cmd)])

    # Prep the tree building command
    tree_fp = '%s/rep_set.tre' % output_dir
    try:
        params_str = get_params_str(params['make_phylogeny'])
    except KeyError:
        params_str = ''
    # Build the tree building command
    make_phylogeny_cmd = 'make_phylogeny.py -i %s -o %s %s' %\
        (filtered_aln_fp, tree_fp, params_str)
    commands.append([('Build phylogenetic tree', make_phylogeny_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)

    return abspath(tree_fp), abspath(otu_table_fp)
run_qiime_data_preparation = run_pick_otus_through_otu_table = run_pick_de_novo_otus


def run_pick_closed_reference_otus(
        input_fp,
        refseqs_fp,
        output_dir,
        taxonomy_fp,
        command_handler,
        params,
        qiime_config,
        assign_taxonomy=False,
        parallel=False,
        logger=None,
        suppress_md5=False,
        status_update_callback=print_to_stdout):
    """ Run the data preparation steps of Qiime

        The steps performed by this function are:
          1) Pick OTUs;
          2) If assignment_taxonomy is True, choose representative sequence
             for OTUs and assign taxonomy using a classifier.
          3) Build an OTU table with optional predefined taxonomy
             (if assign_taxonomy=False) or taxonomic assignments from step 2
             (if assign_taxonomy=True).

    """

    # confirm that a valid otu picking method was supplied before doing
    # any work
    reference_otu_picking_methods = ['blast', 'uclust_ref', 'usearch61_ref',
                                     'usearch_ref', 'sortmerna']

    try:
        otu_picking_method = params['pick_otus']['otu_picking_method']
    except KeyError:
        otu_picking_method = 'uclust_ref'
    assert otu_picking_method in reference_otu_picking_methods,\
        "Invalid OTU picking method supplied: %s. Valid choices are: %s"\
        % (otu_picking_method, ' '.join(reference_otu_picking_methods))

    # Prepare some variables for the later steps
    input_dir, input_filename = split(input_fp)
    input_basename, input_ext = splitext(input_filename)
    create_dir(output_dir)
    commands = []
    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    if not suppress_md5:
        log_input_md5s(logger, [input_fp, refseqs_fp, taxonomy_fp])

    # Prep the OTU picking command
    pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
    otu_fp = '%s/%s_otus.txt' % (pick_otu_dir, input_basename)
    if parallel and (otu_picking_method == 'blast' or
                     otu_picking_method == 'uclust_ref' or
                     otu_picking_method == 'usearch61_ref' or
                     otu_picking_method == 'sortmerna'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''

        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --alignment_method
            # option. This works for now though.
            d = params['pick_otus'].copy()
            if 'otu_picking_method' in d:
                del d['otu_picking_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass
        otu_picking_script = 'parallel_pick_otus_%s.py' % otu_picking_method
        # Build the OTU picking command
        pick_otus_cmd = '%s -i %s -o %s -r %s -T %s' %\
            (otu_picking_script,
             input_fp,
             pick_otu_dir,
             refseqs_fp,
             params_str)
    else:
        try:
            params_str = get_params_str(params['pick_otus'])
        except KeyError:
            params_str = ''
        # Since this is reference-based OTU picking we always want to
        # suppress new clusters -- force it here.
        params_str += ' --suppress_new_clusters'
        logger.write(
            "Forcing --suppress_new_clusters as this is "
            "closed-reference OTU picking.\n\n")
        # Build the OTU picking command
        pick_otus_cmd = 'pick_otus.py -i %s -o %s -r %s -m %s %s' %\
            (input_fp,
             pick_otu_dir,
             refseqs_fp,
             otu_picking_method,
             params_str)

    commands.append([('Pick OTUs', pick_otus_cmd)])

    # Assign taxonomy using a taxonomy classifier, if request by the user.
    # (Alternatively predefined taxonomic assignments will be used, if provided.)
    if assign_taxonomy:
        # Prep the representative set picking command
        rep_set_dir = '%s/rep_set/' % output_dir
        create_dir(rep_set_dir)
        rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir, input_basename)
        rep_set_log_fp = '%s/%s_rep_set.log' % (rep_set_dir, input_basename)

        try:
            params_str = get_params_str(params['pick_rep_set'])
        except KeyError:
            params_str = ''
        # Build the representative set picking command
        pick_rep_set_cmd = 'pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
            (otu_fp, input_fp, rep_set_log_fp, rep_set_fp, params_str)
        commands.append([('Pick representative set', pick_rep_set_cmd)])

        # Prep the taxonomy assignment command
        try:
            assignment_method = params['assign_taxonomy']['assignment_method']
        except KeyError:
            assignment_method = 'uclust'
        assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
            (output_dir, assignment_method)
        taxonomy_fp = '%s/%s_rep_set_tax_assignments.txt' % \
            (assign_taxonomy_dir, input_basename)
        if parallel and (assignment_method == 'rdp' or
                         assignment_method == 'blast' or
                         assignment_method == 'uclust'):
            # Grab the parallel-specific parameters
            try:
                params_str = get_params_str(params['parallel'])
            except KeyError:
                params_str = ''

            # Grab the taxonomy assignment parameters
            try:
                # Want to find a cleaner strategy for this: the parallel script
                # is method-specific, so doesn't take a --assignment_method
                # option. This works for now though.
                d = params['assign_taxonomy'].copy()
                if 'assignment_method' in d:
                    del d['assignment_method']
                params_str += ' %s' % get_params_str(d)
            except KeyError:
                pass

            # Build the parallel taxonomy assignment command
            assign_taxonomy_cmd = \
                'parallel_assign_taxonomy_%s.py -i %s -o %s -T %s' %\
                (assignment_method, rep_set_fp, assign_taxonomy_dir, params_str)
        else:
            try:
                params_str = get_params_str(params['assign_taxonomy'])
            except KeyError:
                params_str = ''
            # Build the taxonomy assignment command
            assign_taxonomy_cmd = 'assign_taxonomy.py -o %s -i %s %s' %\
                (assign_taxonomy_dir, rep_set_fp, params_str)

        commands.append([('Assign taxonomy', assign_taxonomy_cmd)])

    # Prep the OTU table building command
    otu_table_fp = '%s/otu_table.biom' % output_dir
    try:
        params_str = get_params_str(params['make_otu_table'])
    except KeyError:
        params_str = ''
    # If assign_taxonomy is True, this will be the path to the taxonomic
    # assignment results. If assign_taxonomy is False this will be either
    # the precomputed taxonomic assignments that the user passed in,
    # or None.
    if taxonomy_fp:
        taxonomy_str = '-t %s' % taxonomy_fp
    else:
        taxonomy_str = ''
    # Build the OTU table building command
    make_otu_table_cmd = 'make_otu_table.py -i %s %s -o %s %s' %\
        (otu_fp, taxonomy_str, otu_table_fp, params_str)

    commands.append([('Make OTU table', make_otu_table_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)

run_pick_reference_otus_through_otu_table = run_pick_closed_reference_otus
