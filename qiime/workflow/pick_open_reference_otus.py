#!/usr/bin/env python
# File created on 21 Mar 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import split, splitext, getsize, exists, abspath, join
from shutil import copyfile, rmtree
from numpy import inf
from copy import deepcopy
from skbio.util import create_dir, remove_files
from skbio.parse.sequences import parse_fasta
from biom import load_table

from qiime.util import (subsample_fasta, count_seqs_from_file)
from qiime.filter import (filter_otus_from_otu_table,
                          get_seq_ids_from_fasta_file,
                          filter_otus_from_otu_map)
from qiime.workflow.util import (print_to_stdout,
                                 WorkflowLogger,
                                 generate_log_fp,
                                 log_input_md5s,
                                 get_params_str,
                                 WorkflowError)
from qiime.util import write_biom_table
from qiime.workflow.core_diversity_analyses import (format_index_link,
                                                    generate_index_page,
                                                    _index_headers)

def final_repset_from_iteration_repsets(repset_fasta_fs):
    """
        The first observation of each otu is chosen as the representative -
         this ensures that the representative sequence is the centroid of
         the cluster.
    """
    observed = {}
    for repset_fasta_f in repset_fasta_fs:
        for otu_id, seq in parse_fasta(repset_fasta_f):
            o = otu_id.split()[0]
            if not o in observed:
                yield (otu_id, seq)
                observed[o] = None
            else:
                # we already have a representative for this otu id
                pass


def final_repset_from_iteration_repsets_fps(repset_fasta_fps, final_repset_fp):
    final_repset_f = open(final_repset_fp, 'w')
    repset_fasta_fs = map(open, repset_fasta_fps)
    for record in final_repset_from_iteration_repsets(repset_fasta_fs):
        final_repset_f.write('>%s\n%s\n' % record)
    final_repset_f.close()

#####################
# Start functions to port to new Qiime/qiime/workflow/util.py
#####################
# The following functions are currently all tested via
# the wrapper functions in PickSubsampledReferenceOtusThroughOtuTableTests.
# In an up-coming workflow-refactoring, I want to use these in other workflow
# scripts and test directly to simplify WorkflowTests. I split these out when
# writing this code as it became obvious that they're reusable.


def pick_reference_otus(input_fp,
                        output_dir,
                        otu_picking_method,
                        refseqs_fp,
                        parallel,
                        params,
                        logger,
                        similarity_override=None):
    params_copy = deepcopy(params)
    if 'pick_otus' in params_copy and 'refseqs_fp' in params_copy['pick_otus']:
        raise WorkflowError("Cannot pass pick_otus:refseqs_fp in parameters file. This can only be"
                            " passed on the command line or through the API.")
    if similarity_override is not None:
        logger.write(
            'Similiarity of %1.3f being used for pre-filtering.\n' %
            similarity_override)
        if 'pick_otus' in params_copy:
            params_copy['pick_otus']['similarity'] = str(similarity_override)
        else:
            params_copy['pick_otus'] = {'similarity': str(similarity_override)}

    if parallel and (otu_picking_method == 'uclust_ref' or otu_picking_method == "sortmerna"):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params_copy['parallel'])
        except KeyError:
            params_str = ''

        # Grab the OTU picker parameters
        try:
            # Want to find a cleaner strategy for this: the parallel script
            # is method-specific, so doesn't take a --otu_picking_method
            # option. This works for now though.
            if 'otu_picking_method' in params_copy['pick_otus']:
                del params_copy['pick_otus']['otu_picking_method']
        except KeyError:
            pass

        params_str += ' %s' % get_params_str(params_copy['pick_otus'])
        otu_picking_script = 'parallel_pick_otus_%s.py' % otu_picking_method
        # Build the OTU picking command
        pick_otus_cmd = '%s -i %s -o %s -r %s -T %s' %\
            (otu_picking_script,
             input_fp,
             output_dir,
             refseqs_fp,
             params_str)
    else:
        try:
            params_str = get_params_str(params_copy['pick_otus'])
        except KeyError:
            params_str = ''
        # Since this is reference-based OTU picking we always want to
        # suppress new clusters -- force it here.
        params_str += ' --suppress_new_clusters'
        logger.write(
            "Forcing --suppress_new_clusters as this is reference-based OTU picking.\n\n")
        # Build the OTU picking command
        pick_otus_cmd = 'pick_otus.py -i %s -o %s -r %s -m %s %s' %\
            (input_fp,
             output_dir,
             refseqs_fp,
             otu_picking_method,
             params_str)

    return pick_otus_cmd


def pick_denovo_otus(input_fp,
                     output_dir,
                     new_ref_set_id,
                     otu_picking_method,
                     params,
                     logger):
    try:
        d = params['pick_otus'].copy()
        del d['otu_picking_method']
    except KeyError:
        pass

    d['denovo_otu_id_prefix'] = '%s.ReferenceOTU' % new_ref_set_id

    params_str = ' %s' % get_params_str(d)
    # Build the OTU picking command
    result = 'pick_otus.py -i %s -o %s -m %s %s' %\
        (input_fp, output_dir, otu_picking_method, params_str)

    return result


def assign_tax(repset_fasta_fp,
               output_dir,
               command_handler,
               params,
               qiime_config,
               parallel=False,
               logger=None,
               status_update_callback=print_to_stdout):

    input_dir, input_filename = split(repset_fasta_fp)
    input_basename, input_ext = splitext(input_filename)
    commands = []
    if logger is None:
        log_fp = generate_log_fp(output_dir)
        logger = WorkflowLogger(log_fp,
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    # Prep the taxonomy assignment command
    try:
        assignment_method = params['assign_taxonomy']['assignment_method']
    except KeyError:
        assignment_method = 'uclust'
    assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
        (output_dir, assignment_method)
    taxonomy_fp = '%s/%s_tax_assignments.txt' % \
        (assign_taxonomy_dir, input_basename)
    if parallel and (assignment_method == 'rdp' or
                     assignment_method == 'blast' or
                     assignment_method == 'uclust'):
        # Grab the parallel-specific parameters
        try:
            params_str = get_params_str(params['parallel'])
        except KeyError:
            params_str = ''

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
            (assignment_method, repset_fasta_fp,
             assign_taxonomy_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['assign_taxonomy'])
        except KeyError:
            params_str = ''
        # Build the taxonomy assignment command
        assign_taxonomy_cmd = 'assign_taxonomy.py -o %s -i %s %s' %\
            (assign_taxonomy_dir, repset_fasta_fp, params_str)
    if exists(assign_taxonomy_dir):
        rmtree(assign_taxonomy_dir)
    commands.append([('Assign taxonomy', assign_taxonomy_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    return taxonomy_fp


def align_and_tree(repset_fasta_fp,
                   output_dir,
                   command_handler,
                   params,
                   qiime_config,
                   parallel=False,
                   logger=None,
                   status_update_callback=print_to_stdout):

    input_dir, input_filename = split(repset_fasta_fp)
    input_basename, input_ext = splitext(input_filename)
    commands = []
    if logger is None:
        log_fp = generate_log_fp(output_dir)
        logger = WorkflowLogger(log_fp,
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False


    # Prep the pynast alignment command
    alignment_method = 'pynast'
    pynast_dir = '%s/%s_aligned_seqs' % (output_dir, alignment_method)
    aln_fp = '%s/%s_aligned.fasta' % (pynast_dir, input_basename)
    failures_fp = '%s/%s_failures.fasta' % (pynast_dir, input_basename)
    if exists(pynast_dir):
        rmtree(pynast_dir)

    if parallel:
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
            d = params['align_seqs'].copy()
            if 'alignment_method' in d:
                del d['alignment_method']
            params_str += ' %s' % get_params_str(d)
        except KeyError:
            pass

        # Build the parallel pynast alignment command
        align_seqs_cmd = 'parallel_align_seqs_pynast.py -i %s -o %s -T %s' %\
            (repset_fasta_fp, pynast_dir, params_str)
    else:
        try:
            params_str = get_params_str(params['align_seqs'])
        except KeyError:
            params_str = ''
        # Build the pynast alignment command
        align_seqs_cmd = 'align_seqs.py -i %s -o %s %s' %\
            (repset_fasta_fp, pynast_dir, params_str)
    commands.append([('Align sequences', align_seqs_cmd)])

    # Prep the alignment filtering command
    filtered_aln_fp = '%s/%s_aligned_pfiltered.fasta' %\
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
    if exists(tree_fp):
        remove_files([tree_fp])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
    return failures_fp

#####################
# End functions to port to new Qiime/qiime/workflow/util.py
#####################


def iteration_output_exists(
        iteration_output_dir, min_otu_size, remove_partial_output=True):
    """  """
    if not exists(iteration_output_dir):
        return False

    expected_fps = ['%s/new_refseqs.fna' % iteration_output_dir,
                    '%s/rep_set.fna' % iteration_output_dir,
                    '%s/otu_table_mc%d.biom' % (iteration_output_dir, min_otu_size)]

    for fp in expected_fps:
        if not (exists(fp) and getsize(fp) > 0):
            if remove_partial_output:
                # if any of the expected filepaths don't exist or have
                # size == 0, remove the iteration output directory
                rmtree(iteration_output_dir)
            return False

    return True


def iterative_pick_subsampled_open_reference_otus(
        input_fps,
        refseqs_fp,
        output_dir,
        percent_subsample,
        new_ref_set_id,
        command_handler,
        params,
        qiime_config,
        prefilter_refseqs_fp=None,
        prefilter_percent_id=None,
        min_otu_size=2,
        run_assign_tax=True,
        run_align_and_tree=True,
        step1_otu_map_fp=None,
        step1_failures_fasta_fp=None,
        parallel=False,
        suppress_step4=False,
        logger=None,
        suppress_md5=False,
        denovo_otu_picking_method='uclust',
        reference_otu_picking_method='uclust_ref',
        status_update_callback=print_to_stdout,
        minimum_failure_threshold=100000):
    """ Call the pick_subsampled_open_reference_otus workflow on multiple inputs
         and handle processing of the results.
    """
    create_dir(output_dir)
    commands = []

    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False

    # if the user has not passed a different reference collection for the pre-filter,
    # used the input refseqs_fp for all iterations. we want to pre-filter all data against
    # the input data as lower percent identity searches with uclust can be slow, so we
    # want the reference collection to stay at a reasonable size.
    if prefilter_refseqs_fp is None:
        prefilter_refseqs_fp = refseqs_fp

    otu_table_fps = []
    repset_fasta_fps = []
    for i, input_fp in enumerate(input_fps):
        iteration_output_dir = '%s/%d/' % (output_dir, i)
        if iteration_output_exists(iteration_output_dir, min_otu_size):
            # if the output from an iteration already exists, skip that
            # iteration (useful for continuing failed runs)
            log_input_md5s(logger, [input_fp, refseqs_fp])
            logger.write('Iteration %d (input file: %s) output data already exists. '
                         'Skipping and moving to next.\n\n' % (i, input_fp))
        else:
            pick_subsampled_open_reference_otus(input_fp=input_fp,
                                                refseqs_fp=refseqs_fp,
                                                output_dir=iteration_output_dir,
                                                percent_subsample=percent_subsample,
                                                new_ref_set_id='.'.join(
                                                    [new_ref_set_id, str(i)]),
                                                command_handler=command_handler,
                                                params=params,
                                                qiime_config=qiime_config,
                                                run_assign_tax=False,
                                                run_align_and_tree=False,
                                                prefilter_refseqs_fp=prefilter_refseqs_fp,
                                                prefilter_percent_id=prefilter_percent_id,
                                                min_otu_size=min_otu_size,
                                                step1_otu_map_fp=step1_otu_map_fp,
                                                step1_failures_fasta_fp=step1_failures_fasta_fp,
                                                parallel=parallel,
                                                suppress_step4=suppress_step4,
                                                logger=logger,
                                                suppress_md5=suppress_md5,
                                                suppress_index_page=True,
                                                denovo_otu_picking_method=denovo_otu_picking_method,
                                                reference_otu_picking_method=reference_otu_picking_method,
                                                status_update_callback=status_update_callback,
                                                minimum_failure_threshold=minimum_failure_threshold)
        # perform post-iteration file shuffling whether the previous iteration's
        # data previously existed or was just computed.
        # step1 otu map and failures can only be used for the first iteration
        # as subsequent iterations need to use updated refseqs files
        step1_otu_map_fp = step1_failures_fasta_fp = None
        new_refseqs_fp = '%s/new_refseqs.fna' % iteration_output_dir
        refseqs_fp = new_refseqs_fp

        otu_table_fps.append(
            '%s/otu_table_mc%d.biom' %
            (iteration_output_dir, min_otu_size))

        repset_fasta_fps.append('%s/rep_set.fna' % iteration_output_dir)

    # Merge OTU tables - check for existence first as this step has historically
    # been a frequent failure, so is sometimes run manually in failed runs.
    otu_table_fp = '%s/otu_table_mc%d.biom' % (output_dir, min_otu_size)
    if not (exists(otu_table_fp) and getsize(otu_table_fp) > 0):
        merge_cmd = 'merge_otu_tables.py -i %s -o %s' %\
            (','.join(otu_table_fps), otu_table_fp)
        commands.append([("Merge OTU tables", merge_cmd)])

    # Build master rep set
    final_repset_fp = '%s/rep_set.fna' % output_dir
    final_repset_from_iteration_repsets_fps(repset_fasta_fps, final_repset_fp)

    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=False)
    commands = []

    # initialize output file names - these differ based on what combination of
    # taxonomy assignment and alignment/tree building is happening.
    if run_assign_tax and run_align_and_tree:
        tax_input_otu_table_fp = otu_table_fp
        otu_table_w_tax_fp = \
            '%s/otu_table_mc%d_w_tax.biom' % (output_dir, min_otu_size)
        align_and_tree_input_otu_table = otu_table_w_tax_fp
        pynast_failure_filtered_otu_table_fp = \
            '%s/otu_table_mc%d_w_tax_no_pynast_failures.biom' % (output_dir,
                                                                 min_otu_size)
    elif run_assign_tax:
        tax_input_otu_table_fp = otu_table_fp
        otu_table_w_tax_fp = \
            '%s/otu_table_mc%d_w_tax.biom' % (output_dir, min_otu_size)
    elif run_align_and_tree:
        align_and_tree_input_otu_table = otu_table_fp
        pynast_failure_filtered_otu_table_fp = \
            '%s/otu_table_mc%d_no_pynast_failures.biom' % (output_dir,
                                                           min_otu_size)

    if run_assign_tax:
        if exists(otu_table_w_tax_fp) and getsize(otu_table_w_tax_fp) > 0:
            logger.write(
                "Final output file exists (%s). Will not rebuild." %
                otu_table_w_tax_fp)
        else:
            # remove files from partially completed runs
            remove_files([otu_table_w_tax_fp], error_on_missing=False)

            taxonomy_fp = assign_tax(
                repset_fasta_fp=final_repset_fp,
                output_dir=output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                parallel=parallel,
                logger=logger,
                status_update_callback=status_update_callback)

            # Add taxa to otu table
            add_metadata_cmd = 'biom add-metadata -i %s --observation-metadata-fp %s -o %s --sc-separated taxonomy --observation-header OTUID,taxonomy' %\
                (tax_input_otu_table_fp, taxonomy_fp, otu_table_w_tax_fp)
            commands.append([("Add taxa to OTU table", add_metadata_cmd)])

            command_handler(commands,
                            status_update_callback,
                            logger=logger,
                            close_logger_on_success=False)
            commands = []

    if run_align_and_tree:
        if exists(pynast_failure_filtered_otu_table_fp) and\
           getsize(pynast_failure_filtered_otu_table_fp) > 0:
            logger.write("Final output file exists (%s). Will not rebuild." %
                         pynast_failure_filtered_otu_table_fp)
        else:
            # remove files from partially completed runs
            remove_files([pynast_failure_filtered_otu_table_fp],
                         error_on_missing=False)

            pynast_failures_fp = align_and_tree(
                repset_fasta_fp=final_repset_fp,
                output_dir=output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                parallel=parallel,
                logger=logger,
                status_update_callback=status_update_callback)

            # Build OTU table without PyNAST failures
            table = load_table(align_and_tree_input_otu_table)
            filtered_otu_table = filter_otus_from_otu_table(table,
                get_seq_ids_from_fasta_file(open(pynast_failures_fp, 'U')),
                0, inf, 0, inf, negate_ids_to_keep=True)
            write_biom_table(filtered_otu_table,
                             pynast_failure_filtered_otu_table_fp)

            command_handler(commands,
                            status_update_callback,
                            logger=logger,
                            close_logger_on_success=False)
            commands = []

    logger.close()


def pick_subsampled_open_reference_otus(input_fp,
                                        refseqs_fp,
                                        output_dir,
                                        percent_subsample,
                                        new_ref_set_id,
                                        command_handler,
                                        params,
                                        qiime_config,
                                        prefilter_refseqs_fp=None,
                                        run_assign_tax=True,
                                        run_align_and_tree=True,
                                        prefilter_percent_id=None,
                                        min_otu_size=2,
                                        step1_otu_map_fp=None,
                                        step1_failures_fasta_fp=None,
                                        parallel=False,
                                        suppress_step4=False,
                                        logger=None,
                                        suppress_md5=False,
                                        suppress_index_page=False,
                                        denovo_otu_picking_method='uclust',
                                        reference_otu_picking_method='uclust_ref',
                                        status_update_callback=print_to_stdout,
                                        minimum_failure_threshold=100000):
    """ Run the data preparation steps of Qiime

        The steps performed by this function are:
          - Pick reference OTUs against refseqs_fp
          - Subsample the failures to n sequences.
          - Pick OTUs de novo on the n failures.
          - Pick representative sequences for the resulting OTUs.
          - Pick reference OTUs on all failures using the
             representative set from step 4 as the reference set.

    """
    # for now only allowing uclust/usearch/sortmerna+sumaclust for otu picking
    allowed_denovo_otu_picking_methods = ['uclust', 'usearch61', 'sumaclust']
    allowed_reference_otu_picking_methods = ['uclust_ref', 'usearch61_ref',
                                             'sortmerna']
    assert denovo_otu_picking_method in allowed_denovo_otu_picking_methods,\
        "Unknown de novo OTU picking method: %s. Known methods are: %s"\
        % (denovo_otu_picking_method,
           ','.join(allowed_denovo_otu_picking_methods))

    assert reference_otu_picking_method in allowed_reference_otu_picking_methods,\
        "Unknown reference OTU picking method: %s. Known methods are: %s"\
        % (reference_otu_picking_method,
           ','.join(allowed_reference_otu_picking_methods))

    # Prepare some variables for the later steps
    index_links = []
    input_dir, input_filename = split(input_fp)
    input_basename, input_ext = splitext(input_filename)
    create_dir(output_dir)
    commands = []
    if logger is None:
        log_fp = generate_log_fp(output_dir)
        logger = WorkflowLogger(log_fp,
                                params=params,
                                qiime_config=qiime_config)

        close_logger_on_success = True
        index_links.append(
                ('Run summary data',
                log_fp,
                _index_headers['run_summary']))
    else:
        close_logger_on_success = False


    if not suppress_md5:
        log_input_md5s(logger, [input_fp,
                                refseqs_fp,
                                step1_otu_map_fp,
                                step1_failures_fasta_fp])

    # if the user has not passed a different reference collection for the pre-filter,
    # used the main refseqs_fp. this is useful if the user wants to provide a smaller
    # reference collection, or to use the input reference collection when running in
    # iterative mode (rather than an iteration's new refseqs)
    if prefilter_refseqs_fp is None:
        prefilter_refseqs_fp = refseqs_fp

    # Step 1: Closed-reference OTU picking on the input file (if not already
    # complete)
    if step1_otu_map_fp and step1_failures_fasta_fp:
        step1_dir = '%s/step1_otus' % output_dir
        create_dir(step1_dir)
        logger.write("Using pre-existing reference otu map and failures.\n\n")
    else:
        if prefilter_percent_id is not None:
            prefilter_dir = '%s/prefilter_otus/' % output_dir
            prefilter_failures_list_fp = '%s/%s_failures.txt' % \
                (prefilter_dir, input_basename)
            prefilter_pick_otu_cmd = pick_reference_otus(
                input_fp, prefilter_dir, reference_otu_picking_method,
                prefilter_refseqs_fp, parallel, params, logger, prefilter_percent_id)
            commands.append(
                [('Pick Reference OTUs (prefilter)', prefilter_pick_otu_cmd)])

            prefiltered_input_fp = '%s/prefiltered_%s%s' %\
                (prefilter_dir, input_basename, input_ext)
            filter_fasta_cmd = 'filter_fasta.py -f %s -o %s -s %s -n' %\
                (input_fp, prefiltered_input_fp, prefilter_failures_list_fp)
            commands.append(
                [('Filter prefilter failures from input', filter_fasta_cmd)])
            index_links.append(
            ('Pre-filtered sequence identifiers '
             '(failed to hit reference at %1.1f%% identity)' % (float(prefilter_percent_id)*100),
                        prefilter_failures_list_fp,
                        _index_headers['sequences']))


            # Call the command handler on the list of commands
            command_handler(commands,
                            status_update_callback,
                            logger=logger,
                            close_logger_on_success=False)
            commands = []

            input_fp = prefiltered_input_fp
            input_dir, input_filename = split(input_fp)
            input_basename, input_ext = splitext(input_filename)
            if getsize(prefiltered_input_fp) == 0:
                raise ValueError(
                    "All sequences were discarded by the prefilter. "
                    "Are the input sequences in the same orientation "
                    "in your input file and reference file (you can "
                    "add 'pick_otus:enable_rev_strand_match True' to "
                    "your parameters file if not)? Are you using the "
                    "correct reference file?")

        # Build the OTU picking command
        step1_dir = \
            '%s/step1_otus' % output_dir
        step1_otu_map_fp = \
            '%s/%s_otus.txt' % (step1_dir, input_basename)
        step1_pick_otu_cmd = pick_reference_otus(
            input_fp, step1_dir, reference_otu_picking_method,
            refseqs_fp, parallel, params, logger)
        commands.append([('Pick Reference OTUs', step1_pick_otu_cmd)])

        # Build the failures fasta file
        step1_failures_list_fp = '%s/%s_failures.txt' % \
            (step1_dir, input_basename)
        step1_failures_fasta_fp = \
            '%s/failures.fasta' % step1_dir
        step1_filter_fasta_cmd = 'filter_fasta.py -f %s -s %s -o %s' %\
            (input_fp, step1_failures_list_fp, step1_failures_fasta_fp)

        commands.append([('Generate full failures fasta file',
                          step1_filter_fasta_cmd)])

        # Call the command handler on the list of commands
        command_handler(commands,
                        status_update_callback,
                        logger=logger,
                        close_logger_on_success=False)
        commands = []

    step1_repset_fasta_fp = \
        '%s/step1_rep_set.fna' % step1_dir
    step1_pick_rep_set_cmd = 'pick_rep_set.py -i %s -o %s -f %s' %\
        (step1_otu_map_fp, step1_repset_fasta_fp, input_fp)
    commands.append([('Pick rep set', step1_pick_rep_set_cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=False)
    commands = []
    # name the final otu map
    merged_otu_map_fp = '%s/final_otu_map.txt' % output_dir

    # count number of sequences in step 1 failures fasta file
    with open(abspath(step1_failures_fasta_fp), 'U') as step1_failures_fasta_f:
        num_failure_seqs, mean, std = count_seqs_from_file(step1_failures_fasta_f)

    # number of failures sequences is greater than the threshold,
    # continue to step 2,3 and 4
    run_step_2_and_3 = num_failure_seqs > minimum_failure_threshold

    if run_step_2_and_3:

        # Subsample the failures fasta file to retain (roughly) the
        # percent_subsample
        step2_dir = '%s/step2_otus/' % output_dir
        create_dir(step2_dir)
        step2_input_fasta_fp = \
                               '%s/subsampled_failures.fasta' % step2_dir
        subsample_fasta(step1_failures_fasta_fp,
                        step2_input_fasta_fp,
                        percent_subsample)

        logger.write('# Subsample the failures fasta file using API \n' +
                 'python -c "import qiime; qiime.util.subsample_fasta' +
                 '(\'%s\', \'%s\', \'%f\')\n\n"' % (abspath(step1_failures_fasta_fp),
                                                    abspath(
                                                        step2_input_fasta_fp),
                                                    percent_subsample))

        # Prep the OTU picking command for the subsampled failures
        step2_cmd = pick_denovo_otus(step2_input_fasta_fp,
                                     step2_dir,
                                     new_ref_set_id,
                                     denovo_otu_picking_method,
                                     params,
                                     logger)
        step2_otu_map_fp = '%s/subsampled_failures_otus.txt' % step2_dir

        commands.append([('Pick de novo OTUs for new clusters', step2_cmd)])

        # Prep the rep set picking command for the subsampled failures
        step2_repset_fasta_fp = '%s/step2_rep_set.fna' % step2_dir
        step2_rep_set_cmd = 'pick_rep_set.py -i %s -o %s -f %s' %\
            (step2_otu_map_fp, step2_repset_fasta_fp, step2_input_fasta_fp)
        commands.append(
            [('Pick representative set for subsampled failures', step2_rep_set_cmd)])

        step3_dir = '%s/step3_otus/' % output_dir
        step3_otu_map_fp = '%s/failures_otus.txt' % step3_dir
        step3_failures_list_fp = '%s/failures_failures.txt' % step3_dir

        # remove the indexed reference database from the dictionary of
        # parameters as it must be forced to build a new database
        # using the step2_repset_fasta_fp
        if reference_otu_picking_method == 'sortmerna':
            if 'sortmerna_db' in params['pick_otus']:
                del params['pick_otus']['sortmerna_db']

        step3_cmd = pick_reference_otus(
            step1_failures_fasta_fp,
            step3_dir,
            reference_otu_picking_method,
            step2_repset_fasta_fp,
            parallel,
            params,
            logger)

        commands.append([
            ('Pick reference OTUs using de novo rep set', step3_cmd)])

        index_links.append(
            ('Final map of OTU identifier to sequence identifers (i.e., "OTU map")',
             merged_otu_map_fp,
             _index_headers['otu_maps']))

    if not suppress_step4:
        step4_dir = '%s/step4_otus/' % output_dir
        if run_step_2_and_3:
            step3_failures_fasta_fp = '%s/failures_failures.fasta' % step3_dir
            step3_filter_fasta_cmd = 'filter_fasta.py -f %s -s %s -o %s' %\
                (step1_failures_fasta_fp,
                 step3_failures_list_fp, step3_failures_fasta_fp)
            commands.append([('Create fasta file of step3 failures',
                            step3_filter_fasta_cmd)])

            failures_fp = step3_failures_fasta_fp
            failures_otus_fp = 'failures_failures_otus.txt'
            failures_step = 'step3'
        else:
            failures_fp = step1_failures_fasta_fp
            failures_otus_fp = 'failures_otus.txt'
            failures_step = 'step1'
            step3_otu_map_fp = ""

        step4_cmd = pick_denovo_otus(failures_fp,
                                     step4_dir,
                                     '.'.join([new_ref_set_id, 'CleanUp']),
                                     denovo_otu_picking_method,
                                     params,
                                     logger)

        step4_otu_map_fp = '%s/%s' % (step4_dir, failures_otus_fp)
        commands.append([('Pick de novo OTUs on %s failures' % failures_step, step4_cmd)])

        # Merge the otu maps, note that we are explicitly using the '>' operator
        # otherwise passing the --force flag on the script interface would
        # append the newly created maps to the map that was previously created
        cat_otu_tables_cmd = 'cat %s %s %s > %s' %\
            (step1_otu_map_fp, step3_otu_map_fp,
             step4_otu_map_fp, merged_otu_map_fp)
        commands.append([('Merge OTU maps', cat_otu_tables_cmd)])
        step4_repset_fasta_fp = '%s/step4_rep_set.fna' % step4_dir
        step4_rep_set_cmd = 'pick_rep_set.py -i %s -o %s -f %s' %\
            (step4_otu_map_fp, step4_repset_fasta_fp, failures_fp)
        commands.append(
            [('Pick representative set for subsampled failures', step4_rep_set_cmd)])
    else:
        # Merge the otu maps, note that we are explicitly using the '>' operator
        # otherwise passing the --force flag on the script interface would
        # append the newly created maps to the map that was previously created
        if run_step_2_and_3:
            failures_fp = step3_failures_list_fp
        else:
            failures_fp = step1_failures_list_fp
            step3_otu_map_fp = ""

        cat_otu_tables_cmd = 'cat %s %s > %s' %\
            (step1_otu_map_fp, step3_otu_map_fp, merged_otu_map_fp)
        commands.append([('Merge OTU maps', cat_otu_tables_cmd)])

        # Move the step 3 failures file to the top-level directory
        commands.append([('Move final failures file to top-level directory',
                          'mv %s %s/final_failures.txt' % (failures_fp, output_dir))])

    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=False)
    commands = []

    otu_fp = merged_otu_map_fp
    # Filter singletons from the otu map
    otu_no_singletons_fp = '%s/final_otu_map_mc%d.txt' % (output_dir,
                                                          min_otu_size)

    otus_to_keep = filter_otus_from_otu_map(
        otu_fp,
        otu_no_singletons_fp,
        min_otu_size)

    index_links.append(('Final map of OTU identifier to sequence identifers excluding '
                        'OTUs with fewer than %d sequences' % min_otu_size,
                        otu_no_singletons_fp,
                        _index_headers['otu_maps']))

    logger.write('# Filter singletons from the otu map using API \n' +
                 'python -c "import qiime; qiime.filter.filter_otus_from_otu_map' +
                 '(\'%s\', \'%s\', \'%d\')"\n\n' % (abspath(otu_fp),
                                                    abspath(
                                                        otu_no_singletons_fp),
                                                    min_otu_size))

    # make the final representative seqs file and a new refseqs file that
    # could be used in subsequent otu picking runs.
    # this is clunky. first, we need to do this without singletons to match
    # the otu map without singletons. next, there is a difference in what
    # we need the reference set to be and what we need the repseqs to be.
    # the reference set needs to be a superset of the input reference set
    # to this set. the repset needs to be only the sequences that were observed
    # in this data set, and we want reps for the step1 reference otus to be
    # reads from this run so we don't hit issues building a tree using
    # sequences of very different lengths. so...
    final_repset_fp = '%s/rep_set.fna' % output_dir
    index_links.append(
        ('OTU representative sequences',
         final_repset_fp,
         _index_headers['sequences']))
    final_repset_f = open(final_repset_fp, 'w')
    new_refseqs_fp = '%s/new_refseqs.fna' % output_dir
    index_links.append(
        ('New reference sequences (i.e., OTU representative sequences plus input '
         'reference sequences)',
         new_refseqs_fp,
         _index_headers['sequences']))
    # write non-singleton otus representative sequences from step1 to the
    # final rep set file
    for otu_id, seq in parse_fasta(open(step1_repset_fasta_fp, 'U')):
        if otu_id.split()[0] in otus_to_keep:
            final_repset_f.write('>%s\n%s\n' % (otu_id, seq))
    logger.write('# Write non-singleton otus representative sequences ' +
                 'from step1 to the final rep set file: %s\n\n' % final_repset_fp)
    # copy the full input refseqs file to the new refseqs_fp
    copyfile(refseqs_fp, new_refseqs_fp)
    new_refseqs_f = open(new_refseqs_fp, 'a')
    new_refseqs_f.write('\n')
    logger.write('# Copy the full input refseqs file to the new refseq file\n' +
                 'cp %s %s\n\n' % (refseqs_fp, new_refseqs_fp))
    # iterate over all representative sequences from step2 and step4 and write
    # those corresponding to non-singleton otus to the final representative set
    # file and the new reference sequences file.
    if run_step_2_and_3:
        for otu_id, seq in parse_fasta(open(step2_repset_fasta_fp, 'U')):
            if otu_id.split()[0] in otus_to_keep:
                new_refseqs_f.write('>%s\n%s\n' % (otu_id, seq))
                final_repset_f.write('>%s\n%s\n' % (otu_id, seq))
    if not suppress_step4:
        for otu_id, seq in parse_fasta(open(step4_repset_fasta_fp, 'U')):
            if otu_id.split()[0] in otus_to_keep:
                new_refseqs_f.write('>%s\n%s\n' % (otu_id, seq))
                final_repset_f.write('>%s\n%s\n' % (otu_id, seq))
    new_refseqs_f.close()
    final_repset_f.close()

    # steps 1-4 executed
    if run_step_2_and_3:
        logger.write('# Write non-singleton otus representative sequences from ' +
                     'step 2 and step 4 to the final representative set and the new reference' +
                     ' set (%s and %s respectively)\n\n' % (final_repset_fp, new_refseqs_fp))
    # only steps 1 and 4 executed
    else:
        logger.write('# Write non-singleton otus representative sequences from ' +
                     'step 4 to the final representative set and the new reference' +
                     ' set (%s and %s respectively)\n\n' % (final_repset_fp, new_refseqs_fp))

    # Prep the make_otu_table.py command
    otu_table_fp = '%s/otu_table_mc%d.biom' % (output_dir, min_otu_size)

    make_otu_table_cmd = 'make_otu_table.py -i %s -o %s' %\
        (otu_no_singletons_fp, otu_table_fp)
    commands.append([("Make the otu table", make_otu_table_cmd)])
    index_links.append(
        ('OTU table exluding OTUs with fewer than %d sequences' % min_otu_size,
         otu_table_fp,
         _index_headers['otu_tables']))
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=False)

    commands = []

    # initialize output file names - these differ based on what combination of
    # taxonomy assignment and alignment/tree building is happening.
    if run_assign_tax and run_align_and_tree:
        tax_input_otu_table_fp = otu_table_fp
        otu_table_w_tax_fp = \
            '%s/otu_table_mc%d_w_tax.biom' % (output_dir, min_otu_size)

        align_and_tree_input_otu_table = otu_table_w_tax_fp
        index_links.append(
            ('OTU table exluding OTUs with fewer than %d sequences and including OTU '
             'taxonomy assignments' % min_otu_size,
             otu_table_w_tax_fp,
             _index_headers['otu_tables']))

        pynast_failure_filtered_otu_table_fp = \
            '%s/otu_table_mc%d_w_tax_no_pynast_failures.biom' % (output_dir, min_otu_size)
        index_links.append(
            ('OTU table exluding OTUs with fewer than %d sequences and sequences that '
            'fail to align with PyNAST and including OTU taxonomy assignments' % min_otu_size,
             pynast_failure_filtered_otu_table_fp,
             _index_headers['otu_tables']))

    elif run_assign_tax:
        tax_input_otu_table_fp = otu_table_fp
        otu_table_w_tax_fp = \
            '%s/otu_table_mc%d_w_tax.biom' % (output_dir, min_otu_size)
        index_links.append(
            ('OTU table exluding OTUs with fewer than %d sequences and including OTU '
            'taxonomy assignments' % min_otu_size,
             otu_table_w_tax_fp,
             _index_headers['otu_tables']))

    elif run_align_and_tree:
        align_and_tree_input_otu_table = otu_table_fp
        pynast_failure_filtered_otu_table_fp = \
            '%s/otu_table_mc%d_no_pynast_failures.biom' % (output_dir,
                                                           min_otu_size)
        index_links.append(
            ('OTU table exluding OTUs with fewer than %d sequences and sequences that '
             'fail to align with PyNAST' % min_otu_size,
             pynast_failure_filtered_otu_table_fp,
             _index_headers['otu_tables']))

    if run_assign_tax:
        if exists(otu_table_w_tax_fp) and getsize(otu_table_w_tax_fp) > 0:
            logger.write(
                "Final output file exists (%s). Will not rebuild." %
                otu_table_w_tax_fp)
        else:
            # remove files from partially completed runs
            remove_files([otu_table_w_tax_fp], error_on_missing=False)

            taxonomy_fp = assign_tax(
                repset_fasta_fp=final_repset_fp,
                output_dir=output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                parallel=parallel,
                logger=logger,
                status_update_callback=status_update_callback)

            index_links.append(
                    ('OTU taxonomic assignments',
                    taxonomy_fp,
                    _index_headers['taxa_assignments']))

            # Add taxa to otu table
            add_metadata_cmd = 'biom add-metadata -i %s --observation-metadata-fp %s -o %s --sc-separated taxonomy --observation-header OTUID,taxonomy' %\
                (tax_input_otu_table_fp, taxonomy_fp, otu_table_w_tax_fp)
            commands.append([("Add taxa to OTU table", add_metadata_cmd)])

            command_handler(commands,
                            status_update_callback,
                            logger=logger,
                            close_logger_on_success=False)
            commands = []

    if run_align_and_tree:
        rep_set_tree_fp = join(output_dir, 'rep_set.tre')
        index_links.append(
            ('OTU phylogenetic tree',
             rep_set_tree_fp,
             _index_headers['trees']))
        if exists(pynast_failure_filtered_otu_table_fp) and\
           getsize(pynast_failure_filtered_otu_table_fp) > 0:
            logger.write("Final output file exists (%s). Will not rebuild." %
                         pynast_failure_filtered_otu_table_fp)
        else:
            # remove files from partially completed runs
            remove_files([pynast_failure_filtered_otu_table_fp],
                         error_on_missing=False)

            pynast_failures_fp = align_and_tree(
                repset_fasta_fp=final_repset_fp,
                output_dir=output_dir,
                command_handler=command_handler,
                params=params,
                qiime_config=qiime_config,
                parallel=parallel,
                logger=logger,
                status_update_callback=status_update_callback)

            # Build OTU table without PyNAST failures
            table = load_table(align_and_tree_input_otu_table)
            filtered_otu_table = filter_otus_from_otu_table(table,
                get_seq_ids_from_fasta_file(open(pynast_failures_fp, 'U')),
                0, inf, 0, inf, negate_ids_to_keep=True)
            write_biom_table(filtered_otu_table,
                             pynast_failure_filtered_otu_table_fp)

            command_handler(commands,
                            status_update_callback,
                            logger=logger,
                            close_logger_on_success=False)
            commands = []


    if close_logger_on_success:
        logger.close()

    if not suppress_index_page:
        index_fp = '%s/index.html' % output_dir
        generate_index_page(index_links, index_fp)
