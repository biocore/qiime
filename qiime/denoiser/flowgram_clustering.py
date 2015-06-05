#!/usr/bin/env python

"""Several clustering methods to clean up 454 sequencing data"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Nigel Cook", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, popen, makedirs, rename, close
from os.path import exists
from collections import defaultdict
from itertools import izip, imap, ifilter, chain
from asyncore import loop
import datetime
from time import time
from math import fsum, trunc
from tempfile import mkstemp

from burrito.util import ApplicationNotFoundError, ApplicationError
from bfillings.denoiser import (lazy_parse_sff_handle, Flowgram,
                             FlowgramCollection, seq_to_flow)
from skbio.parse.sequences import parse_fasta

from qiime.format import write_Fasta_from_name_seq_pairs
from qiime.util import get_qiime_project_dir, load_qiime_config

from qiime.denoiser.utils import init_flowgram_file, append_to_flowgram_file,\
    FlowgramContainerFile, FlowgramContainerArray, make_stats, store_mapping,\
    store_clusters, read_denoiser_mapping, check_flowgram_ali_exe,\
    sort_seqs_by_clustersize, get_denoiser_data_dir, get_flowgram_ali_exe,\
    write_checkpoint, read_checkpoint, sort_mapping_by_size

from qiime.denoiser.cluster_utils import setup_cluster, adjust_workers,\
    stop_workers, check_workers, ClientHandler,\
    save_send, send_flowgram_to_socket
from qiime.denoiser.utils import write_sff_header
from qiime.denoiser.flowgram_filter import split_sff
from qiime.denoiser.preprocess import preprocess, preprocess_on_cluster,\
    read_preprocessed_data

DENOISER_DATA_DIR = get_denoiser_data_dir()


def compute_workload(num_cores, num_flows, spread):
    """Compute workload for each individual worker

    num_flows: total number of flows to be processed

    num_cores: total number of workers available for processing the flows

    spread: relative performance of the each worker, with 1.0 being nominal processing rate
    """
    # sigma is the sum of the normalized processing velocity scores
    # for each cluster processor. In a perfect world, sigma == num_cores
    # and the normalized processing velocity (in spread) == 1.0

    sigma = fsum(spread[0:num_cores])
    workload = [trunc((num_flows * x / sigma))for x in spread[0:num_cores]]

    while sum(workload) < num_flows:
        finish = workload[0] / spread[0]
        i = 0
        for x in range(1, num_cores):
            t = workload[x] / spread[x]
            if t < finish:
                finish = t
                i = x
        workload[i] += 1
    return workload


def adjust_processing_time(num_cores, workload, timing, epoch):
    """adjust processing time computes the nomalized relative worker throughput,
       with 1.0 being the nominal processing rate

    num_cores: number of processing workers

    workload: the number of flowgrams assigned to each worker

    timing: the time each worker finished processing

    epoch: the iteration start time
    """

    # sigma is the total throughput in flowgrams/sec
    sigma = 0.0
    for i in range(num_cores):
        timing[i] = workload[i] / (timing[i] - epoch)
        sigma += timing[i]
    #
    # spread represents the normalized flowgram/s processing rate
    # with 1.0 being the nominal processing speed
    #
    spread = [None for x in range(num_cores)]
    for i in range(num_cores):
        spread[i] = (timing[i] * num_cores) / sigma
    return spread


def get_flowgram_distances_on_cluster(
        id, flowgram, flowgrams, fc, ids, num_cores,
        num_flows, spread, client_sockets=[]):
    """Computes distance scores of flowgram to all flowgrams in parser.

    id: The flowgram identifier, also used to name intermediate files

    flowgram: This flowgram is used to filter all the other flowgrams

    flowgrams: iterable filehandle of flowgram file

    fc: a sink of flowgrams, which serves as source in the next round

    ids: list of flowgram ids that should be used from flowgrams

    num_cores: number of cpus

    num_flows: Number of flows in parser

    client_sockets: A list of open sockets for client-server communication

    spread: historical distribution of processing runtimes

    """
    epoch = time()

    check_flowgram_ali_exe()

    qiime_config = load_qiime_config()
    min_per_core = int(qiime_config['denoiser_min_per_core'])
    # if using from future import division this has to be checked,
    # as we want true integer division here

    per_core = max(min_per_core, (num_flows / num_cores) + 1)
    names = []
    scores = []

    # Need to call this here, since we iterate over the same iterator repeatedly.
    # Otherwise the call in ifilter will reset the iterator by implicitely  calling __iter__.
    # test if iter does the same
    flowgrams_iter = flowgrams.__iter__()
    # prepare input files and commands
    # synchronous client-server communication

    workload = compute_workload(num_cores, num_flows, spread)

    debug_count = 0
    for i in range(num_cores):
        socket = client_sockets[i]
        # send master flowgram to file first
        send_flowgram_to_socket(id, flowgram, socket)

        if(workload[i] < 1):
            # no data left for this poor guy
            save_send(socket, "--END--")
            continue
        else:
            # Then add all others which are still valid, i.e. in ids
            for (k, f) in (izip(range(workload[i]),
                                ifilter(lambda f: f.Name in ids, flowgrams_iter))):
                fc.add(f)
                send_flowgram_to_socket(k, f, socket, trim=False)
                names.append(f.Name)
                debug_count += 1
            # send the termination signal
            save_send(socket, "--END--")

    # asynchronous client-server communication
    # ClientHandlers write data in results
    results = [None] * num_cores
    timing = [0.0 for x in xrange(num_cores)]
    for i in range(num_cores):
        socket = client_sockets[i]
        ClientHandler(socket, i, results, timing)
    loop()
    # end asynchronous loop

    spread = adjust_processing_time(num_cores, workload, timing, epoch)

    # flatten list
    scores = [item for list in results for item in list]

    if (debug_count != len(scores)):
        raise RuntimeError("Something bad has happened! I received less " +
                           "alignment scores %d than there are flowgrams %d. Most likely this "
                           % (len(scores), debug_count) +
                           "means that the alignment program is not setup correctly or corrupted. " +
                           "Please run the test scripts to figure out the cause of the error.")

    return (scores, names, fc)


def get_flowgram_distances(id, flowgram, flowgrams, fc, ids, outdir,
                           error_profile=DENOISER_DATA_DIR +
                           'FLX_error_profile.dat'):
    """Computes distance scores of flowgram to all flowgrams in parser.

    id: The flowgram identifier, also used to name intermediate files

    flowgram: This flowgram is used to filter all the other flowgrams

    flowgrams: iterable filehandle of flowgram file

    fc: a sink for flowgrams, either a FlowgramContainerArray or
        FlowgramContainerFile object

    ids: dict of ids of flowgrams in flowgrams that should  be aligned

    outdir: directory for intermediate files

    error_profile: path to error profile *.dat file
    """
    check_flowgram_ali_exe()
    # File that serves as input for external alignment program
    (fh, tmpfile) = init_flowgram_file(prefix=outdir)
    append_to_flowgram_file(id, flowgram, fh)

    k = 0
    names = []
    for f in flowgrams:
        if(f.Name in ids):
            fc.add(f)
            append_to_flowgram_file(f.Name, f, fh, trim=False)
            k += 1
            names.append(f.Name)
    fh.close()

    # TODO: capture stderr and warn user
    scores_fh = popen("%s -relscore_pairid %s %s " %
                      (get_flowgram_ali_exe(),
                       error_profile, tmpfile), 'r')
    scores = [map(float, (s.split())) for s in scores_fh if s != "\n"]

    if (k != len(scores)):
        raise RuntimeError("Something bad has happened! I received less " +
                           "alignment scores than there are flowgrams. Most likely this " +
                           "means that the alignment program is not setup or corrupted. " +
                           "Please run the test scripts to figure out the cause of the error.")

    remove(tmpfile)

    return (scores, names, fc)


def filter_with_flowgram(
        id, flowgram, flowgrams, header, ids, num_flows, bestscores, log_fh,
        outdir="/tmp/", threshold=3.75, num_cpus=32,
        fast_method=True, on_cluster=False, mapping=None, spread=[],
        verbose=False, pair_id_thresh=0.97, client_sockets=[],
        error_profile=DENOISER_DATA_DIR + 'FLX_error_profile.dat'):
    """Filter all files in flows_filename with flowgram and split according to threshold.

    id: The flowgram identifier of the master flowgram of this round

    flowgram: This flowgram is used to filter all the other flowgrams

    flowgrams: iterator containing the flowgrams to be filtered

    header: a valid sff.txt header

    ids: this list marks the active flowgrams, i.e. flowgrams that are unclustered

    num_flows: Number of flows remaining in the current round

    bestscores: dictionary that stores for each unclustered flowgram the best
                score it has to to one of the centroids previously seen
                and the id of the centroid. Used in the second denoising phase.

    outdir: directory where intermediate and result files go

    threshold: Filtering threshold

    num_cpus: number of cpus to run on, if on_cluster == True

    fast_method: Boolean value for fast denoising with lots of memory

    on_cluster: Boolean flag for local vs cluster

    mapping: the current cluster mapping

    spread: worker processing throughput

    error_profile: Path to error profile *.dat file


    Implementation detail:
    The iterator behind 'flowgrams' is big and thus we want to keep its traversals
    at a minimum. The naive implementation of this filter function would traverse the
    iterator once to create the input file for the alignment routine, then a second
    time to do the actual filtering. To get rid of the second run through the iterator,
    we keep a list (in fact a dict) of active 'ids' and do the filtering only in the next
    round. A cleaner but still fast solution would be great, as this definitly poses a
    pitfall for future modifications.

    Returns filename of file containing all non-filtered flows and the number of flows
    """
    if verbose:
        log_fh.write("Filtering with %s: %d flowgrams\n" % (id, num_flows))

    # set up the flowgram storage
    if (not fast_method):
        fc = FlowgramContainerFile(header, outdir)
    else:
        fc = FlowgramContainerArray()

    # calculate distance scores
    if on_cluster:
        (scores, names, flowgrams) =\
            get_flowgram_distances_on_cluster(
                id, flowgram, flowgrams, fc, ids, num_cpus,
                num_flows, spread=spread, client_sockets=client_sockets)
    else:
        (scores, names, flowgrams) =\
            get_flowgram_distances(
                id, flowgram, flowgrams, fc, ids, outdir=outdir,
                error_profile=error_profile)

    # shortcut for non-matching flowgrams
    survivors = filter(
        lambda a_b: a_b[
            0] < threshold or a_b[
            1] >= pair_id_thresh,
        scores)
    if(len(survivors) == 0):
        # put it in its own cluster
        # and remove it from any further searches
        if (id in bestscores):
            del(bestscores[id])
        del(ids[id])
        return (flowgrams, num_flows - 1)

    # Do the filtering
    non_clustered_ctr = 0
    for ((score, pair_id), name) in zip(scores, names):
        if (score < threshold or name == id or pair_id >= pair_id_thresh):
            # make sure the original flowgram gets into this cluster
            del(ids[name])
            if (name in bestscores):
                del(bestscores[name])
            if(id != name):
                # update the mapping information
                mapping[id].extend(mapping[name])
                mapping[id].append(name)
                # delete the old cluster from the mapping
                del(mapping[name])
        else:
            non_clustered_ctr += 1
            # keep track of the best match of this guy to any centroid
            if (name not in bestscores or score < bestscores[name][1]):
                bestscores[name] = (id, score)

    # Some extra safety that we are not missing anything
    if (len(ids) != non_clustered_ctr
            or len(bestscores) != non_clustered_ctr):
        raise ApplicationError("filterWithFlowgram failed")

    return (flowgrams, non_clustered_ctr)


def secondary_clustering(sff_file, mapping, bestscores, log_fh,
                         threshold=4.5, verbose=False):
    """Clusters sequences based on their best distance to any of the centroids.

    Does not actually compute distances but uses the results of the first
    phase stored in bestscores.


    sff_file: name of unclustered flowgram file

    mapping: preliminary mapping file, dictionary of ids to list of ids

    bestscores: dictionary that stores for each unclustered flowgram the best
             score it has to to one of the centroid previously seen
             and the id of the centroid. Used in the second denoising phase.

    threshold: Secondary clustering threshold.

    """
    if(len(bestscores) == 0):
        # Either all sequence are already clustered or
        # we had no seq exceeding the bail out limit
        return

    (flowgrams, header) = lazy_parse_sff_handle(open(sff_file))

    counter = 0
    for f in flowgrams:
        (id, score) = bestscores[f.Name]
        if (score < threshold):
            counter += 1
            # update the mapping information
            mapping[id].extend(mapping[f.Name])
            mapping[id].append(f.Name)
            del(mapping[f.Name])
    if verbose:
        log_fh.write("Secondary clustering removed %d flowgrams\n" % counter)


def log_remaining_rounds(ids, cluster_mapping, bail_out, log_fh=None):
    """estimate the worst case number of rounds remaining

    ids: dict of active ids

    cluster_mapping: cluster mapping as dict

    bail_out: minimally required cluster size

    log_fh: log file handle
    """

    # this doesn't look very clever and might run a lot faster if rewritten
    remaining_rounds = len([id for id in ids.keys()
                            if len(cluster_mapping[id]) >= bail_out])
    # Remember, this is an unlikely worst case estimate
    if log_fh:
        log_fh.write("Rounds remaining in worst case: %d\n" % remaining_rounds)
    return remaining_rounds


def greedy_clustering(sff_fp, seqs, cluster_mapping, outdir, num_flows,
                      log_fh, num_cpus=1, on_cluster=False,
                      bail_out=1, pair_id_thresh=0.97, verbose=False,
                      threshold=3.75, fast_method=True,
                      error_profile=DENOISER_DATA_DIR +
                      'FLX_error_profile.dat',
                      max_num_rounds=None, checkpoint_fp=None):
    """second clustering phase of denoiser.

    sff_fp: flowgram file
    seqs: fasta seqs corresponding to sff_fp
    cluster_mapping: preliminary cluster mapping from phase I
    outdir: output directory
    num_flows: number of flowgrams in sff_fp (need to now before parsing sff_fp)
    log_fh: write verbose info to log_fh if set
    num_cpus:number of cpus to use of on_cluster ==True
    on_cluster: run in paralell if True
    bail_out: stop clustering with first cluster having bail_out members
    pair_id_thresh: always cluster flowgrams whose flowgram alignment implies a seq
                     identity of pair_id_thresh or higher
    verbose: be verbose or not
    threshold: low clustering threshold for phase II
    fast_method: use more memory intensive but faster method
    error_profile: path to error profile *.dat file
    max_num_rounds: If set, will stop clustering after this many rounds
    """

    (flowgrams, header) = lazy_parse_sff_handle(open(sff_fp))
    l = num_flows

    spread = [1.0 for x in range(num_cpus)]
    (client_sockets, workers) = (None, None)
    if on_cluster:
        (client_sockets, workers, server_socket) = \
            setup_cluster(num_cpus, outdir, verbose, error_profile)

    if checkpoint_fp:
        (checkpoint_key, round_ctr, cluster_mapping, ids, bestscores, sorted_keys) = \
            read_checkpoint(checkpoint_fp)
        skipping = True
    else:
        # ids stores all the active sequences
        # we initialize it with the ids from  the seqs dict here,
        # as it starts with all active flows.
        ids = dict.fromkeys(seqs)

        sorted_keys = sort_mapping_by_size(cluster_mapping)

        bestscores = {}
        round_ctr = 1

    # this is the main clustering loop, where most of the compute time is spent
    for key in sorted_keys:
        # skip until we reach the checkpoint
        if checkpoint_fp:
            if (checkpoint_key == key):
                if log_fh:
                    log_fh.write("Resume denoising with %s\n" % key)
                skipping = False
            if (skipping):
                continue

        if(key not in cluster_mapping):
            # this guy already has been clustered
            continue

        if (max_num_rounds and round_ctr > max_num_rounds):
            if log_fh:
                log_fh.write("Max number of rounds reached. " +
                             "Aborting clustering phase II and continuing with phase III.\n")
            break

        prefix_clustersize = len(cluster_mapping[key])
        # abort greedy first phase
        if(prefix_clustersize < bail_out):
            break

        # Do not take bad sequences as cluster seeds, as this will break the
        # code
        if('N' in seqs[key]):
            continue

        # check and delete workers if no longer needed
        if on_cluster:
            num_cpus = adjust_workers(l, num_cpus, client_sockets, log_fh)
            # check for dead workers
            check_workers(workers, client_sockets, log_fh)
            if num_cpus != len(spread):
                spread = [1.0 for x in range(num_cpus)]

        # write checkpoint right before expensive computation starts
        # Currently, write checkpint every 50 rounds,
        # could easily be changed here or exposed to command line
        if (round_ctr % 50) == 0:
            write_checkpoint(key, round_ctr, cluster_mapping, ids, bestscores,
                             sorted_keys, outdir)

        if log_fh:
            log_fh.write("Round %d:\n" % round_ctr)
            log_remaining_rounds(ids, cluster_mapping, bail_out, log_fh)

        ideal_flow = seq_to_flow(seqs[key])
        (new_flowgrams, newl) = filter_with_flowgram(key, ideal_flow, flowgrams, header, ids,
                                                     l, bestscores, log_fh, outdir,
                                                     on_cluster=on_cluster,
                                                     num_cpus=num_cpus,
                                                     fast_method=fast_method,
                                                     mapping=cluster_mapping,
                                                     verbose=verbose,
                                                     threshold=threshold,
                                                     pair_id_thresh=pair_id_thresh,
                                                     client_sockets=client_sockets,
                                                     error_profile=error_profile, spread=spread)
        l = newl
        flowgrams = new_flowgrams
        round_ctr += 1
        if(newl == 0):
            # all flowgrams clustered
            break
         # JR: I think this is too much info for the regular user, I leave it in, so
        # we can simply turn it on for debugging
#        if log_fh:
#            log_fh.write("Throughput Spread %s\n" % str(spread))

    if on_cluster:
        stop_workers(client_sockets, log_fh)
        server_socket.close()

    # write all remaining flowgrams into file for next step
    # TODO: might use abstract FlowgramContainer here as well
    fd, non_clustered_filename = mkstemp(dir=outdir, prefix="ff",
                                        suffix=".sff.txt")
    close(fd)
    non_clustered_fh = open(non_clustered_filename, "w")
    write_sff_header(header, non_clustered_fh)
    for f in flowgrams:
        if (f.Name in ids):
            non_clustered_fh.write(f.createFlowHeader() + "\n")

    return(non_clustered_filename, bestscores, cluster_mapping)


def denoise_seqs(
        sff_fps, fasta_fp, tmpoutdir, preprocess_fp=None, cluster=False,
        num_cpus=1, squeeze=True, percent_id=0.97, bail=1, primer="",
        low_cutoff=3.75, high_cutoff=4.5, log_fp="denoiser.log",
        low_memory=False, verbose=False,
        error_profile=DENOISER_DATA_DIR + 'FLX_error_profile.dat',
        max_num_rounds=None, titanium=False, checkpoint_fp=None):
    """The main routine to denoise flowgrams"""

    # abort if binary is missing
    check_flowgram_ali_exe()

    if verbose:
        # switch of buffering for log file
        log_fh = open(tmpoutdir + "/" + log_fp, "w", 0)
    else:
        log_fh = None

    # overwrite settings if titanium is set
    # This flag is only used from qiime. Remove after qiime integration
    if titanium:
        error_profile = DENOISER_DATA_DIR + "Titanium_error_profile.dat"
        low_cutoff = 4
        high_cutoff = 5

    if verbose:
        log_fh.write("Denoiser version: %s\n" % __version__)
        log_fh.write("SFF files: %s\n" % ', '.join(sff_fps))
        log_fh.write("Fasta file: %s\n" % fasta_fp)
        log_fh.write("Preprocess dir: %s\n" % preprocess_fp)
        if checkpoint_fp:
            log_fh.write("Resuming denoiser from %s\n" % checkpoint_fp)
        log_fh.write("Primer sequence: %s\n" % primer)
        log_fh.write("Running on cluster: %s\n" % cluster)
        log_fh.write("Num CPUs: %d\n" % num_cpus)
        log_fh.write("Squeeze Seqs: %s\n" % squeeze)
        log_fh.write("tmpdir: %s\n" % tmpoutdir)
        log_fh.write("percent_id threshold: %.2f\n" % percent_id)
        log_fh.write("Minimal sequence coverage for first phase: %d\n" % bail)
        log_fh.write("Low cut-off: %.2f\n" % low_cutoff)
        log_fh.write("High cut-off: %.2f\n" % high_cutoff)
        log_fh.write("Error profile: %s\n" % error_profile)
        log_fh.write("Maximal number of iteration: %s\n\n" % max_num_rounds)

    # here we go ...
    # Phase I - clean up and truncate input sff
    if(checkpoint_fp):
        if (preprocess_fp):
            # skip preprocessing as we should have data
            # we already have preprocessed data, so use it
            (deprefixed_sff_fp, l, mapping,
             seqs) = read_preprocessed_data(preprocess_fp)
        else:
            raise ApplicationError(
                "Resuming from checkpoint requires --preprocess option")

    else:
        if(preprocess_fp):
            # we already have preprocessed data, so use it
            (deprefixed_sff_fp, l, mapping,
             seqs) = read_preprocessed_data(preprocess_fp)
        elif(cluster):
            preprocess_on_cluster(sff_fps, log_fp, fasta_fp=fasta_fp,
                                  out_fp=tmpoutdir, verbose=verbose,
                                  squeeze=squeeze, primer=primer)
            (deprefixed_sff_fp, l, mapping,
             seqs) = read_preprocessed_data(tmpoutdir)
        else:
            (deprefixed_sff_fp, l, mapping, seqs) = \
                preprocess(
                    sff_fps, log_fh, fasta_fp=fasta_fp, out_fp=tmpoutdir,
                    verbose=verbose, squeeze=squeeze, primer=primer)

        # preprocessor writes into same file, so better jump to end of file
        if verbose:
            log_fh.close()
            log_fh = open(tmpoutdir + "/" + log_fp, "a", 0)

    # phase II:
    # use prefix map based clustering as initial centroids and greedily
    # add flowgrams to clusters with a low threshold

    (new_sff_file, bestscores, mapping) = \
        greedy_clustering(deprefixed_sff_fp, seqs, mapping, tmpoutdir, l,
                          log_fh, num_cpus=num_cpus, on_cluster=cluster,
                          bail_out=bail, pair_id_thresh=percent_id,
                          threshold=low_cutoff, verbose=verbose,
                          fast_method=not low_memory,
                          error_profile=error_profile,
                          max_num_rounds=max_num_rounds,
                          checkpoint_fp=checkpoint_fp)

    # phase III phase:
    # Assign seqs to nearest existing centroid with high threshold
    secondary_clustering(new_sff_file, mapping, bestscores, log_fh,
                         verbose=verbose, threshold=high_cutoff)
    remove(new_sff_file)
    if (verbose):
        log_fh.write("Finished clustering\n")
        log_fh.write("Writing Clusters\n")
        log_fh.write(make_stats(mapping) + "\n")
    store_clusters(mapping, deprefixed_sff_fp, tmpoutdir)
    store_mapping(mapping, tmpoutdir, "denoiser")


def denoise_per_sample(sff_fps, fasta_fp, tmpoutdir, cluster=False,
                       num_cpus=1, squeeze=True, percent_id=0.97, bail=1,
                       primer="", low_cutoff=3.75, high_cutoff=4.5,
                       log_fp="denoiser.log", low_memory=False, verbose=False,
                       error_profile=DENOISER_DATA_DIR +
                       'FLX_error_profile.dat',
                       max_num_rounds=None, titanium=False):
    """Denoise each sample separately"""

    # abort early if binary is missing
    check_flowgram_ali_exe()

    log_fh = None
    if log_fp:
        # switch of buffering for global log file
        log_fh = open(tmpoutdir + "/" + log_fp, "w", 0)

    # overwrite settings if titanium is set
    # This flag is only used from qiime. Remove after qiime integration
    if titanium:
        error_profile = DENOISER_DATA_DIR + "Titanium_error_profile.dat"
        low_cutoff = 4
        high_cutoff = 5

    if verbose:
        log_fh.write("Denoiser version: %s\n" % __version__)
        log_fh.write("SFF files: %s\n" % ', '.join(sff_fps))
        log_fh.write("Fasta file: %s\n" % fasta_fp)
        log_fh.write("Cluster: %s\n" % cluster)
        log_fh.write("Num CPUs: %d\n" % num_cpus)
        log_fh.write("Squeeze Seqs: %s\n" % squeeze)
        log_fh.write("tmpdir: %s\n\n" % tmpoutdir)
        log_fh.write("percent_id threshold: %.2f\n" % percent_id)
        log_fh.write("Minimal sequence coverage for first phase: %d\n" % bail)
        log_fh.write("Error profile: %s\n" % error_profile)
        log_fh.write("Maximal number of iteration: %s\n\n" % max_num_rounds)

    # here we go ...
    sff_files = split_sff(map(open, sff_fps), open(fasta_fp), tmpoutdir)
    combined_mapping = {}
    result_centroids = []
    result_singletons_files = []
    # denoise each sample separately
    for i, sff_file in enumerate(sff_files):
        if not exists(tmpoutdir + ("/%d" % i)):
            makedirs(tmpoutdir + ("/%d" % i))
        out_fp = tmpoutdir + ("/%d/" % i)
        denoise_seqs([sff_file], fasta_fp, out_fp, None, cluster,
                     num_cpus, squeeze, percent_id, bail, primer,
                     low_cutoff, high_cutoff, log_fp, low_memory,
                     verbose, error_profile, max_num_rounds)

        # collect partial results
        this_rounds_mapping = read_denoiser_mapping(
            open(out_fp + "/denoiser_mapping.txt"))
        combined_mapping.update(this_rounds_mapping)
        result_centroids.append(
            parse_fasta(open(out_fp + "/centroids.fasta")))
        result_singletons_files.append(out_fp + "/singletons.fasta")

    # write the combined files
    store_mapping(combined_mapping, tmpoutdir, "denoiser")
    seqs = chain(*result_centroids)
    fasta_fh = open(tmpoutdir + "/denoised.fasta", "w")
    # write centroids sorted by clustersize
    write_Fasta_from_name_seq_pairs(
        sort_seqs_by_clustersize(seqs, combined_mapping),
        fasta_fh)
    for singleton_file in result_singletons_files:
        write_Fasta_from_name_seq_pairs(
            parse_fasta(open(singleton_file, "r")),
            fasta_fh)
    fasta_fh.close()

    # return outdir for tests/test_denoiser
    return tmpoutdir
