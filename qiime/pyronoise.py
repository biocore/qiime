#!/usr/bin/env python
from __future__ import division

"""Pyronoise app controller and pyronoise OTU picker."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from os import system, listdir, remove, rmdir
from os.path import exists, split
from re import search
from itertools import chain

from cogent.app.util import get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files, app_path
from cogent.parse.flowgram_parser import lazy_parse_sff_handle
from cogent.app.util import ApplicationNotFoundError, ApplicationError
from cogent.parse.record import RecordError

from qiime.util import load_qiime_config

# Adapted from align_seqs.py
# Load Denoiser if it's available. If it's not, skip it if not but set up
# to raise errors if the user tries to use it.
try:
    from Denoiser.flowgram_clustering import denoise_seqs

except ImportError:
    def raise_denoiser_not_found_error(*args, **kwargs):
        raise ApplicationNotFoundError,\
         "Denoiser cannot be found.\nIs it installed? Is it in your $PYTHONPATH?"+\
         "\nYou can obtain the Denoiser from http://www.microbio.me/denoiser .\n"
    # set functions which cannot be imported to raise_denoiser_not_found_error
    denoise_seqs = raise_denoiser_not_found_error

def write_pyronoise_file(flowgrams, num_flows, filename=None, prefix = "/tmp/"):
    """Write flowgrams to a (randomly) named file.

    flowgrams: flowgrams to write
    
    num_flows: Number of flowgrams in file

    filename: name of file to write

    prefix: directory prefix, defaults to /tmp/

    Note: Ids are mapped to 0,1,... - use id_mapping to remap 
    Returns filename and mapping of new ids to old ids
    """

    max_flow_len = 400
    if (filename == None ):
        filename = get_tmp_filename(tmp_dir = prefix, suffix=".dat")

    fh = open(filename,"w")
    #Write header information
    #Note: Flow length hard coded here
    #should be safe as long as no longer flowgram in file
    fh.write("%d %d\n" % (num_flows, max_flow_len))

    id_mapping = {}
    for (i,f) in enumerate(flowgrams):
        id_mapping["Id_"+str(i)] = f.Name
        tabbed_flowgram_seq = " ".join(map(lambda a: "%.2f" %a, f.flowgram))
        l = len(f)
        if (l>max_flow_len):
            raise ApplicationError, "Flowgram longer than maximal allowed length"
        fh.write("%d %d %s\n" % (i, len(f), tabbed_flowgram_seq))
    fh.close()

    return filename, id_mapping

def pyroNoise_app(flows, num_flows, num_cpus=2, outdir = "/tmp/", log_fh=None,
                  precision=15.0, cut_off=0.05):
    """Runs pyronoise on flows and return basename of result files.

    flows: List of flowgrams
    
    num_flows: Number of flowgrams

    outdir: directory where intermediate files are stored

    num_cpus: number of cpus requested from mpirun
    """
    
    if(not (app_path("FDist") and app_path("QCluster") and app_path("PCluster"))):
        raise ApplicationNotFoundError,"PyroNoise binaries (FDist,QCluster,PCluster) not found."

    if(num_cpus>1 and not app_path("mpirun")):
        raise ApplicationError,"Can't run in parallel - mpirun not installed.\n"+\
            "Try running on one processor."
    # if mpi is not found, better raise Error and don't fall back to one cpu
    #        num_cpus = 1 #set to a save value
    #        if log_fh:
    #            log_fh.write("Warning: mpirun not found. Falling back to one cpu")

    basename = get_tmp_filename(tmp_dir=outdir, prefix = "", suffix="")
    #copy flowgrams from input sff.txt to pyronoise-formatted file
    filename, id_mapping = write_pyronoise_file(flows, num_flows, filename = basename+".dat")

    # if value is set use it, otherwise fall back to use default
    # hard coded in Pyronoise header files.
    data_fp = load_qiime_config()["pyronoise_data_fp"]
    data_opt = ""
    if data_fp:
        if not exists(data_fp):
            raise ApplicationError("File %s not exists. Check your setting of pyronoise_data_fp in the .qiime_config." %data_fp)
        data_opt = "-l %s" % data_fp

    if(num_cpus >1):
        mpi = "mpirun -np %d "% num_cpus
    else:
        mpi = ""
    cmd = mpi+ "FDist %s -in %s -out %s > /dev/null" % (data_opt, filename, basename)
    
    if log_fh: 
        log_fh.write("Executing: %s\n" % cmd)
    system(cmd)

    #Check if fdist actally produced an output file
    if (not exists(basename+".fdist")):
        remove(filename)
        raise ApplicationError, "Something went wrong with PyroNoise."+\
            " If using mpi, make sure it's setup properly."

    #Qcluster is fast, so no mpi needed
    cmd = "QCluster -in %s.fdist -out %s >/dev/null" % (basename, basename) 
    if log_fh: 
         log_fh.write("Executing: %s\n" % cmd)
    system(cmd)

    cmd = mpi\
        + "PCluster %s -din %s -out %s -lin %s.list -s %f -c %f > %s.pout"\
        % (data_opt, filename, basename, basename, precision, cut_off, basename)
    if log_fh: 
        log_fh.write("Executing: %s\n" % cmd)
    system(cmd)
  
    remove(filename)
    return basename, id_mapping

def remove_pyronoise_intermediates(basename):
    """Removes all intermediate pyronoise files.

    basename: the prefix for all files.
    """

    remove_files(map(lambda a: basename+a, 
                     [ ".otu", ".tau", "_cd.fa", ".fdist", ".pout", 
                       ".tree", "_cd2.fa", ".cen", ".list", ".qual", ".z",
                       "_cf.fa"]), True)
    for filename in listdir(basename):
        remove(basename+"/"+filename)
    rmdir(basename)

def pyroNoise(flows, num_flows, num_cpus=1, outdir = "/tmp/", log_fh=None,
              keep_intermediates=False, precision=15.0, cut_off=0.05):
    """Run Pyronoise and return centroids and cluster sizes.

    flows: an iterable source of flowgrams

    num_flows: the number of flowgrmas in the source
    
    num_cpus: number of cores used by mpirun
    
    outdir: directory to write temporary and result files to
    
    log_fh: an open file handle to a log
    
    keep_intermediates: If True no intermediate files will be deleted

    precision: precision value passed to pyronoise
    
    cut_off: initial cut-off passed to pyroNoise
    """

    basename, id_mapping = pyroNoise_app(flows, num_flows, num_cpus, outdir, log_fh,
                                         precision, cut_off)

    #retrieve the centroid flows
    centroids = {}
    cluster_size = {}
    seqs = MinimalFastaParser(open(basename+"_cd.fa"))
    try:
        for (label,seq) in seqs:
            (name,id,count)=label.split('/')[-1].split('_')
            name = "%s_%s"%(name,id)
            centroids[name] = seq
            cluster_size[name] = int(count)
    except RecordError:
        # If pyronoise doesn't find its Lookup.dat file it still produces output files
        # Usually, the output file contains just one fasta header and the MinimalFastaParser
        # will raise a RecordError
        raise ApplicationError, "Something failed with pyroNoise. Most likely the lookup file wasnt' found"

    if log_fh: 
        log_fh.write("Number of Cluster: %d" % len(centroids))
    if(not keep_intermediates):
         remove_pyronoise_intermediates(basename)
    return centroids, cluster_size

def pyroNoise_otu_picker(sff_fh, outdir="/tmp/", num_cpus=1,
                         log_fh=None, keep_intermediates=False,
                         precision=15.0, cut_off=0.05):
    """Pick OTUS using pyroNoise.

    sff_fh: an open filehandel to a .sff.txt file

    outdir: directory to write temporary and result files to

    num_cpus: number of cores used by mpirun
    
    log_fh: an open file handle to a log
    
    keep_intermediates: If True no intermediate files will be deleted

    precision: precision value passed to pyronoise
    
    cut_off: initial cut-off passed to pyroNoise
    
    Returns centroids and cluster mapping."""
    
    (flowgrams, header) = lazy_parse_sff_handle(sff_fh)
    basename, id_mapping = pyroNoise_app(flowgrams, int(header["# of Reads"]), 
                                         num_cpus, outdir, log_fh, precision, cut_off)
    centroids = {}
    cluster_size = {}
    #get centroids
    seqs = MinimalFastaParser(open(basename+"_cd.fa"))
    #get information from pyronoise header:
    #e.g.    > /long_dir_list/name_0_12   stands for centroid 0 with 12 cluster members
    try:
        for (i, (label,seq)) in enumerate(seqs):
            (name, id, count)=label.split('/')[-1].split('_')
            name = "%s"%(id)
            centroids[name] = seq
            cluster_size[name] = int(count)
    except RecordError:
        # If pyronoise doesn't find its Lookup.dat file it still produces output files
        # Usually, the output file contains just one fasta header and the MinimalFastaParser
        # will raise a RecordError
        raise ApplicationError, "Something failed with pyroNoise. Most likely the lookup file wasn't found."

    otu_map = {}
    #Read in individual clusters
    for filename in listdir(basename):
        seqs = MinimalFastaParser(open(basename+"/"+filename))
        cluster_ids = [a for (a,b) in seqs] 
        if(cluster_ids != []):
            match = search('i_(\d+)\.fa', filename)
            id = match.group(1)
            assert(cluster_size[id] == len(cluster_ids))
            otu_map[id] = [id_mapping[i] for i in cluster_ids]

    if(not keep_intermediates):
        remove_pyronoise_intermediates(basename)
        
    return centroids, otu_map

def fast_denoiser(sff_fp, fasta_fp, tmp_outdir, num_cpus, primer, verbose=True):
    """wrapper function calling methods from the Denoiser pakcage."""
    if num_cpus>1:
        denoise_seqs(sff_fp, fasta_fp, tmp_outdir,
                     primer=primer,  cluster=True, num_cpus=num_cpus,
                     verbose=verbose)
    else:
        denoise_seqs(sff_fp, fasta_fp, tmp_outdir, primer=primer,
                     verbose=verbose)

    #read centroids and singletons
    centroids = MinimalFastaParser(open(tmp_outdir+"/centroids.fasta"))
    singletons = MinimalFastaParser(open(tmp_outdir+"/singletons.fasta"))
    
    seqs = chain(centroids,singletons)

    #read mapping 
    mapping = {}
    cluster_mapping = open(tmp_outdir+"/denoiser_mapping.txt")
    for i,cluster in enumerate(cluster_mapping):
        cluster, members = cluster.split(':')
        members = members.split()
        clust = [cluster]
        clust.extend(members)
        mapping[i] = clust

    return seqs, mapping

def extract_cluster_size(line):
    """ extract the size of a cluster from the header line.

    line is expected to be of this format:
>GCC6FHY01EQVIC | cluster size: 5
    """
    cluster_size = line.split(":")[-1]

    try:
        cluster_size = int(cluster_size)
    except ValueError:
        return 0
    return cluster_size
        
