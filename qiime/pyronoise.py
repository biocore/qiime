#!/usr/bin/env python

"""Pyronoise app controller and pyronoise OTU picker."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jensreeder@gmail.com"
__status__ = "Pre-release"

from optparse import OptionParser
from os import system, makedirs, listdir, remove, rmdir
from os.path import exists, splitext, split
from commands import getoutput
from re import search

from cogent.app.util import get_tmp_filename
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files, app_path
from cogent.core.alignment import SequenceCollection
from cogent.parse.flowgram_parser import lazy_parse_sff_handle
from cogent.app.util import ApplicationNotFoundError


def parse_command_line_parameters(commandline_args=None):
    """ Parses command line arguments """
    usage = """usage: %prog [options] {-i INPUT_SFF_FILEPATH}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

Example usage:

    # Get detailed usage information
    python pyronoise.py -h

    #Denoise flowgrams in file 454Reads.sff.txt 
    python pyronoise.py -i 454Reads.sff.txt 

    #Denoise flowgrams in file 454Reads.sff.txt using 2 cores
    #son your machine in parallel (requires mpirun)
    python pyronoise.py -n 2 -i 454Reads.sff.txt

"""
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\
                          dest='verbose',\
                          help='Print information during execution to log file-- '+\
                          'useful for debugging [default: %default]')
    parser.add_option('-i','--input_file', action='store',\
                          type='string', dest='sff_fp',\
                          help='path to flowgram file (*.sff.txt) '+\
                          '[REQUIRED]')
    parser.add_option('-o','--output_dir', action='store',\
                          type='string', dest='output_dir',\
                          help='path to output directory '+\
                          '[default: %default]')
    parser.add_option('-n','--num_cpus', action='store',\
                          type='int', dest='num_cpus',\
                          help='number of CPUs '+\
                          '[default: %default]')
    parser.add_option('-s','--precision', action='store',\
                          type='float', dest='precision',\
                          help='(passed to pyroNoise) '+\
                          '[default: %default]')
    parser.add_option('-c','--cut-off', action='store',\
                          type='float', dest='cut_off',\
                          help='(passed to pyroNoise) '+\
                          '[default: %default]')
    parser.add_option('-k','--keep_intermediates', action='store_true',\
                          dest='keep',\
                          help='Do not delete intermediate PyroNoise files -- '+\
                          'useful for debugging [default: %default]')

     # Define defaults
    parser.set_defaults(verbose=False, keep=False, sff_fp=None,
                        output_dir="pyronoise_picked_otus/",
                        num_cpus=1, precision=15.0, cut_off=0.05)
    
    opts,args = parser.parse_args(commandline_args)
    if (not opts.sff_fp or (opts.sff_fp and not exists(opts.sff_fp))):
        parser.error(('Flowgram file path does not exist:\n %s \n'+\
                         'Pass a valid one via -i.')% opts.sff_fp)
    return opts,args


def write_pyronoise_file(flowgrams, num_flows, filename=None, prefix = "/tmp/"):
    """Write flowgrams to a (randomly) named file.

    flowgrams: flowgrams to write
    
    num_flows: Number of flowgrams in file

    filename: name of file to write

    prefix: directory prefix, defaults to /tmp/

    Note: Ids are mapped to 0,1,... - use id_mapping to remap 
    Returns filename and mapping of new ids to old ids
    """

    if (filename == None ):
        filename = get_tmp_filename(tmp_dir = prefix, suffix=".dat")

    fh = open(filename,"w")
    #Write header information
    #Note: Flow length hard coded here
    #should be safe as long as no longer flowgram in file
    fh.write("%d %d\n" % (num_flows, 400))

    id_mapping = {}
    for (i,f) in enumerate(flowgrams):
        id_mapping["Id_"+str(i)] = f.Name
        tabbed_flowgram_seq = " ".join(map(lambda a: "%.2f" %a, f.flowgram))
        l = len(f)
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
        num_cpus = 1 #set to a save value
        if log_fh:
            log_fh.write("Warning: mpirun not found. Falling back to one cpu")

    basename = get_tmp_filename(tmp_dir=outdir, prefix = "", suffix="")
    #copy flowgrams from input sff.txt to pyronoise-formatted file
    filename, id_mapping = write_pyronoise_file(flows, num_flows, filename = basename+".dat")
    
    if(num_cpus >1):
        mpi = "mpirun -np %d "% num_cpus
    else:
        mpi = ""
    cmd = mpi+ "FDist -in %s -out %s > /dev/null" % (filename, basename)
    
    if log_fh: 
        log_fh.write("Executing: %s\n" % cmd)
    system(cmd)

    #Qcluster is fast, so no mpi needed
    cmd = "QCluster -in %s.fdist -out %s >/dev/null" % (basename, basename) 
    if log_fh: 
         log_fh.write("Executing: %s\n" % cmd)
    system(cmd)

    cmd = mpi\
        + "PCluster -din %s -out %s -lin %s.list -s %f -c %f > %s.pout"\
        % (filename, basename, basename, precision, cut_off, basename)
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
    for (label,seq) in seqs:
        (name,id,count)=label.split('/')[-1].split('_')
        name = "%s_%s"%(name,id)
        centroids[name] = seq
        cluster_size[name] = int(count)
    
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
    for (i,(label,seq)) in enumerate(seqs):
        (name, id, count)=label.split('/')[-1].split('_')
        name = "%s"%(id)
        centroids[name] = seq
        cluster_size[name] = int(count)

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

def main(commandline_args=None):
    """run PyroNoise on input flowgrams"""
    from sys import argv
    opts, args = parse_command_line_parameters(commandline_args)

    input_seqs_dir, input_seqs_filename = split(opts.sff_fp)
    #split off .txt
    input_seqs_basename, ext = splitext(input_seqs_filename)
    #split off .sff
    input_seqs_basename, ext = splitext(input_seqs_basename)

    outdir = opts.output_dir

    if (not exists(outdir)):
        makedirs(outdir)

    log_fh=None
    if (opts.verbose):
        log_fh = open(outdir+"pyronoise.log", "w")
    
    centroids, cluster_mapping = pyroNoise_otu_picker(open(opts.sff_fp),
                                                      outdir, opts.num_cpus, log_fh, opts.keep,
                                                      opts.precision, opts.cut_off)

    # store mapping file and centroids
    result_otu_path = '%s/%s_otus.txt' % (outdir, input_seqs_basename)
    of = open(result_otu_path,'w')
    for i,cluster in cluster_mapping.iteritems():
        of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
    of.close()
    
    result_fasta_path = '%s/%s.fasta' % (outdir, input_seqs_basename)
    of = open(result_fasta_path,'w')
    of.write(SequenceCollection(centroids).toFasta()+"\n")

if __name__ == "__main__":
    main()
