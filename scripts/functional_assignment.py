#!/usr/bin/env python
# File created on 09 Feb 2010

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger",
               "Jens Reeder", "William Walters", "Jose Carlos Clemente Litran"] 
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import splitext, split, exists, abspath, join
from qiime.util import (make_option, parse_command_line_parameters, create_dir, load_qiime_config)
from qiime.pick_otus  import BlastxOtuPicker
from qiime.pycogent_backports.usearch import Usearch, clusters_from_blast_uc_file
from qiime.format import format_otu_map

qiime_config = load_qiime_config()

def usearch_search(query_fp,
                   refseqs_fp,
                   output_dir,
                   evalue,
                   min_id,
                   queryalnfract,
                   targetalnfract,
                   maxaccepts,
                   maxrejects,
                   HALT_EXEC=False):
    
    output_fp = join(output_dir,'search_results.uc')
    
    params = {}
    
    app = Usearch(params,  
                  WorkingDir=qiime_config['temp_dir'],
                  HALT_EXEC=HALT_EXEC)
    
    app.Parameters['--evalue'].on(evalue)
    app.Parameters['--id'].on(min_id)
    app.Parameters['--queryalnfract'].on(queryalnfract)
    app.Parameters['--targetalnfract'].on(targetalnfract)
    app.Parameters['--maxaccepts'].on(maxaccepts)
    app.Parameters['--maxrejects'].on(maxrejects)
    
    data = {'--query':query_fp,
            '--uc':output_fp,
            '--db':refseqs_fp,
            }
    
    app_result = app(data)

assignment_constructors = {'blastx':BlastxOtuPicker,'usearch':usearch_search}

script_info={}
script_info['brief_description'] = """ Script for performing functional assignment of reads against a reference database """
script_info['script_description'] = """ """
script_info['script_usage'] = [("""""","""Run assignment with BLAST using default parameters""","""%prog -i query_nt.fasta -r refseqs_pr.fasta -o blast_assigned_functions/""")]
script_info['script_usage'].append(("""""","""Run assignment with BLAST using scricter e-value threshold""","""%prog -i query_nt.fasta -r refseqs_pr.fasta -o blast_assigned_function_strict/ -e 1e-70"""))
script_info['output_description'] = """ """

script_info['required_options'] = [
    make_option('-i', '--input_seqs_filepath',type='existing_filepath',
        help='Path to input sequences file'),
    ]

script_info['optional_options'] = [
    make_option('-m', '--assignment_method', type='choice',
        choices=assignment_constructors.keys(), default = "usearch",
        help=('Method for picking OTUs.  Valid choices are: ' +\
              ', '.join(assignment_constructors.keys()) +\
              '. [default: %default]')),

    make_option('-o', '--output_dir',type='new_dirpath',\
        help=('Path to store result file '
              '[default: ./<METHOD>_assigned_functions/]')),
              
    make_option('-r', '--refseqs_fp',type='existing_filepath',
        help=('Path to reference sequences to search against when using -m '
              'blast [default: %default]')),
              
    make_option('-b', '--blast_db',type='string',
        help=('Pre-existing database to blast against when using -m blast '
              '[default: %default]')),
              
    make_option('--min_aligned_percent',
        help=('Minimum percent of query sequence that can be aligned to consider a match'
              ' [default: %default]'),default=0.50,type='float'),
              
    make_option('-e', '--evalue', type='float', default=1e-10,
        help=('Max e-value to consider a match [default: %default]')),
              
    make_option('-s', '--min_percent_id', type='float', default=0.75,
        help=('Min percent id to consider a match [default: %default]')),
              
    make_option('--queryalnfract', type='float', default=0.75,
        help=('Min percent of the query seq that must match to consider a match [default: %default]')),
              
    make_option('--targetalnfract', type='float', default=0.75,
        help=('Min percent of the target/reference seq that must match to consider a match [default: %default]')),

    make_option('--max_accepts',type='int',default=1,
              help="max_accepts value to uclust and "
                   "uclust_ref [default: %default]"),
                   
    make_option('--max_rejects',type='int',default=32,
              help="max_rejects value to uclust and "
                   "uclust_ref [default: %default]"),

]

script_info['version'] = __version__

def main():
    # Parse the command line parameters
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    
    # Create local copies of the options to avoid repetitive lookups
    assignment_method = opts.assignment_method
    refseqs_fp = opts.refseqs_fp
    if refseqs_fp:
        refseqs_fp = abspath(refseqs_fp)
    blast_db = opts.blast_db
    min_percent_similarity = opts.min_percent_id
    min_aligned_percent = opts.min_aligned_percent
    verbose = opts.verbose
    
    if (assignment_method == 'blastx' and
        refseqs_fp == None and
        blast_db == None):
           option_parser.error('blastx requires refseqs_fp or blast_db')
    
    # use the otu_picking_method value to get the otu picker constructor
    assignment_constructor = assignment_constructors[assignment_method]
    
    # split the input filepath into components used for generating
    # the output file name
    input_seqs_filepath = abspath(opts.input_seqs_filepath)
    input_seqs_dir, input_seqs_filename = split(input_seqs_filepath)
    input_seqs_basename, ext = splitext(input_seqs_filename)
    
    # create the output directory name (if not provided) and 
    # create it if it doesn't already exist
    output_dir = abspath(opts.output_dir or assignment_method + '_assigned_functions')
    create_dir(output_dir, fail_on_exist=False)
    
    # Create the output and log file names
    result_path = '%s/%s_clusters.txt' % (output_dir,input_seqs_basename)
    log_path = '%s/%s_clusters.log' % (output_dir,input_seqs_basename)
    failure_path = '%s/%s_failures.txt' % (output_dir,input_seqs_basename)
    
    if assignment_method == 'blastx':
        params = {'max_e_value':evalue,
                  'Similarity': min_percent_similarity,
                  'min_aligned_percent':min_aligned_percent}
        otu_picker = assignment_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, 
                   log_path=log_path,
                   blast_db=opts.blast_db,
                   refseqs_fp=refseqs_fp)
    elif assignment_method == 'usearch':
        usearch_search(query_fp=input_seqs_filepath,
                       refseqs_fp=refseqs_fp,
                       output_dir=output_dir,
                       evalue=opts.evalue,
                       min_id=opts.min_percent_id,
                       queryalnfract=opts.queryalnfract,
                       targetalnfract=opts.targetalnfract,
                       maxaccepts=opts.max_accepts,
                       maxrejects=opts.max_rejects)
        otus, failures = clusters_from_blast_uc_file(\
         open('%s/search_results.uc' % output_dir,'U'),9)
        function_map_f = open('%s/function_map.txt' % output_dir,'w')
        for line in format_otu_map(otus.items(),''):
            function_map_f.write(line)
        function_map_f.close()
        
    else:
        ## other -- shouldn't be able to get here as a KeyError would have
        ## been raised earlier
        raise ValueError, "Unknown functional assignment method: %s" % assignment_method

if __name__ == "__main__":
    main()

