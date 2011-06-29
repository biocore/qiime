#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option
from os import popen, system, makedirs, mkdir
from os.path import split, splitext, join
from subprocess import check_call, CalledProcessError
from shutil import copy

from qiime.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.parse.fasta import MinimalFastaParser

from qiime.util import parse_command_line_parameters,\
    load_qiime_config, get_options_lookup, get_qiime_scripts_dir,\
    write_degapped_fasta_to_file
from qiime.parallel.identify_chimeric_seqs import get_job_commands,\
    get_poller_command
from qiime.parallel.util import split_fasta, get_random_job_prefix,\
    write_jobs_file,submit_jobs, compute_seqs_per_file,\
    build_filepaths_from_filepaths, write_filepaths_to_file,\
    write_merge_map_file_identify_chimeric_seqs
from qiime.identify_chimeric_seqs import make_cidx_file

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Parallel chimera detection"""
script_info['script_description']="""This script works like the identify_chimeric_seqs.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage']=[]
#copied from identify_chimeric_seqs.py
script_info['script_usage'].append(("""blast_fragments example""","""For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in assign_taxonomy.py) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to the directory "identify_chimeras/" and using default parameters (e.g. chimera detection method ("-m blast_fragments"), number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")):""","""%prog -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -o chimeric_seqs.txt"""))

script_info['script_usage'].append(("""ChimeraSlayer Example:""",
                                    """Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference database. The input sequences need to be provided in aligned (Py)Nast format and the reference database needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences!""",
                                    """%prog -m ChimeraSlayer -i repr_set_seqs_aligned.fasta -a ref_seq_set_aligned.fasta -o chimeric_seqs.txt"""))
                           
script_info['output_description']="""The result of parallel_identify_chimeric_seqs.py is a text file that identifies which sequences are chimeric."""

script_info['required_options']=[
    options_lookup['fasta_as_primary_input'],
]

chimera_detection_method_choices = ['blast_fragments','ChimeraSlayer']

script_info['optional_options']=[\
    make_option('-a', '--aligned_reference_seqs_fp',
        help='Path to (Py)Nast aligned reference sequences. '
        'REQUIRED when method ChimeraSlayer [default: %default]'), 

    make_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. [default: %default; REQUIRED when method is blast_fragments]'),

    make_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences (used to build a blast db when method blast_fragments). '
        '[default: %default; REQUIRED when method blast_fragments'+\
         ' if no blast_db is provided;]'),

    make_option('-b', '--blast_db',
        help='Database to blast against. Must provide either --blast_db or '
        '--reference_seqs_fp when method is blast_fragments [default: %default]'),
        
    make_option('-m','--chimera_detection_method',\
          type='choice',help='Chimera detection method. Choices: '+\
                    " or ".join(chimera_detection_method_choices) +\
                    '. [default:%default]',\
          choices=chimera_detection_method_choices, default='ChimeraSlayer'),
          
    make_option('-n','--num_fragments',\
          type='int',help='Number of fragments to split sequences into' +\
          ' (i.e., number of expected breakpoints + 1) [default: %default]',\
          default=3),
          
    make_option('-d','--taxonomy_depth',\
          type='int',help='Number of taxonomic divisions to consider' +\
          ' when comparing taxonomy assignments [default: %default]',\
          default=4),
          
    make_option('-e','--max_e_value',\
                    type='float',help='Max e-value to assign taxonomy' +\
                    ' [default: %default]', default=1e-30),

    make_option('--min_div_ratio',\
                    type='float',help='min divergence ratio '+\
                    '(passed to ChimeraSlayer). If set to None uses ' +\
                    'ChimeraSlayer default value. '+\
                    ' [default: %default]', default=None),       

    make_option('-o', '--output_fp',
        help='Path to store output [default: derived from input_seqs_fp]'),
          
    #Define parallel-script-specific parameters
    make_option('-N','--identify_chimeric_seqs_fp',action='store',\
           type='string',help='full path to '+\
           'scripts/identify_chimeric_seqs.py [default: %default]',\
           default=join(get_qiime_scripts_dir(),'identify_chimeric_seqs.py')),\
        
 options_lookup['jobs_to_start'],\
 options_lookup['poller_fp'],\
 options_lookup['retain_temp_files'],\
 options_lookup['suppress_submit_jobs'],\
 options_lookup['poll_directly'],\
 options_lookup['cluster_jobs_fp'],\
 options_lookup['suppress_polling'],\
 options_lookup['job_prefix'],\
 options_lookup['python_exe_fp'],\
 options_lookup['seconds_to_sleep']\
]

script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

   # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    identify_chimeric_seqs_fp = opts.identify_chimeric_seqs_fp
    reference_seqs_fp = opts.reference_seqs_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    poll_directly = opts.poll_directly
    
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    aligned_reference_seqs_fp = opts.aligned_reference_seqs_fp
    blast_db = opts.blast_db
    chimera_detection_method = opts.chimera_detection_method
    num_fragments = opts.num_fragments
    taxonomy_depth = opts.taxonomy_depth
    max_e_value = opts.max_e_value
    min_div_ratio = opts.min_div_ratio
    output_fp = opts.output_fp

    created_temp_paths = []
    
    #additional option checks (copied from scripts/identify_chimeric_seqs.py)
    if opts.chimera_detection_method == 'blast_fragments':
        if not (opts.blast_db or opts.reference_seqs_fp):
            option_parser.error('Must provide either --blast_db or'+\
                ' --reference_seqs_fp and --id_to_taxonomy_fp when'+\
                ' method is blast_fragments.')
        if not opts.id_to_taxonomy_fp:
            option_parser.error('Must provide --id_to_taxonomy_fp when method'+\
                ' is blast_fragments.')

        if opts.num_fragments < 2:
            option_parser.error('Invalid number of fragments (-n %d) Must be >= 2.' \
                                    % opts.num_fragments)

    elif opts.chimera_detection_method == 'ChimeraSlayer':
        if not opts.aligned_reference_seqs_fp:
            option_parser.error("Must provide --aligned_reference_seqs_fp "+\
                                    "when using method ChimeraSlayer")
            
    #set the output_fp if not set
    if not output_fp:
        input_basename = splitext(split(input_fasta_fp)[1])[0]
        output_fp = '%s_chimeric.txt' % input_basename

    #get the dir path to the output file
    output_dir, _ = split(output_fp)
    if output_dir=="":
        output_dir ="./"
    
    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with CHIM
    job_prefix = opts.job_prefix or get_random_job_prefix('CHIM')
    
    # A temporary output directory is created in output_dir named
    # job_prefix. Output files are then moved from the temporary 
    # directory to the output directory when they are complete, allowing
    # a poller to detect when runs complete by the presence of their
    # output files.
    working_dir = '%s/%s' % (output_dir,job_prefix)
    try:
        makedirs(working_dir)
        created_temp_paths.append(working_dir)
    except OSError:
        # working dir already exists
        pass

    if chimera_detection_method=="blast_fragments":
        #pre-compute the blast db and pass to parallel jobs
        blast_db, db_files_to_remove = \
             build_blast_db_from_fasta_path(reference_seqs_fp,
                                            output_dir=working_dir)
        #unset option because blast_db shall be used instead of ref_set
        reference_seqs_fp = None
        created_temp_paths.extend(db_files_to_remove)

    if chimera_detection_method=="ChimeraSlayer":
        #copy the reference files to working dir
        #ChimeraSlayer creates an index file of the ref and
        #will crash without write permission in the ref seqs dir
        _, new_ref_filename = split(aligned_reference_seqs_fp)
        copy(aligned_reference_seqs_fp, working_dir)
        aligned_reference_seqs_fp = working_dir + "/"+new_ref_filename
        created_temp_paths.append(aligned_reference_seqs_fp)
 
        #if given, also copy the unaligned ref db
        if reference_seqs_fp:
            _, new_ref_filename = split(reference_seqs_fp)
            copy(reference_seqs_fp, working_dir)
            reference_seqs_fp = working_dir + "/" + new_ref_filename

        #otherwise create it
        else:
            reference_seqs_fp = write_degapped_fasta_to_file(MinimalFastaParser( \
                    open(aligned_reference_seqs_fp)), tmp_dir=working_dir)
        #delete it afterwards
        created_temp_paths.append(reference_seqs_fp)

        #build blast db of reference, otherwise ChimeraSlayer will do it
        #and parallel jobs clash
        _, db_files_to_remove = \
             build_blast_db_from_fasta_path(reference_seqs_fp)
        created_temp_paths.extend(db_files_to_remove)

        #make the index file globally
        #Reason: ChimeraSlayer first checks to see if the index file is there.
        #If not it tries to create it. This can lead to race condition if several
        # parallel jobs try to create it at the same time.
        make_cidx_file(aligned_reference_seqs_fp)
        created_temp_paths.append(aligned_reference_seqs_fp+".cidx")

    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,\
      job_prefix,working_dir=output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)
    
    # generate the list of commands to be pushed out to nodes and the list of
    # output files generated by each job
    commands, job_result_filepaths = \
     get_job_commands(python_exe_fp, identify_chimeric_seqs_fp, tmp_fasta_fps,
                      output_dir, reference_seqs_fp, job_prefix, working_dir,
                      aligned_reference_seqs_fp, blast_db, chimera_detection_method,
                      min_div_ratio, num_fragments, taxonomy_depth, max_e_value,
                      id_to_taxonomy_fp)

    created_temp_paths += job_result_filepaths

    # Set up poller apparatus if the user does not suppress polling
    if not suppress_polling:
        # Write the list of files which must exist for the jobs to be 
        # considered complete
        expected_files_filepath = '%s/expected_out_files.txt' % working_dir
        write_filepaths_to_file(job_result_filepaths,expected_files_filepath)
        created_temp_paths.append(expected_files_filepath)
        
        # Write the mapping file which described how the output files from
        # each job should be merged into the final output files
        merge_map_filepath = '%s/merge_map.txt' % working_dir
        process_run_results_f =\
         'qiime.parallel.poller.basic_process_run_results_f'
        write_merge_map_file_identify_chimeric_seqs(job_result_filepaths,
                                                    output_dir,
                                                    merge_map_filepath,
                                                    input_file_basename,
                                                    out_fp=output_fp)
        created_temp_paths.append(merge_map_filepath)
        
        # Create the filepath listing the temporary files to be deleted,
        # but don't write it yet
        deletion_list_filepath = '%s/deletion_list.txt' % working_dir
        created_temp_paths.append(deletion_list_filepath)
        
        # Generate the command to run the poller, and the list of temp files
        # created by the poller
        if not poll_directly:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
             merge_map_filepath,deletion_list_filepath,process_run_results_f,\
             seconds_to_sleep=seconds_to_sleep)
            created_temp_paths += poller_result_filepaths
            # append the poller command to the list of job commands
            commands.append(poller_command)
        else:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
             merge_map_filepath,deletion_list_filepath,process_run_results_f,\
             seconds_to_sleep=seconds_to_sleep,command_prefix='',command_suffix='')
            created_temp_paths += poller_result_filepaths
        
        if not retain_temp_files:
            # If the user wants temp files deleted, now write the list of 
            # temp files to be deleted
            write_filepaths_to_file(created_temp_paths,deletion_list_filepath)
        else:
            # Otherwise just write an empty file
            write_filepaths_to_file([],deletion_list_filepath)
     
    # write the commands to the 'jobs files'
    write_jobs_file(commands,job_prefix=job_prefix,jobs_fp=jobs_fp)
    
    # submit the jobs file using cluster_jobs, if not suppressed by the
    # user
    if not opts.suppress_submit_jobs:
        submit_jobs(cluster_jobs_fp,jobs_fp,job_prefix)
        
    if poll_directly:
        try:
            check_call(poller_command.split())
        except CalledProcessError, e:
            print '**Error occuring when calling the poller directly. '+\
            'Jobs may have been submitted, but are not being polled.'
            print str(e)
            exit(-1)

if __name__ == "__main__":
    main()
