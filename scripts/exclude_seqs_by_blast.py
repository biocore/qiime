#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division
import numpy
from time import time
from os import system, getcwd
from os.path import join
from qiime.exclude_seqs_by_blast import blast_genome,\
                                        check_options,\
                                        find_homologs,\
                                        sequences_to_file,\
                                        no_filter,\
                                        format_options_as_lines,\
                                        check_align_percent_field,\
                                        make_percent_align_filter,\
                                        query_ids_from_blast_result,\
                                        ids_from_fasta_lines,\
                                        id_from_fasta_label_line,\
                                        seqs_from_file,\
                                        ids_to_seq_file,\
                                        compose_logfile_lines

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jesse Zaneveld","Rob Knight", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option


script_info={}
script_info['brief_description']="""Exclude contaminated sequences using BLAST"""
script_info['script_description']="""

This code is designed to allow users of the QIIME workflow to conveniently exclude unwanted sequences from their data. This is mostly useful for excluding human sequences from runs to comply with Internal Review Board (IRB) requirements, but may also have other uses (e.g. perhaps excluding a major bacterial contaminant). Sequences from a run are searched against a user-specified subject database, where BLAST hits are screened by e-value and the percentage of the query that aligns to the sequence.

For human screening THINK CAREFULLY about the data set that you screen against. Are you excluding human non-coding sequences? What about mitochondrial sequences? This point is CRITICAL because submitting human sequences that are not IRB-approved is BAD.

(e.g. you would NOT want to just screen against just the coding sequences of the human genome as found in the KEGG .nuc files, for example)

One valid approach is to screen all putative 16S rRNA sequences against greengenes to ensure they are bacterial rather than human.

WARNING: You cannot use this script if there are spaces in the path to the database of fasta files because formatdb cannot handle these paths (this is a limitation of NCBI's tools and we have no control over it).
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following is a simple example, where the user can take a given FASTA file (i.e. resulting FASTA file from pick_rep_set.py) and blast those sequences against a reference FASTA file containing the set of sequences which are considered contaminated:""","""exclude_seqs_by_blast.py -i repr_set_seqs.fasta -d ref_seq_set.fna -o exclude_seqs/"""))
script_info['script_usage'].append(("""""","""Alternatively, if the user would like to change the percent of aligned sequence coverage ("-p") or the maximum E-value ("-e"), they can use the following command:""","""exclude_seqs_by_blast.py -i repr_set_seqs.fasta -d ref_seq_set.fna -o exclude_seqs/ -p 0.95 -e 1e-10"""))
script_info['output_description']="""Four output files are generated based on the supplied outputpath + unique suffixes:

1. "filename_prefix".screened: A FASTA file of sequences that did pass the screen (i.e. matched the database and passed all filters).

2. "filename_prefix".excluded: A FASTA file of sequences that did not pass the screen.

3. "filename_prefix".raw_blast_results: Contains the raw BLAST results from the screening.

4. "filename_prefix".sequence_exclusion_log: A log file summarizing the options used and results obtained.
"""
script_info['required_options']=[\
 make_option("-i","--querydb",dest='querydb',default = None,\
        help="The path to a FASTA file containing query sequences"),
 make_option("-d","--subjectdb",dest='subjectdb',default = None,\
        help="The path to a FASTA file to BLAST against"),
 make_option("-o","--outputfilename",dest='outputfilename',\
        default = None,\
        help=""" The base path/filename to save results. Sequences passing the screen, failing the screen, raw BLAST results and the log will be saved to your filename + '.screened', '.excluded', '.raw_blast_results', and '.sequence_exclusion_log' respectively.""")
]
script_info['optional_options']=[\
    make_option("-e","--e_value",type='float',dest='e_value',\
        default = 1e-10,\
        help="The e-value cutoff for blast queries [DEFAULT: %default]"),\
    make_option("-p","--percent_aligned",type='float',\
        dest='percent_aligned',default = 0.97,\
        help="The %% alignment cutoff for blast queries [DEFAULT: %default]"),\
    make_option("--blastmatroot",dest='blastmatroot',default = None,\
            help="Path to a folder containing blast matrices. [DEFAULT: %default]"),\
    make_option("--working_dir",dest='working_dir',default = "/tmp",\
        help="Working dir for BLAST [DEFAULT: %default]"),\
    make_option("-M","--max_hits",type='int',dest='max_hits',\
        default = 100,\
        help="""Max hits parameter for BLAST. CAUTION: Because filtering on alignment percentage occurs after BLAST, a max hits value of 1 in combination with an alignment percent filter could miss valid contaminants. [DEFAULT: %default]"""),\
    make_option("-W","--word_size",type='int',dest='wordsize',\
        default = 28,\
        help="Word size to use for BLAST search [DEFAULT: %default]")
]
script_info['version'] = __version__

FORMAT_BAR =   """------------------------------"""*2


def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)
    DEBUG = options.verbose 
    check_options(option_parser, options)
    start_time = time()
    option_lines = format_options_as_lines(options)
    if DEBUG:
        print FORMAT_BAR
        print "Running with options:"
        for line in sorted(option_lines):
            print line
        print FORMAT_BAR

    #because the blast app controller uses absolute paths, make sure subject
    #db path is fully specified
    subject_db = options.subjectdb
    if not subject_db.startswith('/'):
        subject_db = join(getcwd(), subject_db)
    formatdb_cmd = 'formatdb -p F -o T -i %s' % subject_db
    if DEBUG:
        print "Formatting subject db with command: %s" % formatdb_cmd

    system(formatdb_cmd)
    db_format_time = time() - start_time

    if DEBUG:
        print "Formatting subject db took: %2.f seconds" % db_format_time
        print FORMAT_BAR

    blast_results,hit_ids, removed_hit_ids = find_homologs(options.querydb,\
        subject_db, options.e_value,options.max_hits,\
        options.working_dir,options.blastmatroot, options.wordsize,\
                            options.percent_aligned, DEBUG=DEBUG)

    blast_time = (time() - start_time) - db_format_time

    if DEBUG:
        print "BLAST search took: %2.f minute(s)" % (blast_time/60.0)
        print FORMAT_BAR

    #Record raw blast results
    f=open("%s.raw_blast_results" % options.outputfilename,'w')
    f.writelines(blast_results)
    f.close()

    #Record excluded seqs
    ids_to_seq_file(hit_ids,options.querydb,options.outputfilename,".excluded")

    #Record included (screened) seqs
    all_ids = ids_from_fasta_lines(open(options.querydb).readlines())
    included_ids  = set(all_ids) - hit_ids
    ids_to_seq_file(included_ids,options.querydb,options.outputfilename,".screened")

    log_lines = compose_logfile_lines(start_time, db_format_time, blast_time,\
                                                   option_lines,formatdb_cmd,\
                                               blast_results,options,all_ids,\
                                                     hit_ids,removed_hit_ids,\
                                                          included_ids,DEBUG)

    if DEBUG:
        print "Writing summary to: %s" % options.outputfilename +\
                                          ".sequence_exclusion_log"

    f=open(options.outputfilename + ".sequence_exclusion_log",'w')
    f.writelines(log_lines)
    f.close()

if __name__ == "__main__":
    main()
