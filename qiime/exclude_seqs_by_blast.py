#!/usr/bin/env python
# file exclude_seqs_by_blast.py
from __future__ import division

"""
A lightweight script for BLASTing one or more sequences against a number of BLAST databases, and returning FASTA files a) of the results that did match b) of the results that didn't match c) raw blast results and also d) returning a report containing the parameters used, which sequences were excluded and why.
"""

from os.path import join
from time import strftime, time

from skbio.parse.sequences import parse_fasta
from bfillings.blast import blast_seqs, Blastall, BlastResult


__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Zaneveld", "Rob Knight", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"


FORMAT_BAR =   """------------------------------""" * 2


def blast_genome(seqs, blast_db, e_value, max_hits, word_size, working_dir,
                 blast_mat_root, extra_params=[], DEBUG=True):
    """Blast sequences against all genes in a genome

    seqs -- input sequences as strings
    blast_db -- path to blast database
    e_value -- e_value (float)
    max_hits -- maximum sequences detected by BLAST  to show
    word_size -- word size for initial BLAST screen.
    blast_mat_root -- location of BLAST matrix files
    extra_params -- additional paramters to pass to BLAST
    DEBUG -- display verbose debugging outout
    """

    # set up params to use with blastp or
    params = {
        # matrix
        "-M": "BLOSUM62",

        # max procs
        "-a": "1",

        # expectation
        "-e": e_value,

        # max seqs to show
        "-b": max_hits,

        # Word size
        "-W": word_size,

        # max one line descriptions
        "-v": max_hits,

        # tabular output
        "-m": "9",

        # program
        "-p": "blastn"
    }
    params.update(extra_params)

    output = blast_seqs(seqs,
                        Blastall,
                        blast_db=blast_db,
                        params=params,
                        WorkingDir=working_dir,
                        add_seq_names=False,
                        blast_mat_root=blast_mat_root)

    raw_output = [x for x in output['StdOut']]
    return raw_output


def find_homologs(query_file, subject_genome, e_value, max_hits,
                  working_dir, blast_mat_root, wordsize,
                  percent_aligned, extra_params={},
                  require_hit=False, DEBUG=True):
    """BLAST query_file against subject_genome

    query_file -- .nuc file or other FASTA file to BLAST against all files in file_list

    subject_genome -- path to a KEGG .nuc file or other FASTA formated file.

    e-value -- e-value threshold for blasts

    percent_aligned -- minumum percent alignment, between 0.0 and 1.0

    max_hits,blast_mat_root,extra_params -- these are passed along to blastn

    DEBUG -- if True, display debugging output
    """
    start_time = time()
    raw_blast_output = []
    seqs = open(query_file, "U").readlines()

    if DEBUG:
        print "BLASTING %s vs. %s" % (query_file, subject_genome)

    blast_db = subject_genome

    raw_output_data = blast_genome(seqs,
                                   blast_db, e_value,
                                   max_hits, wordsize, working_dir,
                                   blast_mat_root, extra_params,
                                   DEBUG=DEBUG)

    if DEBUG:
        print "Length of raw BLAST results:", len(raw_output_data)

    curr_blast_result = BlastResult(raw_output_data)

    align_filter = make_percent_align_filter(percent_aligned)
    # should a mismatch filter be added?

    filtered_ids, removed_ids = query_ids_from_blast_result(curr_blast_result,
                                                            align_filter, DEBUG=DEBUG)

    return raw_output_data, filtered_ids, removed_ids


def sequences_to_file(results, outfile_name):
    """Translate a generator of label,seq tuples to an output file """

    f = open(outfile_name, 'w+')
    for label, seq in results:
        output_lines = []
        output_lines.append(">%s\n" % label)
        output_lines.append("%s\n" % seq)
        f.writelines(output_lines)
    f.close()


def no_filter(blast_subject_entry):
    """A placeholder filter function which always returns True"""
    return True


def make_percent_align_filter(min_percent):
    """Return a filter function that filters BLAST results on % alignment
    min_percent -- minimum percent match as a float between 0 and 1"""
    min_percent = float(min_percent) * 100

    def align_filter(blast_result):
        if float(blast_result['% IDENTITY']) < min_percent:
            return False
        else:
            return True

    return align_filter


def check_align_percent_field(d):
    """Check for empty percent identity fields in a dict"""
    if d['% IDENTITY']:
        return True
    else:
        return False


def query_ids_from_blast_result(
        blast_result, filter_fn=no_filter, DEBUG=False):
    """Returns a list of blast query ids, filtered by a given function.

        --blast_result:  BLAST result from BLAST app controller

        --filter_fn:  a function that, given a dict representing a BLAST result
                      returns True or False based on whether the result passes
                      some filter.
    """
    ok_ids = []
    removed_ids = []
    for id in blast_result:
        for entry in blast_result[id]:
            for subentry in entry:
                if not check_align_percent_field(subentry):
                    continue
                if not filter_fn(subentry):
                    removed_ids.append(id)
                    continue
                ok_ids.append(subentry['QUERY ID'])
    ok_ids = set(ok_ids)

    # Ensure query seqs with multiple BLAST hits, only some of which
    # are filtered out, don't end up in removed_ids

    removed_ids = set(removed_ids) - ok_ids
    return ok_ids, removed_ids


def ids_from_fasta_lines(lines):
    """Extract ids from label lines"""
    ids = []
    for line in lines:
        if not line.startswith(">"):
            continue
        id = id_from_fasta_label_line(line)
        ids.append(id)

    return ids


def id_from_fasta_label_line(line):
    "Extract id from fasta label line"
    id_field = line.split()[0]
    id = id_field.strip(">")
    return id


def seqs_from_file(ids, file_lines):
    """Extract labels and seqs from file"""

    for label, seq in parse_fasta(file_lines):

        if id_from_fasta_label_line(label) in ids:
            yield label, seq


def compose_logfile_lines(start_time, db_format_time, blast_time, option_lines,
                          formatdb_cmd, blast_results, options, all_ids,
                          hit_ids, removed_hit_ids,
                          included_ids, DEBUG):
    """Compose lines for a logfile from data on analysis"""

    log_lines = []
    log_lines.append("Sequence exclusion analysis run on %s" % strftime("%c"))
    log_lines.append(
        "Formatting subject database took  %2.f seconds" %
        (db_format_time))
    log_lines.append(
        "BLAST search took  %2.f minute(s)" %
        ((blast_time) / 60.0))
    log_lines.append(
        "Total analysis completed in %2.f minute(s)" %
        ((time() - start_time) / 60.0))

    log_lines.append(FORMAT_BAR)
    log_lines.append(
        "|                      Options                             |")
    log_lines.append(FORMAT_BAR)

    log_lines.extend(option_lines)
    log_lines.append("Subject database formatted with command: %s"
                     % formatdb_cmd)

    log_lines.append(FORMAT_BAR)
    log_lines.append(
        "|                      Results                             |")
    log_lines.append(FORMAT_BAR)

    log_lines.append("BLAST results above e-value threshold:")
    log_lines.append(
        "\t".join(["Query id", "Subject id", "percent identity", "alignment length",
                   "mismatches", "gap openings", "q. start", "q. end", "s. start", "s. end", "e-value", "bit score"]))

    for line in blast_results:
        if line.startswith("#"):
            continue
        else:
            log_lines.append(line)

    log_lines.append(
        "Hits matching e-value and percent alignment filter: %s" %
        ','.join(sorted(hit_ids)))

    log_lines.append(FORMAT_BAR)
    log_lines.append(
        "|                      Summary                             |")
    log_lines.append(FORMAT_BAR)

    log_lines.append("Input query sequences: %i" % len(all_ids))
    log_lines.append(
        "Query hits from BLAST: %i" %
        (len(hit_ids) + len(removed_hit_ids)))
    log_lines.append(
        "Query hits from BLAST lacking minimal percent alignment: %i" %
        len(removed_hit_ids))
    log_lines.append("Final hits: %i" % len(hit_ids))
    log_lines.append("Output screened sequences: %i" % len(included_ids))

    log_lines.append(FORMAT_BAR)
    log_lines.append(
        "|                       Output                             |")
    log_lines.append(FORMAT_BAR)

    log_lines.append(
        "Writing excluded sequences (hits matching filters) to: %s" %
        join(options.outputdir, "matching.fna"))
    log_lines.append(
        "Writing screened sequences (excluding hits matching filters) to: %s" %
        join(options.outputdir, "non-matching.fna"))
    log_lines.append(
        "Writing raw BLAST results to: %s" %
        join(options.outputdir, 'raw_blast_results.txt'))

    # format for printing
    revised_log_lines = []
    for line in log_lines:
        line = line + "\n"
        revised_log_lines.append(line)

    if DEBUG:
        for line in log_lines:
            print line

    return revised_log_lines


def check_options(parser, options):
    """Check to insure required options have been supplied"""
    if options.percent_aligned > 1.0:
        parser.error(
            "Please check -p option: should be between 0.0(0%) and 1.0(100%)")

    if options.querydb is None:
        parser.error(
            "Please check -i option: must specify path to a FASTA file")
    try:
        f = open(options.querydb, 'r')
        f.close()
    except IOError:
        parser.error(
            "Please check -i option: cannot read from query FASTA filepath")
    if options.subjectdb is None:
        parser.error(
            "Please check -d option: must specify path to a FASTA file")
    try:
        f = open(options.subjectdb, 'r')
        f.close()
    except IOError:
        parser.error(
            "Please check -d option: cannot read from subject FASTA filepath")
    if options.outputdir is None:
        parser.error(
            "Please check -o option: must specify the output directory path")


def format_options_as_lines(options):
    """Format options as a string for log file"""
    option_lines = []
    option_fields = str(options).split(",")

    for field in option_fields:
        option_lines.append(str(field).strip("{").strip("}"))

    return option_lines


def ids_to_seq_file(ids, infile, outfile, suffix=''):
    """Lookup FASTA recs for ids and record to file
    ids -- list of ids to lookup seqs for in infile

    infile -- path to FASTA file

    outfile -- base path to which to write FASTA entries
               with ids in supplied ids

    suffix  -- will be appended to outfile base path
    """

    seqs = seqs_from_file(ids, open(infile).readlines())
    out_path = outfile + suffix
    sequences_to_file(seqs, out_path)
