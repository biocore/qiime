#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"
 
import os
from qiime.util import parse_command_line_parameters
from optparse import make_option
from os.path import splitext
from qiime.util import get_qiime_project_dir
from qiime.make_sra_submission import (write_xml_generic, 
    make_run_and_experiment, make_submission, make_study, make_sample,
    generate_output_fp)

sra_template_dir = os.path.join(
    get_qiime_project_dir(), 'qiime', 'support_files', 'sra_xml_templates')

script_info={}
script_info['brief_description']="""Makes SRA submission files (xml-format)"""
script_info['script_description']="""This script makes the submission xml files for SRA (study, experiment, etc.).  This script assumes that there is a simple tab-delimited text input (allowing for examples and comments)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Read the sample data from sample.txt, the study data from study.txt, and the submission data from submission.txt, which writes out the corresponding XML files.""","""make_sra_submission.py -a sample.txt -t study.txt -u submission.txt"""))
script_info['script_usage'].append(("""""","""Produces the files study.xml, submission.xml, sample.xml (based on the filenames of the .txt files) using the default xml templates.""",""""""))
script_info['output_description']="""This script produces 3 xml-formatted files."""
script_info['required_options']=[]
script_info['optional_options']=[\
    make_option('-a','--input_sample_fp',\
        help='the tab-delimited text file with info about samples [default: %default]'),
    make_option('--template_sample_fp',
        default=os.path.join(sra_template_dir, 'sample_template.xml'),
        help='the template file for samples [default: %default]'),
    make_option('-t','--input_study_fp',\
        help='the tab-delimited text file with info about the study [default: %default]'),
    make_option('--template_study_fp',
        default=os.path.join(sra_template_dir, 'study_template.xml'),
        help='the template file for the study [default: %default]'),
    make_option('-u','--input_submission_fp',\
        help='the tab-delimited text file with info about the submission [default: %default]'),
    make_option('--template_submission_fp',
        default=os.path.join(sra_template_dir, 'submission_template.xml'),
        help='the template file for the submission [default: %default]'),
    make_option('-e', '--input_experiment_fp', \
        help='the tab-delimited text file with info about the experiment [default: %default]'),
    make_option('--experiment_attribute_fp',
        help='three-column, tab-delimited file of experiment attributes [default: %default]'),
    make_option('--experiment_link_fp',
        help='three-column, tab-delimited file of experiment links [default: %default]'),
    make_option('--multicolumn_input_format', action='store_true', 
        help='enables the use of input files for study and submission where field values are specified in columns, rather than rows [default: %default]'),
    make_option('-s', '--sff_dir', 
        help='the directory containing the demultiplexed sff files: 1 dir per run [default: %default]'),
    make_option('-o', '--output_dir',
        help='Submission output directory')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    docnames = {}

    xml_kwargs = {
        'twocol_input_format': not opts.multicolumn_input_format,
        }

    if opts.input_study_fp:
        docnames['study'] = write_xml_generic(opts.input_study_fp,
            opts.template_study_fp, make_study, xml_kwargs)

    if opts.input_sample_fp:
        docnames['sample'] = write_xml_generic(opts.input_sample_fp,
            opts.template_sample_fp, make_sample, xml_kwargs)

    if opts.input_experiment_fp:
        if not opts.sff_dir:
            option_parser.error("Must specify an sff dir if making an experiment.")

        if opts.experiment_attribute_fp:
            attribute_file = open(opts.experiment_attribute_fp, 'U')
        else:
            attribute_file = None
        if opts.experiment_link_fp:
            link_file = open(opts.experiment_link_fp, 'U')
        else:
            link_file = None

        experiment_xml, run_xml = make_run_and_experiment(
            open(opts.input_experiment_fp, 'U'), opts.sff_dir,
            attribute_file=attribute_file, link_file=link_file)

        output_experiment_fp = generate_output_fp(
            opts.input_experiment_fp, '.xml', opts.output_dir)
        with open(output_experiment_fp, 'w') as f:
            f.write(experiment_xml)
        docnames['experiment'] = os.path.basename(output_experiment_fp)

        output_run_fp = generate_output_fp(
            opts.input_experiment_fp, '_run.xml', opts.output_dir)
        with open(output_run_fp, 'w') as f:
            f.write(run_xml)
        docnames['run'] = os.path.basename(output_run_fp)

    if opts.input_submission_fp:
        input_submission_file = open(opts.input_submission_fp, 'U')
        submission_template = open(opts.template_submission_fp, 'U').read()
        submission_xml = make_submission(
            input_submission_file, submission_template, docnames,
            submission_dir=opts.output_dir,)

        output_submission_fp = generate_output_fp(
            opts.input_submission_fp, '.xml', opts.output_dir)
        with open(output_submission_fp, 'w') as f:
            f.write(submission_xml)


if __name__ == "__main__":
    main()
