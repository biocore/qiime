#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os.path import splitext
from qiime.make_sra_submission import (write_xml_generic, 
    make_run_and_experiment, make_submission)

script_description = """Makes the submission xml files for SRA (study, experiment, etc.).

Assumes simple tab-delimited text input (allowing examples/comments; produces
xml output. """

script_usage = """ Read the sample data from sample.txt, the study data from study.txt,
 and the submission data from submission.txt. Write out the corresponding
 XML files.
 
 make_study_and_sample_xml.py -a sample.txt -t study.txt -u submission.txt

 Produces files study.xml, submission.xml, sample.xml (based on filenames of
 the .txt files) using the default xml templates.
"""

required_options = [\
]

optional_options = [\
    make_option('-a','--input_sample_fp',\
        help='the tab-delimited text file with info about samples [default: %default]'),
    make_option('--template_sample_fp', default='sample_template.xml',\
        help='the template file for samples [default: %default]'),
    make_option('-t','--input_study_fp',\
        help='the tab-delimited text file with info about the study [default: %default]'),
    make_option('--template_study_fp', default='study_template.xml',\
        help='the template file for the study [default: %default]'),
    make_option('-u','--input_submission_fp',\
        help='the tab-delimited text file with info about the submission [default: %default]'),
    make_option('--template_submission_fp', default='submission_template.xml',\
        help='the template file for the submission [default: %default]'),
    make_option('-e', '--input_experiment_fp', \
        help ='the tab-delimited text file with info about the experiment [default: %default]'),
    make_option('-s', '--sff_dir', 
        help = 'the directory containing the demultiplexed sff files: 1 dir per run [default: %default]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    docnames = {}
    if opts.input_study_fp:
        docnames['study'] = write_xml_generic(opts.input_study_fp,
            opts.template_study_fp, make_study)

    if opts.input_sample_fp:
        docnames['sample'] = write_xml_generic(opts.input_sample_fp,
            opts.template_sample_fp, make_sample)

    if opts.input_experiment_fp:
        #in this case, we need to need to also get the sff dir
        if not opts.sff_dir:
            raise IOError, "Must specify an sff dir if making an experiment."
        base_name, ext = splitext(opts.input_experiment_fp)
        base_name = base_name.split('experiment')[-1]
        if base_name:
            base_name += '_'
        run_path = base_name + 'run.xml'
        run_file = open(run_path, 'w')
        experiment_path = base_name + 'experiment.xml'
        experiment_file = open(experiment_path, 'w')
        experiment_xml, run_xml = make_run_and_experiment(
            open(opts.input_experiment_fp, 'U'), opts.sff_dir)
        run_file.write(run_xml)
        experiment_file.write(experiment_xml)
        run_file.close()
        experiment_file.close()
        docnames['run'] = run_path
        docnames['experiment'] = experiment_path

    if opts.input_submission_fp:
        submission_template = open(opts.template_submission_fp, 'U').read()
        base_name, ext = splitext(opts.input_submission_fp)
        outfilename = base_name + '.xml'
        outfile = open(outfilename, 'w')
        outfile.write(make_submission(open(opts.input_submission_fp, 'U'),
            submission_template, docnames))
        outfile.close()


if __name__ == "__main__":
    main()
