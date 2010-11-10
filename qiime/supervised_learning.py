#!/usr/bin/env python

__author__ = "Dan Knights"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Release"

import subprocess
from os import remove, path, devnull
from os.path import join
from tempfile import mkdtemp
from cogent.app.util import get_tmp_filename
from qiime.util import get_qiime_project_dir
from qiime.format import format_otu_table
from cogent.app.util import CommandLineApplication, CommandLineAppResult, \
    FilePath, ResultPath, ApplicationError
from cogent.app.parameters import Parameters
from qiime.parse import parse_otu_table, parse_mapping_file
from numpy import array

def run_R_supervised_learner(
    predictor_fp, response_fp, response_name,
    model_names=['random_forest'],output_dir=None,
    param_file=None):
    """ Wrapper for the RSupervisedLearner App controller"""
 
    learner = RSupervisedLearner()
    results = learner(predictor_fp, response_fp, response_name, 
                        model_names, output_dir,
                        param_file=param_file)

    return results

def R_format_otu_table(otu_filepath, output_dir=None, write_to_tmp_file=True):
    """Formats OTU table for R (remove comments & column 1 header)
       If write_to_tmp_file, writes formatted file to tmp file and returns path
       else, returns lines to go in file
    """
    sample_ids, otu_ids, otu_matrix, lineages = \
        parse_otu_table(open(otu_filepath,'U').readlines())
    # first line is sample ids, no header for first column (how R likes it)
    lines = ['\t'.join(sample_ids)]
    for i in xrange(len(otu_ids)):
        # note: casting array as a string and calling "split" is much faster
        # than mapping "str" onto the array
        array_as_strings = str(otu_matrix[i,:])[1:-1].split()
        lines.append(otu_ids[i] + '\t' + '\t'.join(array_as_strings))
    if write_to_tmp_file:
        if output_dir is None:
            tmp_fp = get_tmp_filename(
                prefix='otus_R_format', suffix='.txt')
        else:
            tmp_fp = join(output_dir, 'otus_R_format.txt')
        fout = open(tmp_fp, 'w')
        fout.write('\n'.join(lines))
        fout.close()
        return tmp_fp
    else:
        return lines

def R_format_map_file(map_filepath, output_dir=None, write_to_tmp_file=True):
    """Formats map file for R (remove comments & column 1 header)
       If write_to_tmp_file, writes formatted file to tmp file and returns path
       else, returns lines to go in file
    """
    map_list = parse_mapping_file(open(map_filepath,'U').readlines())
    rows, headers = map_list[0], map_list[1]

    # first line is column headers, no header for first column (how R likes it)
    lines = ['\t'.join(headers[1:])]
    lines += ['\t'.join(row) for row in rows]
    if write_to_tmp_file:
        if output_dir is None:
            tmp_fp = get_tmp_filename(
                prefix='map_R_format', suffix='.txt')
        else:
            tmp_fp = join(output_dir, 'map_R_format.txt')
        fout = open(tmp_fp, 'w')
        fout.write('\n'.join(lines))
        fout.close()
        return tmp_fp
    else:
        return lines

class RSupervisedLearner(CommandLineApplication):
    """R Supervised Learner application controller
       Runs R with a source script (from qiime/support_files/R), and 
       passes in an OTU table and mapping file. Causes R to run a supervised
       classifier to predict labels from a given category from the mapping file
       using the provided OTUs.
    """
    _input_handler = '_input_as_path'
    _command = "R"
    _options ={}

    _R_parameters = {
        'flags': '--vanilla --slave'
        }

    # The name of the R script (located under qiime/support_files/R/)
    _R_script = 'error.est.r'

    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_R_parameters)

    def getHelp(self):
        """Returns documentation string"""
        help_str =\
        """
        Runs a supervised classifier with an OTU table as predictors and one
        column from a mapping file as the category labels.
        
        Outputs:
            predictions.txt: the labels predicted by the classifier for the given
                samples. Each sample is predicted by a model that was trained 
                without it. 
            probabilities.txt: the label probabilities for each of the given 
                samples. (if available)
            summary.txt: a summary of the results, including the expected
                generalization error of the classifier
            features.txt: a list of discriminative OTUs with their associated
                importance scores (if available)
            params.txt: a list of any non-default parameters used in training
                the model.
        
        For an overview of the application of supervised classification to 
        microbiota, see PubMed ID 21039646.
        """
        return help_str

    def __call__(self, predictor_fp, response_fp, response_name, 
        model_names, output_dir=None, remove_tmp=True,
        param_file=None):
        """Run the application with the specified kwargs on data
        
            data: A file nameinput_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
            
            returns a dict of CommandLineAppResult objects, one for each machine
            learning model, keyed by the model name
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = devnull
        else:
            outfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            outfile = open(outfilepath,'w')
        if suppress_stderr:
            errfile = devnull
        else:
            errfilepath = FilePath(self.getTmpFilename(self.TmpDir))
            errfile = open(errfilepath, 'w')
        predictor_fp = getattr(self,input_handler)(predictor_fp)
        response_fp = getattr(self,input_handler)(response_fp)
        # create random output dir if needed
        if output_dir is None:
            output_dir = mkdtemp(prefix='R_output_')

        rflags = self.RParameters['flags']
        rscript = self._get_R_script_path()
        base_command = self._get_base_command()
        cd_command, base_command = base_command.split(';')
        cd_command += ';'
        R_source_dir = self._get_R_script_dir()
        
        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        pre_command = 'cat'

        if param_file is None:
            param_file = ''
        command = self._commandline_join(
            [   cd_command, pre_command, '%s |' %(rscript), base_command,
                '--args', R_source_dir, predictor_fp, response_fp, response_name,
                output_dir, ','.join(model_names), param_file\
                #~ ,'>',outfile,'2>', errfile \
            ]
            )


        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command

        # run command, wait for output, get exit status
        proc = subprocess.Popen(command, shell=True, stdout=outfile, stderr=errfile)
        proc.wait()
        exit_status = proc.returncode

        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            if exit_status == 2:
                raise ApplicationError, \
                    'R library not installed: \n' + \
                    ''.join(open(errfilepath,'r').readlines()) + '\n'
            else:
                raise ApplicationError, \
                    'Unacceptable application exit status: %s, command: %s'\
                    % (str(exit_status),command) +\
                    ' Program output: \n\n%s\n'\
                     %(''.join(open(errfilepath,'r').readlines()))

        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfilepath,"r")
        err = None        
        if not suppress_stderr:
            err = open(errfilepath,"r")
       
        result = {}
        for i, model in enumerate(model_names):
            subdir = join(output_dir, model)
            # don't attempt to open the out/err files more than once
            if i == 1:
                out = err = None
            result[model] = CommandLineAppResult(
                out, err, exit_status, 
                result_paths=self._get_result_paths(subdir))

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _get_result_paths(self, output_dir):
        """Returns the filepaths for all result files"""
        files = {
            'features': 'features.txt',
            'summary': 'summary.txt',
            'predictions': 'predictions.txt',
            'probabilities': 'probabilities.txt',
            'params': 'params.txt',
        }
        result_paths = {}
        for name, file in files.iteritems():
            result_paths[name] = ResultPath(
                Path=path.join(output_dir, file), IsWritten=True)
        return result_paths

    def _get_R_script_dir(self):
        """Returns the path to the qiime R source directory
        """
        qiime_dir = get_qiime_project_dir()
        script_dir = path.join(qiime_dir,'qiime','support_files','R')
        return script_dir

    def _get_R_script_path(self):
        """Returns the path to the R script to be executed
        """
        return path.join(self._get_R_script_dir(), self._R_script)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self,exit_status):
        """ Return False to raise an error due to exit_status !=0 of application
        """
        if exit_status != 0:
            return False
        return True

    @property
    def RParameters(self):
        return self.__extract_parameters('R')

    def __extract_parameters(self, name):
        """Extracts parameters in self._<name>_parameters from self.Parameters

        Allows the program to conveniently access a subset of user-
        adjusted parameters, which are stored in the Parameters
        attribute.
        
        Relies on the convention of providing dicts named according to
        "_<name>_parameters" and "_<name>_synonyms".  The main
        parameters object is expected to be initialized with the
        contents of these dicts.  This method will throw an exception
        if either convention is not adhered to.
        """
        parameters = getattr(self, '_' + name + '_parameters')
        result = Parameters(parameters)
        for key in result.keys():
            result[key] = self.Parameters[key]
        return result
