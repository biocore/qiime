#!/usr/bin/env python

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Development"

import subprocess
from os import remove, path, devnull
from os.path import join
from sys import stdout
from time import sleep
from tempfile import mkdtemp
from qiime.util import get_tmp_filename
from qiime.util import get_qiime_project_dir
from cogent.app.util import CommandLineApplication, CommandLineAppResult, \
    FilePath, ResultPath, ApplicationError
from cogent.app.parameters import Parameters
from qiime.parse import parse_otu_table, parse_mapping_file
from numpy import array, set_printoptions, nan

def parse_feature_importances(filepath):
    """Returns vector of feature IDs, vector of importance scores
    """
    lines = open(filepath,'U').readlines()
    feature_IDs = []
    scores = []
    for line in lines[1:]:
        words = line.strip().split('\t')
        feature_IDs.append(words[0].strip())
        scores.append(float(words[1].strip()))
    return array(feature_IDs), array(scores)

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
        'flags': '--slave'
        }

    # The name of the R script (located under qiime/support_files/R/)
    _R_script = 'randomforests.r'

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
            cv_probabilities.txt: the label probabilities for each of the given 
                samples. (if available)
            mislabeling.txt: A convenient presentation of cv_probabilities for 
                mislabeling detection.
            confusion_matrix.txt: confusion matrix for hold-out predictions.
            summary.txt: a summary of the results, including the expected
                generalization error of the classifier
            feature_importance_scores.txt: a list of discriminative OTUs with their associated
                importance scores (if available)
        
        For an overview of the application of supervised classification to 
        microbiota, see PubMed ID 21039646.
        """
        return help_str

    def __call__(self, predictor_fp, response_fp, response_name, 
        ntree=1000, errortype='oob', output_dir=None, remove_tmp=True,verbose=False):
        """Run the application with the specified kwargs on data
        
            data: A file nameinput_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
            
            returns a CommandLineAppResult object
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

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        args = ['--qiime_dir', get_qiime_project_dir(), 
                '-i', predictor_fp, 
                '-m', response_fp,
                '-c', response_name,
                '-o', output_dir, 
                '-e', errortype,
                '--ntree', str(ntree)]
        if verbose:
                args += ['-v']

        command = self._commandline_join(
            [   cd_command, base_command,
                '--args'
            ] + args + [' < %s ' %(rscript)]
            )

        if verbose:
            print "Command: ", command
        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command

        # run command, wait for output, get exit status
        proc = subprocess.Popen(command, shell=True, stdout=outfile, stderr=errfile)

        if verbose:
            print '\nR output\n'
            tmpoutfile = open(outfilepath,'U')
            while proc.poll() is None:
                stdout.write(tmpoutfile.readline())
                sleep(0.01)
            tmpoutfile.close()
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
        try:
            result = CommandLineAppResult(
                out, err, exit_status, 
                result_paths=self._get_result_paths(output_dir))
        except ApplicationError, ae:
            msg = str(ae) + \
                '\n\ncommand: %s'\
                % (command) +\
                ' \n\nProgram stdout:\n%s'\
                  %(''.join(open(outfilepath,'r').readlines())) +\
                 ' \n\nProgram stderr:\n%s'\
                 %(''.join(open(errfilepath,'r').readlines()))
            raise ApplicationError, msg

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _get_result_paths(self, output_dir):
        """Returns the filepaths for all result files"""
        files = {
            'features': 'feature_importance_scores.txt',
            'summary': 'summary.txt',
            'cv_probabilities': 'cv_probabilities.txt',
            'mislabeling': 'mislabeling.txt',
            'confusion_matrix': 'confusion_matrix.txt',
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

