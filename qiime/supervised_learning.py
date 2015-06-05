#!/usr/bin/env python

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Dan Knights", "Luke Ursell", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"

from os import remove, path, devnull
from os.path import join, split, splitext, exists

from numpy import array, set_printoptions, nan, sqrt, mean, square
from cogent.app.util import (CommandLineApplication, CommandLineAppResult,
                             ResultPath, ApplicationError)
from burrito.parameters import ValuedParameter, FlagParameter
from biom.parse import convert_biom_to_table

from qiime.util import get_qiime_project_dir


def parse_feature_importances(filepath):
    """Returns vector of feature IDs, vector of importance scores
    """
    lines = open(filepath, 'U').readlines()
    feature_IDs = []
    scores = []
    for line in lines[1:]:
        words = line.strip().split('\t')
        feature_IDs.append(words[0].strip())
        scores.append(float(words[1].strip()))
    return array(feature_IDs), array(scores)


def run_supervised_learning(predictor_fp, response_fp, response_name,
                            ntree=1000, errortype='oob', output_dir='.', verbose=False, HALT_EXEC=False):
    """Run supervised learning (random forests here)

        predictor_fp: path to otu table
        response_fp: path to metadata table
        response_name: Column header for gradient variable in metadata table
        ntree: Number of trees in forest
        errortype: method for estimating generalization error
        output_dir: output directory
        verbose: print verbose output
        output_dir: directory where output should be written (default '.')
        HALT_EXEC: halt just before running the formatdb command and
    """
    # instantiate the object
    rsl = RSupervisedLearner(HALT_EXEC=HALT_EXEC)

    # set options
    rsl.Parameters['-m'].on(response_fp)
    rsl.Parameters['-c'].on(response_name)
    rsl.Parameters['-n'].on(str(ntree))
    rsl.Parameters['-o'].on(output_dir)
    rsl.Parameters['-e'].on(errortype)

    if verbose:
        rsl.Parameters['-v'].on()

    # run the app
    app_result = rsl(predictor_fp)

    # Hack: delete the temporary otu table left behind by hack biom conversion
    remove(join(output_dir, splitext(split(predictor_fp)[1])[0] + '.txt'))

    return app_result


def pooled_standard_deviation(input_variances):
    """Returns the average standard deviation for a list of st dev's

    Input is a list of floats (ideally)
    Output is a single float value
    """
    # compute and return pooled standard deviation
    return sqrt(mean(square([float(i) for i in input_variances])))


def calc_baseline_error_to_observed_error(baseline_error, est_error):
    '''Calculate (baseline error / observed error) for results file.

        Input two values, return single value float
    '''
    return float(baseline_error) / float(est_error)


def assemble_results(
        input_averages, input_variances, baseline_error, errortype,
        ntree):
    '''The summary format below is done on the R backend, so
        if a directory of input tables is used, form a results file
        that matches this format.

        Inputs: list of averages, list of variances, baseline error, and
            error type (oob, loo, cv5, cv10), and ntree
        Output: a list of lines to write out in the format below

        Model    Random Forest
        Error type  leave-one-out cross validation
        Estimated error (mean +/- s.d.) 0.33333 +/- 0.50000
        Baseline error (for random guessing)    0.44444
        Ratio baseline error to observed error  1.33333
        Number of trees 500
'''
    # initiate the results file
    results = ['Model   Random Forests', 'Error type    %s' % errortype]

    # average together the input_averages (regardless of errortype)
    ave_error = float(mean(input_averages))

    # a stdev is return for errortypes cv5, cv10 so get pooled st dev
    if errortype in ['cv5', 'cv10']:
        ave_stdev = pooled_standard_deviation(input_variances)
        est_error = '%s +/- %s' % (ave_error, ave_stdev)
        est_error_line = '\t'.join(['Estimated Error (mean +/- s.d)',
                                    est_error])

    # no stdev is provided for oob or loo, so only return est error
    elif errortype in ['oob', 'loo']:
        est_error_line = '\t'.join(['Estimated Error (mean)', str(ave_error)])

    # write Baseline Error line
    results.append(est_error_line)
    results.append('\t'.join(['Baseline Error (for random guessing',
                              str(baseline_error)]))

    # write Ratio error line
    ratio = calc_baseline_error_to_observed_error(baseline_error, ave_error)
    results.append(
        '\t'.join(['Ratio baseline error to observed error', str(ratio)]))

    # write Number of Trees line
    results.append('\t'.join(['Number of trees', str(ntree)]))

    return results


class RSupervisedLearner(CommandLineApplication):

    """ ApplicationController for detrending ordination coordinates
    """

    _command = 'R'
    _r_script = 'randomforests.r'

    _parameters = {
        # input data table, e.g. otu table
        '-i':
        ValuedParameter(Prefix='-', Name='i', Delimiter=' ', IsPath=True),
        # metadata table filepath
        '-m':
        ValuedParameter(Prefix='-', Name='m', Delimiter=' ', IsPath=True),
        # metadata category header
        '-c': ValuedParameter(Prefix='-', Name='c', Delimiter=' '),
        # error type = 'oob', 'cv5', 'cv10', 'loo'
        '-e': ValuedParameter(Prefix='-', Name='e', Delimiter=' '),
        '-n': ValuedParameter(Prefix='-', Name='n', Delimiter=' '),
        # output dir
        '-o':
        ValuedParameter(Prefix='-', Name='o', Delimiter=' ', IsPath=True),
        '-v': FlagParameter(Prefix='-', Name='v'),
    }
    _input_handler = '_input_as_parameter'
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameter(self, data):
        """ Set the input path and log path based on data (a fasta filepath)
        """
        # temporary hack: this converts a biom file to classic otu table
        # format for use within R
        if self.Parameters['-v'].Value:
            print 'Converting BIOM format to tab-delimited...'
        temp_predictor_fp = join(self.Parameters['-o'].Value,
                                 splitext(split(data)[1])[0] + '.txt')
        temp_predictor_f = open(temp_predictor_fp, 'w')
        temp_predictor_f.write(convert_biom_to_table(data))
        temp_predictor_f.close()
        predictor_fp = temp_predictor_fp

        self.Parameters['-i'].on(predictor_fp)
        # access data through self.Parameters so we know it's been cast
        # to a FilePath
        return ''

    def _get_result_paths(self, data):
        """ Build the dict of result filepaths
        """
        # access data through self.Parameters so we know it's been cast
        # to a FilePath
        wd = self.WorkingDir
        od = self.Parameters['-o'].Value
        result = {}
        result['confusion_matrix'] = ResultPath(
            Path=join(od, 'confusion_matrix.txt'), IsWritten=True)
        result['cv_probabilities'] = ResultPath(
            Path=join(od, 'cv_probabilities.txt'), IsWritten=True)
        result['features'] = ResultPath(
            Path=join(od,
                      'feature_importance_scores.txt'),
            IsWritten=True)
        result['mislabeling'] = ResultPath(
            Path=join(od, 'mislabeling.txt'), IsWritten=True)
        result['summary'] = ResultPath(
            Path=join(od, 'summary.txt'), IsWritten=True)
        return result

    def _get_R_script_dir(self):
        """Returns the path to the qiime R source directory
        """
        qiime_dir = get_qiime_project_dir()
        script_dir = join(qiime_dir, 'qiime', 'support_files', 'R')
        return script_dir

    def _get_R_script_path(self):
        """Returns the path to the R script to be executed
        """
        return join(self._get_R_script_dir(), self._r_script)

    # Overridden to add R-specific command-line arguments
    # This means:
    # R --slave --vanilla --args --source_dir $QIIMEDIR/qiime/support/R/
    # <normal params> < detrend.r
    def _get_base_command(self):
        """Returns the base command plus command-line options.
        """
        cd_command = ''.join(['cd ', str(self.WorkingDir), ';'])
        r_command = self._commandline_join(
            ['R', '--slave', '--no-restore', '--args'])
        source_dir_arg = self._commandline_join(['--source_dir',
                                                 self._get_R_script_dir()])

        script_arguments = self._commandline_join(
            [self.Parameters[k] for k in self._parameters])

        command_parts = [
            cd_command, r_command, source_dir_arg,
            script_arguments, '<', self._get_R_script_path()]
        return self._commandline_join(command_parts).strip()

    BaseCommand = property(_get_base_command)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command

           Taken from RDP classifier app controller
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self, exit_status):
        """ Return True when the exit status was 0
        """
        return exit_status == 0
