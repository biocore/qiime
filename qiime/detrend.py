from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Dan Knights", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"

from os import remove

from burrito.util import CommandLineApplication, ResultPath 
from os.path import join

from burrito.parameters import ValuedParameter, FlagParameter, FilePath

from qiime.util import get_qiime_project_dir
from qiime.parse import parse_mapping_file


class Detrender(CommandLineApplication):

    """ ApplicationController for detrending ordination coordinates
    """

    _command = 'R'
    _r_script = 'detrend.r'
    _parameters = {
        # input PCoA file
        '-i':
        ValuedParameter(Prefix='-', Name='i', Delimiter=' ', IsPath=True),
        # metadata table
        '-m':
        ValuedParameter(Prefix='-', Name='m', Delimiter=' ', IsPath=True),
        # gradient variable column header name
        '-c': ValuedParameter(Prefix='-', Name='c', Delimiter=' '),
        # output directory
        '-o':
        ValuedParameter(Prefix='-', Name='o', Delimiter=' ', IsPath=True),
        # use pre-rotation for optimal alignment with known gradient
        '-r': FlagParameter(Prefix='-', Name='r'),
    }
    _input_handler = '_input_as_parameter'
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameter(self, data):
        """ Set the input path based on data (data is the filepath)
        """
        self.Parameters['-i'].on(data)
        return ''

    def _get_result_paths(self, data):
        """ Build the dict of result filepaths
        """
        # access the output dir through self.Parameters so we know it's been cast
        # to a FilePath
        od = self.Parameters['-o'].Value

        result = {}
        # the before/after coords plot
        result['plot'] = ResultPath(
            Path=join(od,
                      'PCoA_vs_projection.pdf'),
            IsWritten=True)
        # the detrended coords file
        result['coords'] = ResultPath(
            Path=join(od, 'detrended_pcoa.txt'), IsWritten=True)
        # the summary file, only present if metadata was included
        summary_fp = join(od, 'summary.txt')
        result['summary'] = ResultPath(Path=summary_fp,
                                       IsWritten=self.Parameters['-c'].isOn())
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
    # R --slave --vanilla --args --source_dir $QIIMEDIR/qiime/support/R/ ...
    # <normal params> ... < detrend.r
    def _get_base_command(self):
        """Returns the base command plus R command-line options.
           Adapted from RDP Classifier app controller
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

           Taken from RDP Classifier app controller
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    def _accept_exit_status(self, exit_status):
        """ Return True when the exit status was 0
        """
        return exit_status == 0


def detrend_pcoa(input_fp, map_fp=None, gradient_variable=None,
                 suppress_prerotate=False, output_dir='.', HALT_EXEC=False):
    """Detrend PCoA scores in input_fp

        input_fp: path to PCOA file
        map_fp: path to metadata table
        gradient_variable: Column header for gradient variable in metadata table
        suppress_prerotate: True including metadata but don't want prerotation (default: False)
        output_dir: directory where output should be written (default '.')
        HALT_EXEC: halt just before running the formatdb command and
    """

    # instantiate the object
    dt = Detrender(HALT_EXEC=HALT_EXEC)

    # set params
    if map_fp is not None:
        dt.Parameters['-m'].on(map_fp)

    if gradient_variable is not None:
        dt.Parameters['-c'].on(gradient_variable)

    if suppress_prerotate:
        dt.Parameters['-r'].on()

    dt.Parameters['-o'].on(output_dir)

    # run the app
    app_result = dt(input_fp)
    return app_result


def validate_metadata(map_fp, gradient_variable):
    """Ensures that gradient variable is indeed present and real-valued

       If not, raises an exception
    """
    metadata = parse_mapping_file(map_fp)
    if not gradient_variable in metadata[1]:
        raise Exception("Header %s not found in metadata headers: " % (gradient_variable) +
                        ", ".join(metadata[1]))

    column_ix = metadata[1].index(gradient_variable)

    gradient = zip(*metadata[0])[column_ix]
    for i in xrange(len(gradient)):
        try:
            float(gradient[i])
        except ValueError:
            if gradient[i] != 'NA':
                raise Exception(
                    "Non-numeric value %s found in metadata gradient column: " % (gradient[i]) +
                    ", ".join(gradient))
