#!/usr/bin/env python
"""Application controller for rdp_classifier-2.0
"""

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Kyle Bittinger","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"


import re
from os import remove, environ, getenv, path
from os.path import exists
from optparse import OptionParser
from shutil import rmtree
from tempfile import mkdtemp
from cogent.app.parameters import Parameter, ValuedParameter, Parameters
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.rdp_classifier import RdpClassifier
from cogent.app.util import CommandLineApplication, CommandLineAppResult, \
    FilePath, ResultPath, guess_input_handler, system,\
    ApplicationNotFoundError, ApplicationError

class RdpClassifier20(CommandLineApplication):
    """RDP Classifier version 2.0 application controller

    The RDP Classifier program is distributed as a java archive (.jar)
    file.  If the file 'rdp_classifier-2.0.jar' is not found in the
    current directory, the app controller looks in the directory
    specified by the environment variable RDP_JAR_PATH.  If this
    variable is not set, and 'rdp_classifier-2.0.jar' is not found in
    the current directory, the application controller raises an
    ApplicationNotFoundError.

    The RDP Classifier often requires memory in excess of Java's
    default 64M. To correct this situation, the authors recommend
    increasing the maximum heap size for the java virtual machine.  An
    option '-Xmx' (default 1000M) is provided for this purpose.
    Details on this option may be found at
    http://java.sun.com/j2se/1.5.0/docs/tooldocs/solaris/java.html

    The classifier may optionally use a custom training set.  The full
    path to the training set may be provided in the option
    '-training-data'.
    """
    _input_handler = '_input_as_multiline_string'
    _command = "rdp_classifier-2.0.jar"
    _options ={}

    # The following are available in the attributes JvmParameters,
    # JarParameters, and PositionalParameters

    _jvm_synonyms = {}
    _jvm_parameters = {
        # Maximum heap size for JVM.
        '-Xmx': ValuedParameter('-', Name='Xmx', Delimiter='', Value='1000m'),
        }
    _positional_synonyms = {}
    _positional_parameters = {
        '-training-data': ValuedParameter('', Name='', Delimiter='', Value='', IsPath=True),
        }

    _parameters = {}
    _parameters.update(_options)
    _parameters.update(_jvm_parameters)
    _parameters.update(_positional_parameters)

    def getHelp(self):
        """Returns documentation string"""
        # Summary paragraph copied from rdp_classifier-2.0, which is
        # licensed under the GPL 2.0 and Copyright 2008 Michigan State
        # University Board of Trustees
        help_str =\
        """
        Ribosomal Database Project - Classifier
        http://rdp.cme.msu.edu/classifier/

        The RDP Classifier is a naive Bayesian classifier which was
        developed to provide rapid taxonomic placement based on rRNA
        sequence data. The RDP Classifier can rapidly and accurately
        classify bacterial 16s rRNA sequences into the new
        higher-order taxonomy proposed by Bergey's Trust. It provides
        taxonomic assignments from domain to genus, with confidence
        estimates for each assignment. The RDP Classifier is not
        limited to using the bacterial taxonomy proposed by the
        Bergey's editors. It worked equally well when trained on the
        NCBI taxonomy. The RDP Classifier likely can be adapted to
        additional phylogenetically coherent bacterial taxonomies.

        The following paper should be cited if this resource is used:

        Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007.
        Naive Bayesian Classifier for Rapid Assignment of rRNA
        Sequences into the New Bacterial Taxonomy.  Appl Environ
        Microbiol. 73(16):5261-7.
        """
        return help_str

    def __call__(self, data=None, remove_tmp=True):
        """Run the application with the specified kwargs on data
        
            data: anything that can be cast into a string or written out to
                a file. Usually either a list of things or a single string or 
                number. input_handler will be called on this data before it 
                is passed as part of the command-line argument, so by creating
                your own input handlers you can customize what kind of data
                you want your application to accept

            remove_tmp: if True, removes tmp files
        """
        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        assignment_fp = FilePath(self.getTmpFilename(self.TmpDir))
        if suppress_stdout:
            outfile = FilePath('/dev/null')
        else:
            outfile = FilePath(self.getTmpFilename(self.TmpDir))
        if suppress_stderr:
            errfile = FilePath('/dev/null')
        else:
            errfile = FilePath(self.getTmpFilename(self.TmpDir))
        if data is None:
            input_arg = ''
        else:
            input_arg = getattr(self,input_handler)(data)

        training_data = self.PositionalParameters['-training-data']

        # Build up the command, consisting of a BaseCommand followed by
        # input and output (file) specifications
        command = self._commandline_join(
            [self.BaseCommand, input_arg, assignment_fp, training_data, 
             '>', outfile, '2>', errfile,]
            )

        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command
        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        exit_status = system(command) >> 8
      
        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            raise ApplicationError, \
             'Unacceptable application exit status: %s, command: %s'\
                % (str(exit_status),command)
        
        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfile,"r")
        err = None        
        if not suppress_stderr:
            err = open(errfile,"r")

        result_paths = self._get_result_paths(data)
        result_paths['Assignments'] = ResultPath(assignment_fp)
        result = CommandLineAppResult(
            out, err, exit_status, result_paths=result_paths)

        # Clean up the input file if one was created
        if remove_tmp:
            if self._input_filename:
                remove(self._input_filename)
                self._input_filename = None

        return result

    def _accept_exit_status(self, status):
        """Returns false if an error occurred in execution
        """
        return (status == 0)

    def _error_on_missing_application(self,params):
        """Raise an ApplicationNotFoundError if the app is not accessible
        """
        command = self._get_jar_fp()
        if not exists(command):
            raise ApplicationNotFoundError,\
             "Cannot find jar file. Is it installed? Is $RDP_JAR_PATH"+\
             " set correctly?"

    def _get_jar_fp(self):
        """Returns the full path to the JAR file.

        Raises an ApplicationError if the JAR file cannot be
        found in the (1) current directory or (2) the path specified
        in the RDP_JAR_PATH environment variable.
        """
        # handles case where the jar file is in the current working directory
        if exists(self._command):
            return self._command
        # handles the case where the user has specified the location via
        # an environment variable
        elif 'RDP_JAR_PATH' in environ:
            return getenv('RDP_JAR_PATH')
        # error otherwise
        else:
            raise ApplicationError,\
             "$RDP_JAR_PATH is not set -- this must be set to use the"+\
             " RDP classifier application controller."

    # Overridden to pull out JVM-specific command-line arguments.
    def _get_base_command(self):
        """Returns the base command plus command-line options.

        Does not include input file, output file, and training set.
        """
        # Necessary? Preserve for consistency.
        if self._command is None:
            raise ApplicationError, '_command has not been set.'

        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])

        jvm_command = "java"
        jvm_arguments = self._commandline_join(self.JvmParameters.values())
        jar_arguments = '-jar "%s"' % self._get_jar_fp()

        result = self._commandline_join(
            [cd_command, jvm_command, jvm_arguments, jar_arguments]
            )
        return result
    
    BaseCommand = property(_get_base_command)

    def _commandline_join(self, tokens):
        """Formats a list of tokens as a shell command

        This seems to be a repeated pattern; may be useful in
        superclass.
        """
        commands = filter(None, map(str, tokens))
        return self._command_delimiter.join(commands).strip()

    @property
    def JvmParameters(self):
        return self.__extract_parameters('jvm')

    @property
    def PositionalParameters(self):
        return self.__extract_parameters('positional')

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
        synonyms   = getattr(self, '_' + name + '_synonyms')
        result = Parameters(parameters, synonyms)
        for key in result.keys():
            result[key] = self.Parameters[key]
        return result


class RdpTrainer20(RdpClassifier20):
    _input_handler = '_input_as_lines'
    TrainingClass = 'edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker'
    PropertiesFile = 'RdpClassifier.properties'

    def __call__(self, training_seqs_file, taxonomy_file, model_output_dir,
        remove_tmp=True):
        return self._train_with_rdp_files(
            training_seqs_file, taxonomy_file, model_output_dir, remove_tmp)

    def _train_with_mapping_file(self, training_seqs_file, lineage_file,
        model_output_dir, remove_tmp=True):
        """Creates a set of training data for the RDP Classifier

            training_seqs_file: The set of training sequences, in
                fasta format.

            lineage_file: A File-like object that specifies a lineage
                for each sequence. Each line must contain a Sequence
                ID, followed by a tab, then followed by the assigned
                lineage.  The taxa comprising the lineage must be
                separated with a comma.

            model_output_dir: Directory in which to store training data.

            remove_tmp: if True, removes tmp files

        To use the resulting model with the RdpClassifier, set
        '-training_data' to the following path: model_output_dir +
        RdpClassifier.PropertiesFile
        """

    def _train_with_rdp_files(self, training_seqs_file, taxonomy_file, 
        model_output_dir, remove_tmp=True):
        """Creates a set of training data for the RDP Classifier

            training_seqs_file: A pre-classified set of training
                sequences, in fasta-like format.  Each sequence must
                be labelled with an identifier (no spaces) and an
                assigned lineage (taxa separated by ';'). Example of
                a valid label: ">seq1 ROOT;Ph1;Fam1;G1;"

            taxonomy_file: A File-like object that specifies a
                taxonomic heirarchy. Each line in the file must
                contain a '*'-separated list of the following items:
                Taxon ID, Taxon Name, Parent Taxon ID, Depth, and
                Rank.  IDs should have an integer format.  Example of
                a valid line: "1*Bacteria*0*0*domain"

            model_output_dir: Directory in which to store training data.

            remove_tmp: if True, removes tmp files

        To use the resulting model with the RdpClassifier, set
        '-training_data' to the following path: model_output_dir +
        RdpClassifier.PropertiesFile
        """
        # Three extra pieces of information are required to create
        # training data.  Unless we want built-in support for
        # versioned training sets, these may be set to sensible
        # defaults.
        training_set_id = '1'
        taxonomy_version = 'version1'
        modification_info = 'cogent'

        # The properties file specifies the names of the files in the
        # training directory.  We use the example properties file
        # directly from the rdp_classifier distribution, which lists
        # the default set of files created by the application.  We
        # must write this file explicitly after generating the
        # training data.
        properties = (
            "# Sample ResourceBundle properties file\n"
            "bergeyTree=bergeyTrainingTree.xml\n"
            "probabilityList=genus_wordConditionalProbList.txt\n"
            "probabilityIndex=wordConditionalProbIndexArr.txt\n"
            "wordPrior=logWordPrior.txt\n"
            "classifierVersion=Naive Bayesian rRNA Classifier Version 1.0, November 2003\n"
            )

        input_handler = self.InputHandler
        suppress_stdout = self.SuppressStdout
        suppress_stderr = self.SuppressStderr
        if suppress_stdout:
            outfile = FilePath('/dev/null')
        else:
            outfile = self.getTmpFilename(self.TmpDir)
        if suppress_stderr:
            errfile = FilePath('/dev/null')
        else:
            errfile = FilePath(self.getTmpFilename(self.TmpDir))

        input_handler_function = getattr(self, input_handler)
        taxonomy_filename = input_handler_function(taxonomy_file)
        training_seqs_filename = input_handler_function(training_seqs_file)

        # Build up the command, consisting of a BaseCommand followed
        # by input and output (file) specifications 

        # Example from rdp_classifier/sampledata/README: 
        # java -Xmx400m -cp rdp_classifier-2.0.jar
        # edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker
        # mydata/mytaxon.txt mydata/mytrainseq.fasta 1 version1 test
        # mydata
        command = self._commandline_join(
            [self.BaseCommand, taxonomy_filename, training_seqs_filename,
             training_set_id, taxonomy_version, modification_info,
             model_output_dir, '>', outfile, '2>', errfile]
            )

        if self.HaltExec: 
            raise AssertionError, "Halted exec with command:\n" + command
        # The return value of system is a 16-bit number containing the signal 
        # number that killed the process, and then the exit status. 
        # We only want to keep the exit status so do a right bitwise shift to 
        # get rid of the signal number byte
        exit_status = system(command) >> 8

        # Determine if error should be raised due to exit status of 
        # appliciation
        if not self._accept_exit_status(exit_status):
            raise ApplicationError, \
             'Unacceptable application exit status: %s, command: %s'\
                % (str(exit_status),command)

        # must write properties file to output directory manually
        properties_fp = path.join(model_output_dir, self.PropertiesFile)
        properties_file = open(properties_fp, 'w')
        properties_file.write(properties)
        properties_file.close()

        # open the stdout and stderr if not being suppressed
        out = None
        if not suppress_stdout:
            out = open(outfile,"r")
        err = None        
        if not suppress_stderr:
            err = open(errfile,"r")
       
        result = CommandLineAppResult(out, err, exit_status, 
            result_paths=self._get_result_paths(model_output_dir))

        # Clean up the input files
        if remove_tmp:
            remove(taxonomy_filename)
            remove(training_seqs_filename)

        return result

    def _input_as_lines(self, data):
        """ Write a seq of lines to a temp file and return the filename string.

        This method has been overridden for RdpTrainer so that the
        _input_filename attribute is not assigned.

            data: a sequence to be written to a file, each element of the 
                sequence will compose a line in the file
           * Note: the result will be the filename as a FilePath object 
            (which is a string subclass).

           * Note: '\n' will be stripped off the end of each sequence element
                before writing to a file in order to avoid multiple new lines
                accidentally be written to a file
        """
        filename = FilePath(self.getTmpFilename(self.TmpDir))
        data_file = open(filename, 'w')
        # Parent method does not take advantage of laziness, due to
        # temporary variable that contains entire file contents --
        # better to write explicit loop over lines in the data source,
        # storing only each line in turn.
        for line in data:
            line = str(line).strip('\n')
            data_file.write(line)
            data_file.write('\n')
        data_file.close()
        return filename

    def _get_result_paths(self, output_dir):
        files = {
            'bergeyTree': 'bergeyTrainingTree.xml',
            'probabilityList': 'genus_wordConditionalProbList.txt',
            'probabilityIndex': 'wordConditionalProbIndexArr.txt',
            'wordPrior': 'logWordPrior.txt',
            'properties': self.PropertiesFile,
        }
        result_paths = {}
        for name, file in files.iteritems():
            result_paths[name] = ResultPath(
                Path=path.join(output_dir, file), IsWritten=True)
        return result_paths
    
    # Overridden to pull out JVM-specific command-line arguments.
    def _get_base_command(self):
        """Returns the base command plus command-line options.

        Does not include input file, output file, and training set.
        """
        # Necessary? Preserve for consistency.
        if self._command is None:
            raise ApplicationError, '_command has not been set.'

        # Append a change directory to the beginning of the command to change 
        # to self.WorkingDir before running the command
        # WorkingDir should be in quotes -- filenames might contain spaces
        cd_command = ''.join(['cd ',str(self.WorkingDir),';'])

        jvm_command = "java"
        jvm_arguments = self._commandline_join(self.JvmParameters.values())
        jar_arguments = '-cp "%s"' % self._get_jar_fp()

        result = self._commandline_join(
            [cd_command, jvm_command, jvm_arguments, jar_arguments, self.TrainingClass]
            )
        return result

    BaseCommand = property(_get_base_command)


def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] input_sequences_filepath'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)
          
    parser.add_option('-o','--output_fp',action='store',\
          type='string',dest='output_fp',help='Path to store '+\
          'output file [default: generated from input_sequences_filepath]')
          
    parser.add_option('-c','--min_confidence',action='store',\
          type='float',dest='min_confidence',help='minimum confidence '+\
          'level to return a classification [default: %default]')

    parser.set_defaults(verbose=False,min_confidence=0.80)

    opts,args = parser.parse_args()
    num_args = 1
    if len(args) != num_args:
       parser.error('Exactly one argument is required.')


    return opts,args
    
def assign_taxonomy(data, min_confidence=0.80, output_fp=None,
    training_data_fp=None, max_memory=None):
    """ Assign taxonomy to each sequence in data with the RDP classifier 
    
        data: open fasta file object or list of fasta lines
        confidence: minimum support threshold to assign taxonomy to a sequence
        output_fp: path to write output; if not provided, result will be 
         returned in a dict of {seq_id:(taxonomy_assignment,confidence)}
    
    """
    data = list(data)
    
    # build a map of seq identifiers as the RDP classifier doesn't 
    # preserve these perfectly
    identifier_lookup = {}
    for seq_id, seq in MinimalFastaParser(data):
        identifier_lookup[seq_id.split()[0]] = seq_id
    
    # build the classifier object
    app = RdpClassifier20()
    if max_memory is not None:
        app.Parameters['-Xmx'].on(max_memory)
    if training_data_fp is not None:
        app.Parameters['-training-data'].on(training_data_fp)

    # apply the rdp app controller
    rdp_result = app('\n'.join(data))
    # grab assignment output
    result_lines = rdp_result['Assignments']
    
    # start a list to store the assignments
    results = {}

    # ShortSequenceException messages are written to stdout
    # Tag these ID's as unassignable
    stdout_lines = rdp_result['StdOut']
    for line in stdout_lines:
        if line.startswith('ShortSequenceException'):
            matchobj = re.search('recordID=(\S+)', line)
            if matchobj:
                rdp_id = matchobj.group(1)
                orig_id = identifier_lookup[rdp_id]
                results[orig_id] = ('Unassignable', 1.0)
    
    # iterate over the identifier, assignment strings (this is a bit
    # of an abuse of the MinimalFastaParser, as these are not truely
    # fasta lines)
    for identifier, assignment_str in MinimalFastaParser(result_lines):
        # get the original identifier from the one in the rdp result
        identifier = identifier_lookup[\
         identifier[:identifier.index('reverse=')].strip()]
        # build a list to store the assignments we're confident in
        # (i.e., the ones that have a confidence greater than min_confidence)
        confident_assignments = []
        # keep track of the lowest acceptable confidence value that
        # has been encountered
        lowest_confidence = 0.0
        
        # split the taxonomy assignment string
        assignment_fields = assignment_str.split(';')
        # iterate over (assignment, assignment confidence) pairs
        for i in range(0,len(assignment_fields),2):
            assignment = assignment_fields[i]
            try:
                assignment_confidence = float(assignment_fields[i+1])
            except IndexError:
                break
            # check the confidence of the current assignment
            if assignment_confidence >= min_confidence:
                # if the current assignment confidence is greater than
                # the min, store the assignment and confidence value
                confident_assignments.append(assignment.strip())
                lowest_confidence = assignment_confidence 
            else:
                # otherwise, we've made it to the lowest assignment that
                # met the confidence threshold, so bail out of the loop
                break

        # store the identifier, the semi-colon-separated assignments, and the
        # confidence for the last assignment
        results[identifier] = \
             (';'.join(confident_assignments),lowest_confidence)
            
    if output_fp:
        try:
            output_file = open(output_fp,'w')
        except OSError:
            raise OSError, "Can't open output file for writing: %s" % output_fp
            
        for seq_id, values in results.items():
            output_file.write('%s\t%s\t%1.3f\n' % (seq_id,values[0],values[1]))
            
        output_file.close()
        return None
    else:   
        return results


def train_rdp_classifier(training_seqs_file, taxonomy_file, model_output_dir,
    max_memory=None):
    """ Train RDP Classifier, saving to model_output_dir

        training_seqs_file, taxonomy_file: file-like objects used to
            train the RDP Classifier (see RdpTrainer documentation for
            format of training data)

        model_output_dir: directory in which to save the files
            necessary to classify sequences according to the training
            data

    Once the model data has been generated, the RDP Classifier may 
    """
    app = RdpTrainer20()
    if max_memory is not None:
        app.Parameters['-Xmx'].on(max_memory)
    return app(training_seqs_file, taxonomy_file, model_output_dir)


def train_rdp_classifier_and_assign_taxonomy(
    training_seqs_file, taxonomy_file, seqs_to_classify, min_confidence=0.80, 
    model_output_dir=None, classification_output_fp=None, max_memory=None):
    """ Train RDP Classifier and assign taxonomy in one fell swoop

    The file objects training_seqs_file and taxonomy_file are used to
    train the RDP Classifier (see RdpTrainer documentation for
    details).  Model data is stored in model_output_dir.  If
    model_output_dir is not provided, a temporary directory is created
    and removed after classification.

    The sequences in seqs_to_classify are classified according to the
    model and filtered at the desired confidence level (default:
    0.80).

    The results are saved to classification_output_fp if provided,
    otherwise a dict of {seq_id:(taxonomy_assignment,confidence)} is
    returned.
    """
    if model_output_dir is None:
        training_dir = mkdtemp(prefix='RdpTrainer_')
    else:
        training_dir = model_output_dir

    trainer = RdpTrainer20()
    if max_memory is not None:
        trainer.Parameters['-Xmx'].on(max_memory)
    training_results = trainer(
        training_seqs_file, taxonomy_file, training_dir)

    training_data_fp = training_results['properties'].name
    assignment_results = assign_taxonomy(
        seqs_to_classify, min_confidence=min_confidence, 
        output_fp=classification_output_fp, training_data_fp=training_data_fp,
        max_memory=max_memory)

    if model_output_dir is None:
        rmtree(training_dir)

    return assignment_results


if __name__ == "__main__":
    
    opts,args = parse_command_line_parameters()
    assign_taxonomy(open(args[0]),min_confidence=opts.min_confidence,\
     output_fp=opts.output_fp)

