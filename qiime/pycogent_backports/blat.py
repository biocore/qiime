#!/usr/bin/env python

"""Application controller for BLAT v34"""

from cogent.app.parameters import FlagParameter, ValuedParameter, \
    MixedParameter, FilePath
from cogent.app.util import CommandLineApplication, ResultPath, \
    ApplicationError, get_tmp_filename
from cogent import DNA
from cogent.core.genetic_code import GeneticCodes
from cogent.parse.fasta import MinimalFastaParser
from os import remove
from os.path import isabs

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2007-2012, The QIIME Project"
__credits__ = ["Adam Robbins-Pianka", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Prototype"


class Blat(CommandLineApplication):

    """BLAT generic application controller"""

    _command = 'blat'
    _input_handler = "_input_as_list"

    _database_types = ['dna', 'prot', 'dnax']
    _query_types = ['dna', 'rna', 'prot', 'dnax', 'rnax']
    _mask_types = ['lower', 'upper', 'out', 'file.out']
    _out_types = ['psl', 'pslx', 'axt', 'maf', 'sim4', 'wublast', 'blast',
                  'blast8', 'blast9']
    _valid_combinations = [('dna', 'dna'), ('dna', 'rna'), ('prot', 'prot'),
                           ('dnax', 'prot'), ('dnax', 'dnax'),
                           ('dnax', 'rnax')]
    _database = None
    _query = None
    _output = None

    _parameters = {
        # database type (dna, prot, or dnax, where dnax is DNA sequence
        # translated in six frames to protein
        '-t': ValuedParameter('-', Delimiter='=', Name='t'),

        # query type (dna, rna, prot, dnax, rnax, where rnax is DNA sequence
        # translated in three frames to protein
        '-q': ValuedParameter('-', Delimiter='=', Name='q'),

        # Use overused tile file N.ooc, and N should correspond to the tileSize
        '-ooc': ValuedParameter('-', Delimiter='=', Name='ooc', IsPath=True),

        # Sets the size of at match that that triggers an alignment
        '-tileSize': ValuedParameter('-', Delimiter='=', Name='tileSize'),

        # Spacing between tiles.
        '-stepSize': ValuedParameter('-', Delimiter='=', Name='stepSize'),

        # If set to 1, allows one mismatch in the tile and still triggers
        # an alignment.
        '-oneOff': ValuedParameter('-', Delimiter='=', Name='oneOff'),

        # sets the number of tile matches
        '-minMatch': ValuedParameter('-', Delimiter='=', Name='minMatch'),

        # sets the minimum score
        '-minScore': ValuedParameter('-', Delimiter='=', Name='minScore'),

        # sets the minimum sequence identity in percent
        '-minIdentity':
        ValuedParameter('-', Delimiter='=', Name='minIdentity'),

        # sets the size o the maximum gap between tiles in a clump
        '-maxGap': ValuedParameter('-', Delimiter='=', Name='maxGap'),

        # make an overused tile file. Target needs to be complete genome.
        '-makeOoc': ValuedParameter('-', Delimiter='=', Name='makeOoc',
                                    IsPath=True),

        # sets the number of repetitions of a tile allowed before it is marked
        # as overused
        '-repMatch': ValuedParameter('-', Delimiter='=', Name='repMatch'),

        # mask out repeats.  Alignments won't be started in masked region but
        # may extend through it in nucleotide searches.  Masked areas are
        # ignored entirely in protein or translated searches.  Types are:
        # lower, upper, out, file.out (file.out - mask database according to
        # RepeatMasker file.out
        '-mask': ValuedParameter('-', Delimiter='=', Name='mask'),

        # Mask out repeats in query sequence.  similar to -mask but for query
        # rather than target sequence
        '-qMask': ValuedParameter('-', Delimiter='=', Name='qMask'),

        # repeat bases will not be masked in any way, but matches in
        # repeat areas will be reported separately from matches in other
        # areas in the pls output
        '-repeats': ValuedParameter('-', Delimiter='=', Name='repeats'),

        # minimum percent divergence of repeats to allow them to be unmasked
        '-minRepDivergence': ValuedParameter('-', Delimiter='=',
                                             Name='minRepDivergence'),

        # output dot every N sequences to show program's progress
        '-dots': ValuedParameter('-', Delimiter='=', Name='dots'),

        # controls output file format.  One of:
        # psl - Default.  Tab separated format, no sequence
        # pslx - Tab separated format with sequence
        # axt - blastz-associated axt format
        # maf - multiz-associated maf format
        # sim4 - similar to sim4 format
        # wublast - similar to wublast format
        # blast - similar to NCBI blast format
        # blast8- NCBI blast tabular format
        # blast9 - NCBI blast tabular format with comments
        '-out': ValuedParameter('-', Delimiter='=', Name='out'),

        # sets maximum intron size
        '-maxIntron': ValuedParameter('-', Delimiter='=', Name='maxIntron'),

        # suppress column headers in psl output
        '-noHead': FlagParameter('-', Name='noHead'),

        # trim leading poly-T
        '-trimT': FlagParameter('-', Name='trimT'),

        # do not trim trailing poly-A
        '-noTrimA': FlagParameter('-', Name='noTrimA'),

        # Remove poly-A tail from qSize as well as alignments in psl output
        '-trimHardA': FlagParameter('-', Name='trimHardA'),

        # run for fast DNA/DNA remapping - not allowing introns,
        # requiring high %ID
        '-fastMap': FlagParameter('-', Name='fastMap'),

        # for high quality mRNAs, look harder for small initial and terminal
        # exons
        '-fine': FlagParameter('-', Name='fine'),

        # Allows extension of alignment through large blocks of N's
        '-extendThroughN': FlagParameter('-', Name='extendThroughN')
    }

    def _get_result_paths(self, data):
        """Returns the file location for result output
        """

        return {'output': ResultPath(data[2], IsWritten=True)}

    def _get_base_command(self):
        """Gets the command that will be run when the app controller is
        called.
        """
        command_parts = []
        cd_command = ''.join(['cd ', str(self.WorkingDir), ';'])
        if self._command is None:
            raise ApplicationError('_command has not been set.')
        command = self._command
        parameters = sorted([str(x) for x in self.Parameters.values()
                            if str(x)])

        synonyms = self._synonyms

        command_parts.append(cd_command)
        command_parts.append(command)
        command_parts.append(self._database)  # Positional argument
        command_parts.append(self._query)  # Positional argument
        command_parts += parameters
        if self._output:
            command_parts.append(self._output.Path)  # Positional

        return (
            self._command_delimiter.join(filter(None, command_parts)).strip()
        )

    BaseCommand = property(_get_base_command)

    def _input_as_list(self, data):
        '''Takes the positional arguments as input in a list.

        The list input here should be [query_file_path, database_file_path,
        output_file_path]'''
        query, database, output = data
        if (not isabs(database)) \
                or (not isabs(query)) \
                or (not isabs(output)):
            raise ApplicationError("Only absolute paths allowed.\n%s" %
                                   ', '.join(data))

        self._database = FilePath(database)
        self._query = FilePath(query)
        self._output = ResultPath(output, IsWritten=True)

        # check parameters that can only take a particular set of values
        # check combination of databse and query type
        if self.Parameters['-t'].isOn() and self.Parameters['-q'].isOn() and \
                (self.Parameters['-t'].Value, self.Parameters['-q'].Value) not in \
                self._valid_combinations:
            error_message = "Invalid combination of database and query " + \
                            "types ('%s', '%s').\n" % \
                            (self.Paramters['-t'].Value,
                             self.Parameters['-q'].Value)

            error_message += "Must be one of: %s\n" % \
                             repr(self._valid_combinations)

            raise ApplicationError(error_message)

        # check database type
        if self.Parameters['-t'].isOn() and \
                self.Parameters['-t'].Value not in self._database_types:
            error_message = "Invalid database type %s\n" % \
                            self.Parameters['-t'].Value

            error_message += "Allowed values: %s\n" % \
                             ', '.join(self._database_types)

            raise ApplicationError(error_message)

        # check query type
        if self.Parameters['-q'].isOn() and \
                self.Parameters['-q'].Value not in self._query_types:
            error_message = "Invalid query type %s\n" % \
                            self.Parameters['-q'].Value

            error_message += "Allowed values: %s\n" % \
                ', '.join(self._query_types)

            raise ApplicationError(error_message)

        # check mask type
        if self.Parameters['-mask'].isOn() and \
                self.Parameters['-mask'].Value not in self._mask_types:
            error_message = "Invalid mask type %s\n" % \
                            self.Parameters['-mask']

            error_message += "Allowed Values: %s\n" % \
                ', '.join(self._mask_types)

            raise ApplicationError(error_message)

        # check qmask type
        if self.Parameters['-qMask'].isOn() and \
                self.Parameters['-qMask'].Value not in self._mask_types:
            error_message = "Invalid qMask type %s\n" % \
                            self.Parameters['-qMask'].Value

            error_message += "Allowed values: %s\n" % \
                             ', '.join(self._mask_types)

            raise ApplicationError(error_message)

        # check repeat type
        if self.Parameters['-repeats'].isOn() and \
                self.Parameters['-repeats'].Value not in self._mask_types:
            error_message = "Invalid repeat type %s\n" % \
                            self.Parameters['-repeat'].Value

            error_message += "Allowed values: %s\n" % \
                             ', '.join(self._mask_types)

            raise ApplicationError(error_message)

        # check output format
        if self.Parameters['-out'].isOn() and \
                self.Parameters['-out'].Value not in self._out_types:
            error_message = "Invalid output type %s\n" % \
                            self.Parameters['-out']

            error_message += "Allowed values: %s\n" % \
                             ', '.join(self._out_types)

            raise ApplicationError(error_message)

        return ''


def assign_reads_to_database(query_fasta_fp, database_fasta_fp, output_fp,
                             params=None):
    """Assign a set of query sequences to a reference database

    query_fasta_fp : absolute file path to query sequences
    database_fasta_fp : absolute file path to the reference database
    output_fp : absolute file path of the output file to write
    params : dict of BLAT specific parameters.

    This method returns an open file object. The output format
    defaults to blast9 and should be parsable by the PyCogent BLAST parsers.
    """
    if params is None:
        params = {}
    if '-out' not in params:
        params['-out'] = 'blast9'
    blat = Blat(params=params)

    result = blat([query_fasta_fp, database_fasta_fp, output_fp])
    return result['output']


def assign_dna_reads_to_dna_database(query_fasta_fp, database_fasta_fp,
                                     output_fp, params=None):
    """Assign DNA reads to a database fasta of DNA sequences.

    Wraps assign_reads_to_database, setting database and query types. All
    parameters are set to default unless params is passed.

    query_fasta_fp: absolute path to the query fasta file containing DNA
                   sequences.
    database_fasta_fp: absolute path to the database fasta file containing
                      DNA sequences.
    output_fp: absolute path where the output file will be generated.
    params: optional. dict containing parameter settings to be used
                  instead of default values. Cannot change database or query
                  file types from dna and dna, respectively.

    This method returns an open file object. The output format
    defaults to blast9 and should be parsable by the PyCogent BLAST parsers.
    """
    if params is None:
        params = {}

    my_params = {'-t': 'dna',
                 '-q': 'dna'
                 }

    # if the user specified parameters other than default, then use them.
    # However, if they try to change the database or query types, raise an
    # applciation error.
    if '-t' in params or '-q' in params:
        raise ApplicationError("Cannot change database or query types when " +
                               "using assign_dna_reads_to_dna_database. " +
                               "Use assign_reads_to_database instead.\n")

    my_params.update(params)

    result = assign_reads_to_database(query_fasta_fp, database_fasta_fp,
                                      output_fp, my_params)

    return result


def assign_dna_reads_to_protein_database(query_fasta_fp, database_fasta_fp,
                                         output_fp, temp_dir="/tmp", params=None):
    """Assign DNA reads to a database fasta of protein sequences.

    Wraps assign_reads_to_database, setting database and query types. All
    parameters are set to default unless params is passed. A temporary
    file must be written containing the translated sequences from the input
    query fasta file because BLAT cannot do this automatically.

    query_fasta_fp: absolute path to the query fasta file containing DNA
                   sequences.
    database_fasta_fp: absolute path to the database fasta file containing
                      protein sequences.
    output_fp: absolute path where the output file will be generated.
    temp_dir: optional. Change the location where the translated sequences
              will be written before being used as the query. Defaults to
              /tmp.
    params: optional. dict containing parameter settings to be used
                  instead of default values. Cannot change database or query
                  file types from protein and dna, respectively.

    This method returns an open file object. The output format
    defaults to blast9 and should be parsable by the PyCogent BLAST parsers.
    """
    if params is None:
        params = {}

    my_params = {'-t': 'prot', '-q': 'prot'}

    # make sure temp_dir specifies an absolute path
    if not isabs(temp_dir):
        raise ApplicationError("temp_dir must be an absolute path.")

    # if the user specified parameters other than default, then use them.
    # However, if they try to change the database or query types, raise an
    # applciation error.
    if '-t' in params or '-q' in params:
        raise ApplicationError("Cannot change database or query types "
                               "when using assign_dna_reads_to_dna_database. Use "
                               "assign_reads_to_database instead.")

    if 'genetic_code' in params:
        my_genetic_code = GeneticCodes[params['genetic_code']]
        del params['genetic_code']
    else:
        my_genetic_code = GeneticCodes[1]

    my_params.update(params)

    # get six-frame translation of the input DNA sequences and write them to
    # temporary file.
    tmp = get_tmp_filename(tmp_dir=temp_dir, result_constructor=str)
    tmp_out = open(tmp, 'w')

    for label, sequence in MinimalFastaParser(open(query_fasta_fp)):
        seq_id = label.split()[0]

        s = DNA.makeSequence(sequence)
        translations = my_genetic_code.sixframes(s)
        frames = [1, 2, 3, -1, -2, -3]
        translations = dict(zip(frames, translations))

        for frame, translation in sorted(translations.iteritems()):
            entry = '>{seq_id}_frame_{frame}\n{trans}\n'
            entry = entry.format(seq_id=seq_id, frame=frame, trans=translation)
            tmp_out.write(entry)

    tmp_out.close()
    result = assign_reads_to_database(tmp, database_fasta_fp, output_fp,
                                      params=my_params)

    remove(tmp)

    return result
