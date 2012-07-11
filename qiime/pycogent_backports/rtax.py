#!/usr/bin/env python
"""Application controller for RTAX version 1.0

Includes application controller for RTAX.

Modified from uclust.py and rdp_classifier.py on 12-27-11
"""

__author__ = "David Soergel"
__copyright__ = "Copyright 2007-2011, The PyCogent Project"
__credits__ = ["David Soergel", "William Walters", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "David Soergel"
__email__ = "soergel@cs.umass.edu"
__status__ = "Development"

from os import remove, makedirs
from os.path import exists, split, splitext, basename, isdir, abspath, isfile
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename, ApplicationError, ApplicationNotFoundError
from cogent.util.misc import remove_files, app_path
from cogent import DNA
import tempfile
import os.path
import re
from sys import stderr

from qiime.util import (load_qiime_config, get_tmp_filename)
from shutil import rmtree

class RtaxParseError(Exception):
    pass

class Rtax(CommandLineApplication):
    """ Rtax ApplicationController

    """

    _command = 'rtax'
    _input_handler = '_input_as_parameters'
    _parameters = {\

        # -r a reference database in FASTA format
        '-r':ValuedParameter('-',Name='r',Delimiter=' ', IsPath=True),

        # -t a taxonomy file with sequence IDs matching the reference database
        '-t':ValuedParameter('-',Name='t',Delimiter=' ', IsPath=True),

        # -a a FASTA file containing query sequences (single-ended, read 1, or paired-end delimited)
        '-a':ValuedParameter('-',Name='a',Delimiter=' ', IsPath=True),

        # -b a FASTA file containing query sequences (read 2, with matching IDs)
        '-b':ValuedParameter('-',Name='b',Delimiter=' ', IsPath=True),

        # -l a text file containing sequence IDs to process, one per line
        '-l':ValuedParameter('-',Name='l',Delimiter=' ', IsPath=True),

        # -d a delimiter separating the two reads when provided in a single file
        '-d':ValuedParameter('-',Name='d',Delimiter=' ', IsPath=False, Quote="\""),

        # -i a regular expression used to select part of the fasta header to use as the sequence id.
        '-i':ValuedParameter('-',Name='i',Delimiter=' ', IsPath=False, Quote="'"),

        # -o output file name for classification assignment
        '-o': ValuedParameter('-', Name='o', Delimiter=' ', IsPath=True),

        # -m temporary directory
        '-m': ValuedParameter('-', Name='m', Delimiter=' ', IsPath=True),

        # -f allow fallback from paired-end to single-ended classification when one read is missing
        '-f':FlagParameter(Prefix='-',Name='f'),

        # -g do not allow fallback from paired-end to single-ended classification when one read is too generic
        '-g':FlagParameter(Prefix='-',Name='g')
    }

    _suppress_stdout = False
    _suppress_stderr = False

    #def __init__(self):
    #    super().__init__()...
    #    usearch_command = "usearch"
    #    if not (exists(usearch_command) or app_path(usearch_command)):
    #        raise ApplicationNotFoundError,\
 	#        "Cannot find %s. Is it installed? Is it in your path?"\
 	#        % usearch_command


    def _input_as_parameters(self,data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['-r','-t','-a','-b','-l','-d','i','-o','-m','-v','-f', '-g']

        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError,\
             "Unsupported parameter(s) passed when calling rtax: %s" %\
              ' '.join(unsupported_parameters)

        for v in allowed_values:
            # turn the parameter off so subsequent runs are not
            # affected by parameter settings from previous runs
            self.Parameters[v].off()
            if v in data:
                # turn the parameter on if specified by the user
                self.Parameters[v].on(data[v])

        return ''

    def _get_result_paths(self,data):
        """ Return a dict of ResultPath objects representing all possible output
        """
        assignment_fp = str(self.Parameters['-o'].Value).strip('"')
        if not os.path.isabs(assignment_fp):
            assignment_fp = os.path.relpath(assignment_fp, self.WorkingDir)
        return {'Assignments': ResultPath(assignment_fp, IsWritten=True)}



    def _accept_exit_status(self,exit_status):
        """ Test for acceptable exit status

            uclust can seg fault and still generate a parsable .uc file
            so we explicitly check the exit status

        """
        return exit_status == 0

    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        RTAX is hosted at:
        http://dev.davidsoergel.com/rtax/

        The following paper should be cited if this resource is used:

        Soergel D.A.W., Dey N., Knight R., and Brenner S.E.  2012.
        Selection of primers for optimal taxonomic classification
        of environmental 16S rRNA gene sequences.  ISME J (6), 1440-1444
        """
        return help_str

def assign_taxonomy(dataPath, reference_sequences_fp, id_to_taxonomy_fp, read_1_seqs_fp, read_2_seqs_fp, single_ok=False, no_single_ok_generic=False,
                    header_id_regex=None, read_id_regex = "\S+\s+(\S+)", amplicon_id_regex = "(\S+)\s+(\S+?)\/",
                    output_fp=None, log_path=None, HALT_EXEC=False):
    """Assign taxonomy to each sequence in data with the RTAX classifier

        # data: open fasta file object or list of fasta lines
        dataPath: path to a fasta file

        output_fp: path to write output; if not provided, result will be
         returned in a dict of {seq_id:(taxonomy_assignment,confidence)}
    """

    usearch_command = "usearch"
    if not (exists(usearch_command) or app_path(usearch_command)):
        raise ApplicationNotFoundError,\
         "Cannot find %s. Is it installed? Is it in your path?"\
         % usearch_command

    qiime_config = load_qiime_config()
    base_tmp_dir = qiime_config['temp_dir'] or '/tmp/'
    my_tmp_dir = get_tmp_filename(tmp_dir=base_tmp_dir,prefix='rtax_',suffix='',result_constructor=str)
    os.makedirs(my_tmp_dir)


    try:
        # RTAX classifier doesn't necessarily preserve identifiers
        # it reports back only the id extracted as $1 using header_id_regex
        # since rtax takes the original unclustered sequence files as input,
        # the usual case is that the regex extracts the amplicon ID from the second field



        # Use lookup table
        read_1_id_to_orig_id = {}
        readIdExtractor = re.compile(read_id_regex)  # OTU clustering produces ">clusterID read_1_id"
        data = open(dataPath,'r')
        for seq_id, seq in MinimalFastaParser(data):
            # apply the regex
            extract = readIdExtractor.match(seq_id)
            if extract is None:
                stderr.write("Matched no ID with read_id_regex " + read_id_regex +" in '" + seq_id + "' from file " + dataPath + "\n")
            else:
                read_1_id_to_orig_id[extract.group(1)] = seq_id
                #stderr.write(extract.group(1) + " => " +  seq_id + "\n")
            #seq_id_lookup[seq_id.split()[1]] = seq_id
        data.close()



        # make list of amplicon IDs to pass to RTAX

        id_list_fp = open(my_tmp_dir+"/ampliconIdsToClassify", "w")

        # Establish mapping of amplicon IDs to read_1 IDs
        # simultaneously write the amplicon ID file for those IDs found in the input mapping above

        amplicon_to_read_1_id = {}
        ampliconIdExtractor = re.compile(amplicon_id_regex)  # split_libraries produces >read_1_id ampliconID/1 ...  // see also assign_taxonomy 631
        read_1_data = open(read_1_seqs_fp,'r')
        for seq_id, seq in MinimalFastaParser(read_1_data):
            # apply the regex
            extract = ampliconIdExtractor.match(seq_id)
            if extract is None:
                stderr.write("Matched no ID with amplicon_id_regex " + amplicon_id_regex + " in '" + seq_id + "' from file " + read_1_seqs_fp + "\n")
            else:
                read_1_id = extract.group(1)
                amplicon_id = extract.group(2)
                try:
                    amplicon_to_read_1_id[amplicon_id] = read_1_id
                    bogus = read_1_id_to_orig_id[read_1_id]  # verify that the id is valid
                    id_list_fp.write('%s\n' % (amplicon_id))
                except KeyError:
                    pass
        data.close()
        id_list_fp.close()

        app = Rtax(HALT_EXEC=HALT_EXEC)

        temp_output_file = tempfile.NamedTemporaryFile(
            prefix='RtaxAssignments_', suffix='.txt')
        app.Parameters['-o'].on(temp_output_file.name)
        app.Parameters['-r'].on(reference_sequences_fp)
        app.Parameters['-t'].on(id_to_taxonomy_fp)
        # app.Parameters['-d'].on(delimiter)
        app.Parameters['-l'].on(id_list_fp.name)  # these are amplicon IDs
        app.Parameters['-a'].on(read_1_seqs_fp)
        if read_2_seqs_fp is not None:
            app.Parameters['-b'].on(read_2_seqs_fp)
        app.Parameters['-i'].on(header_id_regex)
        app.Parameters['-m'].on(my_tmp_dir)
        if single_ok: app.Parameters['-f'].on();
        if no_single_ok_generic: app.Parameters['-g'].on();
        #app.Parameters['-v'].on()

        app_result = app()

        if log_path:
            f=open(log_path, 'a')
            errString=''.join(app_result['StdErr'].readlines()) + '\n'
            f.write(errString)
            f.close()

        assignments = {}

        # restore original sequence IDs with spaces

        for line in app_result['Assignments']:
            toks = line.strip().split('\t')
            rtax_id = toks.pop(0)
            if len(toks):
                bestpcid = toks.pop(0)  # ignored
            lineage = toks

            # RTAX does not provide a measure of confidence.  We could pass one in,
            # based on the choice of primers, or even look it up on the fly in the tables
            # from the "optimal primers" paper; but it would be the same for every
            # query sequence anyway.
            # we could also return bestpcid, but that's not the same thing as confidence.
            confidence = 1.0

            read_1_id = amplicon_to_read_1_id[rtax_id]
            orig_id = read_1_id_to_orig_id[read_1_id]
            if lineage:
                assignments[orig_id] = (';'.join(lineage), confidence)
            else:
                assignments[orig_id] = ('Unclassified', 1.0)

        if output_fp:
            try:
                output_file = open(output_fp, 'w')
            except OSError:
                raise OSError("Can't open output file for writing: %s" % output_fp)
            for seq_id, assignment in assignments.items():
                lineage, confidence = assignment
                output_file.write(
                    '%s\t%s\t%1.3f\n' % (seq_id, lineage, confidence))
            output_file.close()
            return None
        else:
            return assignments
    finally:
        try:
            rmtree(my_tmp_dir)
        except OSError:
            pass
