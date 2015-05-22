#!/usr/bin/env python
# File created on 05 Oct 2009.

from __future__ import division

from os.path import split, basename, abspath, exists, join
from subprocess import PIPE, Popen

from skbio.util import remove_files

from burrito.parameters import ValuedParameter
from burrito.util import (which, CommandLineApplication, ResultPath,
                            ApplicationError, ApplicationNotFoundError)
from skbio.parse.sequences import parse_fasta

from qiime.util import (FunctionWithParams, write_degapped_fasta_to_file,
                        split_sequence_file_on_sample_ids_to_files)
from qiime.assign_taxonomy import BlastTaxonAssigner

from bfillings.formatdb import build_blast_db_from_fasta_path
from bfillings.usearch import (usearch61_smallmem_cluster,
                            usearch61_chimera_check_denovo,
                            parse_usearch61_clusters,
                            usearch61_chimera_check_ref)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "William Walters",
               "Jai Ram Rideout", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

# NOTE: The next ChimeraSlayer release will have an option to increase the number of
#      iterations performed. Currently the hard coded default is 100, which leads to
#      considerable fluctuations in the classification. As of Brian Haas recommendation
#      we should increase this value to (optimally to 1000, need to test if a lower value
#      is a good compromise)


class ChimeraChecker(FunctionWithParams):

    """A ChimeraChecker takes a sequence collection and returns list of chimeric seqs.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the alignment app is
    destructive to them.
    """
    Name = 'ChimeraChecker'

    def __init__(self, params):
        """Return new ChimeraChecker object with specified params.

        """
        self.Params = params

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns alignment from sequences.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
         dump the result to the desired path as fasta, otherwise will
         return list of sequence identifiers.
        log_path: path to log, which should include dump of params.

        """
        # return list of the chimeric sequences
        return self.getResult(seq_path)

    def getResult(self, seq_path):
        """ Must be overwritten in subclasses """
        raise NotImplementedError("ChimeraChecker is an abstract class")

    def cleanUp(self):
        """ Overwrite if any clean-up is necessasry

            Included here so all ChimeraCheckers can call cleanUp() without
             having to know if it is defined for a certain sub-class.
        """
        pass


class MetaChimeraChecker(ChimeraChecker):

    """ """
    Name = 'MetaChimeraChecker'

    def __init__(self, params):
        """
        """

        try:
            chimera_checker_constructors = \
                params['chimera_checker_constructors']
        except KeyError:
            raise ValueError("Must provide 'chimera_checker_constructors' to define %s."
                             % self.Name)

        chimera_checkers = []
        for constructor, constructor_params in chimera_checker_constructors:
            chimera_checkers.append(constructor(constructor_params))
        self.chimera_checkers = chimera_checkers

    def getResult(self, seqs_fp):
        raise NotImplementedError("%s is not yet written" % self.Name)


class ChimeraSlayerChimeraChecker(ChimeraChecker):

    """ """
    # Most of the work is done in the ChimeraSlayer_app controller,
    # so we have a very lean Class here
    Name = 'ChimeraSlayerChimeraChecker'

    def __init__(self, params):
        """Return new ChimeraSlayerChimeraChecker object with specified params."""

        _params = {}

        _params.update(params)

        ChimeraChecker.__init__(self, _params)

    def cleanUp(self):
        """Remove all intermediate files"""
        # All intermediates should be removed by app controller
        pass

    def __call__(self, seq_path, db_FASTA_fp=None, db_NAST_fp=None,
                 result_path=None, log_path=None, min_div_ratio=None,
                 keep_intermediates=False):
        """Run chimeraSlayer on input seqs."""

        chimeras = get_chimeras_from_Nast_aligned(
            seq_path, ref_db_aligned_fp=db_NAST_fp,
            ref_db_fasta_fp=db_FASTA_fp,
            min_div_ratio=min_div_ratio,
            keep_intermediates=keep_intermediates)
        return chimeras


def chimeraSlayer_identify_chimeras(seqs_fp, db_FASTA_fp=None,
                                    db_NAST_fp=None, output_fp=None,
                                    min_div_ratio=None, keep_intermediates=False):
    """ """
    params = {}
    cc = ChimeraSlayerChimeraChecker(params)

    if output_fp:
        of = open(output_fp, 'w')
        for seq_id, parents in cc(seqs_fp, db_FASTA_fp=db_FASTA_fp,
                                  db_NAST_fp=db_NAST_fp,
                                  min_div_ratio=min_div_ratio,
                                  keep_intermediates=keep_intermediates):
            of.write('\t'.join([seq_id] + parents))
            of.write('\n')
        of.close()
        result = None
    else:
        result = list(cc(seqs_fp, db_FASTA_fp=db_FASTA_fp,
                         db_NAST_fp=db_NAST_fp,
                         min_div_ratio=min_div_ratio,
                         keep_intermediates=keep_intermediates))

    if not keep_intermediates:
        cc.cleanUp()
    return result


class BlastFragmentsChimeraChecker(ChimeraChecker):

    """  """
    Name = 'BlastFragmentsChimeraChecker'

    def __init__(self, params):
        """Return new BlastFragmentsChimeraChecker object with specified params.

        """
        _params = {'max_e_value': 1e-30,
                   'min_pct_id': 0.90,
                   'num_fragments': 3,
                   'taxonomy_depth': 4}
        _params.update(params)

        try:
            id_to_taxonomy_fp = params['id_to_taxonomy_fp']
        except KeyError:
            raise ValueError(
                "id_to_taxonomy_filepath must be provided to %s" %
                self.Name)

        # Create the blast database if it hasn't been provided
        if 'blast_db' not in params or params['blast_db'] is None:
            try:
                reference_seqs_fp = params['reference_seqs_fp']
            except KeyError:
                raise ValueError(
                    "refseqs_fp or blast_db must be provided to  %s" %
                    self.Name)
            blast_db, self._db_files_to_remove = \
                build_blast_db_from_fasta_path(reference_seqs_fp)
        else:
            blast_db = params['blast_db']
            self._db_files_to_remove = []

        self._taxon_assigner = BlastTaxonAssigner(
            {'blast_db': blast_db,
             'id_to_taxonomy_filepath': id_to_taxonomy_fp,
             'Max E value': _params['max_e_value'],
             'Min percent identity': _params['min_pct_id']
             })

        ChimeraChecker.__init__(self, _params)

    def cleanUp(self):
        """ Remove temporary blast database files, if applicable
        """
        remove_files(self._db_files_to_remove, error_on_missing=False)

    def _fragment_seq(self, seq):
        """ Returns list of n roughly equal-sized fragments of seq

            Modified from Python Cookbook entry comment at:
             http://code.activestate.com/recipes/425397/

            ** This method needs to be parameterized, as there are other ways
            we might want to fragment sequences. For example, we might want to
            specify breakpoints. Currently, this will split seq into
            self.num_fragments roughly equal-sized fragments.
        """
        num_fragments = self.Params['num_fragments']
        results = []
        start = 0
        for i in range(num_fragments):
            # My notes:
            # len(seq[i::n]) gives the number of even multiples of
            # num_fragments exist between i (inclusive) and the end of the seq.
            stop = start + len(seq[i::num_fragments])
            results.append(seq[start:stop])
            start = stop
        return results

    def getResult(self, seq_path):
        """ """
        # Iterate over seq_id, seq pairs from seq_path
        for seq_id, seq in parse_fasta(open(seq_path)):
            # Generate a list of fragments
            fragments = self._fragment_seq(seq)
            # Assign the taxonomy for each of the fragments
            taxon_assignments = \
                [self._get_taxonomy(fragment) for fragment in fragments]
            # Test whether the taxon_assignments suggest that the
            # current seq is a chimera
            if self._is_chimeric(taxon_assignments):
                # If so, yield the seq_id and the list of taxon_assignments
                yield seq_id, taxon_assignments

    def _is_chimeric(self, taxon_assignments):
        """From list of taxon assignments, determine if sequence is chimeric

            ** This method needs to be parameterized, as there are other ways
            we might want to test for chimeras based on taxon_assignments.

            This method checks that taxon assignments are identical to a
            depth of depth (the input parameter). Assignments to
            'No blast hit' are ignored, as these cannot provide evidence
            for a sequence being chimeric. This function considers a sequnce
            non-chimeric until it finds evidence that it is chimeric.
            Currently only the first blast hit is considered -- this is
            something we will want to change.

            When tested against a balanced test set containing 84 chimeric
             sequences and 84 non-chimeric sequences
             (precision, recall, f-measure)
             respectively are:
              @ depth=3 (0.67,0.36,0.47)
              @ depth=4 (0.67,0.60,0.63)
              @ depth=5 (0.56,0.64,0.60)

        """
        depth = self.Params['taxonomy_depth']
        # split the assignments, ignoring fragments which had no blast
        # result (as these cannot provide evidence for a sequence being
        # chimeric)
        split_assignments = [ta.split(';')
                             for ta in taxon_assignments if ta != 'No blast hit']
        for i in range(depth):
            try:
                assignments_at_depth_i =\
                    [a[i].strip() for a in split_assignments]
            except IndexError:
                # We return False if there is no information at current
                # depth because there is not enough information to call
                # the sequence chimeric
                return False
            if len(set(assignments_at_depth_i)) > 1:
                return True

        return False

    def _get_taxonomy(self, fragment):
        """ Return the taxonomy of fragment

            This function is awful right now, because the blast app
             controller takes a file. This therefore results in each
             fragment being written to file and subsequently cleaned up.
             This can be _very_ slow. I recently updated BlastTaxonAssigner
             to handle this better, but it will involve some re-writes here
             to get it all going.
        """

        # Pass the temporary file to the taxon assigner, and get the
        # taxonomy and quality score back
        r = self._taxon_assigner(seqs=[('fragment', fragment)])['fragment']
        taxonomy, quality_score = r[0], r[1]

        # Return the taxonomy
        return taxonomy


def blast_fragments_identify_chimeras(seqs_fp, id_to_taxonomy_fp,
                                      reference_seqs_fp=None, blast_db=None, min_pct_id=0.90,
                                      max_e_value=1e-30, num_fragments=3, output_fp=None, taxonomy_depth=4):
    """ """
    params = {'id_to_taxonomy_fp': id_to_taxonomy_fp,
              'reference_seqs_fp': reference_seqs_fp,
              'blast_db': blast_db,
              'min_pct_id': min_pct_id,
              'max_e_value': max_e_value,
              'num_fragments': num_fragments,
              'taxonomy_depth': taxonomy_depth}
    bcc = BlastFragmentsChimeraChecker(params)

    if output_fp:
        of = open(output_fp, 'w')
        for seq_id, tax_assignments in bcc(seqs_fp):
            of.write('\t'.join([seq_id] + tax_assignments))
            of.write('\n')
        of.close()
        result = None
    else:
        result = list(bcc(seqs_fp))

    bcc.cleanUp()
    return result


# App controller for ChimeraSlayer and convenience functions
class ChimeraSlayer(CommandLineApplication):

    """ ChimeraSlayer ApplicationController

    """

    _command = 'ChimeraSlayer.pl'
    _input_handler = '_input_as_parameters'
    _parameters = {

        # multi-fasta file containing query sequences in alignment format
        '--query_NAST': ValuedParameter('--', Name='query_NAST', Delimiter=' ',
                                        IsPath=True),

        '--db_NAST': ValuedParameter('--', Name='db_NAST', Delimiter=' ',
                                     IsPath=True),

        '--db_FASTA': ValuedParameter('--', Name='db_FASTA', Delimiter=' ',
                                      IsPath=True),

        '--exec_dir': ValuedParameter('--', Name='exec_dir', Delimiter=' ',
                                      IsPath=True),

        '-R': ValuedParameter('-', Name='R', Delimiter=' ')
    }

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input paths (a NAST aligned fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['--query_NAST', '--db_NAST', '--db_FASTA', '-R']

        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError(
                "Unsupported parameter(s) passed when calling ChimeraSlayer: %s" %
                ' '.join(unsupported_parameters))

        return ''

    def _get_result_paths(self, data):
        """ Set the result paths """

        result = {}

        inp_file_name = str(self.Parameters['--query_NAST'].Value)
        inp_file_name = inp_file_name.rstrip('"')
        inp_file_name = inp_file_name.lstrip('"')

        exec_dir = self.Parameters['--exec_dir']
        if exec_dir.isOn():
            exec_dir = str(exec_dir.Value)
            exec_dir = exec_dir.lstrip('"')
            exec_dir = exec_dir.rstrip('"')

            if inp_file_name[0] == '/':
                # path is already absolute
                pass
            else:
                inp_file_name = exec_dir + "/" + inp_file_name

        if not exists(inp_file_name + ".CPS.CPC"):
            raise ApplicationError("Calling ChimeraSlayer failed.")

        result['CPS'] = ResultPath(Path=inp_file_name + ".CPS.CPC",
                                   IsWritten=True)
        return result

    def remove_intermediate_files(self):
        """Remove all intermediate files."""

        # tmp files are written in the current dir,
        # app controller always jumps into dir specified via exec_dir
        # Note: blast intermediates are not removed
        exec_dir = str(self.Parameters['--exec_dir'].Value)
        inp_file_name = str(self.Parameters['--query_NAST'].Value)

        exec_dir = exec_dir.rstrip('"')
        exec_dir = exec_dir.lstrip('"')

        inp_file_name = inp_file_name.rstrip('"')
        inp_file_name = inp_file_name.lstrip('"')

        tmp_suffixes = [".CPS", ".CPS.CPC", ".CPS_RENAST", ".CPS_RENAST.cidx",
                        ".CPS.CPC.wTaxons", ".cidx"]
        cs_tmp_files = [
            exec_dir +
            '/' +
            inp_file_name +
            x for x in tmp_suffixes]
        remove_files(cs_tmp_files, error_on_missing=False)

        db_param = self.Parameters['--db_NAST']
        if db_param.isOn():
            nast_db_name = str(db_param.Value)
            nast_db_name = nast_db_name.rstrip('"')
            nast_db_name = nast_db_name.lstrip('"')

            # Better do not remove this file since other ChimeraSlayer
            # instances running on the same ref set might use this file
            # Should be rather deleted in the calling function
#            remove_files([nast_db_name + ".cidx"],
#                         error_on_missing=False)

        fasta_param = self.Parameters['--db_FASTA']
        if fasta_param.isOn():
            fasta_name = str(fasta_param.Value)
            fasta_name = fasta_name.rstrip('"')
            fasta_name = fasta_name.lstrip('"')

            blast_db_files = [
                fasta_name +
                x for x in [
                    ".nsq",
                    ".nin",
                    ".nhr",
                    ".cidx"]]
            remove_files(blast_db_files, error_on_missing=False)

    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
            """##########################################################################################
#
#  Required:
#
#    --query_NAST      multi-fasta file containing query sequences in alignment format
#
#  Common opts:
#
#    --db_NAST        db in NAST format
#    --db_FASTA       db in fasta format (megablast formatted)
#
#
#    -n       number of top matching database sequences to compare to (default 15)
#    -R       min divergence ratio default: 1.007
#    -P       min percent identity among matching sequences (default: 90)
#
#  ## parameters to tune ChimeraParentSelector:
#
#  Scoring parameters:
#   -M match score   (default: +5)
#   -N mismatch penalty  (default: -4)
#   -Q min query coverage by matching database sequence (default: 70)
#   -T maximum traverses of the multiple alignment  (default: 1)

#
#  ## parameters to tune ChimeraPhyloChecker:
#
#
#    --windowSize                default 50
#    --windowStep                default 5
#    --minBS      minimum bootstrap support for calling chimera (default: 90)
#    -S           percent of SNPs to sample on each side of breakpoint for computing bootstrap support (default: 10)
#    --num_parents_test       number of potential parents to test for chimeras (default: 3)
#    --MAX_CHIMERA_PARENT_PER_ID    Chimera/Parent alignments with perID above this are considered non-chimeras (default 100; turned off)
#
#  ## misc opts
#
#   --printFinalAlignments          shows alignment between query sequence and pair of candidate chimera parents
#   --printCSalignments             print ChimeraSlayer alignments in ChimeraSlayer output
#   --exec_dir                      chdir to here before running
#
#########################################################################################
        """
        return help_str


# Start functions for processing output files
def parse_CPS_file(lines):
    """Parse the CPS file from ChimeraSlayer.

One line if this file look like this:

ChimeraSlayer   chimera_AJ007403        7000004131495956        S000469847      1.0360  99.35   100     0.9354  89.70   0       YES    NAST:4595-4596  ECO:941-942

0      ChimeraSlayer
1      chimera_AJ007403            # the accession of the chimera query
2      S000387216                  # reference parent A
3      S000001688                  # reference parent B
4      0.9422                      # divergence ratio of query to chimera (left_A, right_B)
5      90.00                       # percent identity between query and chimera(left_A, right_B)
6      0                           # confidence in query as a chimera related to (left_A, right_B)
7      1.0419                      # divergence ratio of query to chimera (right_A, left_B)
8      99.52                       # percent identity between query and chimera(right_A, left_B)
9      100                         # confidence in query as a chimera related to (right_A, left_B)
10     YES                         # ** verdict as a chimera or not **
11     NAST:4032-4033              # estimated approximate chimera breakpoint in NAST coordinates
12     ECO:767-768                 # estimated approximate chimera breakpoint according to the E. coli unaligned reference seq coordinates
"""

    chimeras = []
    for line in lines:
        record = line.split()
        try:
            id = record[1]
            parent1 = record[2]
            parent2 = record[3]
            verdict = record[10]
        except IndexError:
            raise ValueError("Error parsing ChimeraSlayer CPS file.")
        if verdict == "YES":
            chimeras.append((id, [parent1, parent2]))
    return chimeras

# End functions for processing output files


def make_cidx_file(fp):
    """make a CS index file.

    fp: reference DB to make an index of

    ChimeraSlayer implicetely performs this task when the file is not
    already created and deletes it when done. Especially in a parallel
    environment it mkaes sense to manually make and delete the file.
    """

    if which("cdbfasta"):
        # cdbfasta comes with the microbiometools/CS
        # should always be installed when using CS
        args = ["cdbfasta", fp]
        # cdbfasta write one line to stderr
        stdout, stderr = Popen(args, stderr=PIPE, stdout=PIPE).communicate()
    else:
        raise ApplicationNotFoundError

# Start convenience functions


def get_chimeras_from_Nast_aligned(seqs_fp, ref_db_aligned_fp=None,
                                   ref_db_fasta_fp=None,
                                   HALT_EXEC=False, min_div_ratio=None,
                                   keep_intermediates=False):
    """remove chimeras from seqs_fp using chimeraSlayer.

    seqs_fp:  a filepath with the seqs to check in the file
    ref_db_aligned_fp: fp to (pynast) aligned reference sequences
    ref_db_fasta_fp: same seqs as above, just unaligned. Will be computed on the fly if not provided,
    HALT_EXEC: stop execution if true
    min_div_ratio: passed to ChimeraSlayer App
    """

    files_to_remove = []
    # might come in as FilePath object with quotes
    seqs_fp = str(seqs_fp)
    seqs_fp = seqs_fp.rstrip('"')
    seqs_fp = seqs_fp.lstrip('"')

    seqs_dir, new_seqs_fp = split(seqs_fp)

    # if fp is in current dir, we fake a dir change
    if seqs_dir == "":
        seqs_dir = "./"

    # Chimera Slayer puts some temp files in current dir and some in dir of input file
    # use exe_dir to change to dir of input file, so to have all tmp files in
    # one place
    params = {'--query_NAST': new_seqs_fp,
              '--exec_dir': seqs_dir}

    if ref_db_aligned_fp is None and ref_db_fasta_fp is None:
        # use default db, whose relative position to the
        # ChimeraSlayer binary is hardcoded
        pass

    else:
        if not ref_db_fasta_fp:
            # make degapped reference file
            ref_db_fasta_fp = write_degapped_fasta_to_file(parse_fasta(
                open(ref_db_aligned_fp)))
            files_to_remove.append(ref_db_fasta_fp)
        # use user db
        params.update({'--db_NAST': abspath(ref_db_aligned_fp),
                       '--db_FASTA': abspath(ref_db_fasta_fp)})

    if min_div_ratio is not None:
        params.update({'-R': min_div_ratio})

    app = ChimeraSlayer(params=params, HALT_EXEC=HALT_EXEC)
    app_results = app()

#    this is a FilePath object in case of success.
#    How can we test for failure here?
    #    if not exists(app_results['CPS']):
#         raise ApplicationError, "ChimeraSlayer failed. No output file."

    chimeras = parse_CPS_file((app_results['CPS']))
    if not keep_intermediates:
        app.remove_intermediate_files()
        remove_files(files_to_remove)

    return chimeras


def usearch61_chimera_check(input_seqs_fp,
                            output_dir,
                            reference_seqs_fp=None,
                            suppress_usearch61_intermediates=False,
                            suppress_usearch61_ref=False,
                            suppress_usearch61_denovo=False,
                            split_by_sampleid=False,
                            non_chimeras_retention="union",
                            usearch61_minh=0.28,
                            usearch61_xn=8.0,
                            usearch61_dn=1.4,
                            usearch61_mindiffs=3,
                            usearch61_mindiv=0.8,
                            usearch61_abundance_skew=2.0,
                            percent_id_usearch61=0.97,
                            minlen=64,
                            word_length=8,
                            max_accepts=1,
                            max_rejects=8,
                            verbose=False,
                            threads=1.0,
                            HALT_EXEC=False):
    """ Main convenience function for usearch61 chimera checking

    input_seqs_fp:  filepath of input fasta file.
    output_dir:  output directory
    reference_seqs_fp: fasta filepath for reference chimera detection.
    suppress_usearch61_intermediates:  Suppress retention of .uc and log files.
    suppress_usearch61_ref:  Suppress usearch61 reference chimera detection.
    suppress_usearch61_denovo:  Suppress usearch61 de novo chimera detection.
    split_by_sampleid:  Split by sample ID for de novo chimera detection.
    non_chimeras_retention: Set to "union" or "intersection" to retain
     non-chimeras between de novo and reference based results.
    usearch61_minh: Minimum score (h) to be classified as chimera.
     Increasing this value tends to the number of false positives (and also
     sensitivity).
    usearch61_xn:  Weight of "no" vote.  Increasing this value tends to the
     number of false positives (and also sensitivity).
    usearch61_dn:  Pseudo-count prior for "no" votes. (n). Increasing this
     value tends to the number of false positives (and also sensitivity).
    usearch61_mindiffs:  Minimum number of diffs in a segment. Increasing this
     value tends to reduce the number of false positives while reducing
     sensitivity to very low-divergence chimeras.
    usearch61_mindiv:  Minimum divergence, i.e. 100% - identity between the
     query and closest reference database sequence. Expressed as a percentage,
     so the default is 0.8%, which allows chimeras that are up to 99.2% similar
     to a reference sequence.
    usearch61_abundance_skew: abundance skew for de novo chimera comparisons.
    percent_id_usearch61: identity to cluster sequences at
    minlen: minimum sequence length for use with usearch61
    word_length: length of nucleotide 'words' for usearch61
    max_accepts: max number of accepts for hits with usearch61
    max_rejects: max number of rejects for usearch61, increasing allows more
     sensitivity at a cost of speed
    threads: Specify number of threads used per core per CPU
    HALT_EXEC=application controller option to halt execution and print command
    """

    """
    Need to cluster sequences de novo first to get 1. abundance information
    and 2 consensus sequence for each cluster.  Using dereplication followed
    by clustering does not appear to automatically update complete cluster
    size, will directly cluster raw seqs with the small_mem clustering option.

    This means without additional parsing steps to recalculate
    actual cluster sizes, the sizeorder option can't be used for de novo
    clustering and downstream chimera detection."""

    files_to_remove = []

    # Get absolute paths to avoid issues with calling usearch
    input_seqs_fp = abspath(input_seqs_fp)
    output_dir = abspath(output_dir)
    if reference_seqs_fp:
        reference_seqs_fp = abspath(reference_seqs_fp)
    log_fp = join(output_dir, "identify_chimeric_seqs.log")
    chimeras_fp = join(output_dir, "chimeras.txt")
    non_chimeras_fp = join(output_dir, "non_chimeras.txt")

    non_chimeras = []
    chimeras = []
    log_lines = {'denovo_chimeras': 0,
                 'denovo_non_chimeras': 0,
                 'ref_chimeras': 0,
                 'ref_non_chimeras': 0}

    if split_by_sampleid:
        if verbose:
            print "Splitting fasta according to SampleID..."

        with open(input_seqs_fp, 'U') as full_seqs:
            sep_fastas =split_sequence_file_on_sample_ids_to_files(
                full_seqs,
                'fasta',
                output_dir)

        if suppress_usearch61_intermediates:
            files_to_remove += sep_fastas

        for curr_fasta in sep_fastas:
            curr_chimeras, curr_non_chimeras, files_to_remove, log_lines =\
                identify_chimeras_usearch61(curr_fasta, output_dir,
                                            reference_seqs_fp, suppress_usearch61_intermediates,
                                            suppress_usearch61_ref, suppress_usearch61_denovo,
                                            non_chimeras_retention, usearch61_minh, usearch61_xn,
                                            usearch61_dn, usearch61_mindiffs, usearch61_mindiv,
                                            usearch61_abundance_skew, percent_id_usearch61, minlen,
                                            word_length, max_accepts, max_rejects, files_to_remove, HALT_EXEC,
                                            log_lines, verbose, threads)

            chimeras += curr_chimeras
            non_chimeras += curr_non_chimeras

    else:
        chimeras, non_chimeras, files_to_remove, log_lines =\
            identify_chimeras_usearch61(input_seqs_fp, output_dir,
                                        reference_seqs_fp, suppress_usearch61_intermediates,
                                        suppress_usearch61_ref, suppress_usearch61_denovo,
                                        non_chimeras_retention, usearch61_minh, usearch61_xn,
                                        usearch61_dn, usearch61_mindiffs, usearch61_mindiv,
                                        usearch61_abundance_skew, percent_id_usearch61, minlen,
                                        word_length, max_accepts, max_rejects, files_to_remove, HALT_EXEC,
                                        log_lines, verbose, threads)

    # write log, non chimeras, chimeras.
    write_usearch61_log(log_fp, input_seqs_fp, output_dir,
                        reference_seqs_fp, suppress_usearch61_intermediates,
                        suppress_usearch61_ref, suppress_usearch61_denovo,
                        split_by_sampleid, non_chimeras_retention, usearch61_minh,
                        usearch61_xn, usearch61_dn, usearch61_mindiffs, usearch61_mindiv,
                        usearch61_abundance_skew, percent_id_usearch61, minlen,
                        word_length, max_accepts, max_rejects, HALT_EXEC, log_lines)

    chimeras_f = open(chimeras_fp, "w")
    non_chimeras_f = open(non_chimeras_fp, "w")
    for curr_chimera in chimeras:
        chimeras_f.write("%s\n" % curr_chimera)
    for curr_non_chimera in non_chimeras:
        non_chimeras_f.write("%s\n" % curr_non_chimera)
    chimeras_f.close()
    non_chimeras_f.close()

    remove_files(files_to_remove)


def identify_chimeras_usearch61(input_seqs_fp,
                                output_dir,
                                reference_seqs_fp=None,
                                suppress_usearch61_intermediates=False,
                                suppress_usearch61_ref=False,
                                suppress_usearch61_denovo=False,
                                non_chimeras_retention="union",
                                usearch61_minh=0.28,
                                usearch61_xn=8.0,
                                usearch61_dn=1.4,
                                usearch61_mindiffs=3,
                                usearch61_mindiv=0.8,
                                usearch61_abundance_skew=2.0,
                                percent_id_usearch61=0.97,
                                minlen=64,
                                word_length=8,
                                max_accepts=1,
                                max_rejects=8,
                                files_to_remove=[],
                                HALT_EXEC=False,
                                log_lines={},
                                verbose=False,
                                threads=1.0):
    """ Calls de novo, ref useach61 chimera checking, returns flagged seqs

    input_seqs_fp:  filepath of input fasta file.
    output_dir:  output directory
    reference_seqs_fp: fasta filepath for reference chimera detection.
    suppress_usearch61_intermediates:  Suppress retention of .uc and log files.
    suppress_usearch61_ref:  Suppress usearch61 reference chimera detection.
    suppress_usearch61_denovo:  Suppress usearch61 de novo chimera detection.
    non_chimeras_retention: Set to "union" or "intersection" to retain
     non-chimeras between de novo and reference based results.
    usearch61_minh: Minimum score (h) to be classified as chimera.
     Increasing this value tends to the number of false positives (and also
     sensitivity).
    usearch61_xn:  Weight of "no" vote.  Increasing this value tends to the
     number of false positives (and also sensitivity).
    usearch61_dn:  Pseudo-count prior for "no" votes. (n). Increasing this
     value tends to the number of false positives (and also sensitivity).
    usearch61_mindiffs:  Minimum number of diffs in a segment. Increasing this
     value tends to reduce the number of false positives while reducing
     sensitivity to very low-divergence chimeras.
    usearch61_mindiv:  Minimum divergence, i.e. 100% - identity between the
     query and closest reference database sequence. Expressed as a percentage,
     so the default is 0.8%, which allows chimeras that are up to 99.2% similar
     to a reference sequence.
    usearch61_abundance_skew: abundance skew for de novo chimera comparisons.
    percent_id_usearch61: identity to cluster sequences at
    minlen: minimum sequence length for use with usearch61
    word_length: length of nucleotide 'words' for usearch61
    max_accepts: max number of accepts for hits with usearch61
    max_rejects: max number of rejects for usearch61, increasing allows more
     sensitivity at a cost of speed
    files_to_remove: list of filepaths to remove if intermediate files are to
     be removed.
    HALT_EXEC=application controller option to halt execution and print command
    threads: Specify number of threads used per core per CPU
    """

    output_consensus_fp = join(output_dir,
                               "%s_consensus_with_abundance.fasta" % basename(input_seqs_fp))
    filtered_consensus_fp = join(output_dir,
                                 "%s_consensus_fixed.fasta" % basename(input_seqs_fp))
    output_consensus_uc = join(output_dir,
                               "%s_consensus_with_abundance.uc" % basename(input_seqs_fp))
    smallmem_log_fp = join(output_dir,
                           "%s_smallmem_clustered.log" % basename(input_seqs_fp))
    uchime_denovo_fp = join(output_dir,
                            "%s_chimeras_denovo.uchime" % basename(input_seqs_fp))
    uchime_denovo_log_fp = join(output_dir,
                                "%s_chimeras_denovo.log" % basename(input_seqs_fp))
    uchime_ref_fp = join(output_dir,
                         "%s_chimeras_ref.uchime" % basename(input_seqs_fp))
    uchime_ref_log_fp = join(output_dir,
                             "%s_chimeras_ref.log" % basename(input_seqs_fp))

    if verbose:
        print "Clustering sequence %s..." % basename(input_seqs_fp)

    abundance_fp, app_result = usearch61_smallmem_cluster(input_seqs_fp,
                                                          percent_id=percent_id_usearch61,
                                                          minlen=minlen,
                                                          rev=False,
                                                          output_dir=output_dir,
                                                          remove_usearch_logs=
                                                          suppress_usearch61_intermediates,
                                                          wordlength=word_length,
                                                          usearch61_maxrejects=max_rejects,
                                                          usearch61_maxaccepts=max_accepts,
                                                          HALT_EXEC=HALT_EXEC,
                                                          output_uc_filepath=output_consensus_uc,
                                                          log_name=smallmem_log_fp,
                                                          sizeout=True,
                                                          consout_filepath=output_consensus_fp)

    de_novo_clusters, failures =\
        parse_usearch61_clusters(open(output_consensus_uc, "U"))

    # RE's parsing of size=X part of label is position dependent, and only
    # creates this in correct format with derep and fast cluster options,
    # so must do manual fix of fasta for later steps.
    fix_abundance_labels(output_consensus_fp, filtered_consensus_fp)

    if suppress_usearch61_intermediates:
        files_to_remove.append(output_consensus_fp)
        files_to_remove.append(output_consensus_uc)
        files_to_remove.append(filtered_consensus_fp)

    chimeras = []
    non_chimeras = []

    if not suppress_usearch61_denovo:
        if verbose:
            print "De novo chimera checking %s..." % basename(input_seqs_fp)
        denovo_uchime_fp, app_result =\
            usearch61_chimera_check_denovo(filtered_consensus_fp,
                                           uchime_denovo_fp,
                                           minlen=minlen,
                                           output_dir=output_dir,
                                           remove_usearch_logs=
                                           suppress_usearch61_intermediates,
                                           uchime_denovo_log_fp=
                                           uchime_denovo_log_fp,
                                           usearch61_minh=usearch61_minh,
                                           usearch61_xn=usearch61_xn,
                                           usearch61_dn=usearch61_dn,
                                           usearch61_mindiffs=usearch61_mindiffs,
                                           usearch61_mindiv=usearch61_mindiv,
                                           usearch61_abundance_skew=
                                           usearch61_abundance_skew,
                                           HALT_EXEC=HALT_EXEC)

        if suppress_usearch61_intermediates:
            files_to_remove.append(uchime_denovo_fp)

        de_novo_chimeras = get_chimera_clusters(denovo_uchime_fp)

        expanded_denovo_chimeras, expanded_non_denovo_chimeras =\
            merge_clusters_chimeras(de_novo_clusters, set(de_novo_chimeras))

        log_lines['denovo_chimeras'] += len(expanded_denovo_chimeras)
        log_lines['denovo_non_chimeras'] += len(expanded_non_denovo_chimeras)

    if not suppress_usearch61_ref:
        if verbose:
            print "Reference chimera checking %s..." % basename(input_seqs_fp)

        ref_uchime_fp, app_result =\
            usearch61_chimera_check_ref(filtered_consensus_fp,
                                        uchime_ref_fp,
                                        reference_seqs_fp,
                                        minlen=minlen,
                                        output_dir=output_dir,
                                        remove_usearch_logs=
                                        suppress_usearch61_intermediates,
                                        uchime_ref_log_fp=
                                        uchime_ref_log_fp,
                                        usearch61_minh=usearch61_minh,
                                        usearch61_xn=usearch61_xn,
                                        usearch61_dn=usearch61_dn,
                                        usearch61_mindiffs=usearch61_mindiffs,
                                        usearch61_mindiv=usearch61_mindiv,
                                        threads=threads,
                                        HALT_EXEC=HALT_EXEC)

        if suppress_usearch61_intermediates:
            files_to_remove.append(uchime_ref_fp)

        ref_chimeras = get_chimera_clusters(uchime_ref_fp)

        expanded_ref_chimeras, expanded_non_ref_chimeras =\
            merge_clusters_chimeras(de_novo_clusters, set(ref_chimeras))

        log_lines['ref_chimeras'] += len(expanded_ref_chimeras)
        log_lines['ref_non_chimeras'] += len(expanded_non_ref_chimeras)

    if suppress_usearch61_ref:
        chimeras = expanded_denovo_chimeras
        non_chimeras = expanded_non_denovo_chimeras
    elif suppress_usearch61_denovo:
        chimeras = expanded_ref_chimeras
        non_chimeras = expanded_non_ref_chimeras
    else:
        expanded_non_denovo_chimeras = set(expanded_non_denovo_chimeras)
        expanded_non_ref_chimeras = set(expanded_non_ref_chimeras)
        expanded_denovo_chimeras = set(expanded_denovo_chimeras)
        expanded_ref_chimeras = set(expanded_ref_chimeras)
        if non_chimeras_retention == "union":
            non_chimeras = expanded_non_denovo_chimeras.union(
                expanded_non_ref_chimeras)
            chimeras = expanded_denovo_chimeras.intersection(
                expanded_ref_chimeras)
        else:
            non_chimeras = expanded_non_denovo_chimeras.intersection(
                expanded_non_ref_chimeras)
            chimeras = expanded_denovo_chimeras.union(
                expanded_ref_chimeras)

    return chimeras, non_chimeras, files_to_remove, log_lines


def write_usearch61_log(log_fp,
                        input_seqs_fp,
                        output_dir,
                        reference_seqs_fp,
                        suppress_usearch61_intermediates,
                        suppress_usearch61_ref,
                        suppress_usearch61_denovo,
                        split_by_sampleid,
                        non_chimeras_retention,
                        usearch61_minh,
                        usearch61_xn,
                        usearch61_dn,
                        usearch61_mindiffs,
                        usearch61_mindiv,
                        usearch61_abundance_skew,
                        percent_id_usearch61,
                        minlen,
                        word_length,
                        max_accepts,
                        max_rejects,
                        HALT_EXEC,
                        log_lines):
    """ Writes overall log data for identify_chimeric_seqs.py call

    log_fp: filepath to write log data.
    """

    out_f = open(log_fp, "w")

    param_names = ["input_seqs_fp", "output_dir",
                   "reference_seqs_fp", "suppress_usearch61_intermediates",
                   "suppress_usearch61_ref", "suppress_usearch61_denovo",
                   "split_by_sampleid", "non_chimeras_retention", "usearch61_minh",
                   "usearch61_xn", "usearch61_dn", "usearch61_mindiffs", "usearch61_mindiv",
                   "usearch61_abundance_skew", "percent_id_usearch61", "minlen",
                   "word_length", "max_accepts", "max_rejects", "HALT_EXEC"]

    param_values = [input_seqs_fp, output_dir,
                    reference_seqs_fp, suppress_usearch61_intermediates,
                    suppress_usearch61_ref, suppress_usearch61_denovo,
                    split_by_sampleid, non_chimeras_retention, usearch61_minh,
                    usearch61_xn, usearch61_dn, usearch61_mindiffs, usearch61_mindiv,
                    usearch61_abundance_skew, percent_id_usearch61, minlen,
                    word_length, max_accepts, max_rejects, HALT_EXEC]

    for curr_param in range(len(param_names)):
        out_f.write("%s\t%s\n" % (param_names[curr_param],
                                  param_values[curr_param]))

    out_f.write("\n")

    for curr_line in log_lines.keys():
        out_f.write("%s\t%s\n" % (curr_line, log_lines[curr_line]))

    out_f.close()

    return


def fix_abundance_labels(output_consensus_fp, filtered_consensus_fp):
    """ puts size= part of label as second component after white space

    output_consensus_fp: consensus filepath with abundance data
    filtered_consensus_fp: output filepath name
    """

    consensus_f = open(output_consensus_fp, "U")

    filtered_f = open(filtered_consensus_fp, "w")

    for label, seq in parse_fasta(consensus_f):
        fasta_label = label.split()[0]
        size = "size=" + label.split('size=')[1].replace(';', '')
        final_label = "%s;%s" % (fasta_label, size)
        filtered_f.write(">%s\n%s\n" % (final_label, seq))

    consensus_f.close()
    filtered_f.close()


def get_chimera_clusters(uchime_fp):
    """ Parses .uchime file, returns list of seq IDs flagged as chimeras

    uchime_fp: filepath of uchime file
    """

    uchime_f = open(uchime_fp, "U")

    chimeras = []

    seq_id_ix = 1
    flag_ix = 17

    for line in uchime_f:
        if line.startswith("#"):
            continue
        if len(line.strip()) == 0:
            continue
        curr_line = line.split()
        if curr_line[flag_ix] == "Y":
            chimeras.append(
                curr_line[seq_id_ix].split('centroid=')[1].split(';')[0])

    return chimeras


def merge_clusters_chimeras(de_novo_clusters,
                            chimeras):
    """ Merges results of chimeras/clusters into list of chimeras, nonchimeras

    de_novo_clusters: dict of clusterID: seq IDs
    chimeras: list of chimeric seq IDs
    """

    expanded_chimeras = []
    expanded_non_chimeras = []

    for curr_cluster in de_novo_clusters:
        curr_seq_ids = []
        matches_chimera = False
        for curr_seq in de_novo_clusters[curr_cluster]:
            if curr_seq in chimeras:
                matches_chimera = True
            curr_seq_ids.append(curr_seq)
        if matches_chimera:
            expanded_chimeras += curr_seq_ids
        else:
            expanded_non_chimeras += curr_seq_ids

    return expanded_chimeras, expanded_non_chimeras
