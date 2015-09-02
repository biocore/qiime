from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"


from os.path import join, basename, splitext, exists
from functools import partial


def create_commands_jpe(pairs, base_output_dir, optional_params = "",
        leading_text = "", trailing_text = "", include_input_dir_path=False,
        remove_filepath_in_name=False, match_barcodes = False,
        bc_pairs = {}):
    """ Creates commands for join_paired_ends.py

    pairs: dictionary of forward:reverse read filepaths
    base_output_dir: output directory to write log, stitched reads
    optional_params: added parameters to join_paired_ends.py calls
    leading_text: Text to add before join_paired_ends.py call
    trailing_text: Text to add after join_paired_ends.py call
    include_input_dir_path: If True, include input directory in output
        directory names
    remove_filepath_in_name: If True, the base filename will not be used in the
        output directory names.
    match_barcodes: True to match barcodes.
    bc_pairs: dictionary of read1:bc_read filepaths (empty if not used)
    """

    commands = []
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq']

    for curr_fp in pairs:
        for extension in extensions:
            if extension in curr_fp:
                curr_ext = extension
        if include_input_dir_path:
            added_output_str = curr_fp.split('/')[-2]
        else:
            added_output_str = ""
        if not remove_filepath_in_name:
            added_output_str += basename(curr_fp).split(curr_ext)[0]


        curr_outputdir = join(base_output_dir, added_output_str)
        if match_barcodes:
            command = "%sjoin_paired_ends.py %s -b %s -f %s -r %s -o %s %s" %\
                (_clean_leading_text(leading_text), optional_params, bc_pairs[curr_fp], curr_fp,
                pairs[curr_fp], curr_outputdir, trailing_text)
        else:
            command = "%sjoin_paired_ends.py %s -f %s -r %s -o %s %s" %\
                (_clean_leading_text(leading_text), optional_params, curr_fp, pairs[curr_fp],
                curr_outputdir, trailing_text)

        commands.append([('join_paired_ends.py: %s' % curr_fp, command)])

    return commands

def create_commands_eb(all_files, ispaired, base_output_dir,
        optional_params = "", leading_text = "", trailing_text = "",
        include_input_dir_path=False, remove_filepath_in_name=False):
    """ Creates commands for extract_barcodes.py

    all_files: list of input filelpaths or dict of paired files
    ispaired: True if paired data
    base_output_dir: output directory to write log, stitched reads
    optional_params: added parameters to extract_barcodes.py calls
    leading_text: Text to add before extract_barcodes.py call
    trailing_text: Text to add after extract_barcodes.py call
    include_input_dir_path: If True, include input directory in output
        directory names
    remove_filepath_in_name: If True, the base filename will not be used in the
        output directory names.
    """

    commands = []
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq']

    for curr_fp in all_files:
        if include_input_dir_path:
            added_output_str = curr_fp.split('/')[-2]
        else:
            added_output_str = ""
        if not remove_filepath_in_name:
            for extension in extensions:
                if extension in curr_fp:
                    curr_ext = extension
            added_output_str += basename(curr_fp).split(curr_ext)[0]

        curr_outputdir = join(base_output_dir, added_output_str)
        if ispaired:
            command = "%sextract_barcodes.py %s -f %s -r %s -o %s %s" %\
            (_clean_leading_text(leading_text), optional_params, curr_fp, all_files[curr_fp],
            curr_outputdir, trailing_text)
        else:
            command = "%sextract_barcodes.py %s -f %s -o %s %s" %\
            (_clean_leading_text(leading_text), optional_params, curr_fp,
            curr_outputdir, trailing_text)

        commands.append([('extract_barcodes.py: %s' % curr_fp, command)])

    return commands

def create_commands_slf(all_files, demultiplexing_method, output_dir,
        params = "", leading_text = "", trailing_text = "",
        include_input_dir_path=False, remove_filepath_in_name=False,
        sampleid_indicator = "_"):
    """ Creates command for split_libraries_fastq.py

    all_files: list of input filelpaths or dict of reads:(barcode,mapping)
    demultiplexing_method: Either 'sampleid_by_file' or 'mapping_barcode_files'
    output_dir: output directory to write split_libraries_fastq output, the
        directory has to exist before calling this function.
    params: added parameters to split_libraries_fastq.py calls
    leading_text: Text to add before split_libraries_fastq.py call
    trailing_text: Text to add after split_libraries_fastq.py call
    include_input_dir_path: If True, include input directory in output
        directory names
    remove_filepath_in_name: If True, the base filename will not be used in the
        output directory names.
    sampleid_indicator: Split on this character in input fastq filenames to
        generate output SampleID name.

    Raises
    ------
    IOError
        If the output_dir doesn't exist.
    """

    # we need the folder to exist so we can write the seqs.txt, barcodes.txt,
    # sample_ids.txt and maps.txt files
    if not exists(output_dir):
        raise IOError("%s directory doesn't exist" % output_dir)

    commands = []
    read_files = []
    barcode_files = []
    mapping_files = []
    sample_ids = []

    params_dict = {}
    path_builder = partial(join, output_dir)

    # Using a set in this case to keep consistent order (needed for unit tests)
    all_fps = set(all_files)

    for curr_fp in all_fps:
        read_files.append(curr_fp)
        # Just need to build up a list of SampleID names
        if demultiplexing_method == 'sampleid_by_file':
            if include_input_dir_path:
                sample_id = curr_fp.split('/')[-2]
            else:
                sample_id = ""
            if not remove_filepath_in_name:
                sample_id += basename(curr_fp).split(sampleid_indicator)[0]
            sample_ids.append(sample_id)
        # Need list of barcode filepaths, mapping filepaths
        else:
            barcode_files.append(all_files[curr_fp][0])
            mapping_files.append(all_files[curr_fp][1])

    # gather all the arguments into a dictionary to format them into a command
    params_dict = {'leading_text': _clean_leading_text(leading_text),
                   'trailing_text': trailing_text, 'output_dir': output_dir,
                   'params': params,
                   'read_files': _make_file(read_files,
                                            path_builder('seqs.txt'))}
    cmd = ("{leading_text}split_libraries_fastq.py {params} "
           "-i {read_files} -o {output_dir} --read_arguments_from_file ")

    if demultiplexing_method == 'sampleid_by_file':
        cmd += ("--sample_ids {sample_ids} --barcode_type 'not-barcoded' ")
        params_dict['sample_ids'] = _make_file(sample_ids,
                                               path_builder('sample_ids.txt'))
    else:
        cmd += ("--barcode_read_fps {barcode_files} "
                "--mapping_fps {mapping_files} ")
        params_dict['barcode_files'] = _make_file(barcode_files,
                                                  path_builder('barcodes.txt'))
        params_dict['mapping_files'] = _make_file(mapping_files,
                                                  path_builder('maps.txt'))

    # regardless of the method, trailing_text should go last
    cmd += '{trailing_text}'
    commands.append([('split_libraries_fastq.py', cmd.format(**params_dict))])

    return commands


def _make_file(elements, fp):
    """Create a file with one element per line

    This function is just a helper for create_commands_slf.

    Parameters
    ----------
    elements : list or tuple
        List of elements to be formatted one by line in the output file.
    fp : str
        Filepath where the file should be written to

    Returns
    -------
    str
        Filepath where the file is written.
    """
    with open(fp, 'w') as f:
        f.write('\n'.join(elements))
    return fp


def get_pairs(all_files, read1_indicator, read2_indicator, match_barcodes=False,
        barcode_indicator="_I1_"):
    """ Finds pairs of files from a list of files, optionally matches barcodes

    all_files: list of filepaths
    read1_indicator: string indicating read 1 of a pair
    read2_indicator: string indicating read 2 of a pair
    match_barcodes: If True, will attempt to match up barcodes file
    barcode_indicator: string indicating barcode file.
    """

    pairs = {}
    bc_pairs = {}

    read1_files = []
    read2_files = []
    bc_files = []

    for curr_file in all_files:
        curr_file_string_r1 = curr_file.split(read1_indicator)
        curr_file_string_r2 = curr_file.split(read2_indicator)
        if match_barcodes:
            curr_file_string_bc = curr_file.split(barcode_indicator)

        if len(curr_file_string_r1) == 2:
            read1_files.append(curr_file_string_r1)
        elif len(curr_file_string_r2) == 2:
            read2_files.append(curr_file_string_r2)
        elif match_barcodes and len(curr_file_string_bc) == 2:
            bc_files.append(curr_file_string_bc)
        else:
            raise ValueError,("Invalid filename found for splitting on input "+\
                "for file %s, " % curr_file + "check input read1_indicator "+\
                "and read2_indicator parameters as well.")

    for curr_read1 in read1_files:
        for curr_read2 in read2_files:
            if curr_read1 == curr_read2:
                pairs[read1_indicator.join(curr_read1)] =\
                    read2_indicator.join(curr_read2)

    if match_barcodes:
        for curr_read1 in read1_files:
            for curr_bc in bc_files:
                if curr_read1 == curr_bc:
                    bc_pairs[read1_indicator.join(curr_read1)] =\
                        barcode_indicator.join(curr_bc)
        # Need a specific test if matched barcodes are used-the barcodes should
        # match both the forward and reverse reads.
        forward_reads = set(pairs.keys())
        bc_reads = set(bc_pairs.keys())
        non_matching_f_reads = forward_reads - bc_reads
        if non_matching_f_reads:
            raise ValueError,("Found forward reads without matching barcodes "
                "file: %s" % non_matching_f_reads)

    return pairs, bc_pairs

def get_matching_files(all_fastq, all_mapping,
        read_indicator, barcode_indicator, mapping_indicator):
    """ Matches up read, barcode, and mapping files based on filenames

    all_fastq: list of sequence filepaths
    all_mapping: list of mapping filepaths
    read_indicator: string indicating read file
    barcode_indicator: string indicating barcode file
    mapping_indicator: string indicating mapping file
    """

    read_files = []
    barcode_files = []
    mapping_files = {}
    matching_files = {}

    # Have to assume trailing text will not match extensions, so have to
    # do some splitting at the extension point to match up.
    for curr_file in all_mapping:
        try:
            curr_mapping = curr_file.split(mapping_indicator)
            mapping_files[curr_mapping[0] +
                splitext(curr_mapping[1])[0]] = curr_file
        except IndexError:
            raise IndexError(
                "Found file with a mapping file extension that does not "
                "contain the mapping file indicators (see mapping_indicator): "
                "%s" % curr_file)

    for curr_file in all_fastq:
        curr_file_string_read = curr_file.split(read_indicator)
        curr_file_string_bc = curr_file.split(barcode_indicator)

        if len(curr_file_string_read) == 2:
            read_files.append(curr_file_string_read)
        elif len(curr_file_string_bc) == 2:
            barcode_files.append(curr_file_string_bc)
        else:
            raise ValueError("Invalid filename found for splitting on input "+\
                "for file %s, " % curr_file + "check input read indicator "+\
                "and barcode indicator parameters.")

    for curr_read in read_files:
        for curr_bc in barcode_files:
            if curr_read == curr_bc:
                curr_read_sans_ext = curr_read[0] + curr_read[1].split('.f')[0]
                try:
                    matching_files[read_indicator.join(curr_read)] =\
                        (barcode_indicator.join(curr_bc),
                        mapping_files[curr_read_sans_ext])
                except KeyError:
                    raise KeyError("Found read file with no matching mapping "
                       "file: %s" % read_indicator.join(curr_read))
    return matching_files


def _clean_leading_text(leading_text):
    leading_text = leading_text.strip()
    if leading_text:
        return leading_text + ' '
    else:
        return leading_text
