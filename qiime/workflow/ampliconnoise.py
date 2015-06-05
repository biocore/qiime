#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Justin Kuczynski",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import os
from skbio.sequence import DNASequence
from qiime.parse import parse_mapping_file
from qiime.util import create_dir
from qiime.workflow.util import (print_to_stdout,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 log_input_md5s)


def run_ampliconnoise(mapping_fp,
                      output_dir, command_handler, params, qiime_config,
                      logger=None, status_update_callback=print_to_stdout,
                      chimera_alpha=-3.8228, chimera_beta=0.6200, sff_txt_fp=None, numnodes=2,
                      suppress_perseus=True, output_filepath=None, platform='flx',
                      seqnoise_resolution=None, truncate_len=None):
    """ Run the ampliconnoise pipeline

        The steps performed by this function are:
1. Split input sff.txt file into one file per sample

2. Run scripts required for PyroNoise

3. Run scripts required for SeqNoise

4. Run scripts requred for Perseus (chimera removal)

5. Merge output files into one file similar to the output of split_libraries.py

    output_filepath should be absolute
    seqnoise_resolution should be string
    environment variable PYRO_LOOKUP_FILE must be set correctly. Thus be
    careful passing command handlers that don't spawn child processes, as they
    may not inherit the correct environment variable setting
    """
    map_data, headers, comments = parse_mapping_file(open(mapping_fp, 'U'))
    create_dir(output_dir)

    if seqnoise_resolution is None:
        if platform == 'flx':
            seqnoise_resolution = '30.0'
        elif platform == 'titanium':
            seqnoise_resolution = '25.0'
        else:
            raise RuntimeError('seqnoise_resolution not set, and no' +
                               ' default for platform ' + platform)

    if truncate_len is None:
        if platform == 'flx':
            truncate_len = '220'
        elif platform == 'titanium':
            truncate_len = '400'
        else:
            raise RuntimeError('truncate_len not set, and no' +
                               ' default for platform ' + platform)

    # these are filenames minus extension, and are sample IDs
    sample_names = []
    primer_seqs = []  # same order as sample_names
    bc_seqs = []  # same order as sample_names
    for i in range(len(map_data)):
        sample_names.append(map_data[i][headers.index('SampleID')])
        bc_seqs.append(map_data[i][headers.index('BarcodeSequence')])
        primer = (map_data[i][headers.index('LinkerPrimerSequence')])
        for char, bases in DNASequence.iupac_degeneracies().iteritems():
            primer = primer.replace(char, '[' + ''.join(bases) + ']')
        primer_seqs.append(primer)

    if len(set(primer_seqs)) != 1:
        raise RuntimeError(
            'Error: only one primer per mapping file supported.')
    one_primer = primer_seqs[0]

    commands = []

    if logger is None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    log_input_md5s(logger, [mapping_fp, sff_txt_fp])

    # execute commands in output_dir
    called_dir = os.getcwd()
    os.chdir(output_dir)
    fh = open(os.path.join(output_dir, 'map.csv'), 'w')
    for i in range(len(sample_names)):
        fh.write(sample_names[i] + ',' + bc_seqs[i] + '\n')
    fh.close()

    # these are the fasta results, e.g. PC.636_Good.fa
    # later we merge them and copy to output file
    post_pyro_tail = '_' + truncate_len
    if suppress_perseus:
        fasta_result_names = [sample_name + post_pyro_tail + '_seqnoise_cd.fa'
                              for sample_name in sample_names]
    else:
        fasta_result_names = [sample_name + '_Good.fa'
                              for sample_name in sample_names]

    cmd = 'cd ' + output_dir  # see also os.chdir above
    commands.append([('change to output dir', cmd)])

    cmd = 'echo $PYRO_LOOKUP_FILE > pyro_lookup_filepath.txt'
    commands.append([('confirm pyro lookup filepath environment variable',
                      cmd)])

    cmd = 'SplitKeys.pl ' + one_primer + ' map.csv < ' +\
        os.path.join(called_dir, sff_txt_fp) +\
        ' > splitkeys_log.txt 2> unassigned.fna'
    commands.append([('split sff.txt via barcodes (keys)', cmd)])

    for i, sample_name in enumerate(sample_names):

        # Build the summarize taxonomy command
        if platform == 'flx':
            cmd = 'Clean360.pl ' + one_primer + ' ' + sample_name + ' < ' +\
                sample_name + '.raw'
            commands.append([('clean flows ' + sample_name, cmd)])

            # these run through the whole sff file once per sample, I think
            # cmd = "FlowsFA.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        elif platform == 'titanium':
            cmd = 'CleanMinMax.pl ' + one_primer + ' ' + sample_name + ' < ' +\
                sample_name + '.raw'
            commands.append([('clean flows ' + sample_name, cmd)])

            # cmd = "FlowsMinMax.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        else:
            raise RuntimeError("platform " + platform + " not supported")

        cmd = "mpirun -np " + str(numnodes) + " PyroDist -in " +\
            sample_name + ".dat -out " + \
            sample_name + " > " + sample_name + ".pdout"
        commands.append([('pyrodist ' + sample_name, cmd)])

        cmd = "FCluster -in " + sample_name + ".fdist -out " + sample_name +\
            " > " + sample_name + ".fcout"
        commands.append([('fcluster pyrodist ' + sample_name, cmd)])

# e.g.:
# mpirun -np 2 PyroNoise -din PC.354.dat -out PC.354_pyronoise -lin
# PC.354.list -s 60.0 -c 0.01 > PC.354_pyronoise.pnout
        cmd = "mpirun -np " + str(numnodes) + " PyroNoise -din " +\
            sample_name + ".dat -out " +\
            sample_name + "_pyronoise " + "-lin " +\
            sample_name + ".list -s 60.0 -c 0.01 > " +\
            sample_name + "_pyronoise.pnout"
        commands.append([('pyronoise ' + sample_name, cmd)])

        cmd = 'Parse.pl ' + bc_seqs[i] + one_primer + ' ' + truncate_len + ' < ' +\
            sample_name + '_pyronoise_cd.fa' + ' > ' + sample_name + '_' +\
            truncate_len + '.fa'
        commands.append([('truncate ' + sample_name, cmd)])

        # now start with post_pyro_tail
        cmd = "mpirun -np " + str(numnodes) + " SeqDist -in " +\
            sample_name + post_pyro_tail +\
            ".fa > " + sample_name + post_pyro_tail + ".seqdist"
        commands.append([('seqdist ' + sample_name, cmd)])

        cmd = "FCluster -in " + sample_name + post_pyro_tail + ".seqdist -out " +\
            sample_name + post_pyro_tail + "fcl > " +\
            sample_name + post_pyro_tail + ".fcout"
        commands.append([('fcluster seqdist ' + sample_name, cmd)])

# e.g.:
# mpirun -np 2 SeqNoise -in PC.354_pyronoise_cd.fa -din
# PC.354_pyronoise_cd.seqdist -out PC.354_pyronoise_cd_seqnoise -lin
# PC.354_pyronoise_cdfcl.list -min PC.354_pyronoise.mapping -s 30.0 -c 0.08 >
# PC.354_pyronoise_cd.snout

        cmd = "mpirun -np " + str(numnodes) + " SeqNoise -in " +\
            sample_name + post_pyro_tail +\
            ".fa -din " + sample_name + post_pyro_tail + ".seqdist -out " +\
            sample_name + post_pyro_tail +\
            "_seqnoise -lin " + sample_name + post_pyro_tail + 'fcl.list -min ' +\
            sample_name + '_pyronoise' +\
            '.mapping -s ' + seqnoise_resolution + ' -c 0.08 > ' +\
            sample_name + post_pyro_tail + '.snout'
        commands.append([('seqnoise ' + sample_name, cmd)])

        if not suppress_perseus:

            cmd = 'Perseus -sin ' + sample_name + post_pyro_tail +\
                '_seqnoise_cd.fa > ' +\
                sample_name + '.per'
            commands.append([('Perseus ' + sample_name, cmd)])

            cmd = 'Class.pl ' + sample_name + '.per ' +\
                str(chimera_alpha) + ' ' + str(chimera_beta) +\
                ' > ' + sample_name + '.class'
            commands.append([('Class.pl ' + sample_name, cmd)])

            cmd = 'FilterGoodClass.pl ' + sample_name + post_pyro_tail +\
                '_seqnoise_cd.fa ' +\
                sample_name + '.class 0.5 > ' + sample_name + '_Chi.fa 2> ' +\
                sample_name + '_Good.fa'
            commands.append([('FilterGoodClass ' + sample_name, cmd)])

        cmd = 'unweight_fasta.py -i %s -o %s -l %s' %\
            (fasta_result_names[i], sample_name + '_unw.fna', sample_name)
        commands.append([('unweight fasta ' + sample_name, cmd)])

    cmd = 'cat ' +\
        ' '.join([sample_name + '_unw.fna' for sample_name in sample_names]) +\
        ' > ' + output_filepath  # this should be an abs filepath
    commands.append([('cat into one fasta file', cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)
