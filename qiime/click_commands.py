import os

import click

from .cli import qiime_cli


@qiime_cli.command()
# I/O options
@click.option('--sequence-read-fp', '-i', multiple=True, required=True,
              type=click.Path(exists=True), help='Input sequence reads')
@click.option('--output-dir', '-o', type=click.Path(exists=False),
              required=True)
@click.option('--mapping-fp', '-m', required=True,
              type=click.File('U'), help='Mapping file')
@click.option('--barcode-read-fp', '-b', multiple=True, required=False,
              type=click.Path(exists=True), help='Barcode read files')
@click.option('--rev-comp/--no-rev-comp', default=False,
              help='Reverse complement sequences on output')
@click.option('--start-seq-id', type=int, default=0,
              help='The starting unique ID for sequences')
@click.option('--to-fastq', is_flag=True,
              help='Write out in fastq')
# Iterator options
@click.option('--rev-comp-barcodes/--no-rev-comp-barcodes', default=False,
              help='Reverse complement barcode reads')
@click.option('--phred-offset', default='33', type=click.Choice(['33', '64']),
              help='The ASCII offset used to decode PHRED scores')
# Runtime options
@click.option('--phred-quality-threshold', '-q', default=3, type=int,
              help='Minimum PHRED quality score')
@click.option('--barcode-type', help='The type of barcode used',
              default='golay_12', type=click.Choice(['golay_12', 'hamming_8',
                                                     'not-barcoded']))
@click.option('--max-barcode-error', default=1.5, type=float,
              help='The maximum number of barcode errors allowed')
@click.option('--retain-primer/--no-retain-primer', default=False,
              help='Whether to retain the primers or not (if applicable)')
@click.option('--max-primer-mismatch', type=int, default=0,
              help='Maximum mismatches allowed within the primers')
@click.option('--min-seq-len', type=int,
              help='The minimum sequence length')
@click.option('--max-ambig-count', default=0, type=int,
              help='Maximum ambiguous bases allowed')
# Other options
@click.option('--rev-comp-mapping-barcodes/--no-rev-comp-mapping-barcodes',
              default=False, help='Reverse complement the mapping barcodes')
@click.pass_context
def slib(ctx, **kwargs):
    """Quality filter and demultiplex sequences

    Examples
    --------

    Demultiplex and quality filter (at Phred >= Q20) one lane of Illumina fastq
    data and write results to ./slout_q20:

    $ qiime slib -i $PWD/lane1_read1.fastq.gz -b $PWD/lane1_barcode.fastq.gz \
            -m $PWD/map.txt -o slout_q20 --rev-comp-mapping-barcodes -q 20

    Demultiplex and quality filter (at Phred >= Q20) two lanes of Illumina
    fastq data and write results to ./slout_q20:

    $ qiime slib -i $PWD/lane1_read1.fastq.gz -i $PWD/lane2_read1.fastq.gz \
            -b $PWD/lane1_barcode.fastq.gz -b $PWD/lane2_barcode.fastq.gz \
            -m $PWD/map.txt -o slout_q20 --rev-comp-mapping-barcodes -q 20
    """
    print kwargs
    from skbio import DNA
    from skbio.parse.sequences.factory import load
    from skbio.format.sequences.fastq import format_fastq_record

    from qiime.parse import parse_mapping_file_to_dict
    from qiime.process_seqs import SequenceWorkflow, IterAdapter

    # qiime_config is available under ctx.obj['qiime_config']

    phred_offset = int(kwargs.pop('phred_offset'))

    # reverse complement for reversing mapping barcodes
    def rc(seq):
        return str(DNA(seq).rc())

    # reverse complement for reversing barcode reads
    def rc_it(st):
        st['Sequence'] = rc(st['Sequence'])
        st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None

    # id formatter for writing
    def format_id(idx, state):
        seq_id = "%s_%d" % (state['Sample'], idx)
        ori_id = state['SequenceID']
        ori_bc = "orig_bc=%s" % state['Original barcode']
        new_bc = "new_bc=%s" % state['Final barcode']
        bc_diff = "bc_diffs=%d" % (state['Barcode errors'] or 0)
        return " ".join([seq_id, ori_id, ori_bc, new_bc, bc_diff])

    # should be sourced from skbio but there doesn't appear to be a method that
    # takes a single seq and ID
    def format_fasta(id_, item):
        return ">%s\n%s\n" % (id_, item['Sequence'])

    def make_format_fastq(offset):
        def f(id_, state):
            seq = state['Sequence']
            qual = state['Qual']
            return format_fastq_record(id_, seq, qual, offset)
        return f
    format_fastq = make_format_fastq(phred_offset)

    # setup sequence iterator
    seqs = load(kwargs.pop('sequence_read_fp'), phred_offset=phred_offset)

    # setup barcode iterator
    barcode_read_fp = kwargs.pop('barcode_read_fp')
    if barcode_read_fp:
        transform = rc_it if kwargs.pop('rev_comp_barcodes') else None
        barcodes = load(barcode_read_fp, transform=transform,
                        phred_offset=phred_offset)
    else:
        barcodes = None

    # load mapping, setup barcode and primer maps
    mapping, comments = parse_mapping_file_to_dict(kwargs.pop('mapping_fp'))
    barcode_map = {v['BarcodeSequence']: k for k, v in mapping.items()}
    primer_map = {v['BarcodeSequence']: v['LinkerPrimerSequence'].split(',')
                  for v in mapping.values()}

    # reverse complement barcodes if necessary
    if kwargs.pop('rev_comp_mapping_barcodes'):
        barcode_map = {rc(k): v for k, v in barcode_map.items()}

    # setup outputs and options
    output_dir = kwargs.pop('output_dir')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    to_fastq = kwargs.pop('to_fastq')
    ext = 'fq' if to_fastq else 'fna'
    formatter = format_fastq if to_fastq else format_fasta
    success_fp = os.path.join(output_dir, 'seqs.%s' % ext)
    fail_fp = os.path.join(output_dir, 'unassigned.%s' % ext)

    if os.path.exists(success_fp):
        raise IOError("%s already exists!" % success_fp)

    if os.path.exists(fail_fp):
        raise IOError("%s already exists!" % fail_fp)

    # setup starting sequence ID and whether to RC on write
    seq_id = kwargs.pop('start_seq_id')
    rc = rc if kwargs.pop('rev_comp') else lambda x: x

    # setup sequence/barcode iterator
    iter_ = IterAdapter(seqs, barcodes)
    wf = SequenceWorkflow(options=kwargs, mapping=mapping,
                          barcodes=barcode_map, primers=primer_map)

    with open(success_fp, 'w') as success, open(fail_fp, 'w') as failed:
        for idx, item in enumerate(wf(iter_, fail_callback=lambda x: x.state)):
            id_ = format_id(seq_id + idx, item)
            formatted = formatter(id_, item)
            failed.write(formatted) if wf.failed else success.write(formatted)
