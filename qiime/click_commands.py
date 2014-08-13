import os

import click

from .cli import qiime_cli


@qiime_cli.command()
# I/O options
@click.option('--sequence-read-fp', '-i', multiple=True, required=True,
              type=click.Path(exists=True), help='Input sequence reads')
@click.option('--output-dir', '-o', type=click.Path(exists=False),
              required=True)
@click.option('--mapping_fp', '-m', required=True,
              type=click.File('U'), help='Mapping file')
@click.option('--barcode-read-fp', '-b', multiple=True, required=False,
              type=click.Path(exists=True), help='Barcode read files')
@click.option('--rev-comp/--no-rev-comp', default=False,
              help='Reverse complement sequences on output')
@click.option('--start-seq-id', type=int, default=0,
              help='The starting unique ID for sequences')
# Iterator options
@click.option('--rev-comp-barcodes/--no-rev-comp-barcodes', default=False,
              help='Reverse complement barcode reads')
@click.option('--phred-offset', type=click.Choice(['33', '64']),
              help='The ASCII offset used to decode PHRED scores')
# Runtime options
@click.option('--phred-quality-threshold', '-q', default=3, type=int,
              help='Minimum PHRED quality score')
@click.option('--barcode-type', help='The type of barcode used', default=None,
              type=click.Choice(['golay_12', 'hamming_8', 'not-barcoded']))
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
    """Quality filter and demultiplex sequences"""
    from skbio import DNA
    from skbio.parse.sequences.factory import load

    from qiime.parse import parse_mapping_file_to_dict
    from qiime.process_seqs import SequenceWorkflow, IterAdapter

    # qiime_config is available under ctx.obj['qiime_config']

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

    # should be sourced from skbio
    def format_fasta(id_, seq):
        return ">%s\n%s\n" % (id_, seq)

    # not defining fastq as the method should be sourced from skbio
    # particularly dealing with qual

    # setup sequence iterator
    seqs = load(kwargs.pop('sequence_read_fp'))

    # setup barcode iterator
    barcode_read_fp = kwargs.pop('barcode_read_fp')
    if barcode_read_fp:
        transform = rc_it if kwargs.pop('rev_comp_barcodes') else None
        barcodes = load(barcode_read_fp, transform=transform)
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

    success_fp = os.path.join(output_dir, 'seqs.fna')
    fail_fp = os.path.join(output_dir, 'unassigned.fna')

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
            formatted = format_fasta(id_, item['Sequence'])
            failed.write(formatted) if wf.failed else success.write(formatted)
