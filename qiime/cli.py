import click

import qiime


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo('Version %s' % qiime.__version__)
    ctx.exit()


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
@click.pass_context
def qiime_cli(ctx):
    """QIIME, canonically pronounced 'chime'"""
    pass
