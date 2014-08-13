import inspect
import os
import re
from unittest import TestCase, main

from click.testing import CliRunner

from qiime.click_commands import slib


def get_qiime_test_data(name):
    """Based off of skbio.util.testing.get_data_path"""
    callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
    path = os.path.dirname(os.path.abspath(callers_filename))
    data_path = os.path.join(path, '../qiime_test_data', name)
    return data_path


class ClickMixin(object):
    runner = CliRunner()

    def find_usage_commands(self, text):
        pat = re.compile(r'\$ qiime.*(?:\n|$)')

        for match in pat.finditer(text):
            yield match.group(0)

    def parse_usage_commands(self, command_name, command_str):
        if not command_str.startswith('$ qiime %s' % command_name):
            raise ValueError("Cannot interpret: %s" % command_str)

        _, args = command_str.split(command_name, 1)
        args = args.replace('$PWD', self.testdata_dir)
        return [arg.strip() for arg in args.split()]

    def test_usage(self):
        for usage in self.find_usage_commands(self.command.help):
            args = self.parse_usage_commands(self.command.name, usage)

            with self.runner.isolated_filesystem():
                result = self.runner.invoke(self.command, args=args)
                if result.exit_code != 0:
                    self.fail(result.output)


class SlibTests(ClickMixin, TestCase):
    def setUp(self):
        self.testdata_dir = get_qiime_test_data('split_libraries_fastq')
        self.command = slib


if __name__ == '__main__':
    main()
