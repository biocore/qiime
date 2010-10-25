#!/usr/bin/env python
"""Parser for 454 Flowgram files in native binary format."""

__author__ = 'Kyle Bittinger'
__copyright__ = 'Copyright 2010, The Cogent Project'
__license__ = 'GPL'
__version__ = '1.5.0.dev'
__credits__ = ['Kyle Bittinger']
__maintainer__ = 'Kyle Bittinger'
__email__ = 'kylebittinger@gmail.com'
__status__ = 'Prototype'

from cStringIO import StringIO
import string
import struct

# Sections were inspired by, but not derived from, several other implementations:
# * BioPython (biopython.org)
# * sff_extract (www.melogen.upv.es/sff_extract)
# * Mothur (mothur.org)


class NamedStruct(struct.Struct):
    """Enhanced Struct class that associates names with each item in the struct.
    """

    def __init__(self, format, keys):
        """Create a new NamedStruct with a list of keys.
        """
        self.keys = keys
        super(NamedStruct, self).__init__(format)

    def read_from(self, file):
        """Read the struct from a file object and return the values as a dict.
        """
        buff = file.read(self.size)
        return self.unpack(buff)

    def pack(self, dict_of_vals):
        vals = [dict_of_vals[k] for k in self.keys]
        return super(NamedStruct, self).pack(*vals)

    def unpack(self, buffer):
        vals = super(NamedStruct, self).unpack(buffer)
        return dict(zip(self.keys, vals))


def seek_pad(file, unit=8):
    """Set a file's position to the next multiple of a given number.
    """
    position = file.tell()
    rem = position % unit
    if rem != 0:
        padding = unit - rem
        file.seek(padding, 1)        


def write_pad(file, unit=8):
    """Write zeros until the file's position is a multiple of the given number.
    """
    position = file.tell()
    rem = position % unit
    if rem != 0:
        num_bytes = unit - rem
        padding_bytes = '\x00' * num_bytes
        file.write(padding_bytes)


common_header_fields = [
    'magic_number',
    'version',
    'index_offset',
    'index_length',
    'number_of_reads',
    'header_length',
    'key_length',
    'number_of_flows_per_read',
    'flowgram_format_code',
    ]


common_header_struct = NamedStruct('>IIQIIHHHB', common_header_fields)


def parse_common_header(sff_file):
    """Parse a Common Header section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Common Header section.
    """
    h = common_header_struct.read_from(sff_file)
    h['flow_chars'] = sff_file.read(h['number_of_flows_per_read'])
    h['key_sequence'] = sff_file.read(h['key_length'])
    seek_pad(sff_file)
    return h


def write_common_header(sff_file, header):
    """Write a common header section to a binary SFF file.
    """
    header_bytes = common_header_struct.pack(header)
    sff_file.write(header_bytes)
    sff_file.write(header['flow_chars'])
    sff_file.write(header['key_sequence'])
    write_pad(sff_file)


common_header_formats = [
    '  Magic Number:  0x%X\n',
    '  Version:       %04d\n',
    '  Index Offset:  %d\n',
    '  Index Length:  %d\n',
    '  # of Reads:    %d\n',
    '  Header Length: %d\n',
    '  Key Length:    %d\n',
    '  # of Flows:    %d\n',
    '  Flowgram Code: %d\n',
    ]


def format_common_header(header):
    """Format a dictionary representation of an SFF common header as text.
    """
    out = StringIO()
    out.write('Common Header:\n')
    for key, fmt in zip(common_header_fields, common_header_formats):
        val = header[key]
        out.write(fmt % val)
    out.write('  Flow Chars:    %s\n' % header['flow_chars'])
    out.write('  Key Sequence:  %s\n' % header['key_sequence'])
    return out.getvalue()


class UnsupportedSffError(Exception):
    pass


def validate_common_header(header):
    """Validate the Common Header section of a binary SFF file.

    Raises an UnsupportedSffError if the header is not supported.
    """
    supported_values = {
        'magic_number': 0x2E736666,
        'version': 1,
        'flowgram_format_code': 1,
        }
    for attr_name, expected_value in supported_values.items():
        observed_value = header[attr_name]
        if observed_value != expected_value:
            raise UnsupportedSffError(
                '%s not supported. (Expected %s, observed %s)' % (
                    attr_name, expected_value, observed_value))


read_header_fields = [
    'read_header_length',
    'name_length',
    'number_of_bases',
    'clip_qual_left', 
    'clip_qual_right',
    'clip_adapter_left',
    'clip_adapter_right',
    ]


read_header_struct = NamedStruct('>HHIHHHH', read_header_fields)


def parse_read_header(sff_file):
    """Parse a Read Header section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Read Header section.
    """
    data = read_header_struct.read_from(sff_file)
    data['Name'] = sff_file.read(data['name_length'])
    seek_pad(sff_file)
    return data


def write_read_header(sff_file, read_header):
    """Write a read header section to a binary SFF file.
    """
    header_bytes = read_header_struct.pack(read_header)
    sff_file.write(header_bytes)
    sff_file.write(read_header['Name'])
    write_pad(sff_file)


read_header_formats = [
    '  Read Header Len:  %d\n',
    '  Name Length:      %d\n',
    '  # of Bases:       %d\n',
    '  Clip Qual Left:   %d\n',
    '  Clip Qual Right:  %d\n',
    '  Clip Adap Left:   %d\n',
    '  Clip Adap Right:  %d\n',
    ]


def format_read_header(read_header):
    """Format a dictionary representation of an SFF read header as text.
    """
    out = StringIO()
    out.write('\n>%s\n' % read_header['Name'])
    timestamp, hashchar, region, location = decode_accession(read_header['Name'])
    out.write('  Run Prefix:   R_%d_%02d_%02d_%02d_%02d_%02d_\n' % timestamp)
    out.write('  Region #:     %d\n' % region)
    out.write('  XY Location:  %04d_%04d\n' % location)
    out.write('\n')
    for key, fmt in zip(read_header_fields, read_header_formats):
        val = read_header[key]
        out.write(fmt % val)
    return out.getvalue()


def parse_read_data(sff_file, number_of_bases, number_of_flows=400):
    """Parse a Read Data section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Read Header section.
    """
    data = {}
    flow_fmt = '>' + ('H' * number_of_flows)
    base_fmt = '>' + ('B' * number_of_bases)
    flow_fmt_size = struct.calcsize(flow_fmt)
    base_fmt_size = struct.calcsize(base_fmt)

    buff = sff_file.read(flow_fmt_size)
    data['flowgram_values'] = struct.unpack(flow_fmt, buff)

    buff = sff_file.read(base_fmt_size)
    data['flow_index_per_base'] = struct.unpack(base_fmt, buff)

    data['Bases'] = sff_file.read(number_of_bases)

    buff = sff_file.read(base_fmt_size)
    data['quality_scores'] = struct.unpack(base_fmt, buff)

    seek_pad(sff_file)
    return data


def write_read_data(sff_file, read_data):
    """Write a read data section to a binary SFF file.
    """
    number_of_flows = len(read_data['flowgram_values'])
    number_of_bases = len(read_data['quality_scores'])
    flow_fmt = '>' + ('H' * number_of_flows)
    base_fmt = '>' + ('B' * number_of_bases)

    flow_bytes = struct.pack(flow_fmt, *read_data['flowgram_values'])
    sff_file.write(flow_bytes)

    index_bytes = struct.pack(base_fmt, *read_data['flow_index_per_base'])
    sff_file.write(index_bytes)
    
    sff_file.write(read_data['Bases'])

    qual_bytes = struct.pack(base_fmt, *read_data['quality_scores'])
    sff_file.write(qual_bytes)

    write_pad(sff_file)


def format_read_data(read_data, read_header):
    """Format a dictionary representation of an SFF read data as text.

    The read data is expected to be in native flowgram format.
    """
    out = StringIO()
    out.write('\n')

    out.write('Flowgram:')
    for x in read_data['flowgram_values']:
        out.write('\t%01.2f' % (x * 0.01))
    out.write('\n')

    out.write('Flow Indexes:')
    current_index = 0
    for i in read_data['flow_index_per_base']:
        current_index = current_index + i
        out.write('\t%d' % current_index)
    out.write('\n')

    out.write('Bases:\t')
    # Roche uses 1-based indexing
    left_idx = read_header['clip_qual_left'] - 1
    right_idx = read_header['clip_qual_right'] - 1
    for i, base in enumerate(read_data['Bases']):
        if (i < left_idx) or (i > right_idx):
            out.write(base.lower())
        else:
            out.write(base.upper())
    out.write('\n')

    out.write('Quality Scores:')
    for score in read_data['quality_scores']:
        out.write('\t%d' % score)
    out.write('\n')

    return out.getvalue()


def parse_read(sff_file, number_of_flows=400):
    """Parse a single read from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation for the Read Header and Read Data sections.

    As a side effect, sets the position of the file object to the end
    of the Read Data section.
    """
    header_data = parse_read_header(sff_file)
    read_data = parse_read_data(
        sff_file, header_data['number_of_bases'], number_of_flows)
    read_data.update(header_data)
    return read_data


def write_read(sff_file, read):
    """Write a single read to a binary SFF file.
    """
    write_read_header(sff_file, read)
    write_read_data(sff_file, read)


def format_read(read):
    """Format a dictionary representation of an SFF read as text.
    """
    out = StringIO()
    out.write(format_read_header(read))
    out.write(format_read_data(read, read))
    return out.getvalue()


def parse_binary_sff(sff_file, native_flowgram_values=False):
    """Parse a binary SFF file, returning the header and a sequence of reads.

    In the binary file, flowgram values are stored as integers, 100
    times larger than the normalized floating point value.  Because
    the conversion is relatively expensive, we allow the computation
    to be skipped if the keyword argument native_flowgram_values is
    True.
    """
    header = parse_common_header(sff_file)
    number_of_flows = header['number_of_flows_per_read']
    validate_common_header(header)
    def get_reads():
        for i in range(header['number_of_reads']):

            # Skip the index section
            if sff_file.tell() == header['index_offset']:
                sff_file.seek(header['index_length'], 1)

            read = parse_read(sff_file, number_of_flows)

            if not native_flowgram_values:
                read['flowgram_values'] = [x * 0.01 for x in read['flowgram_values']]

            yield read
    return header, get_reads()


def write_binary_sff(sff_file, header, reads):
    """Write a binary SFF file, using provided header and read dicts.
    """
    sff_file.seek(0)
    sff_file.truncate()
    write_common_header(sff_file, header)
    for read in reads:
        write_read(sff_file, read)


def format_binary_sff(sff_file, output_file=None):
    """Write a text version of a binary SFF file to an output file.

    If no output file is provided, an in-memory file-like buffer is
    used (namely, a StringIO object).
    """
    if output_file is None:
        output_file = StringIO()
    header, reads = parse_binary_sff(sff_file, True)
    output_file.write(format_common_header(header))
    for read in reads:
        output_file.write(format_read(read))
    return output_file


def base36_encode(n):
    """Convert a positive integer to a base36 string.

    Following the conventions outlined in the Roche 454 manual, the
    numbers 0-25 are represented by letters, and the numbers 36-35 are
    represented by digits.

    Based on the code example at http://en.wikipedia.org/wiki/Base_36
    """
    if n < 0:
        raise ValueError('Only poitive numbers are supported.')
    chars = []
    while n != 0:
        n, remainder = divmod(n, 36)
        chars.append(base36_encode.alphabet[remainder])
    return ''.join(chars)


base36_encode.alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'


def base36_decode(base36_str):
    """Convert a base36 string to a positive integer.

    Following the conventions outlined in the Roche 454 manual, the
    numbers 0-25 are represented by letters, and the numbers 36-35 are
    represented by digits.
    """
    base36_str = base36_str.translate(base36_decode.translation)
    return int(base36_str, 36)


base36_decode.translation = string.maketrans(
    'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789',
    '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ',
    )


def decode_location(location_str):
    """Decode a base36-encoded well location, in Roche 454 format.

    Such timestamps are embedded in the final 5 characters of Roche
    \"universal\" accession numbers.
    """
    return divmod(base36_decode(location_str), 4096)


def decode_timestamp(timestamp_str):
    """Decode a base36-encoded timestamp, in Roche 454 format.

    Such timestamps are embedded in the first 6 characters of Roche
    \"universal\" accession numbers and SFF filenames.
    """
    n = base36_decode(timestamp_str)
    year, n = divmod(n, 13 * 32 * 24 * 60 * 60)
    year = year + 2000
    month, n = divmod(n, 32 * 24 * 60 * 60)
    day, n = divmod(n, 24 * 60 * 60)
    hour, n = divmod(n, 60 * 60)
    minute, second = divmod(n, 60)
    return year, month, day, hour, minute, second


def decode_accession(accession):
    """Decode a Roche 454 \"universal\" accession number.
    """
    assert len(accession) == 14
    timestamp = decode_timestamp(accession[:6])
    hashchar = accession[6]
    region = int(accession[7:9])
    location = decode_location(accession[9:14])
    return timestamp, hashchar, region, location


def decode_sff_filename(sff_filename):
    """Decode a Roche 454 SFF filename, returning a timestamp and other info.
    """
    assert len(sff_filename) == 13
    assert sff_filename.endswith('.sff')
    timestamp = decode_timestamp(sff_filename[:6])
    hashchar = sff_filename[6]
    region = int(sff_filename[7:9])
    return timestamp, hashchar, region
