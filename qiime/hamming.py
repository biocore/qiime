#!/usr/bin/env python
"""
Hamming(5, 11) decoder for 8 nt codewords (16 bits)

Requires python and numpy

Run demo via command line: "python hamming.py"

Author: Micah Hamady (hamady@colorado.edu)
"""
from numpy import array
from functools import reduce
__author__ = "Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Micah Hamady", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


# current encoding scheme
INT_TO_BS = {0: "00", 1: "01", 2: "10", 3: "11"}
CUR_ENC_FO = {'A': 3, 'C': 2, 'T': 0, 'G': 1}
CUR_REV_ENC_SI = {"11": "A", "10": "C", "00": "T", "01": "G"}


def calc_parity_vector(parity_vector):
    """ Returns even or odd parit for parity vector """
    return reduce(lambda x, y: x ^ y, parity_vector[1:])


def calc_syndrome(codeword, n):
    """ Calculate syndrome and correct codeword if possible """
    sym = 0
    for i in range(1, n):
        if codeword[i]:
            sym ^= i
    extra_parity = calc_parity_vector(codeword)
    if extra_parity == codeword[0]:
        if sym == 0:
            return 0, sym
        else:
            return 2, sym
    else:
        if sym >= n:
            pass
        else:
            codeword[sym] ^= 1
        return 1, sym


def nt_to_cw(cur_enc, cur_nt):
    """ Convert nt sequence to codeword """
    return array(map(int, ''.join([INT_TO_BS[cur_enc[x]] for x in cur_nt])))


def unpack_bitstr(rev_cur_bit, bitstr):
    """ Unpack bistring into nt sequence """
    bstr_len = len(bitstr)
    return (
        ''.join([rev_cur_bit[bitstr[i:i + 2]] for i in range(0, bstr_len, 2)])
    )


def decode_barcode_8(nt_barcode):
    """ Decode length 8 barcode (16 bits) """
    # check proper length
    if len(nt_barcode) != 8:
        raise ValueError("barcode must be 8 nt long.")

    # check valid characters
    if set(list(nt_barcode)).difference(CUR_ENC_FO.keys()):
        raise ValueError("Only A,T,C,G valid chars.")

    # decode
    decoded = nt_to_cw(CUR_ENC_FO, nt_barcode)
    num_errors, sym = calc_syndrome(decoded, 16)

    # convert corrected codeword back to nt sequence
    if num_errors == 1:
        nt_barcode = unpack_bitstr(CUR_REV_ENC_SI, ''.join(map(str, decoded)))
    elif num_errors > 1:
        nt_barcode = None

    return nt_barcode, num_errors / 2.0
# alt name for function to support consistency with golay
decode_hamming_8 = decode_barcode_8

# mapping from barcode used to original sample id
DEMO_SAMPLE_MAPPING = {
    "AACCATGC": "Sample 1",
    "TCGTAGCA": "Sample 2",
    "ACACCTCT": "Sample 3",
    "CTTCCTAG": "Sample 4",
    "GGTAGCTT": "Sample 5"
}

# dummy barcode tagged sequences (e.g. from 454 run)
DEMO_SEQUENCES = [
    # these all have no bit errors
    "AACCATGCTTTTTTTTTTTTTTTTTTTTTTTTT_0",
    "TCGTAGCATTTTTTTTTTTTTTTTTTTTTTTTT_1",
    "ACACCTCTTTTTTTTTTTTTTTTTTTTTTTTTT_2",
    "CTTCCTAGTTTTTTTTTTTTTTTTTTTTTTTTT_3",
    "GGTAGCTTTTTTTTTTTTTTTTTTTTTTTTTTT_4",

    # these all have single bit errors can correct
    "ACCCATGCAAAAAAAAAAAAAAAAAAAAAAAAA_5",
    "TCGTAGCCAAAAAAAAAAAAAAAAAAAAAAAAA_6",
    "ACACATCTAAAAAAAAAAAAAAAAAAAAAAAAA_7",
    "TTTCCTAGAAAAAAAAAAAAAAAAAAAAAAAAA_8",
    "GGTAGCTTAAAAAAAAAAAAAAAAAAAAAAAAA_9",

    "AACCACGCAAAAAAAAAAAAAAAAAAAAAAAAA_10",
    "TAGTAGCAAAAAAAAAAAAAAAAAAAAAAAAAA_11",
    "ACACCTCTAAAAAAAAAAAAAAAAAAAAAAAAA_12",
    "CGTCCTAGAAAAAAAAAAAAAAAAAAAAAAAAA_13",
    "GTTAGCTTAAAAAAAAAAAAAAAAAAAAAAAAA_14",

    "AGCCATGCAAAAAAAAAAAAAAAAAAAAAAAAA_15",
    "GCGTAGCAAAAAAAAAAAAAAAAAAAAAAAAAA_16",
    "ACACCTCTAAAAAAAAAAAAAAAAAAAAAAAAA_17",
    "ATTCCTAGAAAAAAAAAAAAAAAAAAAAAAAAA_18",
    "GGTAGATTAAAAAAAAAAAAAAAAAAAAAAAAA_19",

    # these all have double bit errors, cannot correct
    "TACCATGCCCCCCCCCCCCCCCCCCCCCCCCCC_20",
    "TCCTAGCACCCCCCCCCCCCCCCCCCCCCCCCC_21",
    "ACAGCTCTCCCCCCCCCCCCCCCCCCCCCCCCC_22",
    "CTTCGTAGCCCCCCCCCCCCCCCCCCCCCCCCC_23",
    "GCTAGCTTCCCCCCCCCCCCCCCCCCCCCCCCC_24",
]


def demo():
    """ Demo decoding a list of tagged reads from several samples """
    print "---------------------------------------"
    print "Processing %d sequences from %d samples" % (
          len(DEMO_SEQUENCES), len(DEMO_SAMPLE_MAPPING))
    print "---------------------------------------"

    for ix, cur_seq in enumerate(DEMO_SEQUENCES):
        barcode = cur_seq[:8]
        seq_read = cur_seq[8:]
        print "---> processing demo sequence", ix
        print "read barcode      :", barcode
        try:
            corrected_barcode = decode_barcode_8(barcode)
            orig_sample_id = DEMO_SAMPLE_MAPPING[corrected_barcode]

            if corrected_barcode != barcode:
                print "*corrected barcode:", corrected_barcode
            else:
                print "-no error  barcode:", corrected_barcode

            print "original sample id:", orig_sample_id
            print "sequence read     :", seq_read

        except ValueError as e:
            print "!", str(e), "skipping..."
            continue
