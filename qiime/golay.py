#!/usr/bin/env python
import numpy

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""
module provides golay encoding/decoding of DNA barcodes.
may not be the same golay code as previously used.
Provides mainly decode_nt()

If you wish to assign a read DNA seq barcode to a list of known originals,
correcting for errors if necessary, I recommend you use the generic version
in barcode.py .  Golay decoding assumes that the sequence can be any valid
golay sequence (4096 options), not just those you used in the study.

this module implements *a* golay (24,12,8) code.  That's 24 bit codewords,
12 bits of information, min 8 bits (hamming distance) between any 2 codewords
there are 2**12 = 4096 codewords, all <= 3 bit errors are corrected, 4 bit
errors are detected, and worse errors can cause erroneously corrected codewords.

Since DNA has 4 nucleotides, each base represents 2 bits (see
DEFAULT_GOLAY_NT_TO_BITS).  The default values represent A <--> C and  G <--> T
as 2 bit errors, all other nucleotide swaps are 1 bit errors.
e.g.:
cw1 = AAACCCGGGTTT (24 bits, 12 nucleotides)
cw2 = AAACCCGGGTTA (bit distance of 2 from cw1)
cw3 = AAACCCGGGTTG (bit distance of 1 from cw1)

A specific golay code was chosen, as there are multiple.

bitvectors as referenced here are often just listlike objects with ints

refs:
http://www-math.mit.edu/phase2/UJM/vol1/MKANEM~1.PDF
http://cwww.ee.nctu.edu.tw/~ftchien/class/ecc08f/topic2.pdf
error correcting coding theory (Rhee)
the art of error correcting coding (Morelos-Zaragoza)


TO DOs:
* types are messy, numpy arrays, lists, tuples, all doing the same stuff
* haven't tested on all 2**24 bitvectors, could do that to be thorough
* test speed performance
"""


def get_invalid_golay_barcodes(seqs):
    result = []
    for e in seqs:
        if len(e) != 12:
            result.append(e)
        elif decode(e)[1] > 0:
            result.append(e)
    return result


def decode(seq, nt_to_bits=None):
    """decodes a nucleotide string of 12 bases, using bitwise error checking

    inputs:
    - seq, a string of nucleotides
    - nt_to_bits, e.g.: { "A":"11",  "C":"00", "T":"10", "G":"01"}
    output:
    corrected_seq (str), num_bit_errors
    corrected_seq is None if 4 bit error detected"""
    if nt_to_bits is None:
        nt_to_bits = DEFAULT_GOLAY_NT_TO_BITS
    received_bits = _seq_to_bits(seq, nt_to_bits)
    corrected_bits, num_errors = decode_bits(received_bits)  # errors in # bits
    if corrected_bits is None:
        return None, num_errors
    else:
        # put match into nucleotide format
        return _bits_to_seq(corrected_bits, nt_to_bits), num_errors
# alt name for the decode function for consistency with hamming decoding
decode_golay_12 = decode


def encode(bits, nt_to_bits=None):
    """ takes any 12 bits, returns the golay 24bit codeword in nucleotide format

    bits is a list/array, 12 long, e.g.: [0,0,0,0,0,0,0,0,0,1,0,0]
    nt_to_bits is e.g.: {"A":"11", "C":"00", "T":"10", "G":"01"},None => default
    output is e.g.: 'AGTCTATTGGCT'
    """
    if nt_to_bits is None:
        nt_to_bits = DEFAULT_GOLAY_NT_TO_BITS

    bits = numpy.array(bits).reshape((12, 1))

    # cheap way to do binary xor in matrix dot
    res = numpy.dot(DEFAULT_G.T, bits)
    codeword = divmod(res.ravel(), 2)[1]

    return _bits_to_seq(codeword, nt_to_bits)


def decode_bits(received_bitvec):
    """ decode a recieved 24 bit vector to a corrected 24 bit vector

    uses golay defaults
    input: received bitvec is 24 bits long, listlike
    output: corrected_vec, num_bit_errors
    corrected_vec is None iff num_errors = 4"""
    rec = received_bitvec
    syn = numpy.dot(DEFAULT_H, rec) % 2
    try:
        err = numpy.array(DEFAULT_SYNDROME_LUT[tuple(syn)])
    except KeyError:
        return None, 4
    corrected = (rec + err) % 2  # best guess for transmitted bitvector

    return corrected, numpy.sum(err)

# begin support fns


def _make_3bit_errors(veclen=24):
    """ return list of all bitvectors with <= 3 bits as 1's, rest 0's

    returns list of lists, each 24 bits long by default.

    not included:
    [0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0]

    included:
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],

    """
    errorvecs = []
    # all zeros
    errorvecs.append([0] * veclen)
    # one 1
    for i in range(veclen):
        vec = [0] * veclen
        vec[i] = 1
        errorvecs.append(vec)

    # two 1s
    for i in range(veclen):
        for j in range(i + 1, veclen):
            vec = [0] * veclen
            vec[i] = 1
            vec[j] = 1
            errorvecs.append(vec)

    # three 1s
    for i in range(veclen):
        for j in range(i + 1, veclen):
            for k in range(j + 1, veclen):
                vec = [0] * veclen
                vec[i] = 1
                vec[j] = 1
                vec[k] = 1
                errorvecs.append(vec)
    return errorvecs


def _seq_to_bits(seq, nt_to_bits):
    """ e.g.: "AAG" -> array([0,0,0,0,1,0])
    output is array of ints, 1's and 0's

    nt_to_bits is e.g.: {"A":"11", "C":"00", "T":"10", "G":"01"}

    """
    bitstring = ""
    for nt in seq:
        bitstring += nt_to_bits[nt]
    bits = numpy.array(map(int, bitstring))
    return bits


def _bits_to_seq(bits, nt_to_bits):
    """ e.g.: array([0,0,0,0,1,0]) -> "AAG"

    nt_to_bits is e.g.: {"A":"11", "C":"00", "T":"10", "G":"01"}

    """
    bits_to_nt = dict(zip(nt_to_bits.values(), nt_to_bits.keys()))
    seq = ""
    for i in range(0, len(bits), 2):  # take bits in twos
        bit1 = str(int(round(bits[i])))
        bit2 = str(int(round(bits[i + 1])))
        seq += bits_to_nt[bit1 + bit2]
    return seq
# end support fns


# BEGIN module level constants
DEFAULT_GOLAY_NT_TO_BITS = {"A": "11", "C": "00", "T": "10", "G": "01"}

# We use this matrix as the parity submatrix P
DEFAULT_P = numpy.array([
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ],
    [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, ],
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, ],
    [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, ],
    [1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, ],
    [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, ],
    [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, ],
    [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, ],
    [1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, ],
    [1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, ],
    [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, ],
    [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, ],
], dtype='int')  # from http://courses.csusm.edu/math540ak/codes.pdf

# generator mtx G, where transmitted codewords (24bits) are
# G.T dot msg or (msg dot G) (msg is 12 bit message)
# 2**12 = 4096 total codewords, one for each msg
# (all mod 2 arithmetic)
# other G matrices give golay (24,12,8) codes, but this one
# matches existing codes from pre 2010 used in knight lab
DEFAULT_G = numpy.concatenate((DEFAULT_P, numpy.eye(12, dtype="int")), axis=1)

# pairity check matrix H satisfies G dot H.T = zeros (mod 2 arithmetic)
# also satisfies syn = H dot rec = H dot err (rec is recieved 24 bits,
# err is 24 bit error string added to transmitted 24 bit vec)
# (all mod 2 arithmetic)
DEFAULT_H = numpy.concatenate(
    (numpy.eye(12, dtype="int"), DEFAULT_P.T), axis=1)

_ALL_3BIT_ERRORS = _make_3bit_errors()
# len = 2325.  (1 (all zeros) + 24 (one 1) + 276 (two 1s) + 2024)

# syndrome lookup table is the key to (fast, syndrome) decoding
# decode() uses syndrome lookup table

DEFAULT_SYNDROME_LUT = {}
# key: syndrome (12 bits).  Val: 24 bit err for that syn
# we include the all zeros error (key = all zeros syndrome)


# build syndrome lookup table
for errvec in _ALL_3BIT_ERRORS:
    syn = tuple(numpy.dot(DEFAULT_H, errvec) % 2)
    DEFAULT_SYNDROME_LUT[syn] = (errvec)

# END module level constants
