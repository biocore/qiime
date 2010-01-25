#!/usr/bin/env python
import numpy

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"


"""
module provides golay encoding of DNA barcodes.  
may not be the same golay code as previously used.
Provides decode() and encode()

If you wish to assign a read DNA seq barcode to a list of known originals,
correcting for errors if necessary, I recommend you use the generic version
in barcode.py .  Golay decoding assumes that the sequence can be any valid
golay sequence, not just those you used in the study.

this module implements *a* golay (24,12,8) code.  That's 24 bit codewords,
12 bits of information, min 8 bits (hamming distance) between any 2 codewords
there are 2**12 = 4096 codewords, all <= 3 bit errors are corrected, 4 bit
errors are detected, and worse errors can cause erroneously corrected codewords.

Since DNA has 4 nucleotides, each base represents 2 bits (see 
DEFAULT_NT_TO_BITS).  The default values represent A <--> T and  C <--> G as
2 bit errors, all other mutations are 1 bit errors.
e.g.: 
cw1 = AAACCCGGGTTT (24 bits, 12 nucleotides)
cw2 = AAACCCGGGTTA (bit distance of 2 from cw1)
cw3 = AAACCCGGGTTG (bit distance of 1 from cw1)

A specific golay code was chosen, as there are at least 2.
those defaults are from http://www-math.mit.edu/phase2/UJM/vol1/MKANEM~1.PDF

other refs:
http://cwww.ee.nctu.edu.tw/~ftchien/class/ecc08f/topic2.pdf
error correcting coding theory (Rhee)
the art of error correcting coding (Morelos-Zaragoza)
"""
DEFAULT_GOLAY_NT_TO_BITS = { "A":"11",  "C":"00", "T":"10", "G":"01"}
DEFAULT_HAMMING_NT_TO_BITS = { "A":"11",  "C":"10", "T":"00", "G":"01"}

def decode(seq, nt_to_bits=DEFAULT_GOLAY_NT_TO_BITS, code=None):
    """ returns (best match sequence, num_biterrors) for a golay (24,12,8) code

    the seq is converted to bits using nt_to_bits,
    which is then matched against the legal golay 24 bit codewords
    to obtain the closest match.  there are at least 2 sets of legal golay 
    24 bit codewords, I chose one arbitrarily as the default.
    
    inputs:
    * seq is e.g.: 'ACATCCGCACAC'
    * nt_to_bits: e.g: {"A":"00", "C":"01", "G":"10", "T":"11"}, or omit
    to use default for golay (different from hamming default)
    * code: a list of valid 24 bit codewords, or omit to use default 

    returns:
    (match_seq, num_biterrors/2).  match_seq is best match, eg 'ACATCCGCACAA', 
    or None if the query seq is equidistant from multiple words (4 bit errors).
    num_biterrors is divided by 2 to maintain compatibility with API for
    general decoder, which works in terms of # bases.
    """
    if code == None:
        code = DEFAULT_GOLAY_CODE

    received_bits = _seq_to_bits(seq, nt_to_bits)
    corrected_bits, num_errors = _decode_bits(received_bits, code)
    if corrected_bits == None:
        return None, num_errors
    else:
        # put match into nucleotide format
        return _bits_to_seq(corrected_bits, nt_to_bits), num_errors/2.0

def encode(bits, nt_to_bits=None, generator_mtx=None):
    """ takes 12 bits, returns the golay 24bit codeword in nucleotide format

    bits is a list/array, 12 long, e.g.: [0,0,0,0,0,0,0,0,0,1,0,0]
    nt_to_bits is e.g.: {"A":"00", "C":"01", "G":"10", "T":"11"},None => default
    generator_mtx specifies which golay code to use.  None => default
    """
    if nt_to_bits == None:
        nt_to_bits = DEFAULT_GOLAY_NT_TO_BITS
    if generator_mtx == None:
        G = DEFAULT_GENERATOR_MTX

    bits = numpy.array(bits)
    bits = bits.reshape((1,12))

    # cheap way to do binary xor in matrix dot
    res = numpy.dot(bits, G)
    codeword = divmod(res.ravel(),2)[1]
 
    return _bits_to_seq(codeword, nt_to_bits)


def _get_default_golay_code():
    """ uses http://www-math.mit.edu/phase2/UJM/vol1/MKANEM~1.PDF
    returns (code, G)
    code is list of all 24 bit valid codewords (numpy 2d)
    G is generator matrix
    """
    B = numpy.array([
        [1,1,0,1,1,1,0,0,0,1,0,1],
        [1,0,1,1,1,0,0,0,1,0,1,1],
        [0,1,1,1,0,0,0,1,0,1,1,1],
        [1,1,1,0,0,0,1,0,1,1,0,1],
        [1,1,0,0,0,1,0,1,1,0,1,1],
        [1,0,0,0,1,0,1,1,0,1,1,1],
        [0,0,0,1,0,1,1,0,1,1,1,1],
        [0,0,1,0,1,1,0,1,1,1,0,1],
        [0,1,0,1,1,0,1,1,1,0,0,1],
        [1,0,1,1,0,1,1,1,0,0,0,1],
        [0,1,1,0,1,1,1,0,0,0,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,0],
        ])

    P = numpy.array([
        [0,1,1,1,1,1,1,1,1,1,1,1,],
        [1,1,1,0,1,1,1,0,0,0,1,0,],
        [1,1,0,1,1,1,0,0,0,1,0,1,],
        [1,0,1,1,1,0,0,0,1,0,1,1,],
        [1,1,1,1,0,0,0,1,0,1,1,0,],
        [1,1,1,0,0,0,1,0,1,1,0,1,],
        [1,1,0,0,0,1,0,1,1,0,1,1,],
        [1,0,0,0,1,0,1,1,0,1,1,1,],
        [1,0,0,1,0,1,1,0,1,1,1,0,],
        [1,0,1,0,1,1,0,1,1,1,0,0,],
        [1,1,0,1,1,0,1,1,1,0,0,0,],
        [1,0,1,1,0,1,1,1,0,0,0,1,],
    ])  #from http://courses.csusm.edu/math540ak/codes.pdf 
    #P is the version used in original code, but can't figure
    #out how to get the same result...

    #~ G = numpy.concatenate((numpy.eye(12,dtype="int") , B),axis=1)
    #~ G = numpy.concatenate((B, numpy.eye(12,dtype="int")),axis=1)
    #~ G = numpy.concatenate((numpy.eye(12,dtype="int") , P),axis=1)
    G = numpy.concatenate((P, numpy.eye(12,dtype="int")),axis=1)

    all_12bits = _make_12bits()
    code = numpy.zeros((2**12,24),dtype="int")
    for i, u in enumerate(all_12bits):
        umtx = u.reshape((1,12))
        res = numpy.dot(umtx, G)
        # cheap way to do binary xor in matrix dot
        code[i] = divmod(res.ravel(),2)[1] 
    return code, G

def _make_12bits():
    n=12
    res = numpy.zeros((2**12,12),dtype="int")
    for i in range(2**12):
        res[i] = tuple((0,1)[i>>j & 1] for j in xrange(n-1,-1,-1)) 
    return res


def _hamming_dist(v1,v2):
    """ two binary arrays, equal len"""
    edits = (v1!=v2)
    return edits.sum()


    
def _bits_to_seq(bits, nt_to_bits):
    """ e.g.: array([1,1,1,1,0,1]) -> "AAG"
    """
    bits_to_nt = dict(zip(nt_to_bits.values(), nt_to_bits.keys()))
    seq = ""
    for i in range(0,len(bits),2): #take bits in twos
        bit1 = str(int(round(bits[i])))
        bit2 = str(int(round(bits[i+1])))
        seq += bits_to_nt[bit1+bit2]
    return seq

def _seq_to_bits(seq, nt_to_bits):
    """ e.g.: "AAG" -> array([1,1,1,1,0,1])
    """
    bitstring = ""
    for nt in seq:
        bitstring += nt_to_bits[nt]
    bits = numpy.array(map(int,bitstring))
    return bits

def _decode_bits(query_bits, code):
    """ query bits is matched against all codewords in code
    both are lists of bits, not strings or DNA sequences"""
    dists = [_hamming_dist(query_bits, codeword) for codeword in code]
    min_dist = min(dists)
    number_mins = dists.count(min_dist)
    if number_mins > 1:
        return None, min_dist
    else:
        best_hit = code[dists.index(min_dist)]
        return best_hit, min_dist

DEFAULT_GOLAY_CODE, DEFAULT_GENERATOR_MTX = _get_default_golay_code()
