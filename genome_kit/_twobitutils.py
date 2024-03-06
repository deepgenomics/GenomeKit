# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from array import array
from itertools import product


# Copied from twobitreader
def _true_long_type():
    """
    OS X uses an 8-byte long, so make sure L (long) is the right size
    and switch to I (int) if needed
    """
    for type_ in ['L', 'I']:
        test_array = array(type_, [0])
        long_size = test_array.itemsize
        if long_size == 4:
            return type_
    raise ImportError("Couldn't determine a valid 4-byte long type to use as \
                       equivalent to _LONG")


_LONG = _true_long_type()

_TWOBIG_SIG_LILEND = 0x1A412743
_TWOBIG_SIG_BIGEND = 0x4327411A
_TWOBIG_VER_LILEND = 0x00000000
_TWOBIG_VER_BIGEND = 0x00000000

_BASES = ['T', 'C', 'A', 'G', 'N']

_BASE_AS_BITS = {
    'T': 0,
    't': 0,
    'C': 1,
    'c': 1,
    'A': 2,
    'a': 2,
    'G': 3,
    'g': 3,
    'N': 0,
    'n': 0,
}


def _quad_as_bits(s):
    """Packs a string of up to four bases (TCAG) into a byte representation."""
    if len(s) < 4:
        s += 'N' * (4 - len(s))
    return _BASE_AS_BITS[s[0]] << 6 | \
           _BASE_AS_BITS[s[1]] << 4 | \
           _BASE_AS_BITS[s[2]] << 2 | \
           _BASE_AS_BITS[s[3]]  # noqa


def _basecombos(reps):
    args = [_BASES] * reps  # List of lists
    for x in product(*args):
        yield "".join(x)


def _find_blocks(seq, includes):
    starts = []
    sizes = []

    # Loop over each position in the string, taking note of each
    # transition from "inside" a block to "outside" a block
    was_inside = False
    for i, c in enumerate(seq):
        is_inside = c in includes
        if is_inside != was_inside:
            if was_inside:
                sizes.append(i - starts[-1])
            else:
                starts.append(i)
            was_inside = is_inside

    if was_inside:
        sizes.append(len(seq) - starts[-1])

    assert len(starts) == len(sizes)
    return {'starts': starts, 'sizes': sizes}


def write2bit(outfile, seqs):
    """Write a 2bit file containing the given DNA sequences.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    seqs : dict
        A dictionary of named DNA sequence strings each of which will be
        represented by name. Sequence names can be arbitrary. Sequences may
        contain a mix of upper/lower case ACGTN characters.

    Examples
    --------
      >>> seqs={
      ...   "chr1" : "AAAANNnnCCAAaattTT",
      ...   "chr2" : "NNNnnnnggggGGGG",
      ... }
      >>> write2bit('mydna.2bit', seqs)
    """

    # Lookup table for all possible substrings of
    QUAD_AS_BITS = {}
    QUAD_AS_BITS.update({q: _quad_as_bits(q) for q in _basecombos(4)})
    QUAD_AS_BITS.update({q: _quad_as_bits(q) for q in _basecombos(3)})
    QUAD_AS_BITS.update({q: _quad_as_bits(q) for q in _basecombos(2)})
    QUAD_AS_BITS.update({q: _quad_as_bits(q) for q in _basecombos(1)})

    # Convert each quad of 4 bases to a byte
    # Pad the resulting byte list with 0s until it's a multiple of 4
    seq_to_bits = lambda seq: \
        [QUAD_AS_BITS[seq[i:i+4]] for i in range(0, len(seq), 4)] \
        + [0]*(((len(seq)+3)//4) % 4)

    # header layout on disk:
    #   signature      (4 bytes)
    #   version        (4 bytes)
    #   sequence_count (4 bytes)
    #   reserved       (4 bytes)
    headsize = 16

    # index layout on disk:
    #   name_size     (1 byte)
    #   name          (1 byte x name_size, no NULL terminator)
    #   offset        (4 bytes)
    idxsize = lambda name: 1 + len(name) + 4

    # record layout on disk:
    #   dna_size      (4 bytes)
    #   nblock_count  (4 bytes)
    #   nblock_starts (4 bytes x number of nblocks)
    #   nblock_sizes  (4 bytes x number of nblocks)
    #   mask_count    (4 bytes)
    #   mask_starts   (4 bytes x number of masks)
    #   mask_sizes    (4 bytes x number of masks)
    #   reserved      (4 bytes)
    #   bits          ((dna_size+3//4) bytes)
    recsize = lambda rec: 4 + \
                          4+8*len(rec['nblocks']['starts']) + \
                          4+8*len(rec['masks']['starts']) + \
                          4 + \
                          len(rec['bits'])  # noqa

    names = sorted(seqs.keys())
    recs = []

    for name in names:
        seq = seqs[name]
        rec = {}
        rec['dnasize'] = len(seq)
        rec['nblocks'] = _find_blocks(seq, 'Nn')
        rec['masks'] = _find_blocks(seq, "acgtn")
        rec['bits'] = seq_to_bits(seq.upper())
        recs.append(rec)

    with open(outfile, 'wb') as f:
        # Write the header
        array(_LONG, [_TWOBIG_SIG_LILEND, _TWOBIG_VER_LILEND, len(seqs), 0]).tofile(f)

        # Write the index, one entry for each sequence
        for i, name in enumerate(names):
            rec_offset = headsize + sum(map(idxsize, names)) \
                                  + sum(map(recsize, recs[:i]))

            # Write the name of this sequence
            array('B', [len(name)]).tofile(f)
            array('B', name.encode('ascii')).tofile(f)
            array(_LONG, [rec_offset]).tofile(f)

        # Write the sequence records themselves
        for rec, name in zip(recs, names):
            # Write the actual DNA sequence length
            array(_LONG, [rec['dnasize']]).tofile(f)

            # Write the nblocks
            array(_LONG, [len(rec['nblocks']['starts'])]).tofile(f)
            array(_LONG, rec['nblocks']['starts']).tofile(f)
            array(_LONG, rec['nblocks']['sizes']).tofile(f)

            # Write the masks
            array(_LONG, [len(rec['masks']['starts'])]).tofile(f)
            array(_LONG, rec['masks']['starts']).tofile(f)
            array(_LONG, rec['masks']['sizes']).tofile(f)

            # Write reserved DWORD
            array(_LONG, [0]).tofile(f)  # reserved

            # Write DNA bits (size rounded up to nearest byte)
            array('B', rec['bits']).tofile(f)
