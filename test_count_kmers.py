# Tomas Kazmar, 2012, BSD 2-clause license, see LICENSE
#

from nose.tools import raises

import numpy as np
from kmers import count_kmers, create_kmers_lookup, count_kmers_lookup
from collections import defaultdict

def reverse_complement(s):
    """
    Reverse complement a DNA sequence.
    """
    tab = { 'A': 'T',
            'a': 't',
            'T': 'A',
            't': 'a',
            'G': 'C',
            'g': 'c',
            'C': 'G',
            'c': 'g',
            'N': 'N' }
    s = ''.join([tab[c] for c in s])
    return s[::-1]

def ocount_kmers(k, seq=None):
    if seq is None:
        return None
    map_to_index = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
    seqi_fw = np.array([map_to_index[c] for c in seq])
    map_to_index = { 'A': 3, 'C': 2, 'G': 1, 'T': 0 }
    seqi_rc = np.array([map_to_index[c] for c in seq[::-1]])
    counts = np.zeros(4**k, dtype=float)
    first_char_factor = 2*(k-1)
    for seqi in [seqi_fw, seqi_rc]:
        index = 0
        for l in range(1, k):
            index = (index + seqi[k-1-l])*4
        for i in range(k-1, len(seqi)):
            c = seqi[i]
            index = (index >> 2) + (c << first_char_factor)
            counts[index] += 1
    return counts

def test_count_7mers():
    seq = 'ACGTCGACCGCGTTA'
    counts = count_kmers(7, seq)
    assert counts.sum() == (len(seq)-6)*2
    counts_rc = count_kmers(7, reverse_complement(seq))
    assert (counts == counts_rc).all()
    map_to_index = defaultdict(lambda: -1)
    map_to_index.update({ 'A': 0, 'C': 1, 'G': 2, 'T': 3 })
    seqi_fw = np.array([map_to_index[c] for c in seq])
    map_to_index.update({ 'A': 3, 'C': 2, 'G': 1, 'T': 0 })
    seqi_rc = np.array([map_to_index[c] for c in seq[::-1]])
    for seqi in [seqi_fw, seqi_rc]:
        for i in range(6, len(seqi)):
            index = (((((seqi[i]*4+seqi[i-1])*4+seqi[i-2])*4+seqi[i-3])*4+seqi[i-4])*4+seqi[i-5])*4+seqi[i-6]
            counts[index] -= 1
    assert counts.sum() == 0

def test_count_2mers():
    seq = 'ACGTCG' #ACCGCGTTA'
    counts = count_kmers(2, seq)
    ocounts = ocount_kmers(2, seq)
    assert (counts == ocounts).all()
    assert counts.sum() == (len(seq)-1)*2
    counts_rc = count_kmers(2, reverse_complement(seq))
    assert (counts == counts_rc).all()
    map_to_index = defaultdict(lambda: -1)
    map_to_index.update({ 'A': 0, 'C': 1, 'G': 2, 'T': 3 })
    seqi_fw = np.array([map_to_index[c] for c in seq])
    map_to_index.update({ 'A': 3, 'C': 2, 'G': 1, 'T': 0 })
    seqi_rc = np.array([map_to_index[c] for c in seq[::-1]])
    for seqi in [seqi_fw, seqi_rc]:
        for i in range(1, len(seqi)):
            index = seqi[i]*4 + seqi[i-1]
            counts[index] -= 1
    assert counts.sum() == 0

def test_kmers_lookup():
    seq = 'ACGTCGACCGCGTTA' * 100
    counts = count_kmers(2, seq)
    lookup = create_kmers_lookup(2, seq)
    lcounts = count_kmers_lookup(lookup, 0, len(seq))
    assert (counts == lcounts).all()


#@raises(ValueError)
#def test_non_acgt():
#    count_kmers('abcdefghijklmnopqrstvwuxyz')
