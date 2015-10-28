#!/usr/bin/env python

from pyfasta import Fasta
from tsh.bio import regions_reader
import sys
from kmers import count_kmers, all_kmers

if __name__ == '__main__':
    genome = Fasta(sys.argv[1])
    k = int(sys.argv[2])
    kmers = 0
    for region in regions_reader(sys.argv[3]):
        seq = genome[region.chrom][region.start:region.stop]
        seq = str(seq).upper()
        for s in seq.split('N'):
            kmers += count_kmers(k, s)
    for kmer, count in zip(all_kmers(k), kmers):
        print '%s\t%d' % (kmer, count)
