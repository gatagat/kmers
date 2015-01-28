import numpy as np
from kmers import count_kmers, count_kmers_lookup, create_kmers_lookup

class KMerIntersectionKernel(object):
    def __init__(self, k):
        self._k = k

    def __call__(self, s, t):
        return self.one2one(s, t)

    def one2one(self, s, t):
        kmers_s = count_kmers(self._k, seq=s)
        kmers_t = count_kmers(self._k, seq=t)
        return np.minimum(kmers_s, kmers_t).sum()

    def one2many(self, s, ts):
        kmers_s = count_kmers(self._k, seq=s)
        return np.array([np.minimum(kmers_s,
                    count_kmers(self._k, seq=t)).sum()
                for t in ts], dtype=float)


def precompute_kmers_lookups(kmer_k, sequences, id_prefix=None):
    if id_prefix is None:
        id_prefix = ''
    lookups = {}
    for i, s in enumerate(sequences):
        lookups['%s%d' % (id_prefix, i)] = create_kmers_lookup(kmer_k, s)
    return lookups


class LookupKMerIntersectionKernel(object):
    def __init__(self, lookups):
        self._lookups = lookups

    def __call__(self, s, t):
        return self.one2one(s, t)

    def one2one(self, s, t):
        kmers_s = count_kmers_lookup(self._lookups[s.index], s.start, s.stop)
        kmers_t = count_kmers_lookup(self._lookups[t.index], t.start, t.stop)
        return np.minimum(kmers_s, kmers_t).sum()

    def one2many(self, s, ts):
        kmers_s = count_kmers_lookup(self._lookups[s.index], s.start, s.stop)
        return np.array([np.minimum(kmers_s,
                    count_kmers_lookup(self._lookups[t.index], t.start, t.stop
                        )).sum()
                for t in ts], dtype=float)
