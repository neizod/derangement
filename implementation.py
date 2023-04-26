#!/usr/bin/env python3

import sys
sys.setrecursionlimit(3000)

from functools import lru_cache

from math import factorial
from random import randrange


## helper fn ##################################

def filtering_out(xs, ys):
    rs = set(ys)
    return [x for x in xs if x not in rs]


## counting ###################################

def choose(n, r):
    if r < 0 or n-r < 0:
        return 0
    return factorial(n) // factorial(r) // factorial(n-r)

def periodic_choose(n, r):
    return factorial(n) // factorial(r)**(n//r) // factorial(n//r)

@lru_cache(maxsize=None)
def stirling1st(n, m, r=1, dp={}):
    if (n,m) == (0,0):
        return 1
    if n <= 0 or m <= 0 or n < r*m:
        return 0
    return ( (n-1) * stirling1st(n-1, m, r)
           + choose(n-1, r-1)*factorial(r-1) * stirling1st(n-r, m-1, r) )

def stirling1st_partialtype(n, m, k, v):
    return ( choose(n, k*v)                 # partition k*v from n
           * periodic_choose(k*v, k)        # split into v groups of size k
           * factorial(k-1)**v              # sattolo cyclic on each group k
           * stirling1st(n-k*v, m-v, k+1) )


## unrank algos ###############################

def unrank_choose(xs, k, i):
    ys = []
    for j, x in enumerate(xs):
        if i >= choose(len(xs)-j-1, k-len(ys)-1):
            i -= choose(len(xs)-j-1, k-len(ys)-1)
        else:
            ys += [x]
    return ys

def unrank_periodic_choose(xs, k, i):
    if not xs:
        return []
    n = len(xs)
    i, j = divmod(i, choose(n-1, k-1))
    ys = index_choose(xs, k, j)
    xs = filtering_out(xs, ys)
    return [ys, *unrank_periodic_choose(xs, k, i)]

def unrank_cyclic_permute(xs, i):
    xs = xs[:]
    for k in reversed(range(1, len(xs))):
        i, j = divmod(i, k)
        xs[k], xs[j] = xs[j], xs[k]
    return xs

def bisect_cycle_length(n, m, r, i):
    lo, hi = r, 1+(n//m)
    while lo < hi:
        mid = (lo + hi) // 2
        if i < stirling1st(n, m, mid):
            lo = mid + 1
        else:
            hi = mid
    k = lo - 1
    return k, i - stirling1st(n, m, k+1)

def find_cycle_amount(n, m, k, i):
    for v in range(1, 1+(n//k)):
        if i < stirling1st_partialtype(n, m, k, v):
            return v, i
        i -= stirling1st_partialtype(n, m, k, v)

def reconstruct_type(n, m, r, i):
    if m == 0:
        return []
    k, i = bisect_cycle_length(n, m, r, i)
    v, i = find_cycle_amount(n, m, k, i)
    ix, j = divmod(i, stirling1st(n-k*v, m-v, k+1))
    return [[k,v,ix], *reconstruct_type(n-k*v, m-v, k+1, j)]

def unrank_derangement(vs, m, i):
    n = len(vs)
    xs = [x for x in range(n)]
    ps = [x for x in range(n)]
    for k, v, ix in reconstruct_type(n, m, 2, i):
        ix, ix1 = divmod(ix, choose(len(xs), k*v))
        ys = index_choose(xs, k*v, ix1)
        ix, ix2 = divmod(ix, periodic_choose(len(ys), k))
        for zs in unrank_periodic_choose(ys, k, ix2):
            ix, ix3 = divmod(ix, factorial(k-1))
            for a, b in zip(zs, unrank_cyclic_permute(zs, ix3)):
                ps[a] = b
        xs = filtering_out(xs, ys)
    return [vs[p] for p in ps]


## random gen #################################

def random_derangement(n, m):
    vs = [x for x in range(n)]
    i = randrange(stirling1st(n, m, 2))
    return unrank_derangement(vs, m, i)
