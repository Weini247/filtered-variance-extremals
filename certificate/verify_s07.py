from base64 import b64decode
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from typing import NamedTuple
from math import gcd
import sys

class K(NamedTuple):
    a: F = F(0)
    b: F = F(0)

    def __add__(x, y):
        return K(x.a + y.a, x.b + y.b)

    def __sub__(x, y):
        return K(x.a - y.a, x.b - y.b)

    def __mul__(x, y):
        return K(x.a*y.a + 2*x.b*y.b, x.a*y.b + x.b*y.a)

    def scale(x, t):
        return K(x.a*t, x.b*t)

    def sign(x):
        a, b = x
        if a >= 0 and b >= 0:
            return int(bool(a or b))
        if a <= 0 and b <= 0:
            return -1
        return 1 if (a > 0) == (a*a > 2*b*b) else -1

def dot(x, y):
    return sum(a*b for a, b in zip(x, y))

def det(rows):
    if len(rows) == 1:
        return rows[0][0]
    return sum((-1)**i * a * det(tuple(r[:i]+r[i+1:] for r in rows[1:]))
               for i, a in enumerate(rows[0]) if a)

def match(indices):
    if not indices:
        yield ()
        return
    i, *rest = indices
    for j in rest:
        for tail in match(tuple(t for t in rest if t != j)):
            yield ((i, j),) + tail

@lru_cache(None)
def A(v):
    return sum(abs(dot(s, v)) for s in signs) / F(32)

@lru_cache(None)
def upper(p):
    nonzero = [x for x in p if x.sign()]
    if len(nonzero) <= 1:
        return nonzero[0] if nonzero else K()
    if len(nonzero) == 2 and nonzero[0] == nonzero[1]:
        return K(0, 1) * nonzero[0]
    target = sum((x*x for x in p), K())
    d = 10**12
    lo, hi = 0, d
    assert (K(1) - target).sign() >= 0
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if (K(F(mid*mid, d*d)) - target).sign() >= 0:
            hi = mid
        else:
            lo = mid
    return K(F(hi, d))

def solve(matrix, right):
    g = [list(row) for row in matrix]
    b = list(right)
    for i in range(len(b)):
        d = g[i][i]
        assert d > 0
        g[i] = [v/d for v in g[i]]
        b[i] = b[i].scale(1/d)
        for j in range(len(b)):
            if j != i and g[j][i]:
                t = g[j][i]
                g[j] = [x-t*y for x, y in zip(g[j], g[i])]
                b[j] = b[j] - b[i].scale(t)
    return b

def verify(vs, pairing, depth=0):
    global visited, leaves, tight, maxdepth
    visited += 1
    maxdepth = max(maxdepth, depth)
    code = next(tape)
    assert 1 <= code <= 78
    if code >= 64:
        i, j = edges[code-64]
        mid = tuple((a+b)/2 for a, b in zip(vs[i], vs[j]))
        for t in (j, i):
            verify(vs[:t] + (mid,) + vs[t+1:], pairing, depth+1)
        return
    leaves += 1
    act = [i for i in range(6) if code >> i & 1]
    moments = [A(v) for v in vs]
    g = [[dot(x,y) - a*b for y,b in zip(vs,moments)]
         for x,a in zip(vs,moments)]
    b = [upper(tuple(K(v[i]-v[j], v[j]) for i,j in pairing)) - K(a)
         for v,a in zip(vs,moments)]
    w = solve([[g[i][j] for j in act] for i in act], [b[i] for i in act])
    for row, target in zip(g, b):
        assert (sum((x.scale(row[j]) for x,j in zip(w,act)), K())
                - target).sign() >= 0
    norm = sum((x*b[j] for x,j in zip(w,act)), K())
    status = (K(4,F(-8,3)) - norm).sign()
    assert norm.sign() >= 0 and status >= 0
    if status == 0:
        tight += 1
        y = [sum((x.scale(vs[i][j]) for x,i in zip(w,act)), K())
             for j in range(6)]
        assert y[0] == y[1] == y[2]*K(F(1,2), F(3,2))
        assert y[2] == y[3] and y[2].sign() > 0 and y[4] == y[5] == K()

signs = [(1,) + tuple(1 if mask >> j & 1 else -1 for j in range(5))
         for mask in range(32)]
prefix = [tuple(sum(s[:i]) for i in range(1,7)) for s in signs]
walls = [h for h in prefix if min(h) < 0 < max(h)]
normals = [tuple(int(i == j) for i in range(6)) for j in range(6)] + walls
rays = set()
for rows in combinations(normals, 5):
    d = tuple((-1)**i * det(tuple(row[:i]+row[i+1:] for row in rows))
              for i in range(6))
    if max(d) <= 0:
        d = tuple(-x for x in d)
    if min(d) >= 0 and max(d) > 0:
        x = tuple(sum(d[i:]) for i in range(6))
        # Normalize only after sorting primitive integer rays.
        rays.add(tuple(t//gcd(*x) for t in x))
rays = sorted(rays)
incidences = [tuple(dot(s,x) for s,p in zip(signs,prefix) if p in walls)
              for x in rays]
cells = {tuple(i for i,h in enumerate(incidences)
               if all(z*t >= 0 for z,t in zip(h, bits)))
         for bits in product((-1,1), repeat=len(walls))}
cells.discard(())
cells = sorted(c for c in cells if not any(set(c) < set(d) for d in cells))
assert (len(walls),len(rays),len(cells)) == (12,18,21)
edges = tuple(combinations(range(6),2))
matchings = tuple(match(tuple(range(6))))
assert len(matchings) == 15

report = Path(sys.argv[1]).read_text(encoding="utf-8")
data = report.split("CERTIFICATE_BEGIN\n")[1].split("\nCERTIFICATE_END")[0]
tape = iter(b64decode("".join(data.split()), validate=True))
visited = leaves = tight = maxdepth = dominant = rootcases = 0
for cell in cells:
    assert len(cell) == 6 and det(tuple(rays[i] for i in cell)) != 0
    xs = [rays[i] for i in cell]
    if any(x[0] > sum(x[1:]) for x in xs):
        assert all(x[0] >= sum(x[1:]) for x in xs)
        dominant += 1
        continue
    vs = tuple(tuple(F(x, sum(v)) for x in v) for v in xs)
    for m in matchings:
        rootcases += 1
        verify(vs, m)
assert next(tape, None) is None
assert (dominant,rootcases,visited,leaves,tight,maxdepth) == (1,300,1854,1077,34,6)
print("PASS: 21 chambers; 15 matchings; 300 subdivided cases; "
      "1,854 nodes; 1,077 exact certificates (34 tight); depth <= 6")
