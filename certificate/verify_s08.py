from base64 import b64decode
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import gcd, isqrt, lcm
from pathlib import Path
import sys

def dot(x, y):
    return sum(a*b for a, b in zip(x, y))

def primitive(x):
    z = [a*lcm(*(v.denominator for v in x)) for a in x]
    z = [int(a) for a in z]
    return tuple(a//gcd(*z) for a in z)

def chamber_partition(n, signs):
    # Keep all feasible vertex pairs when clipping, including non-edges.
    # Extra section points are harmless and are filtered by active rank.
    facets = [tuple(int(i == j)-int(i+1 == j) for j in range(n)) for i in range(n)]
    vertices = tuple(tuple(F(int(j <= i), i+1) for j in range(n)) for i in range(n))
    cells = [(vertices, tuple(facets))]
    for form in signs:
        if not min(sum(form[:j]) for j in range(1,n+1)) < 0 < max(sum(form[:j]) for j in range(1,n+1)):
            continue
        new = []
        for vertices, rows in cells:
            pairs = [(v, dot(form,v)) for v in vertices]
            pos = [(v,t) for v,t in pairs if t > 0]
            neg = [(v,t) for v,t in pairs if t < 0]
            if not pos or not neg:
                # The omitted half only contributes a face of the retained half.
                new.append((vertices, rows))
                continue
            section = {v for v,t in pairs if t == 0}
            for u,a in pos:
                for v,b in neg:
                    section.add(tuple((a*y-b*x)/(a-b) for x,y in zip(u,v)))
            for side, s in ((pos,1),(neg,-1)):
                h = rows + (tuple(s*a for a in form),)
                candidates = section | {v for v,_ in side}
                ext = [v for v in candidates if rank([r for r in h if dot(r,v) == 0]) == n-1]
                new.append((tuple(ext),h))
        cells = new
    return sorted({tuple(sorted(primitive(v) for v in vertices)) for vertices,_ in cells})

def rank(rows):
    basis = {}
    for row in rows:
        v = list(map(F,row))
        for j in range(len(v)):
            if not v[j]:
                continue
            if j not in basis:
                basis[j] = [x/v[j] for x in v]
                break
            v = [x-v[j]*y for x,y in zip(v,basis[j])]
    return len(basis)

def matchings(n):
    if not n:
        yield ()
        return
    i, *rest = n
    for j in rest:
        for m in matchings(tuple(x for x in rest if x != j)):
            yield ((i,j),)+m

# Field elements are plain ordered pairs (a,b) with value a+b*sqrt(2).
def cmp(z):
    a,b = z
    if not b:
        return (a > 0)-(a < 0)
    if b < 0:
        return -cmp((-a,-b))
    return 1 if a >= 0 or 2*b*b > a*a else -1

def plus(z,w): return z[0]+w[0],z[1]+w[1]
def times(z,w): return z[0]*w[0]+2*z[1]*w[1],z[0]*w[1]+z[1]*w[0]
def scale(t,z): return t*z[0],t*z[1]
def sub(z,w): return plus(z,scale(-1,w))
def fsum(xs):
    out = (0,0)
    for x in xs:
        out = plus(out,x)
    return out

@cache
def upper(p):
    p = tuple(z for z in p if z != (0,0))
    assert all(cmp(z)>0 for z in p)
    if len(p) <= 1:
        return p[0] if p else (0,0)
    if len(p) == 2 and p[0] == p[1]:
        return times((0,1),p[0])
    sq = fsum(times(z,z) for z in p)
    d,e = 10**12,10**30
    t = isqrt(2*e*e)
    rough = d*d*(sq[0]+sq[1]*F(t+(sq[1]>0),e))
    j = isqrt(rough.numerator//rough.denominator)+1
    while cmp(sub((F((j-1)**2,d*d),0),sq)) >= 0:
        j -= 1
    assert cmp(sub((F(j*j,d*d),0),sq)) >= 0
    return F(j,d),0

def ldl(g,b):
    # Rational LDL^T factorization, then two triangular solves in the field.
    n = len(b)
    l,d = [[F(int(i==j)) for j in range(n)] for i in range(n)],[]
    for j in range(n):
        d.append(g[j][j]-sum(d[t]*l[j][t]**2 for t in range(j)))
        assert d[j] > 0
        for i in range(j+1,n):
            l[i][j] = (g[i][j]-sum(d[t]*l[i][t]*l[j][t] for t in range(j)))/d[j]
    w = []
    for i in range(n):
        w.append(sub(b[i],fsum(scale(l[i][j],w[j]) for j in range(i))))
    for i in range(n-1,-1,-1):
        w[i] = sub(scale(1/d[i],w[i]),fsum(scale(l[j][i],w[j]) for j in range(i+1,n)))
    return w

@cache
def moment(x):
    return sum(abs(dot(x,s)) for s in all_signs)/len(all_signs)

def audit(path):
    raw = Path(path).read_text(encoding="utf-8").splitlines()
    tape = iter(b64decode("".join(raw[raw.index("CERTIFICATE_BEGIN")+1:
                                     raw.index("CERTIFICATE_END")]),validate=True))
    cells = chamber_partition(6,[s for s in all_signs if s[0]==1])
    assert len(cells)==21 and len(set().union(*map(set,cells)))==18
    assert all(len(c)==6 and rank(c)==6 for c in cells)
    edges = tuple(combinations(range(6),2))
    visited = leaves = tight = 0
    for c in cells:
        for s in all_signs:
            assert not min(dot(s,v) for v in c) < 0 < max(dot(s,v) for v in c)
        if any(2*v[0]>sum(v) for v in c):
            assert all(2*v[0]>=sum(v) for v in c)
            continue
        for pairing in matchings(tuple(range(6))):
            stack = [tuple(tuple(F(t,sum(v)) for t in v) for v in c)]
            while stack:
                v = stack.pop()
                code = next(tape)
                visited += 1
                if 64 <= code <= 78:
                    i,j = edges[code-64]
                    m = tuple((x+y)/2 for x,y in zip(v[i],v[j]))
                    stack.append(v[:i]+(m,)+v[i+1:])
                    stack.append(v[:j]+(m,)+v[j+1:])
                    continue
                assert 1 <= code <= 63
                leaves += 1
                ix = [i for i in range(6) if code & (1<<i)]
                a = [moment(x) for x in v]
                g = [[dot(x,y)-p*q for y,q in zip(v,a)] for x,p in zip(v,a)]
                b = [sub(upper(tuple((x[i]-x[j],x[j]) for i,j in pairing)),(p,0))
                     for x,p in zip(v,a)]
                w = ldl([[g[i][j] for j in ix] for i in ix],[b[i] for i in ix])
                values = [fsum(scale(row[j],t) for j,t in zip(ix,w)) for row in g]
                assert all(values[i]==b[i] for i in ix)
                assert all(cmp(sub(x,y))>=0 for x,y in zip(values,b))
                norm = fsum(times(t,values[i]) for i,t in zip(ix,w))
                h = cmp(sub((4,F(-8,3)),norm))
                assert cmp(norm)>=0 and h>=0
                if h==0:
                    tight += 1
                    y = [fsum(scale(v[i][j],t) for i,t in zip(ix,w)) for j in range(6)]
                    assert y[0]==y[1]==times(y[2],(F(1,2),F(3,2)))
                    assert y[2]==y[3] and cmp(y[2])>0 and y[4]==y[5]==(0,0)
    assert next(tape,None) is None
    assert (visited,leaves,tight)==(1854,1077,34)
    print("AUDIT PASS",len(cells),visited,leaves,tight)

def run_fixtures():
    from collections import defaultdict

    def div(z,w):
        a,b = w
        return scale(1/F(a*a-2*b*b),times(z,(a,-b)))

    def gram(pairs):
        return fsum(times(z,z) for a,b in pairs
                    for z in [(a-b,b)])

    gamma,c2,k = (F(45,97),F(24,97)),(F(1,2),F(1,4)),(-1,1)
    assert cmp(sub(gamma,(F(4,5),0)))>0 and cmp(sub(c2,gamma))>0
    pairs = [(5,2),(4,2),(2,0),(F(1,2),F(1,4))]
    base_signs = tuple(product((-1,1),repeat=8))
    old = tuple(t for a,b in pairs[:3] for t in (a,b))
    old_support = [dot(old,p[:6]) for p in base_signs]
    assert min(map(abs,old_support))==1 and sum(map(abs,old_support))==6*256
    assert dot(old,old)==53
    B = (5,2)
    assert times(B,B)==gram(pairs[:3])
    r3 = div(sub((53,0),scale(6,B)),sub(B,(6,0)))
    assert cmp(sub(r3,(3,0)))>0
    assert cmp(sub((64,0),gram(pairs)))>0
    c4 = gram(pairs[3:])
    assert cmp(sub(c4,(F(5,16),0)))>0
    b0 = (F(5,2),2)  # 9/2+2k
    assert cmp(sub(scale(2,times(b0,b0)),(49,0)))>0
    assert cmp(sub((5,0),plus((F(9,2),0),c4)))>0
    assert cmp(sub(scale(9,c2),scale(8,times(gamma,gamma))))>0
    # Sharp symmetric boundary and exact single-pair merger obstruction.
    a = (1,F(3,2))
    q = (F(23,2),3)
    b = (5,F(-1,2))
    assert sub(q,times(a,a))==(6,0)
    assert cmp(sub(q,times(times(a,a),plus((1,0),times(k,k)))))>0
    d = sub(b,a)
    assert cmp(sub(times(d,d),scale(6,times(k,k))))>0
    assert times(gamma,(F(99,2),-7))==(F(39,2),9)
    assert sub(scale(2,sub((1,0),gamma)),sub(times(gamma,b),a))==(0,0)

    # Identical-pair fixture where every possible deletion raises the optimum.
    rho = F(33,40)
    b3 = F(21,5)
    assert cmp(sub((b3*b3,0),(9,6)))>0
    assert (1-rho)*(15-rho*b3*b3)-(rho*b3-F(51,16))**2==F(39,32000)
    assert sum(abs(dot((2,1)*3,p[:6])) for p in base_signs)==256*F(51,16)
    assert sum(abs(dot((2,1)*4,p)) for p in base_signs)==256*F(57,16)
    assert cmp(sub(scale(rho,(24,16)),(F(153,4),0)))>0

    # Four identical effective pairs: positive-definite metric and feasible optimizer.
    d = F(367,1024)
    u,v = F(3,2),F(9,16)
    assert d == 1-(u*u+v*v)/4
    rr = scale(1/d,sub((F(1,4),0),scale(v/2,k)))
    aa = scale(F(1,4),sub((2,0),scale(u,rr)))
    bb = scale(F(1,4),sub(scale(2,k),scale(v,rr)))
    assert all(cmp(t)>0 for t in [rr,aa,bb,sub(aa,scale(2,bb))])
    inverse = fsum([rr,scale(2,aa),scale(2,times(k,bb))])
    assert inverse==(F(1919,367),F(-1040,367))
    assert cmp(sub((1,0),times(gamma,inverse)))>0

    intervals = [(0,12)] + [(a,b) for j in range(4) for a,b in
        [(3*j,3*j+3),(3*j,3*j+1),(3*j+1,3*j+2)]]
    histories = [[j for j,(_,b) in enumerate(intervals) if b<=a] for a,_ in intervals]
    w = tuple(w for a,b in pairs for w in ((a-b,0),(b,0),(b,0)))
    hist = [defaultdict(int) for _ in histories]
    moment_sum = (0,0)
    for ps in product((-1,1),repeat=9):
        children = tuple(t for j in range(4) for t in (ps[2*j],ps[2*j],ps[2*j+1]))
        s = fsum(scale(t,wj) for t,wj in zip(children,w))
        z = cmp(s) or ps[-1]
        signs = (z,)+children
        for i,h in enumerate(histories):
            hist[i][tuple(signs[j] for j in h)] += signs[i]
        v = plus(s,scale(z,r3))
        moment_sum = plus(moment_sum,times(v,v))
    assert sum(map(len,hist))==341 and all(t==0 for h in hist for t in h.values())
    target = fsum([times(r3,r3),scale(12,r3),(F(853,16),0)])
    assert scale(F(1,512),moment_sum)==target
    print("FIXTURES PASS: strict deletion failures, merger obstruction, identical-pair optimum, "
          "341 conditional cells over 512 draws")

all_signs = tuple(product((-1,1),repeat=6))
audit(sys.argv[1])
run_fixtures()
if "--geometry" in sys.argv[2:]:
    signs = [(1,)+t for t in product((-1,1),repeat=7)]
    cells = chamber_partition(8, signs)
    totals = len(cells),len(set().union(*map(set,cells))),max(map(len,cells))
    assert totals==(2470,206,37)
    print("GEOMETRY PASS",totals)
