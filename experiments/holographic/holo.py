"""Search for holographic 'translations' (Valiant): a 2x2 basis M over F_p such that
   variable side  EQ_r = [1,0,...,0,1]  (a variable used r times)  -> M^{(x)r} EQ_r
   clause side    OR_k = [0,1,1,...,1]  (monotone k-clause)       -> (M^{-T})^{(x)k} OR_k
both become planar-matchgate signatures. Then counting solutions mod p of planar,
read-r, monotone k-CNF is polynomial (FKT / Pfaffian).
Symmetric signature [x0..xn] is matchgate-realizable iff (a) odd-index entries all 0 or
even-index entries all 0, and (b) the remaining entries form a geometric progression."""
import itertools, sys
from math import comb

def sym_transform(sig, M, p):
    """sig symmetric of arity n given by weights; apply M to every input; return symmetric weights."""
    n = len(sig) - 1
    out = []
    for w in range(n + 1):                      # output weight w: first w outputs 1, rest 0
        b = [1] * w + [0] * (n - w)
        tot = 0
        for x in itertools.product((0, 1), repeat=n):
            prod = sig[sum(x)]
            if prod == 0: continue
            for bi, xi in zip(b, x):
                prod = prod * M[bi][xi] % p
                if prod == 0: break
            tot = (tot + prod) % p
        out.append(tot)
    return out

def realizable(s, p):
    if all(v == 0 for v in s): return False
    n = len(s) - 1
    odd = [s[i] for i in range(1, n + 1, 2)]; even = [s[i] for i in range(0, n + 1, 2)]
    seq = even if all(v == 0 for v in odd) else odd if all(v == 0 for v in even) else None
    if seq is None: return False
    for i in range(len(seq) - 2):               # geometric: s_i * s_{i+2} == s_{i+1}^2
        if (seq[i] * seq[i + 2] - seq[i + 1] ** 2) % p: return False
    return True

def inv2(M, p):
    a, b = M[0]; c, d = M[1]
    det = (a * d - b * c) % p
    if det == 0: return None
    di = pow(det, p - 2, p)
    return [[d * di % p, -b * di % p], [-c * di % p, a * di % p]]

def search(p, k, r):
    EQ = [1] + [0] * (r - 1) + [1]
    OR = [0] + [1] * k
    hits = []
    for a, b, c, d in itertools.product(range(p), repeat=4):
        first = next((v for v in (a, b, c, d) if v), 0)
        if first != 1: continue          # scaling M doesn't change realizability
        M = [[a, b], [c, d]]
        Mi = inv2(M, p)
        if Mi is None: continue
        MiT = [[Mi[0][0], Mi[1][0]], [Mi[0][1], Mi[1][1]]]
        g = sym_transform(EQ, M, p)
        if not realizable(g, p): continue
        rr = sym_transform(OR, MiT, p)
        if realizable(rr, p):
            hits.append((M, g, rr))
    return hits

if __name__ == "__main__":
    primes = [q for q in range(3, int(sys.argv[1]) + 1) if all(q % d for d in range(2, q))]
    ks = [int(x) for x in sys.argv[2].split(",")]
    rs = [int(x) for x in sys.argv[3].split(",")]
    for k in ks:
        for r in rs:
            found = []
            for p in primes:
                h = search(p, k, r)
                if h: found.append((p, len(h), h[0]))
            print(f"k={k} (clause size)  r={r} (uses per variable): moduli that work ->",
                  [f[0] for f in found] or "none", flush=True)
            for p, n, (M, g, rr) in found[:6]:
                print(f"    p={p}: {n} bases, e.g. M={M}  variable sig={g}  clause sig={rr}")
