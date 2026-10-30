import itertools, sys
from holo import realizable, inv2
def pmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % p
    return r
def ppow(a, e, p):
    r = [1]
    for _ in range(e): r = pmul(r, a, p)
    return r
def sym_transform(sig, M, p):
    n = len(sig) - 1; out = []
    for w in range(n + 1):
        poly = pmul(ppow([M[1][0], M[1][1]], w, p), ppow([M[0][0], M[0][1]], n - w, p), p)
        out.append(sum(sig[j] * poly[j] for j in range(n + 1)) % p)
    return out
def works(p, k, r=2):
    EQ = [1] + [0] * (r - 1) + [1]; OR = [0] + [1] * k
    for a, b, c, d in itertools.product(range(p), repeat=4):
        if next((v for v in (a, b, c, d) if v), 0) != 1: continue
        M = [[a, b], [c, d]]; Mi = inv2(M, p)
        if Mi is None: continue
        if not realizable(sym_transform(EQ, M, p), p): continue
        if realizable(sym_transform(OR, [[Mi[0][0], Mi[1][0]], [Mi[0][1], Mi[1][1]]], p), p): return True
    return False
if __name__ == "__main__":
  primes = [q for q in range(3, 48) if all(q % d for d in range(2, q))]
  for k in range(3, 13):
        good = [p for p in primes if works(p, k)]
        div = [p for p in primes if (2 ** k - 1) % p == 0]
        print(f"k={k:2d}  works mod: {good}   primes dividing 2^k-1 = {2**k-1}: {div}", flush=True)
