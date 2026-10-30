"""Information-set reduction (upper bound on Q/R for a code).
If I is an information set of C, restriction to I is a group isomorphism C -> F_2^k sending 1 -> all-ones, and
wt_I(c) <= wt(c).  So any solution for C at degree d is a solution for the k-cube at degree d_I(d) = max_{wt(c)<=d} wt_I(c).
Hence tau_C(d) >= tau*_cube(k, d_I(d)).  (tau*_cube = interval value, MILP-verified exact for k <= 8.)
Also the pull-back to F_2^Q gives tau_C(d) >= tau*_cube(Q, d)."""
import numpy as np, random
from gsym import gf2_rank, cube_tau

def is_info(code, I):
    return gf2_rank(code.G[:, I]) == code.k

def dI_profile(code, I, ds):
    mask = np.uint64(sum(1 << t for t in I))
    wI = np.bitwise_count(code.cw & mask).astype(np.int64)
    return {d: int(wI[code.wt <= d].max()) for d in ds}

def best_info_set(code, d, iters=300, seed=0):
    """local search for an information set minimising max_{wt(c)<=d} wt_I(c)."""
    rnd = random.Random(seed); Q, k = code.Q, code.k
    B = code.cw[code.wt <= d]
    def score(I):
        mask = np.uint64(sum(1 << t for t in I))
        w = np.bitwise_count(B & mask)
        return (int(w.max()), int((w == w.max()).sum()))
    best = None
    for start in range(5):
        order = list(range(Q)); rnd.shuffle(order); I = []
        for t in order:
            if gf2_rank(code.G[:, I + [t]]) > len(I): I.append(t)
        s = score(I)
        for _ in range(iters):
            i = rnd.randrange(k); t = rnd.randrange(Q)
            if t in I: continue
            J = I[:i] + [t] + I[i + 1:]
            if not is_info(code, J): continue
            sj = score(J)
            if sj <= s: I, s = J, sj
        if best is None or s < best[0]: best = (s, I)
    return best[1], best[0][0]

def ub_QR(code, thr, ds=None):
    """largest Q/R compatible with tau <= thr under the info-set bound and the Q-cube pull-back."""
    Q, k = code.Q, code.k
    ws = sorted(set(int(w) for w in np.unique(code.wt) if 0 < w < Q))
    for d in ws:
        I, dI = best_info_set(code, d, iters=150)
        if cube_tau(k, dI) <= thr + 1e-12 and cube_tau(Q, d) <= thr + 1e-12:
            return Q / (d + 1), d, dI
    return None
