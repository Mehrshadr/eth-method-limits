"""Cut code of K_{a,b} with symmetric coefficients, solved over symmetry classes (works for large a, b).
Shift z = complement of the cut of S; violations = |cut(S)| = sL(b-sR) + sR(a-sL) <= R-1.
Character T (even vertex set) with tL, tR: eigenvalue (-1)^tL * sum g(sL,sR) K_sL(tL;a) K_sR(tR;b)."""
import json, sys
import numpy as np
from math import comb
from scipy.optimize import milp, LinearConstraint, Bounds

def kraw(s, t, n):
    return sum((-1) ** j * comb(t, j) * comb(n - t, s - j) for j in range(0, min(s, t) + 1))

def solve_cut(a, b, R, tl=60):
    S = [(sl, sr) for sl in range(a + 1) for sr in range(b + 1)
         if sl * (b - sr) + sr * (a - sl) <= R - 1 and not (sl == 0 and sr == 0) and not (sl == a and sr == b)]
    # (0,0) and (a,b) are the all-ones shift itself = M0 term with no violation? no: S=empty gives z = cut(L) = all-ones, allowed (0 violations)
    S = [(0, 0)] + S
    T = [(tl_, tr) for tl_ in range(a + 1) for tr in range(b + 1) if (tl_ + tr) % 2 == 0]
    wT = np.array([comb(a, x) * comb(b, y) for x, y in T], dtype=float)
    K = np.array([[(-1) ** x * kraw(sl, x, a) * kraw(sr, y, b) for (sl, sr) in S] for (x, y) in T], dtype=float)
    scale = np.abs(K).max(0); scale[scale == 0] = 1; K = K / scale
    nS, nT = len(S), len(T)
    B = 1e3; M = B * np.abs(K).sum(1).max() + 2
    A = np.zeros((2 * nT, nS + nT)); hi = np.zeros(2 * nT)
    A[:nT, :nS] = K; A[:nT, nS:] = -M * np.eye(nT); hi[:nT] = 1
    A[nT:, :nS] = -K; A[nT:, nS:] = -M * np.eye(nT); hi[nT:] = -1
    c = np.r_[np.zeros(nS), wT / wT.sum()]
    r = milp(c, constraints=LinearConstraint(A, -np.inf, hi), integrality=np.r_[np.zeros(nS), np.ones(nT)],
             bounds=Bounds(np.r_[np.full(nS, -B), np.zeros(nT)], np.r_[np.full(nS, B), np.ones(nT)]),
             options={"time_limit": tl, "disp": False})
    if r.x is None: return None
    g = r.x[:nS]; m = K @ g
    bad = np.abs(m - 1) > 1e-6
    return float(wT[bad].sum() / wT.sum()), bool(r.status == 0)

def cube_tau(k, d):
    p = [comb(k, j) / 2 ** k for j in range(k + 1)]
    return 1 - max(sum(p[s:s + d + 1]) for s in range(k - d + 1)) if d < k else 0.0

if __name__ == "__main__":
    out = open("cutsym_results.jsonl", "w")
    for a, b in [(3, 3), (4, 4), (5, 5), (4, 8), (6, 6), (7, 7), (5, 10), (8, 8), (10, 10), (12, 12), (13, 13)]:
        k, Q = a + b - 1, a * b
        for R in range(2, Q):
            res = solve_cut(a, b, R)
            if res is None: continue
            tau, opt = res
            rec = {"a": a, "b": b, "k": k, "Q": Q, "R": R, "QR": round(Q / R, 3), "tau": round(tau, 5), "opt": opt}
            out.write(json.dumps(rec) + "\n"); out.flush()
            print(rec, flush=True)
            if tau < 0.005: break
