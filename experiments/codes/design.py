"""Search for one-orbit designs M0 = sum_S g_S * T_{E+e_S} on F_2^k.

Code C = row space of G (k x Q). Codeword c <-> character y.
Eigenvalue of M0 on c: m(c) = sum_S g_S * (-1)^{|c & ~S|}.
Measurements on the orbit = #{c : m(c) != 1}  (tau = that / 2^k).
Group condition (a''): for every point t with some used S not containing t,
there is a codeword c_t with (-1)^{|c_t & ~S|} = -1 for every used S not containing t.
Degree of replacement = max |S| = d  (R = d + 1).
"""
import itertools, sys, json, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix

def codewords(G):
    k, Q = G.shape
    out = []
    for y in itertools.product([0, 1], repeat=k):
        out.append(tuple((np.array(y) @ G) % 2))
    return np.array(out, dtype=np.int8)

def sets_upto(Q, d):
    S = []
    for s in range(d + 1):
        S += [frozenset(c) for c in itertools.combinations(range(Q), s)]
    return S

def solve(G, d, B=200.0, time_limit=60):
    k, Q = G.shape
    C = codewords(G)
    nC = len(C)
    Ss = sets_upto(Q, d)
    nS = len(Ss)
    comp = np.ones((nS, Q), dtype=np.int8)
    for i, S in enumerate(Ss):
        for t in S:
            comp[i, t] = 0
    A = ((-1.0) ** ((C @ comp.T) % 2))           # nC x nS
    # variable layout: g[nS] | w[nS] | z[nC] | x[Q*nC]
    ng, nw, nz, nx = nS, nS, nC, Q * nC
    N = ng + nw + nz + nx
    G0, W0, Z0, X0 = 0, ng, ng + nw, ng + nw + nz
    rows, lo, hi = [], [], []
    M = B * nS + 1
    def row():
        r = lil_matrix((1, N)); rows.append(r); return r
    for i in range(nS):                               # |g| <= B w
        r = row(); r[0, G0 + i] = 1; r[0, W0 + i] = -B; lo.append(-np.inf); hi.append(0)
        r = row(); r[0, G0 + i] = -1; r[0, W0 + i] = -B; lo.append(-np.inf); hi.append(0)
    for c in range(nC):                               # |m(c)-1| <= M z
        r = row(); r[0, G0:G0 + ng] = A[c]; r[0, Z0 + c] = -M; lo.append(-np.inf); hi.append(1)
        r = row(); r[0, G0:G0 + ng] = -A[c]; r[0, Z0 + c] = -M; lo.append(-np.inf); hi.append(-1)
    for t in range(Q):                                # one witness c_t per point
        r = row()
        for c in range(nC):
            r[0, X0 + t * nC + c] = 1
        lo.append(1); hi.append(1)
        avoid = [i for i, S in enumerate(Ss) if t not in S]
        for c in range(nC):                           # used S avoiding t must have A[c,S] = -1
            bad = [i for i in avoid if A[c, i] > 0]
            if not bad:
                continue
            r = row()
            for i in bad:
                r[0, W0 + i] = 1
            r[0, X0 + t * nC + c] = len(bad)
            lo.append(-np.inf); hi.append(len(bad))
    from scipy.sparse import vstack
    Acon = vstack([csr_matrix(r) for r in rows]).tocsr()
    cobj = np.zeros(N); cobj[Z0:Z0 + nz] = 1
    integrality = np.zeros(N); integrality[W0:] = 1
    lb = np.zeros(N); ub = np.ones(N); lb[G0:G0 + ng] = -B; ub[G0:G0 + ng] = B
    res = milp(cobj, constraints=LinearConstraint(Acon, lo, hi), integrality=integrality,
               bounds=Bounds(lb, ub), options={"time_limit": time_limit, "disp": False})
    if res.x is None:
        return None
    g = res.x[G0:G0 + ng]
    w = res.x[W0:W0 + nw] > 0.5
    g = np.where(w, g, 0.0)
    # exact re-check
    m = A @ g
    exc = int(np.sum(np.abs(m - 1) > 1e-6))
    used = [i for i in range(nS) if w[i] and abs(g[i]) > 1e-9]
    ok_groups = True
    for t in range(Q):
        av = [i for i in used if t not in Ss[i]]
        if av and not any(all(A[c, i] < 0 for i in av) for c in range(nC)):
            ok_groups = False
    return {"Q": Q, "k": k, "d": d, "exceptional": exc, "of": nC, "tau": exc / nC,
            "groups_ok": ok_groups, "status": res.message, "optimal": bool(res.status == 0),
            "lower_bound": float(getattr(res, "mip_dual_bound", np.nan) or np.nan),
            "used_sets": [sorted(Ss[i]) for i in used]}
