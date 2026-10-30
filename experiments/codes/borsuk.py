"""Construction family: hybrids = codewords of a binary code C (all-ones word in C).
Night t: hybrid x sits in room A if x_t = 0, B if x_t = 1 (rooms always separate automatically).
M0 = sum_z g_z * (shift by z), z in C with weight >= Q-R+1 ("near-opposite"), violations = zeros of z.
Measurements = #characters y of C with eigenvalue m(y) = sum_z g_z (-1)^<y,z> != 1.
Goal: many nights Q per dimension k, small R, few measurements (tau)."""
import itertools, json, random, sys, time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from math import comb

def span(G):
    k = G.shape[0]
    U = np.array(list(itertools.product([0, 1], repeat=k)), dtype=np.int64)
    return U, (U @ G) % 2

def min_tau(G, R, B=60.0, tl=20):
    k, Q = G.shape
    U, W = span(G)
    wt = W.sum(1)
    Z = [i for i in range(len(U)) if wt[i] >= Q - R + 1]
    if not Z: return 1.0, True
    chi = (-1.0) ** ((U @ U[Z].T) % 2)              # chars y (rows) x shifts z (cols)
    nY, nZ = chi.shape
    M = B * nZ + 1
    A = np.zeros((2 * nY, nZ + nY)); lo = np.full(2 * nY, -np.inf); hi = np.zeros(2 * nY)
    A[:nY, :nZ] = chi; A[:nY, nZ:] = -M * np.eye(nY); hi[:nY] = 1
    A[nY:, :nZ] = -chi; A[nY:, nZ:] = -M * np.eye(nY); hi[nY:] = -1
    c = np.r_[np.zeros(nZ), np.ones(nY)]
    integ = np.r_[np.zeros(nZ), np.ones(nY)]
    lb = np.r_[np.full(nZ, -B), np.zeros(nY)]; ub = np.r_[np.full(nZ, B), np.ones(nY)]
    r = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(lb, ub),
             options={"time_limit": tl, "disp": False})
    g = r.x[:nZ]
    exc = int(np.sum(np.abs(chi @ g - 1) > 1e-6))
    return exc / nY, r.status == 0

def cube(k):
    return np.eye(k, dtype=np.int64)

def rank2(G):
    G = G.copy() % 2; r = 0; k, Q = G.shape
    for c in range(Q):
        p = next((i for i in range(r, k) if G[i, c]), None)
        if p is None: continue
        G[[r, p]] = G[[p, r]]
        for i in range(k):
            if i != r and G[i, c]: G[i] ^= G[r]
        r += 1
    return r

def random_code(k, Q, w, rnd):
    """all-ones row + k-1 sparse rows of weight w; columns distinct and nonzero (no wasted nights)."""
    for _ in range(2000):
        rows = [np.ones(Q, dtype=np.int64)]
        for _ in range(k - 1):
            v = np.zeros(Q, dtype=np.int64); v[rnd.sample(range(Q), w)] = 1; rows.append(v)
        G = np.array(rows)
        cols = {tuple(G[:, t]) for t in range(Q)}
        if len(cols) == Q and rank2(G) == k:
            return G
    return None

def bipartite_cut_code(a, b):
    """nights = edges of K_{a,b}; codewords = cuts; all-ones = the cut between the two sides."""
    edges = [(i, a + j) for i in range(a) for j in range(b)]
    V = a + b
    rows = []
    for v in range(V - 1):
        rows.append([1 if v in e else 0 for e in edges])
    return np.array(rows, dtype=np.int64)
