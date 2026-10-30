"""Generic symmetry-reduced solver for the code version (k up to ~24, Q up to 64).
g is taken constant on orbits (under given coordinate permutations) of the low-weight codewords c
(z = 1 + c, wt(c) <= R-1).  Characters y are grouped by their exact value vector (all orbit sums + target),
which is exact for the reduced LP.  Every result is re-verified on all 2^k characters with the full g."""
import numpy as np, time, warnings
np.seterr(all='ignore'); warnings.filterwarnings('ignore')
from scipy.optimize import milp, LinearConstraint, Bounds
from math import comb

def gf2_rank(G):
    G = G.copy() % 2; r = 0; k, Q = G.shape
    for c in range(Q):
        p = next((i for i in range(r, k) if G[i, c]), None)
        if p is None: continue
        G[[r, p]] = G[[p, r]]
        for i in range(k):
            if i != r and G[i, c]: G[i] ^= G[r]
        r += 1
    return r

def basis(G):
    """row-reduce to an independent basis that contains the all-ones word as row 0."""
    G = np.array(G, dtype=np.int64) % 2
    Q = G.shape[1]
    rows = [np.ones(Q, dtype=np.int64)]
    for g in G:
        if gf2_rank(np.array(rows + [g])) > len(rows): rows.append(g)
    return np.array(rows)

class Code:
    def __init__(self, G, perms=()):
        self.G = basis(G); self.k, self.Q = self.G.shape
        assert self.Q <= 64
        k, Q = self.k, self.Q
        rowmask = [int(sum(int(b) << t for t, b in enumerate(r))) for r in self.G]
        cw = np.zeros(1, dtype=np.uint64)
        for m in rowmask:
            cw = np.concatenate([cw, cw ^ np.uint64(m)])
        self.cw = cw                                  # codeword of message u (bit i of u = row i)
        self.wt = np.bitwise_count(cw).astype(np.int64)
        self.u1 = 1                                   # row 0 = all-ones
        self.perms = [np.array(p) for p in perms]
        # decoding: info set
        I = []; M = np.zeros((0, k), dtype=np.int64)
        for t in range(Q):
            if gf2_rank(np.vstack([M, self.G[:, t]])) > len(I):
                I.append(t); M = np.vstack([M, self.G[:, t]])
        self.I = I
        A = self.G[:, I].T % 2                         # A u = c[I]
        self.Ainv = self._inv(A)

    def _inv(self, A):
        n = A.shape[0]; M = np.hstack([A % 2, np.eye(n, dtype=np.int64)])
        r = 0
        for c in range(n):
            p = next(i for i in range(r, n) if M[i, c]); M[[r, p]] = M[[p, r]]
            for i in range(n):
                if i != r and M[i, c]: M[i] ^= M[r]
            r += 1
        return M[:, n:]

    def decode(self, masks):
        v = np.zeros(len(masks), dtype=np.int64)
        for j, t in enumerate(self.I):
            v |= (((masks >> np.uint64(t)) & np.uint64(1)).astype(np.int64)) << j
        u = np.zeros(len(masks), dtype=np.int64)
        for i in range(self.k):                        # u_i = sum_j Ainv[i,j] v_j
            colmask = int(sum(int(self.Ainv[i, j]) << j for j in range(self.k)))
            u |= (np.bitwise_count(v & colmask) & 1).astype(np.int64) << i
        return u

    def permute(self, masks, p):
        out = np.zeros_like(masks)
        for t in range(self.Q):
            out |= ((masks >> np.uint64(t)) & np.uint64(1)) << np.uint64(p[t])
        return out

    def orbits(self, B):
        """B: sorted array of message indices closed under the perms. returns label per element."""
        lab = np.arange(len(B))
        imgs = []
        for p in self.perms:
            u2 = self.decode(self.permute(self.cw[B], p))
            idx = np.searchsorted(B, u2); assert np.all(B[idx] == u2), "B not closed under perm"
            imgs.append(idx)
        while True:
            old = lab.copy()
            for idx in imgs:
                lab = np.minimum(lab, lab[idx])
            while True:
                l2 = lab[lab]
                if np.array_equal(l2, lab): break
                lab = l2
            if np.array_equal(lab, old): break
        _, inv = np.unique(lab, return_inverse=True)
        return inv

def wht(a):
    a = a.astype(np.float64).copy(); n = a.size; h = 1
    while h < n:
        a = a.reshape(-1, 2, h); a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1); h *= 2
    return a

def solve(code, R, tl=60, B_g=1e3, milp_max=600):
    k, N = code.k, 1 << code.k
    d = R - 1
    B = np.nonzero(code.wt <= d)[0]
    inv = code.orbits(B) if code.perms else np.arange(len(B))
    nO = inv.max() + 1
    par = (np.bitwise_count(np.arange(N, dtype=np.int64) & code.u1) & 1)
    target = 1.0 - 2 * par                            # chi_y(u1)
    rng = np.random.default_rng(1)
    h = (par.astype(np.int64) + 1) * int(rng.integers(1, 1 << 40))
    mults = rng.integers(1, 1 << 40, size=nO)
    for o in range(nO):
        ind = np.zeros(N); ind[B[inv == o]] = 1
        h = h + np.rint(wht(ind)).astype(np.int64) * int(mults[o])
    _, rep, cnt = np.unique(h, return_index=True, return_counts=True)
    del h
    K = np.zeros((len(rep), nO))
    for o in range(nO):
        ind = np.zeros(N); ind[B[inv == o]] = 1
        K[:, o] = wht(ind)[rep]
    tt = target[rep]
    nT = len(rep)
    scale = np.abs(K).max(0); scale[scale == 0] = 1; Ks = K / scale
    cands = []
    gl = l1_heur(Ks, tt, cnt)
    if gl is not None: cands.append((gl / scale, False))
    cands.append((greedy_set(Ks, tt, cnt) / scale, False))
    if nT <= milp_max:
        M = B_g * np.abs(Ks).sum(1).max() + 2
        A = np.zeros((2 * nT, nO + nT)); hi = np.zeros(2 * nT)
        A[:nT, :nO] = Ks; A[:nT, nO:] = -M * np.eye(nT); hi[:nT] = tt
        A[nT:, :nO] = -Ks; A[nT:, nO:] = -M * np.eye(nT); hi[nT:] = -tt
        c = np.r_[np.zeros(nO), cnt / N]
        r = milp(c, constraints=LinearConstraint(A, -np.inf, hi), integrality=np.r_[np.zeros(nO), np.ones(nT)],
                 bounds=Bounds(np.r_[np.full(nO, -B_g), np.zeros(nT)], np.r_[np.full(nO, B_g), np.ones(nT)]),
                 options={"time_limit": tl, "disp": False})
        if r.x is not None: cands.append((r.x[:nO] / scale, bool(r.status == 0)))
    best = None
    for gO, opt in cands:
        # exact re-verification over all characters with the full coefficient vector
        gfull = np.zeros(N); gfull[B] = gO[inv]
        m = target * wht(gfull)                       # m(y) = chi_y(u1) * sum_c g_c chi_y(c)
        tau = float(np.mean(~(np.abs(m - 1) <= 1e-6)))
        if best is None or tau < best[0] or (tau == best[0] and opt): best = (tau, opt)
    return dict(tau=best[0], opt=best[1], nO=int(nO), nT=int(nT), nB=int(len(B)))

def l1_heur(Ks, tt, cnt, iters=12):
    from scipy.optimize import linprog
    nT, nO = Ks.shape
    w = cnt / cnt.sum(); bestg, bestbad = None, None
    for _ in range(iters):
        # min sum w s  s.t. -s <= K g - t <= s
        A = np.block([[Ks, -np.eye(nT)], [-Ks, -np.eye(nT)]]); b = np.r_[tt, -tt]
        r = linprog(np.r_[np.zeros(nO), w], A_ub=A, b_ub=b, bounds=[(-1e4, 1e4)] * nO + [(0, None)] * nT, method="highs")
        if r.x is None: break
        g = r.x[:nO]; res = np.abs(Ks @ g - tt)
        bad = cnt[res > 1e-7].sum()
        if bestbad is None or bad < bestbad: bestg, bestbad = g, bad
        w = cnt / cnt.sum() / (res + 1e-3)
    return bestg

def cube_tau(k, d):
    p = [comb(k, j) / 2 ** k for j in range(k + 1)]
    return 1 - max(sum(p[s:s + d + 1]) for s in range(k - d + 1)) if d < k else 0.0

def cube_QR(k, thr):
    for d in range(k + 1):
        if cube_tau(k, d) <= thr + 1e-12: return k / (d + 1)

def greedy_set(Ks, tt, cnt, restarts=8, seed=0):
    """Greedy: accept character types (heaviest first, randomly perturbed) while K_S g = t_S stays consistent.
    Incremental Gram-Schmidt carrying the target column along."""
    rng = np.random.default_rng(seed)
    nT, nO = Ks.shape
    rn = np.linalg.norm(Ks, axis=1); rn[rn == 0] = 1
    Kn = Ks / rn[:, None]; tn = tt / rn
    best = None
    for rs in range(restarts):
        pert = cnt * (1 if rs == 0 else np.exp(rng.normal(0, 0.7, nT)))
        order = np.argsort(-pert)
        U = np.zeros((0, nO)); Sv = np.zeros(0); acc = []
        for i in order:
            kv, tv = Kn[i], tn[i]
            c = U @ kv
            rk = kv - c @ U; rt = tv - c @ Sv
            nr = np.linalg.norm(rk)
            if nr < 1e-9:
                if abs(rt) < 1e-9: acc.append(i)
                continue
            U = np.vstack([U, rk / nr]); Sv = np.r_[Sv, rt / nr]; acc.append(i)
        acc = np.array(acc)
        g = np.linalg.lstsq(Ks[acc], tt[acc], rcond=None)[0]
        bad = cnt[~(np.abs(Ks @ g - tt) <= 1e-7)].sum()
        if best is None or bad < best[0]: best = (bad, g)
    return best[1]
