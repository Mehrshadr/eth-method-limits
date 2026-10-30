"""Cut code of a blow-up H[s]: each vertex l of a bipartite 'layer graph' H becomes an independent set of size s_l,
adjacent layers fully joined.  Nights = edges, k = sum s_l - 1, Q = sum_{lm in E(H)} s_l s_m.
g constant on S_{s_1} x ... x S_{s_L} orbits (and optionally on orbits of layer automorphisms).
Shift class a = (|A_l|), violations = |cut(A)| = sum_{lm} a_l (s_m - a_m) + a_m (s_l - a_l).
Character class t = (|T cap V_l|), T even;  m(t) = (-1)^{sum_{l in side0} t_l} sum_a g_a prod_l K_{a_l}(t_l; s_l)."""
import itertools, json, sys, time
import numpy as np
np.seterr(all="ignore")
import warnings; warnings.filterwarnings("ignore")
from math import comb
from scipy.optimize import milp, LinearConstraint, Bounds
from gsym import l1_heur, cube_tau, greedy_set

def kraw(s, t, n):
    return sum((-1) ** j * comb(t, j) * comb(n - t, s - j) for j in range(0, min(s, t) + 1))

def solve_blowup(sizes, Hedges, side0, R, layer_perms=(), tl=60, milp_max=600):
    L = len(sizes)
    cutw = lambda a: sum(a[l] * (sizes[m] - a[m]) + a[m] * (sizes[l] - a[l]) for l, m in Hedges)
    S = [a for a in itertools.product(*[range(s + 1) for s in sizes]) if cutw(a) <= R - 1]
    T = [t for t in itertools.product(*[range(s + 1) for s in sizes]) if sum(t) % 2 == 0]
    def canon(v):
        return min(tuple(v[p[i]] for i in range(L)) for p in [list(range(L))] + list(layer_perms))
    Sc = sorted({canon(a) for a in S}); Sidx = {a: i for i, a in enumerate(Sc)}
    tval = lambda t: (-1) ** sum(t[l] for l in side0)
    # character classes: group T by canonical form (layer perms preserve target, checked below)
    Tc = {}
    for t in T:
        Tc.setdefault(canon(t), []).append(t)
    reps = sorted(Tc)
    wT = np.array([sum(np.prod([comb(sizes[l], t[l]) for l in range(L)]) for t in Tc[r]) for r in reps], dtype=float)
    kr = [np.array([[kraw(x, y, sizes[l]) for y in range(sizes[l] + 1)] for x in range(sizes[l] + 1)], dtype=float) for l in range(L)]
    Sarr = np.array(S); Tarr = np.array(T)
    KF = np.ones((len(T), len(S)))
    for l in range(L): KF *= kr[l][Sarr[:, l]][:, Tarr[:, l]].T          # full chars x full shifts
    sc = np.array([Sidx[canon(a)] for a in S])
    P = np.zeros((len(S), len(Sc))); P[np.arange(len(S)), sc] = 1
    KFc = KF @ P                                                          # full chars x shift classes
    tpos = {t: i for i, t in enumerate(T)}
    K = np.array([KFc[tpos[r]] for r in reps]); tt = np.array([float(tval(r)) for r in reps])
    tvF = np.array([float(tval(t)) for t in T])
    wF = np.array([np.prod([comb(sizes[l], t[l]) for l in range(L)]) for t in T], dtype=float)
    scale = np.abs(K).max(0); scale[scale == 0] = 1; Ks = K / scale
    nS, nT = len(Sc), len(reps)
    cands = []
    g = l1_heur(Ks, tt, wT)
    if g is not None: cands.append((g, False))
    cands.append((greedy_set(Ks, tt, wT), False))
    if nT <= milp_max:
        B = 1e3; M = B * np.abs(Ks).sum(1).max() + 2
        A = np.zeros((2 * nT, nS + nT)); hi = np.zeros(2 * nT)
        A[:nT, :nS] = Ks; A[:nT, nS:] = -M * np.eye(nT); hi[:nT] = tt
        A[nT:, :nS] = -Ks; A[nT:, nS:] = -M * np.eye(nT); hi[nT:] = -tt
        r = milp(np.r_[np.zeros(nS), wT / wT.sum()], constraints=LinearConstraint(A, -np.inf, hi),
                 integrality=np.r_[np.zeros(nS), np.ones(nT)],
                 bounds=Bounds(np.r_[np.full(nS, -B), np.zeros(nT)], np.r_[np.full(nS, B), np.ones(nT)]),
                 options={"time_limit": tl, "disp": False})
        if r.x is not None: cands.append((r.x[:nS], r.status == 0))
    best = None
    for g, opt in cands:
        # re-verify on every (unreduced) character class with exact Krawtchouk products
        gS = g / scale
        m = tvF * (KFc @ gS)
        bad = wF[~(np.abs(m - 1) <= 1e-6)].sum()
        tau = bad / 2 ** (sum(sizes) - 1)
        if best is None or tau < best[0]: best = (tau, bool(opt))
    return best[0], best[1], nS, nT

def cycle_blowup(L, s):
    sizes = [s] * L; Hedges = [(l, (l + 1) % L) for l in range(L)]
    rot = [(i + 1) % L for i in range(L)]; ref = [(-i) % L for i in range(L)]
    perms = []
    for r in range(0, L, 2):
        perms.append([(i + r) % L for i in range(L)]); perms.append([(r - i) % L for i in range(L)])
    return sizes, Hedges, list(range(0, L, 2)), perms

if __name__ == "__main__":
    out = open("blowup_results.jsonl", "a")
    L, s = int(sys.argv[1]), int(sys.argv[2])
    sizes, He, side0, perms = cycle_blowup(L, s)
    k, Q = sum(sizes) - 1, sum(sizes[a] * sizes[b] for a, b in He)
    t0 = time.time()
    cutw = lambda a: sum(a[l] * (sizes[m] - a[m]) + a[m] * (sizes[l] - a[l]) for l, m in He)
    ws = sorted({cutw(a) for a in itertools.product(*[range(x + 1) for x in sizes])})
    for w in ws:
        if w == 0 or w > Q / 2: continue
        if time.time() - t0 > 540: break
        tau, opt, nS, nT = solve_blowup(sizes, He, side0, w + 1, perms, tl=30)
        rec = dict(family=f"cut(C{L}[{s}])", k=k, Q=Q, R=w + 1, QR=round(Q / (w + 1), 3), tau=round(tau, 5), opt=opt, nS=nS, nT=nT)
        print(rec, flush=True); out.write(json.dumps(rec) + "\n"); out.flush()
        if tau <= 0.005: break
