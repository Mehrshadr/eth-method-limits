import itertools, json, time, sys
import numpy as np
from math import comb
from design import solve

def interval_tau(Q, d):
    if d >= Q: return 0.0
    p = [comb(Q, j) / 2**Q for j in range(Q + 1)]
    return 1 - max(sum(p[a:a + d + 1]) for a in range(Q - d + 1))

def kernel_f2(H):
    H = H.copy() % 2; m, Q = H.shape
    piv, r = [], 0
    for c in range(Q):
        p = next((i for i in range(r, m) if H[i, c]), None)
        if p is None: continue
        H[[r, p]] = H[[p, r]]
        for i in range(m):
            if i != r and H[i, c]: H[i] ^= H[r]
        piv.append(c); r += 1
        if r == m: break
    free = [c for c in range(Q) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(Q, dtype=np.int8); v[f] = 1
        for i, c in enumerate(piv):
            if H[i, f]: v[c] = 1
        basis.append(v)
    return np.array(basis, dtype=np.int8)

def duals(Q, m):
    """Parity-check matrices H (m x Q) of rank m, one per column-multiset class up to GL(m)."""
    types = [np.array(v, dtype=np.int8) for v in itertools.product([0, 1], repeat=m)]
    GL = []
    for M in itertools.product([0, 1], repeat=m * m):
        M = np.array(M, dtype=np.int8).reshape(m, m)
        if np.linalg.matrix_rank(M) == m and round(abs(np.linalg.det(M))) % 2 == 1:
            GL.append(M)
    seen = set()
    for counts in itertools.product(range(Q + 1), repeat=len(types)):
        if sum(counts) != Q: continue
        cols = [types[i] for i, c in enumerate(counts) for _ in range(c)]
        H = np.array(cols, dtype=np.int8).T
        if np.linalg.matrix_rank(H) < m: continue
        # rank over F2
        if len(kernel_f2(H.T)) != 0 and False: pass
        keys = []
        for M in GL:
            HM = (M @ H) % 2
            keys.append(tuple(sorted(tuple(c) for c in HM.T)))
        key = min(keys)
        if key in seen: continue
        seen.add(key)
        yield H

out = open("codes_results.jsonl", "w")
for Q in range(5, 8):
    for m in (1, 2):
        for H in duals(Q, m):
            G = kernel_f2(H)
            k = G.shape[0]
            if k != Q - m: continue
            for d in (1, 2, 3):
                if d >= Q: continue
                r = solve(G, d, time_limit=20)
                if r is None: continue
                r["H_cols"] = ["".join(map(str, c)) for c in H.T]
                r["prod_same_k"] = interval_tau(k, d)
                r["prod_same_Q"] = interval_tau(Q, d)
                out.write(json.dumps(r) + "\n"); out.flush()
                flag = "  <-- beats product at same dimension" if r["groups_ok"] and r["tau"] < r["prod_same_k"] - 1e-9 else ""
                print(f"Q={Q} k={k} d={d} H={r['H_cols']} tau={r['tau']:.3f} groups_ok={r['groups_ok']} opt={r['optimal']} | product k={k}: {r['prod_same_k']:.3f}{flag}", flush=True)
