import numpy as np, itertools, json, sys, time
from gsym import Code, solve, cube_tau, cube_QR

def cut_code(V, edges):
    G = np.zeros((V, len(edges)), dtype=np.int64)
    for t, (a, b) in enumerate(edges): G[a, t] = G[b, t] = 1
    return G

def edge_perm(edges, vp):
    idx = {frozenset(e): t for t, e in enumerate(edges)}
    return [idx[frozenset((vp[a], vp[b]))] for a, b in edges]

def torus(a, b, sym=True):
    vid = lambda i, j: (i % a) * b + (j % b)
    edges = [(vid(i, j), vid(i + 1, j)) for i in range(a) for j in range(b)] + \
            [(vid(i, j), vid(i, j + 1)) for i in range(a) for j in range(b)]
    vps = [[vid(i + 2, j) for i in range(a) for j in range(b)], [vid(i, j + 2) for i in range(a) for j in range(b)],
           [vid(i + 1, j + 1) for i in range(a) for j in range(b)],
           [vid(-i, j) for i in range(a) for j in range(b)], [vid(i, -j) for i in range(a) for j in range(b)]]
    if a == b: vps.append([vid(j, i) for i in range(a) for j in range(b)])
    return cut_code(a * b, edges), [edge_perm(edges, p) for p in vps] if sym else []

def hypercube_graph(h):
    V = 1 << h
    edges = [(x, x ^ (1 << i)) for x in range(V) for i in range(h) if not (x >> i) & 1]
    def cperm(p): return [sum(((x >> i) & 1) << p[i] for i in range(h)) for x in range(V)]
    vps = [[x ^ 3 for x in range(V)], cperm(list(range(1, h)) + [0]), cperm([1, 0] + list(range(2, h)))]
    return cut_code(V, edges), [edge_perm(edges, p) for p in vps]

def rm(r, m):
    pts = list(itertools.product([0, 1], repeat=m))
    rows = []
    for deg in range(r + 1):
        for S in itertools.combinations(range(m), deg):
            rows.append([int(all(x[i] for i in S)) for x in pts])
    pid = {x: t for t, x in enumerate(pts)}
    def pp(f): return [pid[tuple(f(list(x)))] for x in pts]
    perms = [pp(lambda x: [x[0] ^ 1] + x[1:]), pp(lambda x: x[1:] + x[:1]), pp(lambda x: [x[1], x[0]] + x[2:]),
             pp(lambda x: [x[0], x[1] ^ x[0]] + x[2:])]
    return np.array(rows), perms

def run(name, G, perms, Rs, tl=60, out=None, budget=600):
    t0 = time.time()
    code = Code(G, perms)
    k, Q = code.k, code.Q
    if Rs is None:
        Rs = sorted({int(w) + 1 for w in np.unique(code.wt) if 0 < w < Q / 2 + 1})
    for R in Rs:
        if time.time() - t0 > budget: break
        res = solve(code, R, tl=30)
        if res is None: continue
        rec = dict(family=name, k=k, Q=Q, R=R, QR=round(Q / R, 3), **res, sec=round(time.time() - t0, 1))
        print(rec, flush=True)
        if out: out.write(json.dumps(rec) + "\n"); out.flush()
        if res["tau"] <= 0.005: break
    return code
