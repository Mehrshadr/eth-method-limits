import random, json, sys
import numpy as np
from borsuk import *
out = open("borsuk_results.jsonl", "w")
def report(name, G, Rs, tl=20):
    k, Q = G.shape
    for R in Rs:
        if R >= Q: continue
        t, opt = min_tau(G, R, tl=tl)
        rec = {"name": name, "k": k, "Q": Q, "R": R, "Q_over_R": round(Q / R, 2), "tau": round(t, 3), "optimal": opt}
        out.write(json.dumps(rec) + "\n"); out.flush()
        print(rec, flush=True)
for k in range(4, 8):
    report(f"cube{k}", cube(k), range(2, k))
for a, b in [(2, 3), (3, 3), (3, 4), (2, 5), (4, 4)]:
    report(f"K{a},{b}-cuts", bipartite_cut_code(a, b), range(2, 9))
rnd = random.Random(7)
for k in (5, 6, 7):
    for Q in (k + 3, k + 5, 2 * k):
        for w in (2, 3):
            best = None
            for trial in range(6):
                G = random_code(k, Q, w, rnd)
                if G is None: continue
                report(f"rand k{k} Q{Q} w{w} #{trial}", G, range(2, min(Q, 7)), tl=10)
