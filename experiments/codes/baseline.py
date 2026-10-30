import numpy as np, time
from math import comb
from design import solve
def interval_tau(Q, d):
    p = [comb(Q, j) / 2**Q for j in range(Q + 1)]
    best = max(sum(p[a:a + d + 1]) for a in range(Q - d + 1)) if d < Q else 1
    return 1 - best
for Q in range(3, 9):
    for d in range(1, 4):
        if d >= Q: continue
        t = time.time()
        r = solve(np.eye(Q, dtype=np.int8), d, time_limit=120)
        print(f"Q={Q} k={Q} d={d}  search tau={r['tau']:.3f} ({r['exceptional']}/{r['of']}) optimal={r['optimal']} groups_ok={r['groups_ok']}  paper-style tau={interval_tau(Q,d):.3f}  {time.time()-t:.1f}s", flush=True)
