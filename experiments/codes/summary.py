import json, glob, ast, collections
from gsym import cube_QR
recs=collections.defaultdict(list)
def add(r):
    recs[(r['family'],r['k'],r['Q'])].append((r['R'],float(r['tau']),bool(r['opt'])))
for f in glob.glob('log_*.txt')+glob.glob('old/log_rm*.txt'):
    for l in open(f):
        if l.startswith("{'family'"):
            add(ast.literal_eval(l.replace('np.float64(','(')))
for l in open('cutsym_results.jsonl'):
    r=json.loads(l); r['family']=f"cut(K{r['a']},{r['b']})"; add(r)
print(f"{'family':18s} {'k':>3} {'Q':>4} | thr  | {'R':>3} {'Q/R':>6} {'tau':>7} opt | cube Q/R | ratio")
for key in sorted(recs, key=lambda x:(x[0][:5],x[1])):
    fam,k,Q=key; L=sorted(set(recs[key]))
    for thr in (0.05,0.01):
        best=None; cur=(9,False)
        for R,t,o in L:
            if t<cur[0]: cur=(t,o)
            if cur[0]<=thr: best=(R,cur); break
        cq=cube_QR(k,thr)
        if best: print(f"{fam:18s} {k:3d} {Q:4d} | {thr:.2f} | {best[0]:3d} {Q/best[0]:6.3f} {best[1][0]:7.4f} {'Y' if best[1][1] else 'n'}   | {cq:6.3f}   | {Q/best[0]/cq:5.2f}")
        else: print(f"{fam:18s} {k:3d} {Q:4d} | {thr:.2f} | not reached (max R tried {L[-1][0]}, best tau {min(t for _,t,_ in L):.3f}) | {cq:6.3f}")
