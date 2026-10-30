import sys, json, numpy as np, itertools
from geom import ag_lines, pg2_lines, gf
from families import run
which=sys.argv[1] if __name__=="__main__" else ""
out=open('fam_results.jsonl','a')
def ag_perms(q, pts, mul):
    pid={p:i for i,p in enumerate(pts)}
    gen=2 if q>2 else 1
    fs=[lambda x,y:(x^1,y), lambda x,y:(x,y^1), lambda x,y:(y,x), lambda x,y:(x,y^x),
        lambda x,y:(mul(gen,x),y), lambda x,y:(mul(x,x),mul(y,y))]
    return [[pid[f(*p)] for p in pts] for f in fs]
if which.startswith('ag'):
    q=int(which[2:]); G,pts,mul=ag_lines(q)
    G=np.vstack([np.ones(G.shape[1],int),G])
    run(f'lines AG(2,{q})', G, ag_perms(q,pts,mul), None, out=out, tl=60, budget=560)
if which.startswith('pg'):
    q=int(which[2:]); G,pts=pg2_lines(q)
    run(f'lines PG(2,{q})', G, [], None, out=out, tl=60, budget=560)
