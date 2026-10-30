import numpy as np, itertools
from gsym import *
from borsuk import bipartite_cut_code
k=8; cyc=list(range(1,k))+[0]; sw=[1,0]+list(range(2,k))
c=Code(np.eye(k,dtype=int),[cyc,sw])
for R in range(2,8):
    print('cube',k,R,solve(c,R), round(cube_tau(k,R-1),4))
G=bipartite_cut_code(4,4)
edges=[(i,4+j) for i in range(4) for j in range(4)]
def vperm(p): return [edges.index(tuple(sorted((p[a],p[b])))) for a,b in edges]
gens=[[1,2,3,0,4,5,6,7],[1,0,2,3,4,5,6,7],[0,1,2,3,5,6,7,4],[0,1,2,3,5,4,6,7],[4,5,6,7,0,1,2,3]]
c=Code(G,[vperm(p) for p in gens])
for R in range(2,9): print('K44',R,solve(c,R))
