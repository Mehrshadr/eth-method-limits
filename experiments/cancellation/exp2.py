import random, itertools, numpy as np
from exp import rank_mod, P, sat
random.seed(7)
def rcnf(n,m): return [tuple((v,random.random()<.5) for v in random.sample(range(n),3)) for _ in range(m)]
def models(cnf,n): return [x for x in itertools.product([0,1],repeat=n) if sat(cnf,x)]
print("\n=== E2: holographic (basis-change) rewrite: terms wrong, total right, cut-rank unchanged ===")
for (n,m) in [(4,5),(5,6)]:
    cnf=rcnf(n,m)
    occ={v:[] for v in range(n)}; edges=[]
    for ci,c in enumerate(cnf):
        for j,(v,s) in enumerate(c): occ[v].append(len(edges)); edges.append((ci,v,s))
    E=len(edges); T=np.random.randn(2,2); Ti=np.linalg.inv(T)
    # variable tensor EQ_d then apply T on each leg; clause tensor apply Ti on each leg
    def vt(v,y):  # transformed EQ: sum_b prod_e T[b,y_e]
        return sum(np.prod([T[b,y[e]] for e in occ[v]]) for b in (0,1))
    def ct(ci,y):
        es=[i for i,e in enumerate(edges) if e[0]==ci]
        tot=0
        for bits in itertools.product([0,1],repeat=3):
            if any(bits[j]==cnf[ci][j][1] for j in range(3)):
                tot+=np.prod([Ti[y[es[j]],bits[j]] for j in range(3)])
        return tot
    total=0; nz=0; wrong=0
    for y in itertools.product([0,1],repeat=E):
        g=np.prod([vt(v,y) for v in range(n)])*np.prod([ct(c,y) for c in range(m)])
        total+=g; nz+=abs(g)>1e-12
        # original term: 1 iff y consistent per variable and satisfies all clauses
        cons=all(len({y[e] for e in occ[v]})<=1 for v in range(n))
        orig=1.0 if cons and all(any(y[es]==edges[es][2] for es in [i for i,e in enumerate(edges) if e[0]==ci]) for ci in range(m)) else 0.0
        wrong+=abs(g-orig)>1e-9
    print(f"n={n} m={m} edges={E}: #F={len(models(cnf,n))}  holographic total={total:.6f}  nonzero terms={nz}/{2**E}  terms differing from original={wrong}")
# cut rank invariance: matrix of partial contraction across variable split
print("cut-rank invariance: rank of left-half partial tensor over cut edges, original vs transformed")
n,m=8,34
cnf=rcnf(n,m); A=set(range(4))
cut=[(ci,j) for ci,c in enumerate(cnf) for j,(v,s) in enumerate(c) if (v in A)!=(cnf[ci][0][0] in A)]
# simpler: pointwise matrix 1_F[xA,xB] rank and its transform by random invertible maps on each side
M=np.array([[1.0 if sat(cnf,a+b) else 0 for b in itertools.product([0,1],repeat=4)] for a in itertools.product([0,1],repeat=4)])
X=np.random.randn(16,16); Y=np.random.randn(16,16)
print(" rank 1_F =",np.linalg.matrix_rank(M)," rank X M Y =",np.linalg.matrix_rank(X@M@Y), " sum(M)=",M.sum()," (1^T X^-1)(XMY)(Y^-1 1)=",round(float(np.ones(16)@np.linalg.inv(X)@(X@M@Y)@np.linalg.inv(Y)@np.ones(16)),6))

print("\n=== E3: cut-rank (Nisan partial-derivative) of det_n vs perm_n, rows split h|n-h ===")
from itertools import permutations
def sign(p):
    s=1;p=list(p)
    for i in range(len(p)):
        for j in range(i+1,len(p)):
            if p[i]>p[j]: s=-s
    return s
for n in range(2,8):
    h=n//2
    L=list(permutations(range(n),h)); R=list(permutations(range(n),n-h))
    li={x:i for i,x in enumerate(L)}; ri={x:i for i,x in enumerate(R)}
    Md=[[0]*len(R) for _ in L]; Mp=[[0]*len(R) for _ in L]
    for p in permutations(range(n)):
        a,b=li[p[:h]],ri[p[h:]]; Md[a][b]=sign(p); Mp[a][b]=1
    print(f"n={n}: rank(det)={rank_mod(Md,P)} rank(perm)={rank_mod(Mp,P)} binom={__import__('math').comb(n,h)}")
