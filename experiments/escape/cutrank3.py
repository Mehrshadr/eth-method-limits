# heuristic min over balanced cuts (local search, many restarts) of log2 cut-rank of 1_SAT for random 3-XOR, m=n/2
import random, sys
from cutrank import make, code_basis_dim, formula_rank
def ls(n,ch,rnd,iters):
    A=set(rnd.sample(range(n),n//2)); cur=formula_rank(n,ch,list(A))
    for _ in range(iters):
        a=rnd.choice(list(A)); b=rnd.choice([v for v in range(n) if v not in A])
        A2=(A-{a})|{b}; f=formula_rank(n,ch,list(A2))
        if f<=cur: A,cur=A2,f
    return cur
for n in [20,30,40,60,80]:
    m=n//2; out=[]
    for s in range(2):
        ch=make(n,m,1000+s+n); rnd=random.Random(s)
        best=min(ls(n,ch,rnd,1500) for _ in range(6))
        out.append((code_basis_dim(n,ch),best))
    print(f'n={n} m={m} (dimC, heuristic min log2 cut-rank):',out,flush=True)
