import random, itertools, math, numpy as np, time
random.seed(1)
P=1000000007
def rank_mod(M,p):
    M=[[x%p for x in r] for r in M]; r=0
    if not M: return 0
    cols=len(M[0])
    for c in range(cols):
        piv=next((i for i in range(r,len(M)) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p)
        M[r]=[x*inv%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        r+=1
    return r
def rand_cnf(vars_,m):
    return [tuple((v,random.random()<.5) for v in random.sample(vars_,3)) for _ in range(m)]
def sat(cnf,x): return all(any(x[v]==s for v,s in c) for c in cnf)
def count_vec(cnf,S,L):
    # a_F(s) = # assignments of L satisfying F given s
    out=[]
    for s in itertools.product([0,1],repeat=len(S)):
        x=dict(zip(S,s)); t=0
        for l in itertools.product([0,1],repeat=len(L)):
            x.update(zip(L,l)); t+=sat(cnf,x)
        out.append(t)
    return out
print("=== E1: summation-only rank of count matrix C[F,G]=#(F^G) vs pointwise 2^k ===")
print("k | contexts | rank_Q | GF2 | GF3 | GF5 | GF7 | pointwise 2^k")
for k in range(1,7):
    S=[f"s{i}" for i in range(k)]; L=[f"l{i}" for i in range(4)]; R=[f"r{i}" for i in range(4)]
    N=max(3*2**k,40)
    A=[count_vec(rand_cnf(S+L,random.randint(1,int(4.26*(k+4)))),S,L) for _ in range(N)]
    B=[count_vec(rand_cnf(S+R,random.randint(1,int(4.26*(k+4)))),S,R) for _ in range(N)]
    C=[[sum(a*b for a,b in zip(ra,rb)) for rb in B] for ra in A]
    print(k,"|",N,"|",rank_mod(C,P),"|"," | ".join(str(rank_mod(C,p)) for p in (2,3,5,7)),"|",2**k)
