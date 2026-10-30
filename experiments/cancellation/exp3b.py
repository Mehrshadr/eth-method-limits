import random, itertools, numpy as np, math
from exp import rank_mod, P, sat
random.seed(11)
def rcnf(n,m): return [tuple((v,random.random()<.5) for v in random.sample(range(n),3)) for _ in range(m)]
print("=== cut-rank, single formula, balanced cut, density ~4.0 ===")
for n in (8,10,12):
    cnf=rcnf(n,int(2.0*n)); h=n//2
    M=[[int(sat(cnf,a+b)) for b in itertools.product([0,1],repeat=n-h)] for a in itertools.product([0,1],repeat=h)]
    X=np.linalg.qr(np.random.randn(2**h,2**h))[0]
    Mm=np.array(M,float)
    print(f"n={n} #F={int(Mm.sum())} rank(1_F)={rank_mod(M,P)} rank(X·1_F·X^T)={np.linalg.matrix_rank(X@Mm@X.T)}  sum preserved: {round(float(np.ones(2**h)@X.T@(X@Mm@X.T)@X@np.ones(2**h)),6)}")
print("\n=== E5: determinantal complexity of generating poly G_F(z)=sum_{x sat} prod z_i^{x_i} ===")
print("n | m | #F | Hessian-rank LB on dc (Mignon-Ressayre) | best OBDD nodes (UB on dc via ABP) over 300 orders")
def obdd_size(cnf,n,order):
    models=set(x for x in itertools.product([0,1],repeat=n) if sat(cnf,x))
    tot=0
    for lvl in range(n+1):
        fixed=order[:lvl]; rest=order[lvl:]
        subs=set()
        for a in itertools.product([0,1],repeat=lvl):
            key=[]
            for b in itertools.product([0,1],repeat=n-lvl):
                x=[0]*n
                for v,val in zip(fixed,a): x[v]=val
                for v,val in zip(rest,b): x[v]=val
                key.append(tuple(x) in models)
            if any(key): subs.add(tuple(key))
        tot+=len(subs)
    return tot
for n in range(4,12):
    while True:
        cnf=rcnf(n,round(2.0*n)); mods=[x for x in itertools.product([0,1],repeat=n) if sat(cnf,x)]
        if len(mods)>=1: break
    if len(mods)<2:
        cnf=rcnf(n,round(3.5*n)); mods=[x for x in itertools.product([0,1],repeat=n) if sat(cnf,x)]
    def G(z): return sum(math.prod(z[i] for i in range(n) if x[i]) for x in mods)%P
    best=0
    for _ in range(20):
        z=[random.randrange(1,P) for _ in range(n)]
        # make G(z)=0 by solving for z0 (multilinear): G=a+b z0
        z0=z[:]; z0[0]=0; a=G(z0); z0[0]=1; b=(G(z0)-a)%P
        if b==0: continue
        z[0]=(-a*pow(b,P-2,P))%P
        H=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i!=j: H[i][j]=sum(math.prod(z[k] for k in range(n) if x[k] and k not in (i,j)) for x in mods if x[i] and x[j])%P
        best=max(best,rank_mod(H,P))
    ob=min(obdd_size(cnf,n,random.sample(range(n),n)) for _ in range(200 if n<=8 else 25))
    print(n,"|",len(cnf),"|",len(mods),"|",math.ceil(best/2),"|",ob)
