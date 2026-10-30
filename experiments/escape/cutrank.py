# Fixed-partition obstruction: rank over Q of 1_SAT(x_A,x_B) for a 3-XOR system (=3-CNF, 4 clauses per check)
import itertools, random, sys
import numpy as np
def gf2rank(rows):
    rows=[r for r in rows if r]; rk=0
    while rows:
        p=max(rows); rk+=1; hb=p.bit_length()-1
        rows=[r^p if (r>>hb)&1 else r for r in rows if r!=p]
        rows=[r for r in rows if r]
    return rk
def make(n,m,seed):
    rnd=random.Random(seed); checks=[]
    occ=[0]*n
    for _ in range(m):
        vs=rnd.sample(range(n),3); checks.append(sum(1<<v for v in vs))
    return checks
def code_basis_dim(n,checks): return n-gf2rank(checks)
def shortened_dim(n,checks,S):  # dim of codewords supported inside S
    # restrict parity checks to columns S
    mask=sum(1<<v for v in S); return len(S)-gf2rank([c&mask for c in checks])
def formula_rank(n,checks,A):
    B=[v for v in range(n) if v not in A]
    return code_basis_dim(n,checks)-shortened_dim(n,checks,A)-shortened_dim(n,checks,B)
def exact_rank(n,checks,A):
    B=[v for v in range(n) if v not in A]
    M=np.zeros((1<<len(A),1<<len(B)))
    for x in range(1<<n):
        if all(bin(x&c).count('1')%2==0 for c in checks):
            a=sum(((x>>v)&1)<<i for i,v in enumerate(A)); b=sum(((x>>v)&1)<<i for i,v in enumerate(B))
            M[a,b]=1
    return np.linalg.matrix_rank(M)
if __name__=="__main__":
    # 1) verify formula on small cases
    for s in range(5):
        n=12; ch=make(n,5,s); A=random.Random(s).sample(range(n),6)
        print('verify', 2**formula_rank(n,ch,A), exact_rank(n,ch,A))
    # 2) min over ALL balanced cuts
    for n in [12,16,20]:
        m=n//3
        for s in range(2):
            ch=make(n,m,100+s); k=code_basis_dim(n,ch); best=10**9
            for A in itertools.combinations(range(n),n//2):
                fr=formula_rank(n,ch,list(A)); best=min(best,fr)
            print(f'n={n} dimC={k} min over balanced cuts log2 rank={best}')
