import itertools, random
from cutrank import make, code_basis_dim, formula_rank
for n in [12,16,20,22]:
  for frac in [0.5,0.6,0.7]:
    m=int(n*frac); res=[]
    for s in range(3):
        ch=make(n,m,7*s+n); k=code_basis_dim(n,ch)
        best=min(formula_rank(n,ch,list(A)) for A in itertools.combinations(range(n),n//2))
        res.append((k,best))
    print(f'n={n} m={m} (dimC, min_cut log2rank):',res, flush=True)
