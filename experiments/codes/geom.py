import numpy as np, itertools
from gsym import gf2_rank
def gf(q):
    # GF(2^h) with polynomial
    polys={2:0b111,4:0b10011,8:0b1011,16:0b10011,32:0b100101}
    h=q.bit_length()-1; P={2:0b111,3:0b1011,4:0b10011,5:0b100101}[h] if h>1 else None
    def mul(a,b):
        r=0
        while b:
            if b&1: r^=a
            b>>=1; a<<=1
            if h>1 and a>>h&1: a^=P
            if h==1: a&=1
        return r
    return mul
def ag_lines(q, m=2):
    mul=gf(q); pts=list(itertools.product(range(q),repeat=m)); pid={p:i for i,p in enumerate(pts)}
    lines=set()
    for p in pts:
        for d in pts:
            if not any(d): continue
            L=frozenset(pid[tuple(p[i]^mul(t,d[i]) for i in range(m))] for t in range(q))
            lines.add(L)
    rows=np.zeros((len(lines),len(pts)),dtype=np.int64)
    for i,L in enumerate(lines): rows[i,list(L)]=1
    return rows, pts, mul
def pg2_lines(q):
    mul=gf(q)
    pts=[]
    for v in itertools.product(range(q),repeat=3):
        if not any(v): continue
        f=next(x for x in v if x)
        # normalise: first nonzero = 1 -> find inverse
        inv=next(t for t in range(1,q) if mul(t,f)==1)
        w=tuple(mul(inv,x) for x in v)
        if w not in pts: pts.append(w)
    lines=[]
    for a in pts:
        L=[i for i,p in enumerate(pts) if (mul(a[0],p[0])^mul(a[1],p[1])^mul(a[2],p[2]))==0]
        lines.append(L)
    rows=np.zeros((len(lines),len(pts)),dtype=np.int64)
    for i,L in enumerate(lines): rows[i,L]=1
    return rows, pts
if __name__=="__main__":
    for q in (4,8):
        G,_,_=ag_lines(q); print('AG(2,%d)'%q, G.shape[1], gf2_rank(np.vstack([G,np.ones(G.shape[1],int)])))
    G,_,_=ag_lines(4,3); print('AG(3,4) lines', G.shape[1], gf2_rank(np.vstack([G,np.ones(G.shape[1],int)])))
    for q in (2,4,8):
        G,_=pg2_lines(q); print('PG(2,%d)'%q, G.shape[1], gf2_rank(np.vstack([G,np.ones(G.shape[1],int)])))
