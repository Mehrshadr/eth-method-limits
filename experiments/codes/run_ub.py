import sys, numpy as np, time
from gsym import Code, cube_QR
from infoset import ub_QR
from families import rm, torus, hypercube_graph, cut_code
from geom import ag_lines
from borsuk import bipartite_cut_code
fams = {
 'K4,8': lambda: (bipartite_cut_code(4,8),), 'K6,6': lambda: (bipartite_cut_code(6,6),), 'K8,8': lambda: (bipartite_cut_code(8,8),),
 'RM(2,5)': lambda: (rm(2,5)[0],), 'RM(2,6)': lambda: (rm(2,6)[0],),
 'C4xC4': lambda: (torus(4,4)[0],), 'C4xC6': lambda: (torus(4,6)[0],),
 'AG(2,4)': lambda: (np.vstack([np.ones(16,int), ag_lines(4)[0]]),),
 'AG(2,8)': lambda: (np.vstack([np.ones(64,int), ag_lines(8)[0]]),),
}
for name in sys.argv[1:]:
    t=time.time(); code=Code(fams[name]()[0])
    for thr in (0.05, 0.01):
        r=ub_QR(code, thr)
        print(name, 'k',code.k,'Q',code.Q,'thr',thr,'UB Q/R',None if r is None else round(r[0],3),'(d,dI)',None if r is None else r[1:], 'cube',round(cube_QR(code.k,thr),3), 'Qcube', round(cube_QR(code.Q,thr),3), f'{time.time()-t:.0f}s', flush=True)
