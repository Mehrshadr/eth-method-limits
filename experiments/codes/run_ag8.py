import numpy as np, json, sys
from geom import ag_lines
from run_geom import ag_perms
from families import run
from gsym import Code
q=8; G,pts,mul=ag_lines(q); G=np.vstack([np.ones(64,int),G])
out=open('fam_results.jsonl','a')
run('lines AG(2,8)', G, ag_perms(q,pts,mul), [int(x) for x in sys.argv[1:]], out=out, budget=560)
