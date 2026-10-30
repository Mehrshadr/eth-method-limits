import sys, json
from families import *
out=open('fam_results.jsonl','a')
which=sys.argv[1]
if which=='rm25': G,p=rm(2,5); run('RM(2,5)',G,p,None,out=out)
if which=='rm26': G,p=rm(2,6); run('RM(2,6)',G,p,None,out=out)
if which=='q4': G,p=hypercube_graph(4); run('cut(Q4 graph)',G,p,None,out=out)
if which=='t44': G,p=torus(4,4); run('cut(C4xC4)',G,p,None,out=out)
if which=='t46': G,p=torus(4,6); run('cut(C4xC6)',G,p,None,out=out)
