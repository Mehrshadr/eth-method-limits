import random, math
random.seed(1)
for trial in range(5):
    n=random.randint(2,30); d=[random.uniform(0,10) for _ in range(n)]
    beta=random.uniform(0.1,2); c=random.uniform(0.1,3); T=200000; ok=0
    for _ in range(T):
        z=sorted((random.expovariate(beta)-di for di in d), reverse=True)
        ok+= (z[0]-z[1]>c)
    print(n, round(ok/T,4), ">=", round(math.exp(-beta*c),4))
