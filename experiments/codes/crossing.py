"""Step 4-5 sanity check: on random 3-CNFs where every variable occurs exactly 3 times,
even a partition into r equal groups that is optimized to keep clauses inside groups
leaves most clauses cut, and almost every group receives a crossing variable."""
import random
def formula(n, seed):
    rnd = random.Random(seed)
    while True:
        slots = [v for v in range(n) for _ in range(3)]; rnd.shuffle(slots)
        cl = [slots[i:i + 3] for i in range(0, 3 * n, 3)]
        if all(len(set(c)) == 3 for c in cl): return cl
def optimize(n, cl, r, seed, rounds=30):
    rnd = random.Random(seed); cap = (n + r - 1) // r
    order = list(range(n)); rnd.shuffle(order)
    grp = [0] * n
    for i, v in enumerate(order): grp[v] = i * r // n
    occ = [[] for _ in range(n)]
    for j, c in enumerate(cl):
        for v in c: occ[v].append(j)
    def inside(j): return len({grp[v] for v in cl[j]}) == 1
    for _ in range(rounds):                       # swap pairs of variables when it keeps more clauses whole
        for v in rnd.sample(range(n), n):
            for w in rnd.sample(range(n), 8):
                if grp[v] == grp[w]: continue
                js = set(occ[v] + occ[w]); before = sum(inside(j) for j in js)
                grp[v], grp[w] = grp[w], grp[v]
                if sum(inside(j) for j in js) < before: grp[v], grp[w] = grp[w], grp[v]
    return grp
for n, r in [(600, 10), (600, 30), (2000, 20), (2000, 60)]:
    cl = formula(n, 1); grp = optimize(n, cl, r, 2)
    whole = sum(len({grp[v] for v in c}) == 1 for c in cl)
    # clause group = majority group of its variables (the cheapest choice); crossing variables = the others
    touched = set()
    for c in cl:
        gs = [grp[v] for v in c]; home = max(set(gs), key=gs.count)
        touched |= {g for g in gs if g != home}
    print(f"n={n} r={r} group size={n//r}: clauses kept whole {whole}/{len(cl)} = {whole/len(cl):.1%}, groups receiving a crossing variable {len(touched)}/{r}")
