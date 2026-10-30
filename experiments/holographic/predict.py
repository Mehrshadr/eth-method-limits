from fast import works
# prediction: works iff p | 2^k-1 and p = +-1 mod 8 (2 is a square mod p)
for k, p in [(14, 43), (7, 127), (20, 41), (21, 127)]:
    pred = (2 ** k - 1) % p == 0 and p % 8 in (1, 7)
    print(f"k={k} p={p}: predicted {'works' if pred else 'fails'} -> actual {'works' if works(p, k) else 'fails'}", flush=True)
