"""Identify the growth exponent from age variation alone.

Model:  u = (L-35)/10 ,  du/dt = (k/10) u^p ,  k fixed per cell, variable across cells.
Integrated:  s(u0 -> u ; p) = k * t.
The q-th quantile cell of a donor has the q-th quantile k, the same number in every
donor.  So  s(u0_d -> L_q,d ; p)  must be proportional to that donor's age.
The exponent is whichever p makes  s/age  constant across donors of different ages.
No assumption about the shape of the rate distribution.
"""
import os, sys
import gzip, csv, math, collections
import numpy as np

U = (sys.argv[1].rstrip('/') + '/' if len(sys.argv) > 1
     else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')
meta = {}
for r in csv.DictReader(open(U + 'donor_metadata.txt'), delimiter='\t'):
    try: meta[r['DONOR']] = (float(r['AGE']), float(r['HD_CAG_LENGTH']))
    except Exception: pass

by = collections.defaultdict(list)
for r in csv.DictReader(gzip.open(U + 'per_nucleus_repeat_lengths.csv.gz', 'rt')):
    try: L = float(r['CAG_LENGTH'])
    except Exception: continue
    if r['DONOR'] not in meta: continue
    by[('striatum' if 'caudate' in r['SOURCE'] else 'cortex', r['DONOR'])].append(L)

QS = [0.5, 0.75, 0.9, 0.95, 0.99]

def sfun(u0, u, p):
    if abs(p-1) < 1e-9: return 10.0*np.log(u/u0)
    return 10.0*(u0**(1-p) - u**(1-p))/(p-1)

# the quantiles do not depend on p, so compute them once
Q = {}   # (region, donor) -> (age, u0, {q: uq or None})
for (rg, d), v in by.items():
    if len(v) < 200: continue
    age, L0 = meta[d]
    u0 = (L0-35)/10.0
    if u0 <= 0: continue
    qq = {}
    for q in QS:
        Lq = float(np.quantile(v, q))
        qq[q] = (Lq-35)/10.0 if Lq > L0 + 1 else None
    Q[(rg, d)] = (age, u0, qq)

def resid(reg, p, donors=None):
    D = [(d, Q[(rg, d)]) for (rg, d) in Q if rg == reg and (donors is None or d in donors)]
    if len(D) < 4: return None, None
    tot, detail = 0.0, {}
    for q in QS:
        ls = [math.log(sfun(u0, qq[q], p)/age) for d, (age, u0, qq) in D if qq[q] is not None]
        if len(ls) < 4: continue
        tot += float(np.var(ls)); detail[q] = (float(np.std(ls)), math.exp(float(np.mean(ls))))
    return tot, detail

GRID = np.round(np.arange(0.0, 6.001, 0.05), 3)

print('exponent p in  d(repeat)/dt  proportional to ((repeat-35)/10)^p')
print('identified by requiring the implied rate to be age-independent\n')
out = {}
for reg in ('striatum', 'cortex'):
    vals = [(resid(reg, p)[0], p) for p in GRID]
    vals = [(v, p) for v, p in vals if v is not None]
    v0, p0 = min(vals)
    # bootstrap over donors
    D = sorted(d for (rg, d) in Q if rg == reg)
    rng = np.random.default_rng(7)
    bs = []
    for _ in range(2000):
        samp = set(rng.choice(D, size=len(D), replace=True))
        vv = [(resid(reg, p, donors=samp)[0], p) for p in GRID]
        vv = [(a, b) for a, b in vv if a is not None]
        if vv: bs.append(min(vv)[1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    tot, det = resid(reg, p0)
    out[reg] = (p0, det, len(D))
    print('%-9s %d donors   p = %.2f   95%% bootstrap interval %.2f to %.2f' % (reg, len(D), p0, lo, hi))
    for q in QS:
        if q in det:
            sd, k = det[q]
            print('     %4.0fth percentile cell: rate k = %.4f /yr   residual spread %.2f (log units)'
                  % (100*q, k, sd))
    print()

np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fit4.npy'), np.array([out['striatum'][0], out['cortex'][0]]))

print('implied years from the inherited length to 150 CAGs')
for reg in ('striatum', 'cortex'):
    p, det, nd = out[reg]
    print(' %s  (p = %.2f)' % (reg, p))
    for q in (0.5, 0.9, 0.99):
        if q not in det: continue
        k = det[q][1]
        row = []
        for L0 in (40, 42, 45, 50):
            t = sfun((L0-35)/10.0, (150-35)/10.0, p)/k
            row.append('%d:%5.0f' % (L0, t))
        print('   %4.0fth pct cell (k=%.4f):  ' % (100*q, k) + '   '.join(row))
