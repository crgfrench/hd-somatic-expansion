"""Why age gives no signal: in a post-mortem HD cohort, age at death is an
outcome of the expansion, not an exposure time."""
import os, sys
import gzip, csv, math, collections
import numpy as np
from scipy.stats import spearmanr

U = (sys.argv[1].rstrip('/') + '/' if len(sys.argv) > 1
     else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')
meta, allm = {}, []
for r in csv.DictReader(open(U + 'donor_metadata.txt'), delimiter='\t'):
    try:
        rec = (r['DONOR'], r['STATUS'], float(r['AGE']), float(r['HD_CAG_LENGTH']), float(r['CAP_SCORE']))
    except Exception: continue
    allm.append(rec); meta[r['DONOR']] = rec[2:]

hd = [m for m in allm if m[1] == 'pwHD']
print('all %d HD donors in the metadata' % len(hd))
rho, p = spearmanr([m[2] for m in hd], [m[3] for m in hd])
print('  age at death vs inherited CAG length:  rho = %+.2f  p = %.2g' % (rho, p))
print('  i.e. the donors who lived longest are the ones with the shortest repeats.')

by = collections.defaultdict(list)
for r in csv.DictReader(gzip.open(U + 'per_nucleus_repeat_lengths.csv.gz', 'rt')):
    try: L = float(r['CAG_LENGTH'])
    except Exception: continue
    if r['DONOR'] not in meta: continue
    by[('striatum' if 'caudate' in r['SOURCE'] else 'cortex', r['DONOR'])].append(L)

print('\nthe 11 cortex donors with per-nucleus lengths')
D = [(d, np.array(v)) for (rg, d), v in by.items() if rg == 'cortex' and len(v) >= 200]
age = np.array([meta[d][0] for d, v in D])
inh = np.array([meta[d][1] for d, v in D])
cap = np.array([meta[d][2] for d, v in D])
gain = np.array([np.quantile(v, .9) - meta[d][1] for d, v in D])
share = np.array([np.mean(v >= 150) for d, v in D])

for lab, x in (('age', age), ('inherited CAG', inh), ('CAP score', cap)):
    for ylab, y in (('p90 expansion', gain), ('share past 150', share)):
        r_, p_ = spearmanr(x, y)
        print('  %-14s vs %-16s rho = %+.2f  p = %.3f' % (lab, ylab, r_, p_))
r_, p_ = spearmanr(age, inh)
print('  %-14s vs %-16s rho = %+.2f  p = %.3f' % ('age', 'inherited CAG', r_, p_))

# partial: does age predict expansion once inherited length is accounted for?
def partial(x, y, z):
    rxy, _ = spearmanr(x, y); rxz, _ = spearmanr(x, z); ryz, _ = spearmanr(y, z)
    num = rxy - rxz*ryz
    den = math.sqrt((1-rxz**2)*(1-ryz**2))
    r = num/den if den else float('nan')
    n = len(x); dof = n-3
    from scipy.stats import t as tdist
    tt = r*math.sqrt(dof/(1-r*r)) if abs(r) < 1 else float('inf')
    return r, float(2*tdist.sf(abs(tt), dof))
r_, p_ = partial(age, np.log(gain), inh)
print('\n  age vs log p90 expansion, holding inherited length fixed: rho = %+.2f  p = %.3f  (n=11, dof=8)'
      % (r_, p_))

print("""
Consequence.  Age at death in this cohort is not a clock.  It is set by the
disease, so the donors with the most time are exactly the donors whose repeats
grow slowest.  The exponent therefore stays unidentified even with ages in hand:
the one thing that would pin it down -- the same inherited length observed after
different amounts of time -- is what dying of Huntington's disease prevents.""")

# ---------- model-free version of question 4
print("""
What survives is model-free.  In any growth law du/dt = k*g(u), rate and time
appear only as a product, so cutting the rate by a fraction f multiplies every
remaining interval by 1/(1-f).  Anchoring on observed clinical onset instead of a
fitted curve (expected onset = 649/(CAG-30), the standard CAP calibration):""")
print('\n%-9s %-9s' % ('inherited', 'onset') + ''.join('%14s' % ('slow %d%%' % (100*f)) for f in (0.2,0.3,0.5,0.75)))
for L0 in (40, 42, 45, 50):
    T = 649.0/(L0-30)
    cells = []
    for f in (0.2, 0.3, 0.5, 0.75):
        A = 30.0   # treatment started at 30
        cells.append('%14s' % ('+%.0f yr' % ((T-A)*f/(1-f))))
    print('%-9d %-9.0f' % (L0, T) + ''.join(cells))
print('(treatment begun at age 30; delay = (onset - start) * f/(1-f))')

print('\nslowing required to buy a given delay, inherited 42 (onset 54)')
T = 649.0/12
for A in (20, 30, 40, 50):
    row = ['%.0f%%' % (100*want/((T-A)+want)) for want in (5, 10, 20)]
    print('   start at %2d: 5 yr needs %s   10 yr needs %s   20 yr needs %s' % (A, *row))
