"""Diagnostics on the exponent fit, then the slowing-to-delay table from measured rates."""
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

Q = {}
for (rg, d), v in by.items():
    if len(v) < 200: continue
    age, L0 = meta[d]; u0 = (L0-35)/10.0
    if u0 <= 0: continue
    Q[(rg, d)] = (age, u0, {q: ((float(np.quantile(v, q))-35)/10.0
                                if float(np.quantile(v, q)) > L0+1 else None) for q in QS})

def obj(reg, p, use_age=True):
    D = [(d, Q[(rg, d)]) for (rg, d) in Q if rg == reg]
    tot = 0.0
    for q in QS:
        ls = [math.log(sfun(u0, qq[q], p)/(age if use_age else 1.0))
              for d, (age, u0, qq) in D if qq[q] is not None]
        if len(ls) >= 4: tot += float(np.var(ls))
    return tot

print('objective = across-donor variance of log(implied rate), summed over 5 quantiles')
print('lower is better; a real exponent shows an interior minimum\n')
print('%-6s %12s %12s' % ('p', 'striatum', 'cortex'))
for p in [0.0,0.5,1.0,1.25,1.5,1.55,1.75,2.0,2.35,2.75,3.5,4.5,6.0]:
    print('%-6.2f %12.3f %12.3f' % (p, obj('striatum', p), obj('cortex', p)))

print('\nreference: same objective with age ignored (i.e. no clock at all)')
for reg in ('striatum','cortex'):
    best_age = min(obj(reg, p) for p in np.arange(0,6.01,0.05))
    best_noage = min(obj(reg, p, use_age=False) for p in np.arange(0,6.01,0.05))
    print('  %-9s best with age %.3f   best ignoring age %.3f   variance explained by age %.0f%%'
          % (reg, best_age, best_noage, 100*(1-best_age/best_noage)))

# ---- fitted parameters (from exp4)
FIT = {'striatum': (2.35, 0.1872), 'cortex': (1.55, 0.0972)}
print('\ninstantaneous expansion speed implied by the fit, CAG per year')
print('%-9s %8s %8s %8s %8s %8s' % ('region','at 42','at 60','at 80','at 120','at 150'))
for reg,(p,k) in FIT.items():
    sp = [k*(( (L-35)/10.0 )**p) for L in (42,60,80,120,150)]
    print('%-9s %8.2f %8.2f %8.2f %8.2f %8.2f' % (reg,*sp))

print("""
Rate and time enter the model only as their product, so slowing expansion by a
fraction f multiplies every remaining interval by 1/(1-f).  Starting treatment at
age A when the untreated milestone would fall at age T buys (T-A)*f/(1-f) years.""")
p, k = FIT['striatum']
print('\nstriatum, median cell, inherited length 42: reaches 150 CAGs at age %.0f'
      % (sfun((42-35)/10.0, (150-35)/10.0, p)/k))
T = sfun((42-35)/10.0, (150-35)/10.0, p)/k
print('\n%-14s %-12s %-12s %-12s' % ('slowing', 'start at 25', 'start at 40', 'start at 50'))
for f in (0.1, 0.2, 0.3, 0.5, 0.75):
    row = ['%.0f yr' % ((T-A)*f/(1-f)) for A in (25, 40, 50)]
    print('%-14s %-12s %-12s %-12s' % ('%d%%' % (100*f), *row))
