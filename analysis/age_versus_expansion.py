"""Does age predict how far the repeat has expanded?  Tested per region, then the
delay table built on the region where it does."""
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

def rank(v):
    idx = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
    while i < len(idx):
        j = i
        while j < len(idx) and v[idx[j]] == v[idx[i]]: j += 1
        for k in range(i, j): r[idx[k]] = (i+j+1)/2
        i = j
    return r

def spearman(x, y):
    rx, ry = rank(x), rank(y)
    rx = np.array(rx); ry = np.array(ry)
    rho = float(np.corrcoef(rx, ry)[0,1]); n = len(x)
    if abs(rho) >= 1 or n < 4: return rho, float('nan')
    t = rho*math.sqrt((n-2)/(1-rho*rho))
    from scipy.stats import t as tdist
    return rho, float(2*tdist.sf(abs(t), n-2))

print('expansion achieved versus age at death, across donors')
print('(gain is the excess over the donor\'s own inherited length)\n')
for reg in ('striatum', 'cortex'):
    D = [(d, v) for (rg, d), v in by.items() if rg == reg and len(v) >= 200]
    D.sort(key=lambda t: meta[t[0]][0])
    print('  %s, %d donors' % (reg, len(D)))
    print('    %-8s %5s %6s %7s %8s %8s %8s' % ('donor','age','inher','med gain','p90 gain','p99 gain','%>150'))
    for d, v in D:
        age, L0 = meta[d]; v = np.array(v)
        print('    %-8s %5.0f %6.0f %7.1f %8.1f %8.1f %8.1f'
              % (d, age, L0, np.median(v)-L0, np.quantile(v,.9)-L0, np.quantile(v,.99)-L0,
                 100*np.mean(v >= 150)))
    ages = [meta[d][0] for d, v in D]
    for lab, f in (('median gain', lambda v,L0: np.median(v)-L0),
                   ('p90 gain',    lambda v,L0: np.quantile(v,.9)-L0),
                   ('p99 gain',    lambda v,L0: np.quantile(v,.99)-L0),
                   ('share >150',  lambda v,L0: np.mean(np.array(v)>=150))):
        y = [f(np.array(v), meta[d][1]) for d, v in D]
        rho, p = spearman(ages, y)
        print('    age vs %-12s rho = %+.2f   p = %.3f' % (lab, rho, p))
    # inherited length is the competing predictor
    inh = [meta[d][1] for d, v in D]
    rho, p = spearman(inh, [np.quantile(np.array(v),.9)-meta[d][1] for d, v in D])
    print('    inherited length vs p90 gain   rho = %+.2f   p = %.3f' % (rho, p))
    print()

# ---------- delay table, cortex only (the region where the clock is visible)
P = 1.55
def sfun(u0, u, p=P):
    return 10.0*(u0**(1-p) - u**(1-p))/(p-1)

KQ = {0.5: 0.0972, 0.75: 0.1512, 0.9: 0.2117, 0.95: 0.2423, 0.99: 0.2928}
print('cortex: age at which a given share of neurons passes 150 CAGs')
print('%-10s' % 'inherited' + ''.join('%12s' % ('top %g%%' % (100*(1-q))) for q in KQ))
AGES = {}
for L0 in (40, 42, 45, 50):
    row = []
    for q, k in KQ.items():
        t = sfun((L0-35)/10.0, (150-35)/10.0)/k
        AGES[(L0, q)] = t
        row.append('%12s' % ('%.0f' % t if t < 110 else '>110'))
    print('%-10d' % L0 + ''.join(row))

print('\nyears of delay bought by slowing expansion, cortex, inherited length 42')
print('milestone: the top 1%% of cortical neurons passing 150, untreated at age %.0f' % AGES[(42,0.99)])
T = AGES[(42, 0.99)]
print('%-10s %-12s %-12s %-12s %-12s' % ('slowing','start at 20','start at 30','start at 40','start at 50'))
for f in (0.1, 0.2, 0.3, 0.5, 0.75, 0.9):
    row = ['%.0f yr' % ((T-A)*f/(1-f)) for A in (20, 30, 40, 50)]
    print('%-10s %-12s %-12s %-12s %-12s' % ('%d%%' % (100*f), *row))

print('\nthe same, read the other way: slowing needed to buy a target delay from age 40')
for want in (5, 10, 20):
    f = want/((T-40)+want)
    print('   %2d years of delay needs expansion slowed by %.0f%%' % (want, 100*f))
