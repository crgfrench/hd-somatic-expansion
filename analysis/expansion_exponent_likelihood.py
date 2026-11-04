import os, sys
import gzip, csv, math, collections
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

U = (sys.argv[1].rstrip('/') + '/' if len(sys.argv) > 1
     else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')

# ---- donor metadata: age at death and the inherited (blood) repeat length
meta = {}
for r in csv.DictReader(open(U + 'donor_metadata.txt'), delimiter='\t'):
    try:
        meta[r['DONOR']] = (float(r['AGE']), float(r['HD_CAG_LENGTH']), float(r['CAP_SCORE']))
    except Exception:
        pass

rows = []
for r in csv.DictReader(gzip.open(U + 'per_nucleus_repeat_lengths.csv.gz', 'rt')):
    try: L = float(r['CAG_LENGTH'])
    except Exception: continue
    d = r['DONOR']
    if d not in meta: continue
    rows.append((d, 'striatum' if 'caudate' in r['SOURCE'] else 'cortex', L))

print('nuclei with a repeat length and donor metadata:', len(rows))

by = collections.defaultdict(list)
for d, reg, L in rows: by[(reg, d)].append(L)

print('\ninherited length: measured in blood vs the 2nd percentile of the nuclei')
print('%-9s %-8s %5s %5s %7s %8s %8s %7s' % ('region','donor','n','age','blood','2nd pct','median','%>150'))
for (reg, d), v in sorted(by.items()):
    if len(v) < 100: continue
    v = sorted(v); age, L0, cap = meta[d]
    print('%-9s %-8s %5d %5.0f %7.0f %8.0f %8.0f %7.1f'
          % (reg, d, len(v), age, L0, v[int(.02*(len(v)-1))], np.median(v),
             100*sum(1 for x in v if x >= 150)/len(v)))

# ---- growth model.  u = (L-35)/10 ,  du/dt = (k/10) u^p
# time to go from u0 to u:   s(u) = k*t = 10*(u0^(1-p) - u^(1-p))/(p-1)
# each nucleus: k_i = s_i / age_d .  k is lognormal(mu, sg) across cells, shared.
# ages differ between donors, so only one p makes the implied k's agree.
def sfun(u0, u, p):
    return 10.0*(u0**(1-p) - u**(1-p))/(p-1)

def prep(reg, p, pad=3.0):
    K, LJ, KMIN, NT = [], [], [], []
    for (rg, d), v in by.items():
        if rg != reg or len(v) < 100: continue
        age, L0, cap = meta[d]
        u0 = (L0-35)/10.0
        if u0 <= 0: continue
        uth = (L0+pad-35)/10.0
        kmin = sfun(u0, uth, p)/age
        sel = [L for L in v if L >= L0+pad]
        if len(sel) < 30: continue
        u = (np.array(sel)-35)/10.0
        s = sfun(u0, u, p)
        ok = s > 0
        k = s[ok]/age
        K.append(k)
        # log|dk/dL| = log( (1/age) * (1/10) * d s/d u * ... ) ; ds/du = 10*u^-p
        LJ.append(-p*np.log(u[ok]) - math.log(age))
        KMIN.append(np.full(k.shape, kmin))
        NT.append(len(sel))
    if not K: return None
    return np.concatenate(K), np.concatenate(LJ), np.concatenate(KMIN)

def nll(par, k, lj, kmin):
    mu, ls = par[0], par[1]
    sg = math.exp(ls)
    z = (np.log(k)-mu)/sg
    lp = -np.log(k) - ls - 0.5*math.log(2*math.pi) - 0.5*z*z
    # left truncation: condition on k > kmin (cells below the pad are excluded)
    tail = norm.sf((np.log(kmin)-mu)/sg)
    tail = np.clip(tail, 1e-12, 1.0)
    return -float(np.sum(lp + lj - np.log(tail)))

print('\nexponent in  d(repeat)/dt  proportional to  ((repeat-35)/10)^p')
print('%-9s %6s %8s %10s %8s %9s' % ('region','p','logLik','median k','spread','n cells'))
best = {}
for reg in ('striatum','cortex'):
    curve = []
    for p in np.arange(1.1, 4.01, 0.05):
        d = prep(reg, p)
        if d is None: continue
        k, lj, kmin = d
        r = minimize(nll, [math.log(np.median(k)), math.log(0.8)], args=(k,lj,kmin), method='Nelder-Mead',
                     options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
        curve.append((-r.fun, p, r.x[0], math.exp(r.x[1]), len(k)))
    curve.sort(reverse=True)
    ll, p, mu, sg, n = curve[0]
    best[reg] = (p, mu, sg, n, curve)
    print('%-9s %6.2f %8.0f %10.4f %8.2f %9d' % (reg, p, ll, math.exp(mu), sg, n))
    # 95% interval on p by a 1.92 drop in log-likelihood
    lo = min(c[1] for c in curve if c[0] > ll-1.92)
    hi = max(c[1] for c in curve if c[0] > ll-1.92)
    print('          p 95%% interval %.2f to %.2f' % (lo, hi))

np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'best.npy'), np.array([best['striatum'][0], best['striatum'][1], best['striatum'][2],
                                              best['cortex'][0],  best['cortex'][1],  best['cortex'][2]]))

# ---- what the fitted model implies: years from inherited length to 150
print('\nyears to reach 150 CAGs, from the fitted median rate')
for reg in ('striatum','cortex'):
    p, mu, sg, n, curve = best[reg]
    k = math.exp(mu)
    print(' %s (p=%.2f, median rate k=%.4f)' % (reg, p, k))
    for L0 in (40, 42, 45, 50):
        u0 = (L0-35)/10.0; u = (150-35)/10.0
        t = sfun(u0, u, p)/k
        # fraction of cells fast enough to get there by 60
        kneed = sfun(u0, u, p)/60.0
        frac = norm.sf((math.log(kneed)-mu)/sg)
        print('   inherited %2d: median cell needs %5.0f years; %.1f%% of cells make it by age 60'
              % (L0, t, 100*frac))
