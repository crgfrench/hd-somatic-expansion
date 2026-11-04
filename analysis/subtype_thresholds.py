import csv, math, os, sys, statistics as st

DATA = (sys.argv[1].rstrip('/') + '/' if len(sys.argv) > 1
     else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')

CAND = ["MBNL1","MBNL2","CREBBP","EP300","EZH2","EED","SUZ12","TCERG1","MED15","HTT","MSH3"]
rows = list(csv.DictReader(open(DATA + 'subtype_thresholds.csv')))

subs, order = {}, []
for r in rows:
    s = r['SUBTYPE']
    if s not in subs:
        subs[s] = {k: float(r[k]) for k in ('n','n_high','median_UMIs','onset_up','up_lo','up_hi',
                                            'up_per_100','onset_dn','dn_lo','dn_hi','dn_per_100')}
        subs[s]['base'] = []
        order.append(s)
    subs[s]['base'].append(float(r['V13']))

bad = [s for s in order if len(subs[s]['base']) != 11]
print('subtypes:', len(order), '| rows per subtype not 11:', bad)
for s in order:
    for g, v in zip(CAND, subs[s]['base']):
        subs[s][g] = v

def rank(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0]*len(v); i = 0
    while i < len(idx):
        j = i
        while j < len(idx) and v[idx[j]] == v[idx[i]]: j += 1
        for k in range(i, j): r[idx[k]] = (i+j+1)/2
        i = j
    return r

def spearman(x, y):
    rx, ry = rank(x), rank(y)
    mx, my = st.mean(rx), st.mean(ry)
    sx = math.sqrt(sum((a-mx)**2 for a in rx)); sy = math.sqrt(sum((b-my)**2 for b in ry))
    if sx == 0 or sy == 0: return float('nan'), float('nan')
    rho = sum((a-mx)*(b-my) for a, b in zip(rx, ry))/(sx*sy)
    n = len(x)
    if abs(rho) >= 1: return rho, 0.0
    t = rho*math.sqrt((n-2)/(1-rho*rho))
    # two-sided t with n-2 df, via incomplete beta
    df = n-2; xb = df/(df+t*t)
    def betacf(a, b, x):
        MAXIT, EPS, FPMIN = 200, 3e-16, 1e-300
        qab, qap, qam = a+b, a+1, a-1
        c = 1.0; d = 1-qab*x/qap
        if abs(d) < FPMIN: d = FPMIN
        d = 1/d; h = d
        for m in range(1, MAXIT+1):
            m2 = 2*m
            aa = m*(b-m)*x/((qam+m2)*(a+m2))
            d = 1+aa*d; c = 1+aa/c
            if abs(d) < FPMIN: d = FPMIN
            if abs(c) < FPMIN: c = FPMIN
            d = 1/d; h *= d*c
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
            d = 1+aa*d; c = 1+aa/c
            if abs(d) < FPMIN: d = FPMIN
            if abs(c) < FPMIN: c = FPMIN
            d = 1/d; de = d*c; h *= de
            if abs(de-1) < EPS: break
        return h
    def betai(a, b, x):
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        lb = math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log(1-x)
        bt = math.exp(lb)
        return bt*betacf(a, b, x)/a if x < (a+1)/(a+b+2) else 1-bt*betacf(b, a, 1-x)/b
    return rho, betai(df/2, 0.5, xb)

print('\nsubtype      n  high  onset_up (95%)      onset_dn (95%)')
for s in sorted(order, key=lambda s: subs[s]['onset_up']):
    d = subs[s]
    print('%-16s %5d %5d   %3d (%3d-%3d)   %3d (%3d-%3d)' % (s, d['n'], d['n_high'],
          d['onset_up'], d['up_lo'], d['up_hi'], d['onset_dn'], d['dn_lo'], d['dn_hi']))

print('\ncorrelation across subtypes (Spearman)')
res = []
for g in CAND + ['median_UMIs', 'n_high']:
    for o in ('onset_up', 'onset_dn'):
        x = [subs[s][g] for s in order]; y = [subs[s][o] for s in order]
        rho, p = spearman(x, y)
        res.append((p, g, o, rho))
for p, g, o, rho in sorted(res):
    print('  %-12s vs %-9s rho = %+.2f   p = %.3f' % (g, o, rho, p))
n = len(order)
print('\n%d subtypes; a two-sided Spearman test needs |rho| of about %.2f to reach p = 0.05 here'
      % (n, 1.96/math.sqrt(n-1)))

print('\nextra checks')
frac = {s: subs[s]['n_high']/subs[s]['n'] for s in order}
for o in ('onset_up','onset_dn'):
    rho,p = spearman([frac[s] for s in order],[subs[s][o] for s in order])
    print('  share of nuclei above 100 CAGs vs %-9s rho = %+.2f  p = %.3f' % (o,rho,p))
for pre in ('ULIT','DLIT'):
    g=[s for s in order if s.startswith(pre)]
    print('  %s: %d subtypes, median onset_up %.0f, median onset_dn %.0f, median share above 100 CAGs %.2f'
          % (pre,len(g),st.median([subs[s]['onset_up'] for s in g]),
             st.median([subs[s]['onset_dn'] for s in g]), st.median([frac[s] for s in g])))
w=[subs[s]['up_hi']-subs[s]['up_lo'] for s in order]
print('  width of the 95%% interval on onset_up: median %.0f CAGs, range %.0f-%.0f' % (st.median(w),min(w),max(w)))
