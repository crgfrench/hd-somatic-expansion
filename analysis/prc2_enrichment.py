import os, sys
import csv, math, statistics as st, sys
import pandas as pd

U = (sys.argv[1].rstrip('/') + '/' if len(sys.argv) > 1
     else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data') + '/')

# --- reference gene sets from von Schimmelmann 2016 (mouse symbols -> uppercase)
t1 = pd.read_excel(U + '41593_2016_BFnn4360_MOESM3_ESM.xls', header=None, skiprows=3)
K27 = set(str(s).upper() for s in t1[6].dropna())

t3 = pd.read_excel(U + '41593_2016_BFnn4360_MOESM5_ESM.xlsx', sheet_name='Gene list', header=None)
blocks, cur = {}, None
for r in t3.values.tolist():
    s = str(r[0])
    if 'genes at' in s:
        cur = s; blocks[cur] = []; continue
    if cur and isinstance(r[5], str) and r[5] != 'Gene Symbol':
        blocks[cur].append((r[5].upper(), r[4]))
UP6 = set(g for g, reg in blocks[[k for k in blocks if '6 months' in k][0]] if reg == 'up')
UP3 = set(g for g, reg in blocks[[k for k in blocks if '3 months' in k][0]] if reg == 'up')
print('reference sets: K27-marked %d, knockout-up 6mo %d, 3mo %d' % (len(K27), len(UP6), len(UP3)))

def lf(n): return math.lgamma(n + 1)

def fisher(a, b, c, d):
    n = a + b + c + d; r1 = a + b; c1 = a + c
    def p(x): return math.exp(lf(r1) + lf(n - r1) + lf(c1) + lf(n - c1) - lf(n) - lf(x) - lf(r1 - x) - lf(c1 - x) - lf(n - r1 - c1 + x))
    lo = max(0, r1 + c1 - n); hi = min(r1, c1); po = p(a)
    return min(1.0, sum(p(x) for x in range(lo, hi + 1) if p(x) <= po * (1 + 1e-9)))

def num(x):
    try: return float(x)
    except Exception: return None

results = {}
for reg, f in (('striatum', 'caudate'), ('cortex', 'cortex')):
    rows = list(csv.DictReader(open(U + f + '.stats.txt'), delimiter='\t'))
    bg = set(r['GENE'].upper() for r in rows)
    sets = {}
    for r in rows:
        g = r['GENE'].upper()
        pc, ec = num(r['PHASE_C_PVALUE']), num(r['PHASE_C_EFFECT'])
        pd_, ed = num(r['PHASE_D_PVALUE']), num(r['PHASE_D_EFFECT'])
        if pc is not None and pc < 1e-10:
            sets.setdefault('C+' if ec > 0 else 'C-', set()).add(g)
        if pd_ is not None and pd_ < 1e-10:
            sets.setdefault('D+' if ed > 0 else 'D-', set()).add(g)
    results[reg] = (bg, sets)
    print('\n== %s | background %d genes' % (reg, len(bg)))
    for refname, ref in (('knockout switched-on (190)', UP6), ('K27-marked (2057)', K27)):
        R = ref & bg
        print('  %s -> %d present in these data' % (refname, len(R)))
        for s in ('C+', 'C-', 'D+', 'D-'):
            S = sets.get(s, set())
            if not S: continue
            a = len(S & R); b = len(S) - a; c = len(R) - a; d = len(bg) - a - b - c
            exp = len(S) * len(R) / len(bg)
            orr = (a * d) / (b * c) if b * c else float('inf')
            print('    %-3s n=%-4d overlap=%-4d expected=%5.1f  OR=%6.2f  p=%.2g' % (s, len(S), a, exp, orr, fisher(a, b, c, d)))
            if s in ('C+', 'C-') and refname.startswith('knockout') and a:
                print('        ', sorted(S & R))

# --- onsets of K27-marked phase C genes vs the rest
def mannwhitney(x, y):
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    vals = [v for v, _ in allv]; rk = [0] * len(allv); i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and vals[j] == vals[i]: j += 1
        for k in range(i, j): rk[k] = (i + j + 1) / 2
        i = j
    R1 = sum(r for r, (v, g) in zip(rk, allv) if g == 0)
    n1, n2 = len(x), len(y)
    U_ = R1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2; sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (U_ - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))

on = list(csv.DictReader(open(U + 'phaseC_gene_onsets.csv')))
print('\n== onset of K27-marked phase C genes vs other phase C genes')
for reg in ('striatum', 'cortex'):
    for s in ('C+', 'C-'):
        sel = [r for r in on if r['region'] == reg and r['SET'] == s]
        a = [float(r['onset']) for r in sel if r['GENE'].upper() in K27]
        b = [float(r['onset']) for r in sel if r['GENE'].upper() not in K27]
        if len(a) > 2:
            print('  %s %s: marked n=%d median=%.0f mean=%.0f | others n=%d median=%.0f mean=%.0f | p=%.3g'
                  % (reg, s, len(a), st.median(a), st.mean(a), len(b), st.median(b), st.mean(b), mannwhitney(a, b)))

print('\n== control: are the cortex phase-C knockout genes also phase D genes?')
bgc, setc = results['cortex']
ten = sorted(setc['C+'] & (UP6 & bgc))
onm = {(r['region'], r['GENE'].upper()): r for r in on}
for g in ten:
    r = onm.get(('cortex', g))
    print('  %-9s alsoD+=%-5s onset=%-5s ci=%s-%s' % (g, g in setc.get('D+', set()),
          r['onset'] if r else 'na', r['ci_low'] if r else '', r['ci_high'] if r else ''))
print('  of %d: %d also phase D+' % (len(ten), sum(g in setc.get('D+', set()) for g in ten)))
bgs, sets_s = results['striatum']
print('  striatum C+ knockout genes also D+:', sorted(sets_s['C+'] & UP6 & bgs & sets_s.get('D+', set())))

print('\n== repeat, with every phase D gene removed from the phase C sets')
for reg in ('striatum','cortex'):
    bg, sets = results[reg]
    D = sets.get('D+', set()) | sets.get('D-', set())
    for refname, ref in (('knockout switched-on', UP6), ('K27-marked', K27)):
        R = (ref & bg) - D
        bg2 = bg - D
        for s in ('C+','C-'):
            S = (sets.get(s,set()) - D)
            a=len(S&R); b=len(S)-a; c=len(R)-a; d=len(bg2)-a-b-c
            exp=len(S)*len(R)/len(bg2); orr=(a*d)/(b*c) if b*c else float('inf')
            print('  %-8s %-20s %-3s n=%-4d overlap=%-3d expected=%5.1f OR=%6.2f p=%.2g'
                  % (reg, refname, s, len(S), a, exp, orr, fisher(a,b,c,d)))
print('\n== onset test, phase D genes removed')
for reg in ('striatum','cortex'):
    bg, sets = results[reg]
    D = sets.get('D+', set()) | sets.get('D-', set())
    for s in ('C+','C-'):
        sel=[r for r in on if r['region']==reg and r['SET']==s and r['GENE'].upper() not in D]
        a=[float(r['onset']) for r in sel if r['GENE'].upper() in K27]
        b=[float(r['onset']) for r in sel if r['GENE'].upper() not in K27]
        if len(a)>2:
            print('  %s %s: marked n=%d median=%.0f | others n=%d median=%.0f | p=%.3g'
                  % (reg,s,len(a),st.median(a),len(b),st.median(b),mannwhitney(a,b)))

print('\n== numbers for the figure')
for reg in ('striatum','cortex'):
    bg, sets = results[reg]
    R6 = UP6 & bg
    for s in ('C+','D+'):
        S = sets.get(s,set()); a=len(S&R6)
        print('  %s %s: %d of %d = %.2f%% | background %d of %d = %.2f%%' % (reg,s,a,len(S),100*a/len(S),len(R6),len(bg),100*len(R6)/len(bg)))
bg, sets = results['cortex']; D = sets.get('D+',set())|sets.get('D-',set())
sel=[r for r in on if r['region']=='cortex' and r['SET']=='C+' and r['GENE'].upper() not in D]
a=sorted(float(r['onset']) for r in sel if r['GENE'].upper() in K27)
b=sorted(float(r['onset']) for r in sel if r['GENE'].upper() not in K27)
def q(v,p):
    i=(len(v)-1)*p; lo=int(i); return v[lo]+(v[min(lo+1,len(v)-1)]-v[lo])*(i-lo)
print('  cortex C+ onsets, PRC2-marked n=%d: q1=%.0f med=%.0f q3=%.0f  %s' % (len(a),q(a,.25),q(a,.5),q(a,.75),a))
print('  cortex C+ onsets, others   n=%d: q1=%.0f med=%.0f q3=%.0f' % (len(b),q(b,.25),q(b,.5),q(b,.75)))
