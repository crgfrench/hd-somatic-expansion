"""Fold change per HOX cluster with a two-level bootstrap over donors then nuclei.

The pooled Poisson test treats every transcript as an independent draw, which
ignores that nuclei differ and that donors differ far more.  Here donors are
resampled with replacement, then nuclei within each drawn donor, and the fold
is recomputed from the resampled band and baseline.
"""
import pandas as pd, numpy as np, json
CL=['HOXA','HOXB','HOXC','HOXD']
d=pd.read_csv('cortex_hox_per_nucleus.csv').dropna(subset=['CAG_LENGTH','umis'])
rng=np.random.default_rng(20261104)
B=4000

def fold(base, band, c):
    bu, nu = base['umis'].sum(), band['umis'].sum()
    if bu<=0 or nu<=0: return np.nan
    r = base[c].sum()/bu
    if r<=0: return np.nan
    return (band[c].sum()/nu)/r

def run(lo, hi, name):
    base_all = d[d.CAG_LENGTH<100]
    band_all = d[(d.CAG_LENGTH>=lo)&(d.CAG_LENGTH<hi)]
    donors = sorted(set(base_all.DONOR) | set(band_all.DONOR))
    bg = {k: v for k, v in base_all.groupby('DONOR')}
    ng = {k: v for k, v in band_all.groupby('DONOR')}
    print('\n=== %s : %d nuclei, %d donors with nuclei in the band'
          % (name, len(band_all), band_all.DONOR.nunique()))
    pt = {c: fold(base_all, band_all, c) for c in CL}
    draws = {c: [] for c in CL}
    for _ in range(B):
        pick = rng.choice(donors, size=len(donors), replace=True)
        bb, nn = [], []
        for p in pick:
            if p in bg:
                s = bg[p]; bb.append(s.iloc[rng.integers(0, len(s), len(s))])
            if p in ng:
                s = ng[p]; nn.append(s.iloc[rng.integers(0, len(s), len(s))])
        if not bb or not nn: continue
        bb = pd.concat(bb); nn = pd.concat(nn)
        for c in CL:
            draws[c].append(fold(bb, nn, c))
    out={}
    print('  %-6s %8s %10s %22s %10s' % ('','fold','transcripts','95% bootstrap interval','P(fold>1)'))
    for c in CL:
        v=np.array([x for x in draws[c] if np.isfinite(x)])
        lo_,hi_=np.percentile(v,[2.5,97.5])
        p=float((v>1).mean())
        obs=int(band_all[c].sum())
        nd=int((band_all.groupby('DONOR')[c].sum()>0).sum())
        out[c]=dict(fold=float(pt[c]), lo=float(lo_), hi=float(hi_), p_gt1=p,
                    obs=obs, donors_with_any=nd, draws=len(v))
        print('  %-6s %8.1f %10d %10.1f to %-9.1f %9.2f   (%d of %d donors contribute any)'
              % (c, pt[c], obs, lo_, hi_, p, nd, band_all.DONOR.nunique()))
    return out

res={'band_140_199': run(140,200,'140 to 199 copies'),
     'band_100_139': run(100,140,'100 to 139 copies'),
     'band_200_plus': run(200,1e9,'200 copies and over')}
res['method']=('two-level nonparametric bootstrap, %d resamples: donors drawn with '
               'replacement, then nuclei within each drawn donor' % B)
json.dump(res, open('hox_bootstrap.json','w'), indent=1)
