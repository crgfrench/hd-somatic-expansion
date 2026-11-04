"""Four tests of when the polycomb silencing fails, and what sets the timing.

Each nucleus is scored by its own repeat length rather than by the donor it
came from, which is what makes these different from the donor-level work.

  A  Where along the repeat does each HOX cluster come on? One shared threshold
     across the four clusters is a switch thrown at the level of the silencing
     system; four different ones mean each domain fails on its own.

  B  At a fixed repeat length, is HOX derepression tied to how far that nucleus
     has gone in phase C? The released score file carries a control gene set,
     which is what separates a real coupling from a drift in transcriptome
     quality.

  C  Does the MutSbeta-to-FAN1 balance predict expansion between cell types of
     one person, as it appears to between people? The balance is read only in
     nuclei whose repeat has barely moved, so it cannot be an effect of the
     expansion it is predicting, and subtypes are compared inside a donor.

  D  Do the two regions derepress the same HOX genes, and at the same repeat
     length?

Writes silencing_data.json, which figures/src/board_silencing.py draws.
"""
import argparse, json, re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, poisson

CLUSTERS = 'ABCD'
A_EDGES = [36, 100, 120, 140, 160, 180, 200, 230, 260, 300, 10 ** 9]
A_MID = [68, 110, 130, 150, 170, 190, 215, 245, 280, 340]
B_EDGES = [36, 170, 200, 230, 260, 300, 10 ** 9]
B_MID = [100, 185, 215, 245, 280, 340]
LOW_GAIN, MIN_LOW, MIN_ALL, MIN_SUBTYPES = 10, 40, 60, 5


def cluster_cols(df):
    out = {}
    for cl in CLUSTERS:
        out[cl] = [c for c in df.columns if re.match(r'^HOX%s\d' % cl, c)]
    return out


def add_hox(df):
    cc = cluster_cols(df)
    for cl in CLUSTERS:
        df['HOX' + cl] = df[cc[cl]].sum(axis=1)
    df['HOX'] = df[['HOX' + c for c in CLUSTERS]].sum(axis=1)
    return cc


def resid(y, cols):
    X = np.column_stack([np.ones(len(y))] + [np.asarray(c, float) for c in cols])
    b, *_ = np.linalg.lstsq(X, np.asarray(y, float), rcond=None)
    return np.asarray(y, float) - X @ b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cortex-rows', required=True, help='CSV from extract_gene_rows.py, cortex')
    ap.add_argument('--striatum-rows', required=True, help='CSV from extract_gene_rows.py, striatum')
    ap.add_argument('--cortex-summary', required=True, help='the cortex summary, in matrix column order')
    ap.add_argument('--lengths', required=True, help='per_nucleus_repeat_lengths.csv.gz')
    ap.add_argument('--phasec-scores', required=True, help='per_nucleus_phaseC_scores.csv.gz')
    ap.add_argument('--donor-metadata', required=True)
    ap.add_argument('--out', default='silencing_data.json')
    a = ap.parse_args()
    out = {}

    cx = pd.read_csv(a.cortex_rows)
    sm = pd.read_csv(a.cortex_summary,
                     usecols=['CELL_BARCODE', 'DONOR', 'CAG_LENGTH', 'NUM_TRANSCRIPTS', 'subcluster_name'])
    assert len(sm) == len(cx), 'cortex summary and matrix disagree on nuclei'
    drift = np.abs(cx.umis.values - sm.NUM_TRANSCRIPTS.values).max()
    assert drift == 0, 'cortex columns are not in summary order (largest drift %g)' % drift
    cx = pd.concat([cx.drop(columns=['CELL_BARCODE']), sm[['DONOR', 'CAG_LENGTH', 'subcluster_name']]], axis=1)
    add_hox(cx)

    # the striatal summary has fewer rows than the matrix has columns, so match by barcode
    st = pd.read_csv(a.striatum_rows)
    ln = pd.read_csv(a.lengths)
    spn = ln[ln.SOURCE.str.contains('caudate/SPN', na=False)][['CELL_BARCODE', 'DONOR', 'CAG_LENGTH']]
    st = st.merge(spn, on='CELL_BARCODE', how='inner')
    print('striatum: %d of the matrix columns carry a measured repeat length' % len(st))
    add_hox(st)

    # ---- A
    cx['bin'] = pd.cut(cx.CAG_LENGTH, A_EDGES, right=False, labels=False)
    rows = []
    for b, g in cx.groupby('bin', observed=True):
        r = {'mid': A_MID[int(b)], 'n': int(len(g))}
        for cl in CLUSTERS:
            r['HOX' + cl] = float(g['HOX' + cl].sum() / g.umis.sum() * 1e6)
        r['all'] = float(g.HOX.sum() / g.umis.sum() * 1e6)
        rows.append(r)
    out['A'] = rows
    print('\nA: HOX transcripts per million by repeat length, cortex')
    print(pd.DataFrame(rows).round(2).to_string(index=False))

    # ---- B
    sc = pd.read_csv(a.phasec_scores)
    assert len(sc) == len(cx), 'the score file and the matrix disagree on nuclei'
    same = (sc.DONOR.values == cx.DONOR.values) & (np.abs(sc.UMIs.values - cx.umis.values) < 1e-9)
    assert same.all(), 'the score file is not in matrix order'
    d = pd.concat([cx, sc[['Cplus', 'Cminus', 'ctrl']]], axis=1)
    d['hox'] = d.HOX / d.umis * 1e6
    d['bin'] = pd.cut(d.CAG_LENGTH, B_EDGES, right=False, labels=False)
    out['B'] = [dict(mid=B_MID[int(b)], n=int(len(g)),
                     Cminus=float(spearmanr(g.hox, g.Cminus).statistic),
                     Cplus=float(spearmanr(g.hox, g.Cplus).statistic),
                     ctrl=float(spearmanr(g.hox, g.ctrl).statistic))
                for b, g in d.groupby('bin', observed=True)]
    print('\nB: rank correlation of HOX with the phase C scores, inside repeat-length bands')
    print(pd.DataFrame(out['B']).round(3).to_string(index=False))
    g = d[d.CAG_LENGTH >= 230]
    part = {}
    for target in ('Cminus', 'Cplus'):
        ry = resid(rankdata(g.hox), [rankdata(g.ctrl), rankdata(g.umis), rankdata(g.CAG_LENGTH)])
        rx = resid(rankdata(g[target]), [rankdata(g.ctrl), rankdata(g.umis), rankdata(g.CAG_LENGTH)])
        r, p = spearmanr(ry, rx)
        part[target] = float(r)
        print('   past 230 copies, holding the control set, depth and length: %-7s rho %+.3f  p = %.2g'
              % (target, r, p))
    out['B_partial'] = dict(n=int(len(g)), **part)

    # ---- C
    meta = pd.read_csv(a.donor_metadata, sep='\t')
    meta = meta[meta.STATUS == 'pwHD'][['DONOR', 'HD_CAG_LENGTH']]
    c2 = cx.merge(meta, on='DONOR', how='left')
    c2['gain'] = c2.CAG_LENGTH - c2.HD_CAG_LENGTH
    c2 = c2[c2.gain.notna() & c2.subcluster_name.notna()]
    rows = []
    for (dn, sub), g in c2.groupby(['DONOR', 'subcluster_name']):
        low = g[g.gain < LOW_GAIN]
        if len(g) < MIN_ALL or len(low) < MIN_LOW:
            continue
        rows.append(dict(donor=dn, n=len(g),
                         expansion=float(np.quantile(g.gain, 0.9)),
                         ratio=(low.MSH2.sum() + low.MSH3.sum()) / max(low.FAN1.sum(), 1)))
    t = pd.DataFrame(rows)
    per = [dict(d=dn, n=int(len(x)), rho=float(spearmanr(x.ratio, x.expansion).statistic))
           for dn, x in t.groupby('donor') if len(x) >= MIN_SUBTYPES]
    big = t[t.donor == max(per, key=lambda r: r['n'])['d']]
    r, p = spearmanr(big.ratio, big.expansion)
    out['C'] = dict(points=[dict(x=float(x), y=float(y)) for x, y in zip(big.ratio, big.expansion)],
                    rho=float(r), p=float(p), n=int(len(big)), donors=per)
    print('\nC: MutSbeta-to-FAN1 against expansion, across subtypes inside each donor')
    for q in per:
        print('   %-8s %3d subtypes  rho %+.2f' % (q['d'], q['n'], q['rho']))
    print('   mean rho %+.3f over %d donors, positive in %d'
          % (np.mean([q['rho'] for q in per]), len(per), sum(q['rho'] > 0 for q in per)))

    # ---- D
    cc = cluster_cols(cx)
    hox_genes = [g for cl in CLUSTERS for g in cc[cl]]
    tab = {}
    for name, df in (('cortex', cx), ('striatum', st)):
        lo, hi = df[df.CAG_LENGTH < 100], df[df.CAG_LENGTH >= 260]
        tab[name] = {g: float(hi[g].sum() / hi.umis.sum() * 1e6) for g in hox_genes if g in df.columns}
    shared = sorted(set(tab['cortex']) & set(tab['striatum']))
    r, p = spearmanr([tab['cortex'][g] for g in shared], [tab['striatum'][g] for g in shared])
    out['D_genes'] = dict(points=[dict(g=g, x=tab['cortex'][g], y=tab['striatum'][g]) for g in shared],
                          rho=float(r), n=len(shared))
    print('\nD: %d HOX genes in both regions, levels in expanded nuclei correlate at rho %+.2f, p = %.2g'
          % (len(shared), r, p))
    bars = []
    for name, df in (('cortex', cx), ('striatum', st)):
        lo = df[df.CAG_LENGTH < 100]
        base = lo.HOX.sum() / lo.umis.sum()
        e = {}
        for lab, sel in (('140-199', (df.CAG_LENGTH >= 140) & (df.CAG_LENGTH < 200)),
                         ('260+', df.CAG_LENGTH >= 260)):
            q = df[sel]
            obs = float(q.HOX.sum())
            exp = base * q.umis.sum()
            e[lab] = dict(fold=obs / exp, n=int(len(q)), obs=int(obs), exp=float(exp),
                          p=float(poisson.sf(obs - 1, exp) if obs > exp else poisson.cdf(obs, exp)))
            print('   %-9s %-8s %5d nuclei, %5d HOX against %6.1f expected, %.2fx, p = %.2g'
                  % (name, lab, e[lab]['n'], e[lab]['obs'], e[lab]['exp'], e[lab]['fold'], e[lab]['p']))
        bars.append(dict(region=name, base_tpm=float(base * 1e6), **e))
    out['D_bars'] = bars

    json.dump(out, open(a.out, 'w'))
    print('\nwrote %s' % a.out)


if __name__ == '__main__':
    main()
