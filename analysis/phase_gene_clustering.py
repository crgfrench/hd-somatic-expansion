"""Do the phase C and phase D genes sit in physical clusters along the genome?

Polycomb silences whole domains, so if the phase D crisis is a domain-level
failure the derepressed genes should be bunched together on the chromosome. If
it is gene-by-gene they should be scattered across the same set of expressed
genes.

Three tests, in the order they should be read:

1. GENE FAMILIES. Membership in the families whose names already say they are
   contiguous — HOX, protocadherins, histone clusters, olfactory receptors.
   This is the test the HOX result rests on.

2. SPATIAL, UNMATCHED NULL. How many of the set have another member of the set
   within 1 Mb, against random sets of the same size drawn from every gene in
   the statistics file. THIS NULL IS WRONG and the script reports it only so
   the mistake is visible: a set selected by p-value is heavily biased toward
   well-expressed genes (median near 190 TPM against 17 for all tested genes),
   while clustered paralog families are quiet, so the comparison is not
   like-for-like and it understates clustering.

3. SPATIAL, ABUNDANCE-MATCHED NULL. The same statistic with each random set
   drawn to match the abundance distribution of the set being tested. This
   fixes the bias for phase C, where cortical clustering goes from p = 0.66
   unmatched to p = 0.002 matched.

   It CANNOT be applied to the HOX result. Every one of the 11 phase D HOX
   genes sits below 1 transcript per million, which is the floor of the
   abundance table, so matching drops all of them. Read test 1 for HOX, not
   test 3.

Inputs: the released per-gene statistics for each region, a GRCh38 gene table
with columns symbol, chr and start, and, for the matched null, a per-gene
abundance table measured in that same region. An abundance table from one
region cannot be used to match the other, so the matched test runs only where
its own table is given. Paths are on the command line; see
data/PUT_INPUT_FILES_HERE.txt.

Also writes hox_clusters.json, which figures/src/panel_hox.py draws.
"""
import argparse, json, re, os
import numpy as np
import pandas as pd

CHR = [str(i) for i in range(1, 23)] + ['X', 'Y']
WIN = 1_000_000          # two genes count as neighbours within this distance
NPERM = 10_000
NBIN = 10                # abundance bins for the matched null
SEED = 20261104

FAMS = {
    'HOX clusters':        re.compile(r'^HOX[ABCD]\d'),
    'protocadherins':      re.compile(r'^PCDH[ABG]'),
    'histone clusters':    re.compile(r'^(H[1-4]-\d|H2[AB][A-Z]\d|H[34]C\d|HIST\d)'),
    'olfactory receptors': re.compile(r'^OR\d'),
    'keratin-associated':  re.compile(r'^KRTAP'),
    'late cornified':      re.compile(r'^(LCE\d|SPRR\d)'),
}


def load_annotation(path):
    ann = pd.read_csv(path)
    missing = {'symbol', 'chr', 'start'} - set(ann.columns)
    if missing:
        raise SystemExit('the gene table needs columns %s' % sorted(missing))
    ann = ann[ann.chr.astype(str).isin(CHR)]
    return ann.drop_duplicates('symbol', keep='first').set_index('symbol')[['chr', 'start']]


def clustered(chrcode, start, ii):
    """how many of these genes have another of them within WIN on the same chromosome"""
    n = 0
    for c in np.unique(chrcode[ii]):
        s = np.sort(start[ii[chrcode[ii] == c]])
        if len(s) < 2:
            continue
        left = np.r_[np.inf, np.diff(s)]
        right = np.r_[np.diff(s), np.inf]
        n += int(((left <= WIN) | (right <= WIN)).sum())
    return n


def spatial(universe, genes, ann, rng, weights=None):
    """observed count, null mean, and the chance of reaching the observed count"""
    univ = [g for g in universe if g in ann.index]
    if weights is not None:
        univ = [g for g in univ if g in weights.index]
    pos = ann.loc[univ]
    chrcode = pd.factorize(pos.chr.values)[0]
    start = pos.start.values
    idx = {g: i for i, g in enumerate(univ)}
    ii = np.array([idx[g] for g in genes if g in idx])
    if len(ii) < 2:
        return len(ii), 0, 0.0, 1.0
    obs = clustered(chrcode, start, ii)

    if weights is None:
        pool = np.arange(len(univ))
        null = np.array([clustered(chrcode, start, rng.choice(pool, len(ii), replace=False))
                         for _ in range(NPERM)])
    else:
        ab = np.log10(weights.loc[univ].values)
        binid = pd.qcut(ab, NBIN, labels=False, duplicates='drop')
        bybin = {b: np.flatnonzero(binid == b) for b in np.unique(binid)}
        want = pd.Series(binid[ii]).value_counts().to_dict()
        null = np.empty(NPERM, dtype=int)
        for k in range(NPERM):
            pick = np.concatenate([rng.choice(bybin[b], n, replace=False) for b, n in want.items()])
            null[k] = clustered(chrcode, start, pick)
    return len(ii), obs, float(null.mean()), float((1 + (null >= obs).sum()) / (NPERM + 1))


def families(genes, universe):
    rows = []
    gs, us = set(genes), list(universe)
    for name, pat in FAMS.items():
        bg = [g for g in us if pat.match(g)]
        if not bg:
            continue
        hit = sorted(g for g in bg if g in gs)
        rows.append(dict(family=name, observed=len(hit), tested=len(bg),
                         expected=round(len(bg) * len(genes) / len(us), 2),
                         genes=','.join(hit)))
    return [r for r in rows if r['observed']]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotation', required=True, help='GRCh38 gene table: symbol, chr, start')
    ap.add_argument('--cortex', required=True, help='cortex.stats.txt')
    ap.add_argument('--striatum', required=True, help='caudate.stats.txt')
    ap.add_argument('--abundance-cortex', help='CSV with columns gene and tpm, cortex')
    ap.add_argument('--abundance-striatum', help='CSV with columns gene and tpm, striatum')
    ap.add_argument('--cut', type=float, default=1e-10, help='p-value defining a phase gene set')
    ap.add_argument('--out', default='hox_clusters.json')
    a = ap.parse_args()

    ann = load_annotation(a.annotation)
    ab = {r: (pd.read_csv(f).set_index('gene').tpm if f else None)
          for r, f in (('cortex', a.abundance_cortex), ('striatum', a.abundance_striatum))}
    rng = np.random.default_rng(SEED)
    hox = {}

    for region, path in (('cortex', a.cortex), ('striatum', a.striatum)):
        tpm = ab[region]
        st = pd.read_csv(path, sep='\t')
        universe = st.GENE.tolist()
        print('\n%s\n%s — %d genes tested' % ('=' * 78, region, len(universe)))
        for phase in ('C', 'D'):
            genes = st.loc[st['PHASE_%s_PVALUE' % phase] < a.cut, 'GENE'].tolist()
            print('\nphase %s: %d genes' % (phase, len(genes)))
            for r in families(genes, universe):
                print('  %-20s %2d of %3d tested, %.2f expected  %s'
                      % (r['family'], r['observed'], r['tested'], r['expected'], r['genes']))
            n, obs, exp, p = spatial(universe, genes, ann, rng)
            print('  spatial, unmatched null (biased, see the header): %d of %d within %.0f Mb, '
                  '%.1f expected, p = %.4f' % (obs, n, WIN / 1e6, exp, p))
            if tpm is not None:
                n, obs, exp, p = spatial(universe, genes, ann, rng, weights=tpm)
                covered = sum(g in tpm.index for g in genes)
                print('  spatial, abundance-matched null: %d of %d within %.0f Mb, %.1f expected, '
                      'p = %.4f  [%d of %d genes have an abundance value]'
                      % (obs, n, WIN / 1e6, exp, p, covered, len(genes)))

        dset = set(st.loc[st.PHASE_D_PVALUE < a.cut, 'GENE'])
        cset = set(st.loc[st.PHASE_C_PVALUE < a.cut, 'GENE'])
        hox[region] = [dict(g=g, cl=g[3], pos=int(ann.start[g]), d=g in dset, c=g in cset)
                       for g in universe if FAMS['HOX clusters'].match(g) and g in ann.index]

    json.dump(hox, open(a.out, 'w'))
    print('\nwrote %s' % os.path.abspath(a.out))


if __name__ == '__main__':
    main()
