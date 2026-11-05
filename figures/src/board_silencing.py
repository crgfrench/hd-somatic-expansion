"""Build the seventh board: four tests of when the silencing fails.

Inputs are silencing_data.json, written by the four analyses;
hox_cluster_bands.json, the per-cluster transcript counts by repeat-length band
that panel A plots as a fold change over each cluster's own baseline below 100
copies; and hox_bootstrap.json, the donor-level bootstrap intervals that panel A
states the claims against. Output is Silencing.dc.html.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'silencing_data.json')))
PA = json.load(open(os.path.join(HERE, 'hox_cluster_bands.json')))
BOOT = json.load(open(os.path.join(HERE, 'hox_bootstrap.json')))
PROJ = HERE + os.sep

BG, INK, MUT, GRID = '#fbfaf7', '#1d1d1b', '#5c5a55', '#e4dfd4'
BLUE, ORANGE, BAND, GREEN, RED = '#2f6db5', '#d9822b', '#efe8d8', '#1d6e33', '#a3202a'
VIOLET, SLATE = '#63309c', '#8a8679'
CLUS = {'A': '#a3202a', 'B': '#c25a00', 'C': '#1d6e33', 'D': '#1d4f93'}
E = lambda s: s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def tx(x, y, s, sz=12, fill=MUT, anchor='start', w=None, it=False):
    a = ' text-anchor="%s"' % anchor if anchor != 'start' else ''
    b = ' font-weight="%s"' % w if w else ''
    i = ' font-style="italic"' if it else ''
    return '<text x="%.1f" y="%.1f" font-size="%d" fill="%s"%s%s%s>%s</text>' % (x, y, sz, fill, a, b, i, E(s))


def para(x, y, text, size=13, fill=MUT, maxw=1392, lh=None, w=None):
    lh = lh or size + 7
    per = max(10, int(maxw / (0.505 * size)))
    out, line = [], ''
    for word in text.split():
        if len(line) + len(word) + 1 > per:
            out.append(line); line = word
        else:
            line = (line + ' ' + word).strip()
    if line:
        out.append(line)
    return [tx(x, y + i * lh, ln, size, fill, w=w) for i, ln in enumerate(out)], y + len(out) * lh


def ols(X, Y):
    n = len(X); mx = sum(X) / n; my = sum(Y) / n
    sl = sum((a - mx) * (b - my) for a, b in zip(X, Y)) / sum((a - mx) ** 2 for a in X)
    ic = my - sl * mx
    ssr = sum((b - (ic + sl * a)) ** 2 for a, b in zip(X, Y))
    sst = sum((b - my) ** 2 for b in Y)
    return sl, ic, 1 - ssr / sst


def eqn(ylab, xlab, sl, ic):
    return '%s = (%s%.2f × %s %s %.2f)' % (
        ylab, '−' if sl < 0 else '', abs(sl), xlab,
        '−' if ic < 0 else '+', abs(ic))


def svg(W, H, label, body):
    return ('<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="%s" '
            'style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">\n%s\n</svg>'
            % (W, H, W, H, E(label), '\n'.join(body)))


# ---------------------------------------------------------------- panel A
def panel_a():
    W, H = 1392, 700
    L, R, T, B = 110, 856, 100, 372
    A = PA['bands']
    PO = PA['pooled']
    xs = lambda v: L + (v - 36) / (330 - 36) * (R - L)
    lo, hi = 0.35, 220.0
    ys = lambda v: B - (math.log10(max(v, lo)) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * (B - T)
    o = []
    o.append(tx(0, 24, 'Past 200 copies all four HOX clusters run far above their own baseline; across 140 to 199, two of the four do', 19, INK, w=600))
    ls, _ = para(0, 46, 'Each nucleus scored by its own repeat length rather than by its donor, across 78,855 cortical neurons. Each cluster is divided by its own rate below '
                        '100 copies, so the four curves start at 1 and height on the chart is a fold change rather than an abundance.', 14, maxw=1240)
    o += ls
    for v in (0.5, 1, 2, 5, 10, 20, 50, 100, 200):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(v), R, ys(v), GRID))
        o.append(tx(L - 10, ys(v) + 4, ('%g×' % v), 12, MUT, 'end'))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6"/>'
             % (L, ys(1), R, ys(1), SLATE))
    o.append(tx(R - 2, ys(1) - 8, 'baseline below 100 copies', 12, SLATE, 'end'))
    for v in (50, 100, 150, 200, 250, 300):
        o.append(tx(xs(v), B + 24, str(v), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'CAG repeat length in that nucleus', 14, INK, 'middle'))
    o.append(tx(40, (T + B) / 2, 'HOX rate, relative to that cluster below 100 copies', 13, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 40 %.1f)" ' % ((T + B) / 2)))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.5" stroke-dasharray="5 4"/>'
             % (xs(150), T - 22, xs(150), B, SLATE))
    o.append(tx(xs(150) + 8, T - 8, '150 CAGs', 13, INK, w=600))
    for cl in 'ABCD':
        k = 'HOX' + cl
        seg, zer = [], []
        for r in A:
            (zer if r[k + '_u'] == 0 else seg).append(r)
        runs, cur = [], []
        for r in A:
            if r[k + '_u'] == 0:
                if cur:
                    runs.append(cur); cur = []
            else:
                cur.append(r)
        if cur:
            runs.append(cur)
        for run in runs:
            if len(run) > 1:
                pts = ' '.join('%.1f %.1f' % (xs(r['mid']), ys(r[k])) for r in run)
                o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.5"/>' % (pts, CLUS[cl]))
        for r in seg:
            o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' % (xs(r['mid']), ys(r[k]), CLUS[cl]))
        for r in zer:
            o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s" stroke="%s" stroke-width="1.6"/>'
                     % (xs(r['mid']), ys(lo), BG, CLUS[cl]))
    ends = sorted((('HOX' + cl, ys(A[-1]['HOX' + cl]), CLUS[cl]) for cl in 'ABCD'), key=lambda t: t[1])
    placed = []
    for nm, yy, col in ends:
        if placed and yy - placed[-1] < 17:
            yy = placed[-1] + 17
        placed.append(yy)
        o.append(tx(xs(A[-1]['mid']) + 12, yy + 4, nm, 13, col, w=600))
    o.append(tx(L + 16, ys(0.44), 'hollow circles: no HOXB transcript in that band,', 12, CLUS['B']))
    o.append(tx(L + 16, ys(0.44) + 16, 'against the 1.5 its baseline predicts', 12, CLUS['B']))

    # right block: fold change with donor-level uncertainty, which is what the claim rests on
    X0 = 940
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s"/>' % (X0 - 40, T - 34, X0 - 40, B + 110, GRID))
    o.append(tx(X0, T - 12, 'Fold change with donor-level uncertainty', 15, INK, w=600))
    ls, _ = para(X0, T + 8, 'Donors drawn with replacement, then nuclei within each drawn donor, 4,000 times. The '
                            'interval is the middle 95% of those draws.', 12, maxw=410, lh=16)
    o += ls
    bl, br = X0 + 96, 1360
    FLO, FHI = 0.25, 260.0
    fs = lambda v: bl + (math.log10(min(max(v, FLO), FHI)) - math.log10(FLO)) / (math.log10(FHI) - math.log10(FLO)) * (br - bl)
    TOP, BOT = T + 56, B + 58
    for v in (1, 10, 100):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (fs(v), TOP, fs(v), BOT, GRID))
        o.append(tx(fs(v), BOT + 20, '%g×' % v, 12, MUT, 'middle'))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6"/>'
             % (fs(1), TOP, fs(1), BOT, SLATE))
    row = TOP + 4
    for key, lab in (('band_140_199', '140 to 199 copies'), ('band_200_plus', '200 copies and over')):
        o.append(tx(X0, row + 12, lab, 13, INK, w=600))
        row += 24
        for cl in 'ABCD':
            p = BOOT[key]['HOX' + cl]
            o.append(tx(X0, row + 4, 'HOX' + cl, 13, CLUS[cl], w=600))
            o.append(tx(X0 + 86, row + 4, '%d' % p['obs'], 11, MUT, 'end'))
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>'
                     % (fs(p['lo']), row, fs(p['hi']), row, CLUS[cl]))
            o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s"/>' % (fs(p['fold']), row, CLUS[cl]))
            o.append(tx(fs(p['hi']) + 8, row + 4, '%.0f×' % p['fold'] if p['fold'] >= 10 else '%.1f×' % p['fold'],
                        12, CLUS[cl], w=600))
            row += 24
        row += 10
    o.append(tx(X0, BOT + 44, 'Left of each bar: transcripts of that cluster in the window. A bar that runs to', 12, MUT))
    o.append(tx(X0, BOT + 60, 'the axis edge has a lower bound of zero.', 12, MUT))

    ls, y = para(0, B + 134, 'Past 200 copies every cluster is up fiftyfold or more and every interval is clear of baseline. Four polycomb domains lifting together is what a '
                             'failure of the shared silencing system would produce, and the failure is beyond doubt by 200 copies.', 14, INK)
    o += ls
    ls, y = para(0, y + 8, 'In the 140 to 199 window, which is the one containing the 150 copies the original papers identify, the evidence is thinner. HOXA at 3.7 times and '
                           'HOXD at 6.5 are above baseline in 98% and 99% of the resamples, with 6 of the 10 donors contributing transcripts to each. HOXC at 5.3 reaches 87%, on '
                           '3 donors. HOXB rests on 9 transcripts from a single donor and its interval runs from 0 to 28, so whatever HOXB does in this window remains unmeasured.', 13)
    o += ls
    ls, y = para(0, y + 8, 'Single-band movements mean less than they appear to. Below 140 copies the largest is HOXC at 1.8 times baseline, on 11 transcripts from 3 donors, with '
                           'a bootstrap interval of 0.0 to 5.2 and 62% of resamples above 1. The curves show the shape; the intervals to the right are what the claims rest on.', 13)
    o += ls
    return svg(W, H, 'HOX rate against repeat length for the four clusters, as a fold change over each cluster own baseline, cortex', o)


# ---------------------------------------------------------------- panel B
def panel_b():
    W, H = 1392, 496
    L, R, T, B = 110, 900, 96, 320
    Bd = D['B']
    xs = lambda v: L + (v - 100) / (360 - 100) * (R - L)
    ys = lambda v: B - (v + 0.6) / 1.2 * (B - T)
    o = []
    o.append(tx(0, 24, 'At the same repeat length, the neurons further into identity loss are the ones with HOX open', 19, INK, w=600))
    o.append(tx(0, 46, 'Correlation between a nucleus’s HOX level and its phase C scores, computed inside repeat-length bands so length is held constant.', 14))
    for v in (-0.6, -0.3, 0, 0.3, 0.6):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" %s/>'
                 % (L, ys(v), R, ys(v), GRID if v else '#b4ae9d', '' if v else 'stroke-width="1.5"'))
        o.append(tx(L - 10, ys(v) + 4, ('%+.1f' % v).replace('-', '−') if v else '0', 12, MUT, 'end'))
    for v in (100, 150, 200, 250, 300, 350):
        o.append(tx(xs(v), B + 24, str(v), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'CAG repeat length', 14, INK, 'middle'))
    o.append(tx(42, (T + B) / 2, 'rank correlation with HOX', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 42 %.1f)" ' % ((T + B) / 2)))
    for i, (key, col, lab) in enumerate((('Cplus', BLUE, 'genes that rise in phase C'),
                                        ('ctrl', SLATE, 'control gene set'),
                                        ('Cminus', ORANGE, 'identity genes, which fall'))):
        pts = ' '.join('%.1f %.1f' % (xs(r['mid']), ys(r[key])) for r in Bd)
        o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.5"/>' % (pts, col))
        for r in Bd:
            o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' % (xs(r['mid']), ys(r[key]), col))
        ly = T + 18 + i * 22
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="3"/>' % (L + 24, ly, L + 54, ly, col))
        o.append(tx(L + 62, ly + 4, lab, 13, col, w=600))

    p = D['B_partial']
    o.append('<rect x="960" y="88" width="432" height="236" rx="10" fill="#ffffff" stroke="%s"/>' % GRID)
    o.append(tx(984, 118, 'Past 230 copies, holding constant', 14, INK, w=600))
    o.append(tx(984, 138, 'the control set, read depth and', 13, MUT))
    o.append(tx(984, 156, 'repeat length itself', 13, MUT))
    o.append(tx(984, 192, 'genes that rise', 13, INK))
    o.append(tx(1368, 192, '%+.2f' % p['Cplus'], 17, BLUE, 'end', w=700))
    o.append(tx(984, 224, 'identity genes', 13, INK))
    o.append(tx(1368, 224, '%+.2f' % p['Cminus'], 17, ORANGE, 'end', w=700))
    o.append(tx(984, 252, 'both p below 10⁻¹³, on %s nuclei' % f"{p['n']:,}", 13, MUT))
    o.append(tx(984, 290, 'The control line staying flat is what', 13, MUT))
    o.append(tx(984, 308, 'rules out a general drift in quality.', 13, MUT))

    ls, y = para(0, 406, 'Two neurons with the same number of repeats can sit at different points. The one that has lost more of its identity is the one whose HOX clusters have opened, '
                         'and the coupling strengthens as the repeat lengthens.', 14, INK)
    o += ls
    ls, _ = para(0, y + 6, 'So phase C and phase D advance together. A within-band correlation establishes that the two are coupled, and it leaves open whether one process drives '
                           'both. Repeat length sets the pressure; something else sets how far a given cell has traveled under it.', 13)
    o += ls
    return svg(W, H, 'Correlation of HOX expression with the phase C scores, within repeat-length bands', o)


# ---------------------------------------------------------------- panel C
def panel_c():
    W, H = 1392, 556
    C = D['C']
    L, R, T, B = 110, 720, 92, 320
    xs = lambda v: L + (v - 3) / (11 - 3) * (R - L)
    ys = lambda v: B - min(v, 190) / 190 * (B - T)
    o = []
    o.append(tx(0, 24, 'Between cell types inside one person, the repair balance does not predict expansion', 19, INK, w=600))
    ls, _ = para(0, 46, 'The test that was meant to answer the sample size of the FAN1 result, and did not. The balance is measured only in nuclei whose repeat has '
                        'moved less than 10 copies, so it cannot be an effect of the expansion it is predicting.', 14)
    o += ls
    for v in (0, 50, 100, 150):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(v), R, ys(v), GRID))
        o.append(tx(L - 10, ys(v) + 4, str(v), 12, MUT, 'end'))
    for v in (4, 6, 8, 10):
        o.append(tx(xs(v), B + 24, str(v), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'MutSβ to FAN1, measured in barely expanded nuclei', 14, INK, 'middle'))
    o.append(tx(42, (T + B) / 2, 'CAGs added', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 42 %.1f)" ' % ((T + B) / 2)))
    X = [p['x'] for p in C['points']]; Y = [p['y'] for p in C['points']]
    sl, ic, r2 = ols(X, Y)
    for p in C['points']:
        o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s" fill-opacity="0.5" stroke="%s" stroke-width="1"/>'
                 % (xs(p['x']), ys(p['y']), VIOLET, VIOLET))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-dasharray="7 5"/>'
             % (xs(3.3), ys(ic + sl * 3.3), xs(10.3), ys(ic + sl * 10.3), ORANGE))
    o.append(tx(R, T + 22, 'donor S02205, %d subtypes' % C['n'], 14, INK, 'end', w=600))
    o.append(tx(R, T + 44, 'rho %+.2f, p = %.0e' % (C['rho'], C['p']), 15, ORANGE, 'end', w=700))
    o.append(tx(R, T + 66, eqn('CAGs added', 'ratio', sl, ic), 14, INK, 'end'))
    o.append(tx(R, T + 84, 'R² = %.2f' % r2, 13, MUT, 'end'))
    o.append(tx(R, T + 104, 'opposite to the predicted sign', 13, MUT, 'end'))

    # the six donors as a strip
    SL, SR, SY = 860, 1340, 200
    o.append(tx(860, 118, 'Each donor’s correlation, across its own subtypes', 14, INK, w=600))
    o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s"/>' % (SL, SY, SR, SY, GRID))
    zx = lambda r: SL + (r + 1) / 2 * (SR - SL)
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#b4ae9d" stroke-width="1.5"/>' % (zx(0), SY - 26, zx(0), SY + 26))
    for v in (-1, -0.5, 0, 0.5, 1):
        o.append(tx(zx(v), SY + 46, ('%+.1f' % v).replace('-', '−') if v else '0', 12, MUT, 'middle'))
    for dd in C['donors']:
        o.append('<circle cx="%.1f" cy="%d" r="%.1f" fill="%s" fill-opacity="0.55" stroke="%s" stroke-width="1.5"/>'
                 % (zx(dd['rho']), SY, 5 + 4 * math.log10(max(dd['n'], 5)), VIOLET, VIOLET))
    o.append(tx(SL, SY + 74, 'Circle size is the number of subtypes. Positive is the direction the', 13, MUT))
    o.append(tx(SL, SY + 92, 'mechanism predicts; two of six land there, mean %+.2f.'
                % (sum(d['rho'] for d in C['donors']) / len(C['donors'])), 13, MUT))

    ls, y = para(0, 410, 'More MutSβ relative to FAN1 should mean a repeat that travels further. Across cell types of one person it means the opposite, or nothing. '
                         'The association between donors on board 6 therefore holds without a cell-autonomous dose of these enzymes behind it.', 14, INK)
    o += ls
    ls, y = para(0, y + 8, 'Two features of the design limit what a null result here can rule out. Restricting the measurement to barely expanded nuclei samples a different '
                           'fraction of each subtype: in a subtype that expands quickly those nuclei are its stragglers, so the ratio is read off an unrepresentative slice of exactly '
                           'the cell types the comparison is meant to rank. And cell types differ in chromatin state, transcription rate, repeat-containing allele expression and much '
                           'else, so a between-subtype comparison attributes to MutSβ and FAN1 whatever those other differences do.', 13)
    o += ls
    return svg(W, H, 'MutS beta to FAN1 ratio against expansion across cortical subtypes', o)


# ---------------------------------------------------------------- panel D
def panel_d():
    W, H = 1392, 608
    G = D['D_genes']; BARS = D['D_bars']
    L, R, T, B = 110, 620, 152, 386
    lo, hi = 0.02, 20.0
    f = lambda v: (math.log10(max(v, lo)) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    xs = lambda v: L + f(v) * (R - L)
    ys = lambda v: B - f(v) * (B - T)
    o = []
    o.append(tx(0, 24, 'The two regions switch on the same genes, and they do it at different repeat lengths', 19, INK, w=600))
    ls, _ = para(0, 46, 'Left: each of the 42 HOX genes present in both expression matrices, at its level in nuclei at 260 copies or more. These are the 39 genes named '
                        'HOXA to HOXD with a paralog number, plus HOXA10-AS, HOXA11-AS and HOXC13-AS; board 4 instead counts the 23 of them that the cortex paper’s own per-gene '
                        'statistics cover. Right: how much each region has moved from its own baseline.', 14)
    o += ls
    for v in (0.1, 1, 10):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(v), R, ys(v), GRID))
        o.append(tx(L - 10, ys(v) + 4, '%g' % v, 12, MUT, 'end'))
        o.append(tx(xs(v), B + 24, '%g' % v, 12, MUT, 'middle'))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#b4ae9d" stroke-dasharray="5 4"/>'
             % (xs(lo), ys(lo), xs(hi), ys(hi)))
    o.append(tx(xs(0.045) + 10, ys(0.045) + 6, 'equal levels', 12, '#8a8679'))
    LX = [math.log10(max(p['x'], lo)) for p in G['points']]
    LY = [math.log10(max(p['y'], lo)) for p in G['points']]
    gsl, gic, gr2 = ols(LX, LY)
    x1, x2 = math.log10(lo), math.log10(hi)
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-dasharray="7 5"/>'
             % (xs(10 ** x1), ys(10 ** (gic + gsl * x1)), xs(10 ** x2), ys(10 ** (gic + gsl * x2)), ORANGE))
    o.append(tx((L + R) / 2, B + 52, 'cortex, transcripts per million', 14, INK, 'middle'))
    o.append(tx(42, (T + B) / 2, 'striatum, transcripts per million', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 42 %.1f)" ' % ((T + B) / 2)))
    for p in G['points']:
        o.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s" fill-opacity="0.5" stroke="%s" stroke-width="1"/>'
                 % (xs(p['x']), ys(p['y']), GREEN, GREEN))
    for nm, dx, dy, anc in (('HOXD3', 10, -10, 'start'), ('HOXD10', -10, 16, 'end'), ('HOXC4', -10, -8, 'end')):
        p = [q for q in G['points'] if q['g'] == nm][0]
        o.append(tx(xs(p['x']) + dx, ys(p['y']) + dy, nm, 12, INK, anc, w=600))
    o.append(tx(R - 6, B - 56, '%d genes, rho %+.2f' % (G['n'], G['rho']), 15, GREEN, 'end', w=700))
    o.append(tx(R - 6, B - 36, eqn('log striatum', 'log cortex', gsl, gic), 13, INK, 'end'))
    o.append(tx(R - 6, B - 18, 'R² = %.2f, fitted on the logarithms' % gr2, 12, MUT, 'end'))

    # bars
    BL, BR, BB, BT = 760, 1340, 386, 166
    o.append(tx(BL, 132, 'Movement from each region’s own baseline', 15, INK, w=600))
    yb = lambda v: BB - (math.log10(max(v, 0.3)) - math.log10(0.3)) / (math.log10(200) - math.log10(0.3)) * (BB - BT)
    for v in (1, 10, 100):
        o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s"/>' % (BL, yb(v), BR, yb(v), GRID))
        o.append(tx(BL - 10, yb(v) + 4, '%g×' % v, 12, MUT, 'end'))
    o.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#b4ae9d" stroke-width="1.5"/>' % (BL, yb(1), BR, yb(1)))
    gw, x0 = 110, BL + 70
    for i, lab in enumerate(('140-199', '260+')):
        for k, bar in enumerate(BARS):
            x = x0 + i * 280 + k * (gw + 16)
            v = bar[lab]['fold']
            col = BLUE if bar['region'] == 'cortex' else ORANGE
            y1, y2 = (yb(v), yb(1)) if v >= 1 else (yb(1), yb(v))
            o.append('<rect x="%.1f" y="%.1f" width="%d" height="%.1f" rx="3" fill="%s"/>' % (x, y1, gw, max(y2 - y1, 2), col))
            o.append(tx(x + gw / 2, y1 - 8 if v >= 1 else y2 + 16,
                        ('%.0f×' % v) if v >= 10 else ('%.1f×' % v if v >= 1 else '%.2f×' % v),
                        13, INK, 'middle', w=700))
            o.append(tx(x + gw / 2, BB + 22, bar['region'], 12, col, 'middle', w=600))
        o.append(tx(x0 + i * 280 + gw + 8, BB + 44, lab + ' CAGs', 13, INK, 'middle', w=600))

    ls, y = para(0, 474, 'The genes agree: the same clusters lead in both regions, HOXD3, HOXD10 and HOXC4 at the top of each, and the fitted line runs at a slope of 1.01, so close to '
                         'equal levels that the two are hard to separate. The timing does not agree. A cortical neuron opens its HOX clusters as soon as it passes about 150 copies; '
                         'a striatal neuron at the same length has not moved.', 14, INK)
    o += ls
    ls, y = para(0, y + 6, 'The 140 to 199 band holds 18 million striatal transcripts, enough that a rise the size of the cortex’s would have shown about 75 HOX counts where 7 were seen, so this is a '
                           'real difference rather than a thin sample. One reading of it, and a reading rather than a result: the region whose repeats travel fastest tolerates a given length best.')
    o += ls
    ls, _ = para(0, y + 6, 'Both regions reach the same place eventually, 102 times baseline in the cortex and 61 in the striatum past 260 copies.', 13)
    o += ls
    return svg(W, H, 'HOX genes compared between regions, and the thresholds at which each region moves', o)


HEIGHT = 2794
page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>When the silencing fails</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:ital,wght@0,400;0,600;0,700;1,400&amp;display=swap">
<style>
body{{margin:0}}
</style>
</helmet>
<div style="width: 1520px; height: {HEIGHT}px; box-sizing: border-box; padding: 64px; background: {BG}; color: {INK}; font-family: 'IBM Plex Sans', system-ui, sans-serif; display: flex; flex-direction: column; gap: 30px">

<div style="display: flex; flex-direction: column; gap: 10px">
<div style="font-size: 13px; line-height: 18px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase; color: {SLATE}">Board 7 of 7</div>
<h1 style="margin: 0; font-size: 34px; line-height: 42px; font-weight: 700">When the silencing fails, and what sets the timing</h1>
<p style="margin: 0; font-size: 16px; line-height: 24px; color: {MUT}">Four tests on the single-nucleus data released with the two papers, scoring each nucleus by its own repeat length rather than by the donor it came from. Three concern the HOX clusters, the polycomb-silenced genes that board 4 shows derepressing in phase D; the fourth returns to the repair genes of board 6.</p>
</div>

{panel_a()}

{panel_b()}

{panel_c()}

{panel_d()}

<p style="margin: 0; font-size: 12px; line-height: 18px; color: {MUT}">Sources: single-nucleus expression and per-nucleus repeat lengths, Kashin, Lee et al., <i>Cell</i> 2026 (Zenodo record 18202664), cortical excitatory neurons and striatal projection neurons; phase C scores from the same archive. Striatal nuclei are matched to the expression matrix by barcode, since that summary has 6,359 rows against 9,426 matrix columns. Code at github.com/crgfrench/hd-somatic-expansion.</p>

</div>
</x-dc>
<script type="module">
class Component extends DCLogic {{
  renderVals() {{
    return {{}};
  }}
}}
</script>
</body>
</html>
"""
open(PROJ + 'Silencing.dc.html', 'w').write(page)
print('wrote Silencing.dc.html', len(page), 'bytes')
