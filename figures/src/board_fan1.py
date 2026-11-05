import os
HERE = os.path.dirname(os.path.abspath(__file__))
import json, math
D = json.load(open(os.path.join(HERE, 'repair_genes.json')))
PROJ = HERE + os.sep

BG, INK, MUT, GRID = '#fbfaf7', '#1d1d1b', '#5c5a55', '#e4dfd4'
BLUE, ORANGE, BAND, GREEN = '#2f6db5', '#d9822b', '#efe8d8', '#4a8a5c'
RED, SLATE = '#a3202a', '#8a8679'
E = lambda s: (s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))

# Same palette and the same donor-to-color assignment as the timing board, so a
# donor keeps its color across the two boards.
DONOR_HUES = ['#a3202a', '#f08076',
              '#c25a00', '#e89a2e',
              '#d9c400',
              '#1d6e33', '#62c566',
              '#1d4f93', '#5fb0ee',
              '#63309c', '#b074dd']
CX = sorted(D['cortex'], key=lambda z: (z['inh'], z['d']))
COLOR = {d['d']: DONOR_HUES[i] for i, d in enumerate(CX)}

# Regional pseudobulk expression, transcripts per million, from the same matrices.
TPM = {'cortex':   {'MSH2': 80.51, 'MSH3': 114.66, 'FAN1': 31.26},
       'striatum': {'MSH2': 143.81, 'MSH3': 197.63, 'FAN1': 29.05}}
RATIO = {k: (v['MSH2'] + v['MSH3']) / v['FAN1'] for k, v in TPM.items()}

FAN1_RHO = D['fan1_rho']
MSH2_RHO = D['msh2_rho']
POOL_MAX = D['pool_max']
POOL_N = sum(D['matched'])
DEEP = '#7d7768'      # the abundance-matched subset, stacked inside each bar
NGENE = D['ngene']


def ols(X, Y):
    n = len(X); mx = sum(X) / n; my = sum(Y) / n
    sl = sum((a - mx) * (b - my) for a, b in zip(X, Y)) / sum((a - mx) ** 2 for a in X)
    ic = my - sl * mx
    ssr = sum((b - (ic + sl * a)) ** 2 for a, b in zip(X, Y))
    sst = sum((b - my) ** 2 for b in Y)
    return sl, ic, 1 - ssr / sst


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(X, Y):
    rx, ry = ranks(X), ranks(Y)
    return ols(rx, ry)[0] * (
        math.sqrt(sum((a - sum(rx) / len(rx)) ** 2 for a in rx) /
                  sum((b - sum(ry) / len(ry)) ** 2 for b in ry)))


def eqn(ylab, xlab, sl, ic):
    return '%s = (%s%.2f × %s %s %.1f)' % (
        ylab, '−' if sl < 0 else '', abs(sl), xlab,
        '−' if ic < 0 else '+', abs(ic))


def tx(x, y, s, sz=12, fill=MUT, anchor='start', w=None, it=False):
    a = ' text-anchor="%s"' % anchor if anchor != 'start' else ''
    b = ' font-weight="%s"' % w if w else ''
    i = ' font-style="italic"' if it else ''
    return '<text x="%.1f" y="%.1f" font-size="%d" fill="%s"%s%s%s>%s</text>' % (
        x, y, sz, fill, a, b, i, E(s))


def para(x, y, text, size=14, fill=MUT, maxw=1392, lh=None, w=None):
    """lay a paragraph out in lines that fit maxw, returning (lines, next y)"""
    lh = lh or size + 10
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


def svg(W, H, label, body):
    return ('<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="%s" '
            'style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">\n%s\n</svg>'
            % (W, H, W, H, E(label), '\n'.join(body)))


# ---------------------------------------------------------------- panel 1
def tract(x, y, n, col, h=15, w=9, gap=3):
    """a CAG tract drawn as n triplet blocks"""
    o = []
    for i in range(n):
        o.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="1.5" fill="%s"/>'
                 % (x + i * (w + gap), y, w, h, col))
    return o, x + n * (w + gap) - gap


def panel_mechanism():
    W, H = 1392, 500
    o = []
    o.append(tx(0, 24, 'Two enzymes act on the same intermediate, and they pull in opposite directions', 19, INK, w=600))
    ls, _ = para(0, 46, 'A repeat tract being copied can slip, leaving a stretch of CAGs looped out of the duplex. '
                        'What happens to that loop decides whether the repeat ends longer or unchanged.')
    o += ls

    # --- the shared intermediate
    cx = 696
    o.append('<rect x="%.1f" y="88" width="300" height="52" rx="6" fill="%s"/>' % (cx - 150, BAND))
    blocks, _ = tract(cx - 128, 106, 8, SLATE)
    o += blocks
    o.append('<path d="M %.1f 106 q 10 -44 34 -44 q 24 0 34 44" fill="none" stroke="%s" stroke-width="3"/>'
             % (cx + 10, SLATE))
    lb, _ = tract(cx + 20, 62, 3, '#b4ae9d', h=11, w=7, gap=3)
    o += lb
    o.append(tx(cx, 160, 'CAGs looped out during repair', 13, INK, 'middle', w=600))

    o.append('<path d="M %.1f 174 L %.1f 206 M %.1f 174 L %.1f 206" stroke="%s" stroke-width="2" fill="none"/>'
             % (cx - 10, 390, cx + 10, 1002, MUT))
    o.append('<path d="M 390 206 l -6 -11 l 12 0 z" fill="%s"/>' % MUT)
    o.append('<path d="M 1002 206 l -6 -11 l 12 0 z" fill="%s"/>' % MUT)

    # --- left branch: MutS beta
    for x0, col, head, sub, n_end, tail in [
            (60, RED, 'MutSβ nicks the opposite strand',
             'MSH2 and MSH3 together, with MutLγ. The gap is filled using the loop as the template, '
             'so the extra CAGs are written into the sequence.', 11,
             'The repeat ends longer than it started.'),
            (732, GREEN, 'FAN1 cuts the looped strand',
             'A structure-specific nuclease, acting with PCNA and RFC. The loop is removed and the '
             'strand resynthesized, so the extra CAGs never enter the sequence.', 8,
             'The repeat returns to its original length.')]:
        o.append('<rect x="%.1f" y="216" width="600" height="212" rx="10" fill="#ffffff" stroke="%s"/>' % (x0, GRID))
        o.append('<rect x="%.1f" y="216" width="600" height="4" rx="2" fill="%s"/>' % (x0, col))
        o.append(tx(x0 + 24, 250, head, 16, col, w=600))
        words = sub.split()
        line, lines = '', []
        for wd in words:
            if len(line) + len(wd) > 72:
                lines.append(line); line = wd
            else:
                line = (line + ' ' + wd).strip()
        lines.append(line)
        for i, ln in enumerate(lines):
            o.append(tx(x0 + 24, 274 + i * 19, ln, 13, MUT))
        blocks, end = tract(x0 + 24, 352, n_end, col)
        o += blocks
        o.append(tx(x0 + 24, 404, tail, 14, INK, w=600))
        o.append(tx(end + 14, 365, '%d triplets' % n_end, 12, MUT))

    ls, y = para(0, 462, 'The two interact: FAN1 also inhibits MutLγ, the nuclease MutSβ needs to make the nick, so it restrains the '
                         'expansion step as well as removing the loop. The balance between them sets the direction a cell drifts, and neither enzyme sets it alone.', 14, INK)
    o += ls
    ls, _ = para(0, y + 10, 'Drawn from work on purified proteins and synthetic substrates. Triplet counts are illustrative.', 13)
    o += ls
    return svg(W, H, 'MutS beta and FAN1 acting on a slipped-out CAG loop, with opposite outcomes', o)


# ---------------------------------------------------------------- panel 2
def panel_cortex_scatter():
    W, H = 1392, 520
    L, R, T, B = 110, 1150, 100, 380
    xs = lambda v: L + v / 190.0 * (R - L)
    ys = lambda v: B - (v - 24) / (42 - 24) * (B - T)
    o = []
    o.append(tx(0, 24, 'Across the 11 cortex donors, the ones whose repeats traveled furthest have the least FAN1', 19, INK, w=600))
    ls, _ = para(0, 46, 'FAN1 transcripts per million in cortical excitatory neurons, standardized so every donor is compared at the same '
                        'distribution of repeat lengths. Expansion is the 90th percentile repeat length in that donor, less the length inherited.')
    o += ls
    for v in range(24, 43, 3):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(v), R, ys(v), GRID))
        o.append(tx(L - 10, ys(v) + 4, str(v), 12, MUT, 'end'))
    for v in range(0, 191, 30):
        o.append(tx(xs(v), B + 24, str(v), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'CAGs added in cortical neurons', 14, INK, 'middle'))
    o.append(tx(44, (T + B) / 2, 'FAN1, transcripts per million', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 44 %.1f)" ' % ((T + B) / 2)))
    X = [d['exp'] for d in CX]; Y = [d['fan1'] for d in CX]
    sl, ic, r2 = ols(X, Y)
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-dasharray="7 5"/>'
             % (xs(0), ys(ic), xs(190), ys(ic + sl * 190), ORANGE))
    for d in sorted(CX, key=lambda z: -z['exp']):
        o.append('<circle cx="%.1f" cy="%.1f" r="7.5" fill="%s" stroke="%s" stroke-width="2"/>'
                 % (xs(d['exp']), ys(d['fan1']), COLOR[d['d']], BG))
    o.append(tx(R - 16, T + 34, 'rho = %+.2f,  p = 7 × 10⁻⁵' % FAN1_RHO, 15, INK, 'end', w=600))
    o.append(tx(R - 16, T + 56, eqn('FAN1', 'CAGs added', sl, ic), 14, INK, 'end'))
    o.append(tx(R - 16, T + 76, 'R² = %.2f' % r2, 13, MUT, 'end'))

    o.append('<rect x="1190" y="88" width="202" height="292" rx="10" fill="#ffffff" stroke="%s"/>' % GRID)
    o.append(tx(1206, 116, 'Holding the inherited', 13, INK, w=600))
    o.append(tx(1206, 134, 'length constant', 13, INK, w=600))
    o.append(tx(1206, 160, 'rho = −0.82', 18, ORANGE, w=700))
    o.append(tx(1206, 182, 'p = 0.004', 13, MUT))
    o.append(tx(1206, 216, 'Longer inherited repeats', 13, MUT))
    o.append(tx(1206, 234, 'both expand faster and', 13, MUT))
    o.append(tx(1206, 252, 'sit lower in FAN1, so the', 13, MUT))
    o.append(tx(1206, 270, 'partial figure matters.', 13, MUT))
    o.append(tx(1206, 304, 'It is the rank correlation', 13, MUT))
    o.append(tx(1206, 322, 'within each inherited', 13, MUT))
    o.append(tx(1206, 340, 'length, pooled.', 13, MUT))

    ls, y = para(0, 456, 'The direction is the one the biology predicts: less of the enzyme that removes loops, further the repeat gets. '
                         'Which way the arrow runs is left open by a figure of this kind.', 14, INK)
    o += ls
    ls, y = para(0, y + 6, 'If the growth of the repeat, or the damage that follows it, lowers FAN1 in the neurons furthest along, the same pattern appears '
                           'with the arrow running the other way. Only material sampled before the expansion could tell the two apart.')
    o += ls
    ls, _ = para(0, y + 6, 'Each donor keeps the color it has on the timing board.', 13)
    o += ls
    return svg(W, H, 'FAN1 expression against expansion achieved in eleven cortex donors', o)


# ---------------------------------------------------------------- panel 3
def panel_null():
    W, H = 1392, 512
    L, R, T, B = 110, 1340, 96, 330
    xs = lambda v: L + (v + 1) / 2.0 * (R - L)
    # With 11 donors a Spearman correlation is built from an even sum of squared
    # rank differences, so it can only land on 221 values spaced 1/110 apart. Bars
    # of an arbitrary width hold 5 of those values and then 6, which alternates the
    # heights; these bars hold 5 each.
    PER, STEP = 5, 1.0 / 110
    nb = 220 // PER + 1
    hist, sub = [0] * nb, [0] * nb
    for r, m in zip(D['rho'], D['matched']):
        j = min(nb - 1, (int(round(r * 110)) + 110) // PER)
        hist[j] += 1
        sub[j] += m
    top = max(hist)
    ys = lambda c: B - c / top * (B - T)
    o = []
    o.append(tx(0, 24, 'Every gene in the matrix was standardized; no gene tracks expansion more closely than FAN1', 19, INK, w=600))
    ls, _ = para(0, 46, '%s genes reaching 1 transcript per million, each correlated with expansion across the same 11 donors. Expansion orders '
                        'the donors, so any gene varying along that order produces a correlation: this is the distribution of what that alone yields.' % f'{NGENE:,}')
    o += ls
    for c in (0, top // 2, top):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(c), R, ys(c), GRID))
        o.append(tx(L - 10, ys(c) + 4, f'{c:,}', 12, MUT, 'end'))
    for j, c in enumerate(hist):
        x0, x1 = xs((-110 + j * PER) * STEP), xs((-110 + (j + 1) * PER) * STEP)
        w = max(1.0, x1 - x0 - 1.2)
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                 % (x0 + 0.6, ys(c), w, B - ys(c), '#c9c4b8'))
        if sub[j]:
            o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                     % (x0 + 0.6, ys(sub[j]), w, B - ys(sub[j]), DEEP))
    for v in (-1, -0.5, 0, 0.5, 1):
        o.append(tx(xs(v), B + 24, ('%+.1f' % v).replace('-', '−') if v else '0', 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'Spearman rho between the gene and expansion achieved', 14, INK, 'middle'))
    o.append(tx(44, (T + B) / 2, 'genes', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 44 %.1f)" ' % ((T + B) / 2)))
    fx = xs(FAN1_RHO)
    # FAN1 is one gene, so its bar is one gene tall. Mark it as a point on the
    # axis with a leader to the label; a full-height rule would read as a bar.
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" '
             'stroke-dasharray="3 4" opacity="0.6"/>' % (fx, B - 10, fx, T + 34, ORANGE))
    o.append('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s" stroke="%s" stroke-width="2"/>'
             % (fx, B, ORANGE, BG))
    o.append(tx(fx + 12, T + 26, 'FAN1,  rho = %+.2f' % FAN1_RHO, 15, ORANGE, w=700))
    o.append(tx(fx + 12, T + 48, 'rank 1 of %s' % f'{NGENE:,}', 13, INK))
    o.append(tx(fx + 12, T + 68, 'no gene reaches it, and its own bar is one gene tall', 13, INK))
    for i, (col, lab) in enumerate([('#c9c4b8', 'all %s genes' % f'{NGENE:,}'),
                                    (DEEP, 'the %s at FAN1’s own abundance, %.0f to %.0f TPM'
                                     % (f'{POOL_N:,}', D['fan1_tpm'] / 1.3, D['fan1_tpm'] * 1.3))]):
        o.append('<rect x="%.1f" y="%.1f" width="14" height="14" fill="%s"/>' % (fx + 12, T + 98 + i * 24, col))
        o.append(tx(fx + 34, T + 110 + i * 24, lab, 13, INK))
    # MSH2, the other half of MutSβ, marked the same way on the opposite side
    mx = xs(MSH2_RHO)
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" '
             'stroke-dasharray="3 4" opacity="0.6"/>' % (mx, B - 10, mx, T + 34, RED))
    o.append('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s" stroke="%s" stroke-width="2"/>'
             % (mx, B, RED, BG))
    o.append(tx(mx - 12, T + 26, 'MSH2,  rho = %+.2f' % MSH2_RHO, 15, RED, 'end', w=700))
    o.append(tx(mx - 12, T + 48, 'the partner that opposes FAN1,', 13, INK, 'end'))
    o.append(tx(mx - 12, T + 68, 'running the other way', 13, INK, 'end'))
    ls, y = para(0, 410, 'The darker part of each bar is the %s genes expressed within a factor of 1.3 of FAN1 itself. They are spread as widely as '
                         'the rest, a standard deviation of 0.33 against 0.34, and the nearest of them reaches %.2f.' % (f'{POOL_N:,}', POOL_MAX), 14, INK)
    o += ls
    ls, y = para(0, y + 6, 'That figure is a reference point rather than a p-value: genes move together, so the tail holds fewer than %s independent draws.'
                 % f'{NGENE:,}')
    o += ls
    ls, y = para(0, y + 6, 'MSH2 runs in the direction the mechanism implies, from the opposite side of the competition.', 13)
    o += ls
    ls, _ = para(0, y + 4, 'Eleven donors allow a correlation only 221 distinct values, 1/110 apart. Each bar holds five of them.', 13)
    o += ls
    return svg(W, H, 'Distribution of gene-by-gene correlations with expansion, FAN1 marked at the extreme', o)


# ---------------------------------------------------------------- panel 4
def bars(o, x0, label, ratio, mutsb, fan1, median_gain, col):
    """MutSbeta and FAN1 as two bars on a shared axis, with the ratio stated"""
    TOP, BASE, MAX = 112, 300, 360
    yb = lambda v: BASE - v / MAX * (BASE - TOP)
    for i, (v, nm, c) in enumerate(((mutsb, 'MutSβ', RED), (fan1, 'FAN1', GREEN))):
        x = x0 + i * 114
        o.append('<rect x="%.1f" y="%.1f" width="90" height="%.1f" rx="3" fill="%s"/>'
                 % (x, yb(v), BASE - yb(v), c))
        o.append(tx(x + 45, yb(v) - 10, '%.0f' % v, 15, INK, 'middle', w=700))
        o.append(tx(x + 45, BASE + 20, nm, 13, c, 'middle', w=600))
    mid = x0 + 102
    o.append(tx(mid, BASE + 48, label, 15, INK, 'middle', w=600))
    o.append(tx(mid, BASE + 78, '%.1f to 1' % ratio, 20, INK, 'middle', w=700))
    o.append(tx(mid, BASE + 100, 'median %d CAGs added' % median_gain, 13, MUT, 'middle'))


def panel_two_regions():
    W, H = 1392, 772
    o = []
    o.append(tx(0, 24, 'The striatum, where the repeats run furthest, is where FAN1 is most outnumbered — and where 6 donors leave the question open', 19, INK, w=600))
    o.append(tx(0, 46, 'Regional pseudobulk expression from the same matrices, and the same standardized correlation run on striatal projection neurons.', 14))

    TOP, BASE, MAX = 112, 300, 360
    yb = lambda v: BASE - v / MAX * (BASE - TOP)
    for v in (0, 100, 200, 300):
        o.append('<line x1="96" y1="%.1f" x2="760" y2="%.1f" stroke="%s"/>' % (yb(v), yb(v), GRID))
        o.append(tx(88, yb(v) + 4, str(v), 12, MUT, 'end'))
    o.append(tx(44, (TOP + BASE) / 2, 'transcripts per million', 13, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 44 %.1f)" ' % ((TOP + BASE) / 2)))
    bars(o, 130, 'Cortical excitatory neurons', RATIO['cortex'],
         TPM['cortex']['MSH2'] + TPM['cortex']['MSH3'], TPM['cortex']['FAN1'], 26, BLUE)
    bars(o, 480, 'Striatal projection neurons', RATIO['striatum'],
         TPM['striatum']['MSH2'] + TPM['striatum']['MSH3'], TPM['striatum']['FAN1'], 145, ORANGE)
    o.append('<line x1="818" y1="104" x2="818" y2="360" stroke="%s"/>' % GRID)

    # striatal scatter
    L, R, T, B = 920, 1340, 124, 296
    S = sorted(D['striatum'], key=lambda z: z['exp'])
    xs = lambda v: L + (v - 60) / (220 - 60) * (R - L)
    ys = lambda v: B - (v - 22) / (34 - 22) * (B - T)
    o.append(tx(L - 44, 102, 'Striatum: 6 donors, no relationship', 15, INK, w=600))
    for v in (22, 26, 30, 34):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(v), R, ys(v), GRID))
        o.append(tx(L - 10, ys(v) + 4, str(v), 12, MUT, 'end'))
    for v in (60, 120, 180):
        o.append(tx(xs(v), B + 22, str(v), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 48, 'CAGs added', 13, INK, 'middle'))
    o.append(tx(L - 44, (T + B) / 2, 'FAN1, TPM', 13, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 %.1f %.1f)" ' % (L - 44, (T + B) / 2)))
    for d in S:
        c = COLOR.get(d['d'], SLATE)
        o.append('<circle cx="%.1f" cy="%.1f" r="7" fill="%s" stroke="%s" stroke-width="2"/>'
                 % (xs(d['exp']), ys(d['fan1']), c, BG))
    rho_s = spearman([d['exp'] for d in S], [d['fan1'] for d in S])
    o.append(tx(R, T - 22, 'rho = %+.2f' % rho_s, 15, INK, 'end', w=600))
    ls, _ = para(L - 44, 368, 'Two of the six are also cortex donors and keep their colors; the other four are gray.', 12, MUT, maxw=466)
    o += ls

    ls, y = para(0, 440, 'All six striatal donors sit far out along the expansion axis, between 72 and 205 CAGs added, and their FAN1 values span '
                         'barely 8 transcripts per million.', 14, INK)
    o += ls
    ls, y = para(0, y + 6, 'With 6 donors and that little spread in either quantity, a correlation near zero is what a real effect of the cortical '
                           'size would most often produce here, so the striatum leaves the question open.')
    o += ls
    o.append(tx(0, y + 34, 'Taken together', 15, INK, w=600))
    ls, y = para(0, y + 58, 'The cortex result is one cohort of 11 and an association, and the direction of cause is open. It is worth a test in material '
                            'where repeat length is known before the neuron is lost: FAN1 measured in presymptomatic tissue, or in carriers matched on '
                            'inherited length and differing in onset. The mechanism already has independent support and the modifier genetics already '
                            'point at FAN1; what is new here is that the ordering holds between individuals, within one tissue.', 14, INK)
    o += ls
    ls, y = para(0, y + 12, 'The standardization, the barcode matching and the streaming reader of the expression matrices are at '
                            'github.com/crgfrench/hd-somatic-expansion; the figures above come from that code.', 13)
    o += ls
    ls, _ = para(0, y + 14, 'Sources: single-nucleus expression and per-nucleus repeat lengths, Kashin, Lee et al., Cell 2026 (Zenodo record 18202664), '
                            'cortical excitatory neurons and striatal projection neurons. Repair mechanism: Senoussi et al., Nature Communications 16:9445, 2025.', 12)
    o += ls
    return svg(W, H, 'MutS beta against FAN1 in two regions, drawn as a balance, beside the striatal scatter', o)


HEIGHT = 2712
page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>FAN1 and the pace of expansion</title>
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
<div style="width: 1520px; height: {HEIGHT}px; box-sizing: border-box; padding: 64px; background: {BG}; color: {INK}; font-family: 'IBM Plex Sans', system-ui, sans-serif; display: flex; flex-direction: column; gap: 32px">

<div style="display: flex; flex-direction: column; gap: 10px">
<div style="font-size: 13px; line-height: 18px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase; color: #8a8679">Board 6 of 7</div>
<h1 style="margin: 0; font-size: 34px; line-height: 42px; font-weight: 700">FAN1, the gene that ranks first across the 11 cortex donors</h1>
<p style="margin: 0; font-size: 16px; line-height: 24px; color: {MUT}">A re-analysis of the single-nucleus data released with the 2026 cortex paper (Zenodo record 18202664). Every gene in the matrix was given a per-donor expression rate standardized on the distribution of repeat lengths, then correlated with how far each donor&#8217;s repeats had traveled. The gene that comes out first is <i>FAN1</i>, the nuclease that removes the looped-out repeat.</p>
</div>

{panel_mechanism()}

{panel_cortex_scatter()}

{panel_null()}

{panel_two_regions()}

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
open(os.path.join(PROJ, 'FAN1.dc.html'), 'w').write(page)
print('wrote FAN1.dc.html', len(page), 'bytes')
print('cortex ratio %.1f  striatum ratio %.1f' % (RATIO['cortex'], RATIO['striatum']))
print('striatum rho %.3f' % spearman([d['exp'] for d in D['striatum']], [d['fan1'] for d in D['striatum']]))
print('cortex OLS', ols([d['exp'] for d in CX], [d['fan1'] for d in CX]))
print('median gain cortex %.0f striatum %.0f' % (
    sorted(d['exp'] for d in CX)[len(CX)//2], sorted(d['exp'] for d in D['striatum'])[len(D['striatum'])//2]))
