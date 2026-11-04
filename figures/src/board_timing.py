import os
HERE = os.path.dirname(os.path.abspath(__file__))
import json, math
D = json.load(open(os.path.join(HERE, 'donor_table.json')))

BG, INK, MUT, GRID = '#fbfaf7', '#1d1d1b', '#5c5a55', '#e4dfd4'
BLUE, ORANGE, BAND, GREEN = '#2f6db5', '#d9822b', '#efe8d8', '#4a8a5c'
E = lambda s: (s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;'))
# Red, orange, yellow, green, blue, violet, with a deep and a light shade of each
# but yellow. Hues run in spectral order and are handed out by inherited repeat
# length, so the sequence reads left to right on the right-hand panel.
DONOR_HUES = ['#a3202a', '#f08076',   # red
              '#c25a00', '#e89a2e',   # orange
              '#d9c400',              # yellow
              '#1d6e33', '#62c566',   # green
              '#1d4f93', '#5fb0ee',   # blue
              '#63309c', '#b074dd']   # violet
COLOR = {d['d']: DONOR_HUES[i] for i, d in
         enumerate(sorted(D['cortex'], key=lambda z: (z['cag'], z['d'])))}

def ols(X, Y):
    """slope, intercept and R squared"""
    n = len(X); mx = sum(X) / n; my = sum(Y) / n
    sl = sum((a - mx) * (b - my) for a, b in zip(X, Y)) / sum((a - mx) ** 2 for a in X)
    ic = my - sl * mx
    ssr = sum((b - (ic + sl * a)) ** 2 for a, b in zip(X, Y))
    sst = sum((b - my) ** 2 for b in Y)
    return sl, ic, 1 - ssr / sst

def eqn(ylab, xlab, sl, ic):
    """slope-intercept form, y = mx + b"""
    return '%s = (%s%.2f × %s %s %.1f)' % (
        ylab, '−' if sl < 0 else '', abs(sl), xlab,
        '−' if ic < 0 else '+', abs(ic))

def tx(x, y, s, sz=12, fill=MUT, anchor='start', w=None, it=False):
    a = ' text-anchor="%s"' % anchor if anchor != 'start' else ''
    b = ' font-weight="%s"' % w if w else ''
    i = ' font-style="italic"' if it else ''
    return '<text x="%.1f" y="%.1f" font-size="%d" fill="%s"%s%s%s>%s</text>' % (x, y, sz, fill, a, b, i, E(s))

# ---------------------------------------------------------------- panel 1
def panel_age_vs_cag():
    W, H = 1392, 500
    L, R, T, B = 86, 1340, 64, 366
    xs = lambda c: L + (c - 34) / (52 - 34) * (R - L)
    ys = lambda a: B - (a - 30) / (110 - 30) * (B - T)
    o = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="Age at death against inherited repeat length for 54 Huntington\'s disease donors" style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">' % (W, H, W, H)]
    o.append(tx(0, 28, 'Age at death falls as the inherited repeat lengthens', 19, INK, w=600))
    o.append(tx(0, 50, '54 donors in the 2026 archive. Spearman rho = −0.79, p = 2 × 10⁻¹²', 14))
    for a in range(30, 111, 20):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(a), R, ys(a), GRID))
        o.append(tx(L - 10, ys(a) + 4, str(a), 12, MUT, 'end'))
    for c in range(36, 53, 2):
        o.append(tx(xs(c), B + 24, str(c), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'Inherited CAG repeat length', 14, INK, 'middle'))
    o.append(tx(28, (T + B) / 2, 'Age at death', 14, INK, 'middle') .replace('<text ', '<text transform="rotate(-90 28 %.1f)" ' % ((T + B) / 2)))
    # least-squares line on ranks is not drawn; draw the ordinary least-squares line
    X = [d['cag'] for d in D['hd']]; Y = [d['age'] for d in D['hd']]
    sl, ic, r2 = ols(X, Y)
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-dasharray="7 5"/>'
             % (xs(36), ys(ic + sl * 36), xs(51), ys(ic + sl * 51), ORANGE))
    o.append(tx(R, 86, eqn('Age', 'CAG', sl, ic), 15, INK, 'end', w=600))
    o.append(tx(R, 108, 'R² = %.2f' % r2, 14, MUT, 'end'))
    seen = {}
    for d in D['hd']:
        k = (round(d['cag']), round(d['age'] / 3))
        seen[k] = seen.get(k, 0) + 1
        jig = (seen[k] - 1) * 7
        o.append('<circle cx="%.1f" cy="%.1f" r="6" fill="%s" fill-opacity="0.55" stroke="%s" stroke-width="1.2"/>'
                 % (xs(d['cag']) + jig, ys(d['age']), BLUE, BLUE))
    o.append(tx(0, 452, 'Age at death in this group is set by the disease, so the donors who had the most years are the ones whose repeats grow slowest.', 14, INK))
    o.append(tx(0, 476, 'A longer stretch of time and a faster expansion never occur together, and that pairing is the comparison which would give the growth rate a number.', 14, MUT))
    o.append('</svg>')
    return '\n'.join(o)

# ---------------------------------------------------------------- panel 2
def panel_two_scatters():
    W, H = 1392, 534
    o = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="Expansion achieved in cortical neurons against age at death and against inherited repeat length, 11 donors" style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">' % (W, H, W, H)]
    o.append(tx(0, 24, 'In the 11 donors with per-nucleus lengths, the inherited length predicts the expansion and the age does not', 19, INK, w=600))
    o.append(tx(0, 46, 'Expansion is the 90th percentile repeat length in that donor’s cortical neurons, less the length inherited', 14))
    T, B = 114, 374
    ymax = 190
    ys = lambda g: B - g / ymax * (B - T)
    for k, (x0, x1, lab, lo, hi, step, key, rho, pv, col, line) in enumerate([
            (60, 640, 'Age at death', 30, 90, 10, 'age', '−0.11', '0.74', '#a8a49a', 'age'),
            (752, 1332, 'Inherited CAG repeat length', 36, 52, 4, 'cag', '+0.85', '0.001', BLUE, 'CAG')]):
        xs = lambda v, x0=x0, x1=x1, lo=lo, hi=hi: x0 + (v - lo) / (hi - lo) * (x1 - x0)
        for g in (0, 50, 100, 150):
            o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (x0, ys(g), x1, ys(g), GRID))
            o.append(tx(x0 - 10, ys(g) + 4, str(g), 12, MUT, 'end'))
        v = lo
        while v <= hi:
            o.append(tx(xs(v), B + 22, str(v), 12, MUT, 'middle')); v += step
        o.append(tx((x0 + x1) / 2, B + 50, lab, 14, INK, 'middle'))
        o.append(tx(x0 - 46, (T + B) / 2, 'CAGs added', 13, INK, 'middle')
                 .replace('<text ', '<text transform="rotate(-90 %.1f %.1f)" ' % (x0 - 46, (T + B) / 2)))
        X = [d[key] for d in D['cortex']]; Y = [d['gain'] for d in D['cortex']]
        sl, ic, r2 = ols(X, Y)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-dasharray="7 5"/>'
                 % (xs(lo), ys(ic + sl * lo), xs(hi), ys(ic + sl * hi), ORANGE))
        for d in sorted(D['cortex'], key=lambda z: -z['gain']):
            o.append('<circle cx="%.1f" cy="%.1f" r="7.5" fill="%s" stroke="%s" stroke-width="2"/>'
                     % (xs(d[key]), ys(min(d['gain'], ymax)), COLOR[d['d']], BG))
        o.append(tx(x1, T - 42, 'rho = %s,  p = %s' % (rho, pv), 15, INK, 'end', w=600))
        o.append(tx(x1, T - 20, eqn('CAGs added', 'Age' if line == 'age' else 'CAG', sl, ic), 14, INK, 'end'))
        o.append(tx(x1, T - 2, 'R² = %.2f' % r2, 13, MUT, 'end'))
    o.append(tx(0, 470, 'A growth rate is distance divided by time. The distance is the vertical axis of both panels, and the time is measured by the donor’s age.', 14, INK))
    o.append(tx(0, 494, 'The left panel shows age predicting nothing about how far the repeat moved, so this group yields no rate. The length a donor inherited predicts it instead, on the right.', 14, MUT))
    o.append(tx(0, 516, 'The left line falls by less than one CAG per year of age and accounts for 5% of the difference between donors; the right line accounts for 35%.', 14, MUT))
    o.append('</svg>')
    return '\n'.join(o)

# ---------------------------------------------------------------- panel 3
ONSET = 649.0 / 12          # inherited length 42
def panel_delay():
    W, H = 1392, 520
    L, R = 150, 1300
    xs = lambda a: L + (a - 10) / (110 - 10) * (R - L)
    o = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="Years of delay in expected onset bought by slowing expansion, for four treatment start ages" style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">' % (W, H, W, H)]
    o.append(tx(0, 24, 'Slowing expansion stretches the time that is left, so an early start is worth more than a strong effect', 19, INK, w=600))
    o.append(tx(0, 46, 'A person born with 42 repeats, expected onset at age 54. Rate and time enter any growth law only as their product, so cutting the rate by a fraction f multiplies every remaining interval by 1/(1−f).', 14))
    rows = [20, 30, 40, 50]
    FS = [(0.2, '#9fc0e0'), (0.3, BLUE), (0.5, '#1d4f8a')]
    y = 112
    for A in rows:
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, y, R, y, GRID))
        o.append(tx(L - 18, y + 5, 'start at %d' % A, 14, INK, 'end', w=600))
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="12" rx="2" fill="%s"/>' % (xs(10), y - 6, xs(A) - xs(10), BAND))
        # untreated interval left to run
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="12" rx="2" fill="#d5cfc2"/>' % (xs(A), y - 6, xs(ONSET) - xs(A)))
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>' % (xs(A), y - 18, xs(A), y + 18, MUT))
        prev = ONSET
        for j, (f, col) in enumerate(FS):
            nd = ONSET + (ONSET - A) * f / (1 - f)
            yy = y + 0 - 0
            o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="12" rx="2" fill="%s"/>' % (xs(prev), yy - 6, max(xs(nd) - xs(prev), 0.8), col))
            prev = nd
        o.append('<circle cx="%.1f" cy="%.1f" r="5.5" fill="%s"/>' % (xs(ONSET), y, ORANGE))
        o.append(tx(xs(prev) + 10, y + 5, '+%d years at half speed' % round((ONSET - A) * 1.0), 13, INK))
        y += 76
    # axis
    ay = y - 24
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ay, R, ay, GRID))
    for a in range(10, 111, 10):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (xs(a), ay, xs(a), ay + 6, MUT))
        o.append(tx(xs(a), ay + 24, str(a), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, ay + 50, 'Age in years', 14, INK, 'middle'))
    o.append('<line x1="%.1f" y1="92" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="6 4"/>' % (xs(ONSET), xs(ONSET), ay, ORANGE))
    o.append(tx(xs(ONSET) + 8, 86, 'expected onset, untreated: 54', 13, ORANGE, w=600))
    # legend
    lx = L
    o.append(tx(lx, ay + 84, 'Added years of delay:', 14, INK, w=600))
    lx += 170
    for f, col in FS:
        o.append('<rect x="%.1f" y="%.1f" width="26" height="12" rx="2" fill="%s"/>' % (lx, ay + 74, col))
        o.append(tx(lx + 34, ay + 84, 'expansion slowed %d%%' % (100 * f), 14, MUT))
        lx += 230
    o.append('</svg>')
    return '\n'.join(o)

# ---------------------------------------------------------------- panel 4
def panel_required():
    W, H = 1392, 432
    o = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="Slowing required to buy ten years of delay, by age at which treatment starts" style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">' % (W, H, W, H)]
    o.append(tx(0, 24, 'The same relationship read backwards: how much slowing ten years of delay would take', 19, INK, w=600))
    o.append(tx(0, 46, 'Born with 42 repeats, expected onset at age 54', 14))
    L, R, T, B = 110, 1340, 86, 266
    xs = lambda a: L + (a - 15) / (55 - 15) * (R - L)
    ys = lambda f: B - f * (B - T)
    for f in (0, .25, .5, .75, 1.0):
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (L, ys(f), R, ys(f), GRID))
        o.append(tx(L - 10, ys(f) + 4, '%d%%' % (100 * f), 12, MUT, 'end'))
    pts = []
    for A in range(18, 54):
        f = 10.0 / ((ONSET - A) + 10.0)
        pts.append('%.1f %.1f' % (xs(A), ys(f)))
    o.append('<path d="M%s" fill="none" stroke="%s" stroke-width="3"/>' % (' L'.join(pts), BLUE))
    for A, lab in ((20, '23%'), (30, '29%'), (40, '42%'), (50, '71%')):
        f = 10.0 / ((ONSET - A) + 10.0)
        o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 4"/>' % (xs(A), ys(f), xs(A), B, GRID))
        o.append('<circle cx="%.1f" cy="%.1f" r="6" fill="%s"/>' % (xs(A), ys(f), BLUE))
        o.append(tx(xs(A), ys(f) - 14, lab, 15, INK, 'middle', w=600))
    for a in range(20, 55, 5):
        o.append(tx(xs(a), B + 24, str(a), 12, MUT, 'middle'))
    o.append(tx((L + R) / 2, B + 52, 'Age at which treatment starts', 14, INK, 'middle'))
    o.append(tx(L - 72, (T + B) / 2, 'Slowing required', 14, INK, 'middle')
             .replace('<text ', '<text transform="rotate(-90 %.1f %.1f)" ' % (L - 72, (T + B) / 2)))
    o.append(tx(0, 356, 'Silencing MSH3 by about 60% blocked measurable expansion in the striatum of BAC-CAG mice over the four months followed. A partial', 14, INK))
    o.append(tx(0, 378, 'knockdown produced a near-complete stop there, so the percentages on this curve measure how much the expansion slows rather than how much enzyme is removed.', 14, INK))
    o.append(tx(0, 400, 'Read against the curve: ten years of delay takes 29% slowing if it begins at age 30, 42% at 40 and 71% at 50.', 14, MUT))
    o.append('</svg>')
    return '\n'.join(o)

HEIGHT = 2460
page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Timing of the expansion</title>
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
<div style="width: 1520px; height: {HEIGHT}px; box-sizing: border-box; padding: 64px; background: {BG}; color: {INK}; font-family: 'IBM Plex Sans', system-ui, sans-serif; display: flex; flex-direction: column; gap: 36px">

<div style="display: flex; flex-direction: column; gap: 10px">
<div style="font-size: 13px; line-height: 18px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase; color: #8a8679">Board 5 of 7</div>
<h1 style="margin: 0; font-size: 34px; line-height: 42px; font-weight: 700">How fast the repeat grows, and what slowing it would buy</h1>
<p style="margin: 0; font-size: 16px; line-height: 24px; color: {MUT}">Donor records released with the 2026 cortex paper (Zenodo record 18202664): age at death, the repeat length inherited on the affected copy of <i>HTT</i>, and per-nucleus repeat lengths in cortical neurons for 11 of the donors. Expected onset ages are 649/(CAG&#8722;30) years, the Enroll-HD calibration of the CAG-Age Product.</p>
</div>

{panel_age_vs_cag()}

{panel_two_scatters()}

{panel_delay()}

{panel_required()}

<p style="margin: 0; font-size: 12px; line-height: 18px; color: {MUT}">Sources: donor records and per-nucleus repeat lengths, Kashin, Lee et al., <i>Cell</i> 2026 (Zenodo record 18202664). Expected onset from the CAG-Age Product, Warner et al. 2020. MSH3 silencing in mice: O&#8217;Reilly et al., <i>Molecular Therapy</i> 31:1661, 2023; the four-month figures are read from the bioRxiv preprint of that paper. The delay and slowing curves are arithmetic on the onset formula and were never fitted to data.</p>

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
open(os.path.join(HERE, 'Timing.dc.html'), 'w').write(page)
print('wrote Timing.dc.html', len(page), 'bytes')
for A in (20, 30, 40, 50):
    print('  start %d: +%.0f at 20%%, +%.0f at 30%%, +%.0f at 50%%' % (A, (ONSET-A)*.25, (ONSET-A)*3/7, (ONSET-A)))
