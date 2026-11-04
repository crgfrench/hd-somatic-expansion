"""Rebuild out/Repeat_Disease_Thresholds.html from the five board files.

Main.dc.html draws its disease table through the canvas runtime's <sc-for> and
{{ }} markup, so that board is expanded here with the values printed by
expand.js. The other four boards are already plain markup.
"""
import json, re, subprocess, html as H

P = '/tmp/claude-0/-home-claude/c597d61d-b908-51d6-86bf-6a3ff3524319/scratchpad/canvas/project/'
OUT = '/home/claude/out/Repeat_Disease_Thresholds.html'
ORDER = ['Main.dc.html', 'Inferences.dc.html', 'Histogram.dc.html',
         'Reanalysis.dc.html', 'Timing.dc.html']

def board_div(text):
    i = text.find('<div style="width: 1520px')
    assert i >= 0
    depth = 0
    for m in re.finditer(r'<div\b|</div>', text[i:]):
        depth += 1 if m.group().startswith('<div') else -1
        if depth == 0:
            return text[i:i + m.end()]
    raise AssertionError('unbalanced divs')

def lookup(scope, path):
    v = scope
    for part in path.strip().split('.'):
        v = v[part] if isinstance(v, dict) else getattr(v, part)
    return v

def fmt(v):
    if isinstance(v, bool): return 'true' if v else 'false'
    if isinstance(v, float) and v == int(v): return str(int(v))
    return str(v)

def holes(text, scope):
    return re.sub(r'\{\{\s*([A-Za-z_][\w.]*)\s*\}\}',
                  lambda m: H.escape(fmt(lookup(scope, m.group(1))), quote=True), text)

OPEN = re.compile(r'<sc-for\s+list="\{\{\s*([\w.]+)\s*\}\}"\s+as="(\w+)"[^>]*>')

def expand(text, scope):
    """resolve each <sc-for> outermost first, so an outer name is bound before
    an inner loop refers to it, then fill the {{ }} holes"""
    out = []
    pos = 0
    while True:
        m = OPEN.search(text, pos)
        if not m:
            out.append(holes(text[pos:], scope))
            return ''.join(out)
        out.append(holes(text[pos:m.start()], scope))
        depth, k = 1, m.end()
        while depth:
            nxt = re.compile(r'<sc-for\b|</sc-for>').search(text, k)
            assert nxt, 'unclosed sc-for'
            depth += 1 if nxt.group() == '<sc-for' else -1
            k = nxt.end()
        body = text[m.end():k - len('</sc-for>')]
        for item in lookup(scope, m.group(1)):
            out.append(expand(body, dict(scope, **{m.group(2): item})))
        pos = k

vals = json.loads(subprocess.run(['node', '/home/claude/gen3/expand.js'],
                                 capture_output=True, text=True, check=True).stdout)
main = expand(board_div(open(P + 'Main.dc.html').read()), vals)
assert '{{' not in main and '<sc-for' not in main
assert 'Baratela-Scott' in main, 'the disease rows did not expand'

parts = [main]
for fn in ORDER[1:]:
    d = board_div(open(P + fn).read())
    assert '{{' not in d, fn + ' still holds template markup'
    parts.append(d)

open(OUT, 'w').write("""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Repeat Disease Thresholds</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:ital,wght@0,400;0,600;0,700;1,400&amp;display=swap">
<style>
body{margin:0;background:#fbfaf7}
a{color:#2f6db5}a:hover{color:#1d4f8a}
</style>
</head>
<body>
""" + '\n'.join(parts) + """
</body>
</html>
""")
rows = sum(len(g['rows']) for g in vals['groups'])
print('boards %d | disease rows expanded %d | bytes %d' % (len(parts), rows, len(open(OUT).read())))
