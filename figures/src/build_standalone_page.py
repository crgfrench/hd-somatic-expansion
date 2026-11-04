"""Assemble the seven boards into the single self-contained page.

Four of the boards are written by hand and live beside this script as their own
source of record. Three are generated: `board_timing.py`, `board_fan1.py` and
`board_silencing.py` each write their board into this directory, and
`panel_hox.py` rewrites one panel inside the re-analysis board. Run those first
if their inputs have changed; this script only assembles whatever is here.

Board 1 draws its 47-disease table through the canvas runtime's <sc-for> and
{{ }} markup rather than as plain HTML, so that one board is expanded here using
the values printed by `expand_template.js` (Node). The other six are already
plain markup and are copied through unchanged.

    python3 build_standalone_page.py
    python3 build_standalone_page.py --boards some/dir --out some/page.html
"""
import argparse
import html as H
import json
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ['Main.dc.html', 'Inferences.dc.html', 'Histogram.dc.html',
         'Reanalysis.dc.html', 'Timing.dc.html', 'FAN1.dc.html', 'Silencing.dc.html']
GENERATED = {'Timing.dc.html': 'board_timing.py',
             'FAN1.dc.html': 'board_fan1.py',
             'Silencing.dc.html': 'board_silencing.py'}

HEAD = """<!doctype html>
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
"""
TAIL = "\n</body>\n</html>\n"


def board_div(text, name):
    """the board itself, without the canvas wrapper around it"""
    i = text.find('<div style="width: 1520px')
    if i < 0:
        raise SystemExit('%s does not look like a board file' % name)
    depth = 0
    for m in re.finditer(r'<div\b|</div>', text[i:]):
        depth += 1 if m.group().startswith('<div') else -1
        if depth == 0:
            return text[i:i + m.end()]
    raise SystemExit('%s has unbalanced divs' % name)


def lookup(scope, path):
    v = scope
    for part in path.strip().split('.'):
        v = v[part] if isinstance(v, dict) else getattr(v, part)
    return v


def fmt(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, float) and v == int(v):
        return str(int(v))
    return str(v)


def holes(text, scope):
    return re.sub(r'\{\{\s*([A-Za-z_][\w.]*)\s*\}\}',
                  lambda m: H.escape(fmt(lookup(scope, m.group(1))), quote=True), text)


OPEN = re.compile(r'<sc-for\s+list="\{\{\s*([\w.]+)\s*\}\}"\s+as="(\w+)"[^>]*>')


def expand(text, scope):
    """resolve each <sc-for> outermost first, so an outer name is bound before
    an inner loop refers to it, then fill the {{ }} holes"""
    out, pos = [], 0
    while True:
        m = OPEN.search(text, pos)
        if not m:
            out.append(holes(text[pos:], scope))
            return ''.join(out)
        out.append(holes(text[pos:m.start()], scope))
        depth, k = 1, m.end()
        while depth:
            nxt = re.compile(r'<sc-for\b|</sc-for>').search(text, k)
            if not nxt:
                raise SystemExit('unclosed sc-for in the first board')
            depth += 1 if nxt.group() == '<sc-for' else -1
            k = nxt.end()
        body = text[m.end():k - len('</sc-for>')]
        for item in lookup(scope, m.group(1)):
            out.append(expand(body, dict(scope, **{m.group(2): item})))
        pos = k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--boards', default=HERE, help='directory holding the seven board files')
    ap.add_argument('--out', default=os.path.join(HERE, '..', 'Repeat_Disease_Thresholds.html'))
    a = ap.parse_args()

    missing = [f for f in ORDER if not os.path.exists(os.path.join(a.boards, f))]
    if missing:
        lines = ['these board files are not in %s:' % a.boards]
        for f in missing:
            how = ('run %s first' % GENERATED[f]) if f in GENERATED \
                else 'written by hand; it should be in the repository'
            lines.append('  %-20s %s' % (f, how))
        raise SystemExit('\n'.join(lines))

    vals = json.loads(subprocess.run(
        ['node', os.path.join(HERE, 'expand_template.js')],
        capture_output=True, text=True, check=True).stdout)

    first = expand(board_div(open(os.path.join(a.boards, ORDER[0])).read(), ORDER[0]), vals)
    if '{{' in first or '<sc-for' in first:
        raise SystemExit('the first board still holds unexpanded template markup')
    rows = sum(len(g['rows']) for g in vals['groups'])
    parts = [first]
    for f in ORDER[1:]:
        d = board_div(open(os.path.join(a.boards, f)).read(), f)
        if '{{' in d:
            raise SystemExit('%s holds template markup, which only the first board should' % f)
        parts.append(d)

    out = os.path.abspath(a.out)
    open(out, 'w').write(HEAD + '\n'.join(parts) + TAIL)
    print('%d boards | %d disease rows expanded | %d bytes | %s'
          % (len(parts), rows, os.path.getsize(out), out))


if __name__ == '__main__':
    main()
