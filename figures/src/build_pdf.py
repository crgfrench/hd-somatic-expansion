"""Render Repeat_Disease_Thresholds.pdf as one page per board.

Each board is a fixed-size block, so a single page size cuts through whichever
boards do not match it. Here each board is rendered on its own, on a page sized
to that board, and the pages are joined.
"""
import re, json, subprocess, os, tempfile
from pypdf import PdfWriter

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..')
SRC = os.path.join(FIG, 'Repeat_Disease_Thresholds.html')
OUT = os.path.join(FIG, 'Repeat_Disease_Thresholds.pdf')
TMP = os.path.join(tempfile.gettempdir(), 'rdt_pdf_pages')
HEAD = ('<!doctype html><html><head><meta charset="utf-8">'
        '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
        'family=IBM+Plex+Sans:ital,wght@0,400;0,600;0,700;1,400&display=swap">'
        '<style>html,body{margin:0;padding:0;background:#fbfaf7}'
        '@page{margin:0}</style></head><body>')

def boards(text):
    out, pos = [], 0
    while True:
        i = text.find('<div style="width: 1520px', pos)
        if i < 0:
            return out
        depth = 0
        for m in re.finditer(r'<div\b|</div>', text[i:]):
            depth += 1 if m.group().startswith('<div') else -1
            if depth == 0:
                out.append(text[i:i + m.end()])
                pos = i + m.end()
                break
        else:
            raise AssertionError('unbalanced divs')

os.makedirs(TMP, exist_ok=True)
divs = boards(open(SRC).read())
pages = []
for k, d in enumerate(divs):
    h = int(re.search(r'height: (\d+)px', d).group(1))
    p = os.path.join(TMP, 'page%d.html' % k)
    open(p, 'w').write(HEAD + d + '</body></html>')
    pages.append((p, h))
json.dump(pages, open(os.path.join(TMP, 'pages.json'), 'w'))

subprocess.run(['node', os.path.join(HERE, 'build_pdf.js'), TMP], check=True)

w = PdfWriter()
for k, (_, h) in enumerate(pages):
    w.append(os.path.join(TMP, 'page%d.pdf' % k))
w.write(OUT)
print('%d boards | %d pages | %d bytes' % (len(pages), len(w.pages), os.path.getsize(OUT)))
