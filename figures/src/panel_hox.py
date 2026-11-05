"""Draw the HOX-cluster panel and write it into the re-analysis board.

The phase D genes that sit next to each other on a chromosome are, almost
entirely, the four HOX clusters. This panel draws each cluster in paralog
order with its phase D members filled.

Paralog order is the genomic order: a HOX cluster is collinear, so the genes run
in numerical order along the chromosome, ascending for HOXA and HOXB and
descending for HOXC and HOXD. The panel sorts on the paralog number rather than
on the stored annotation start, because the annotation keeps one row per symbol
and for HOXC6 that row is an upstream isoform whose start sits between HOXC10
and HOXC9.

Two cautions are written into the panel because the result is easy to misread:

1. The counts come from the gene-family test rather than from the spatial test below. A
   spatial null drawn from every gene in the statistics file is wrong here,
   because the phase sets have a median abundance of about 190 transcripts per
   million against 17 for all tested genes, and clustered paralog families are
   quiet. Matching the null on abundance fixes that for phase C, where the
   clustering becomes significant (136 against 102 expected, p = 0.002).

2. Matching cannot be applied to the HOX result at all: none of the 11 phase D
   HOX genes reaches 1 transcript per million, so every one of them falls out
   of the abundance-matched comparison. The family counts here are the evidence.

The script rewrites the region between the hox:start and hox:end markers in
Reanalysis.dc.html, so it can be run again after the data change.
"""
import json, re, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BOARD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'Reanalysis.dc.html')
D = json.load(open(os.path.join(HERE, 'hox_clusters.json')))

INK, MUT, GRID, BLUE, PALE = '#1d1d1b', '#5c5a55', '#e4dfd4', '#2f6db5', '#9bb9dc'
E = lambda s: s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def tx(x, y, s, sz=12, fill=MUT, anchor='start', w=None):
    a = ' text-anchor="%s"' % anchor if anchor != 'start' else ''
    b = ' font-weight="%s"' % w if w else ''
    return '<text x="%.1f" y="%.1f" font-size="%d" fill="%s"%s%s>%s</text>' % (x, y, sz, fill, a, b, E(s))


def para(x, y, text, size=13, fill=MUT, maxw=1392, lh=20):
    per = max(10, int(maxw / (0.505 * size)))
    out, line = [], ''
    for word in text.split():
        if len(line) + len(word) + 1 > per:
            out.append(line); line = word
        else:
            line = (line + ' ' + word).strip()
    if line:
        out.append(line)
    return [tx(x, y + i * lh, ln, size, fill) for i, ln in enumerate(out)], y + len(out) * lh


def short(g):
    """HOXA11-AS reads as 11-AS once the cluster letter is the row label"""
    return g[4:]


W, H = 1392, 500
BW, BH, GAP = 84, 30, 8
o = ['<svg width="%d" height="%d" viewBox="0 0 %d %d" role="img" aria-label="The four HOX clusters in genomic order, with the phase D genes filled" style="flex: none; font-family: \'IBM Plex Sans\', system-ui, sans-serif">' % (W, H, W, H)]
o.append(tx(0, 20, 'The phase D genes that sit next to each other are the HOX clusters', 18, INK, w=600))
ls, _ = para(0, 42, 'The 22 genes named HOXA to HOXD with a paralog number that appear in the cortex paper\u2019s per-gene statistics, together with HOXA11-AS, drawn in paralog order within each '
                    'cluster, which is the order they occupy on the chromosome. Filled boxes are phase D genes. The other long non-coding RNAs embedded in these clusters have '
                    'names of their own, so they are counted separately below.', 13, lh=19)
o += ls

cx = {r['g']: r for r in D['cortex']}
y = 108
for cl in 'ABCD':
    row = sorted([r for r in D['cortex'] if r['cl'] == cl],
                 key=lambda r: r['par'], reverse=(cl in 'CD'))
    o.append(tx(0, y + 20, 'HOX' + cl, 14, INK, w=600))
    for i, r in enumerate(row):
        x = 72 + i * (BW + GAP)
        fill, stroke, col = (BLUE, BLUE, '#ffffff') if r['d'] else ('#ffffff', GRID, MUT)
        o.append('<rect x="%.1f" y="%.1f" width="%d" height="%d" rx="4" fill="%s" stroke="%s"/>'
                 % (x, y, BW, BH, fill, stroke))
        o.append(tx(x + BW / 2, y + 20, short(r['g']), 13, col, 'middle', w=600 if r['d'] else None))
    y += BH + 16

o.append('<line x1="%d" y1="100" x2="%d" y2="%d" stroke="%s"/>' % (940, 940, y - 16, GRID))
o.append(tx(980, 122, 'Cortex', 15, INK, w=600))
o.append(tx(980, 146, '11 of these 23 genes are phase D genes.', 13, INK))
o.append(tx(980, 164, 'Chance would put 0.22 of them there.', 13, MUT))
o.append(tx(980, 200, 'Striatum', 15, INK, w=600))
o.append(tx(980, 224, '4 of 17, against 0.54 expected.', 13, INK))
o.append(tx(980, 260, 'What else runs together', 15, INK, w=600))
o.append(tx(980, 284, 'MEG3 with MEG8 on chromosome 14, the imprinted', 13, INK))
o.append(tx(980, 302, 'DLK1-DIO3 locus, whose silencing is held by DNA', 13, INK))
o.append(tx(980, 320, 'methylation at a control region. Nothing else.', 13, INK))

cy = 352
ls, cy = para(0, cy, 'Set the clusters aside and the remaining phase D genes are no more bunched than chance, so what they share is being a polycomb target rather than '
                     'being neighbors. Polycomb targets lie all over the genome; HOX is where they sit in a row.', 13, INK)
o += ls
ls, cy = para(0, cy + 8, 'The clusters bring their own long non-coding RNAs with them: HOTAIRM1, HOTTIP, HOXA-AS2, HOXD-AS2 and HOXB-AS1 all reach phase D. Scored in the '
                         'nuclei themselves rather than per donor, HOX expression runs at 0.33 transcripts per million below 100 CAGs and 33 above 260, present in 1.2% of nuclei '
                         'and then in 23% of them.')
o += ls
ls, cy = para(0, cy + 8, 'Counted by gene family rather than by position, because the phase sets run near 190 transcripts per million against 17 for all tested genes, '
                         'so a null drawn from every gene is mismatched, and matching on abundance drops all 11 HOX genes, which sit below the 1 TPM floor of the abundance table.')
o += ls
o.append('</svg>')
svg = '\n'.join(o)

board = open(BOARD).read()
block = '<!-- hox:start -->\n' + svg + '\n<!-- hox:end -->'
if '<!-- hox:start -->' in board:
    board = re.sub(r'<!-- hox:start -->.*?<!-- hox:end -->', lambda m: block, board, flags=re.S)
else:
    anchor = '\n</div>\n</x-dc>'
    assert anchor in board, 'could not find the end of the board'
    board = board.replace(anchor, '\n' + block + anchor, 1)
height = int(re.search(r'width: 1520px; height: (\d+)px', board).group(1))
open(BOARD, 'w').write(board)
print('panel written into %s, board height now %d' % (os.path.basename(BOARD), height))
