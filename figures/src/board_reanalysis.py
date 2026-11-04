"""The two per-gene onset panels of the re-analysis board.

Takes the per-gene onset summary as JSON on the command line; it carries a
"cand" list built by notebook 01. See data/PUT_INPUT_FILES_HERE.txt.

This writes an SVG fragment, not a whole board. The fragment was pasted into
the hand-written board beside it, which is the source of record the page is
assembled from; the script is kept so the drawing can be regenerated.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import json, sys
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'reanalysis_panels.json')
if not os.path.exists(SRC):
    raise SystemExit('need the onset summary JSON: pass it as the first argument (looked for %s)' % SRC)
D = json.load(open(SRC))
INK='#1d1d1b';Q='#5c5a55';GRID='#e4dfd4';RISE='#2f6db5';FALL='#d9822b';STR='#7b5ea7';CTX='#2a8a7e'
o=[]
def t(x,y,s,size=14,anchor='start',fill=INK,w=None,extra=''):
    fw=f' font-weight="{w}"' if w else ''
    o.append(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}"{fw} fill="{fill}"{extra}>{s}</text>')
# ---------- panel A: two histograms
def hist(reg,x0,label):
    W=640;H=300;y0=70;ymax=75
    xs=lambda c: x0+76+(c-40)/270*(W-96)
    yv=lambda v: y0+H-v/ymax*H
    o.append(f'<g>')
    t(x0+60,40,label,18,w=600)
    for v in (0,25,50,75):
        o.append(f'<line x1="{x0+60}" y1="{yv(v)}" x2="{x0+W}" y2="{yv(v)}" stroke="{GRID}"/>')
        t(x0+50,yv(v)+4,v,12,'end',Q)
    bw=(W-96)/27/2-1
    for s,c,off in (('C+',RISE,-bw),('C-',FALL,1)):
        for k,v in D[reg+s].items():
            k=int(k);x=xs(k)+off
            o.append(f'<rect x="{x:.1f}" y="{yv(v):.1f}" width="{bw:.1f}" height="{yv(0)-yv(v):.1f}" fill="{c}"/>')
    o.append(f'<line x1="{xs(150)}" y1="{y0-6}" x2="{xs(150)}" y2="{yv(0)}" stroke="{INK}" stroke-width="2" stroke-dasharray="4 3"/>')
    t(xs(150)+6,y0+2,'150 CAGs: where both papers',12)
    t(xs(150)+6,y0+18,'place the start of phase C',12,fill=Q)
    o.append(f'<line x1="{x0+60}" y1="{yv(0)}" x2="{x0+W}" y2="{yv(0)}" stroke="{Q}" stroke-width="1.5"/>')
    for c in (50,100,150,200,250,300):
        t(xs(c),yv(0)+20,c,13,'middle')
    t(x0+60+(W-60)/2,yv(0)+46,'CAG length at which the gene starts to change',14,'middle')
    o.append('</g>')
    return xs,yv
SVG_A=[]
o=[];xs,yv=hist('striatum',0,'Striatum: 557 genes, most start between 100 and 170')
# annotations
t(xs(70),yv(44),'rising genes,',12,'middle',RISE,600)
t(xs(70),yv(44)+15,'median 120',12,'middle',RISE,600)
t(xs(185),yv(44),'falling genes,',12,'start',FALL,600)
t(xs(185),yv(44)+15,'median 140',12,'start',FALL,600)
xs2,yv2=hist('cortex',720,'Cortex: 446 genes, half change from the shortest repeats')
t(xs2(40)+14,yv2(71)+4,'71 rising genes with no start point:',12,'start',INK,600)
t(xs2(40)+14,yv2(71)+19,'they change steadily',12,'start',Q)
t(xs2(40)+14,yv2(71)+34,'from the shortest repeats',12,'start',Q)
A=''.join(o)
# ---------- panel B: four small line charts
bins=['[36,100)','[100,150)','[150,175)','[175,200)','[200,250)','[250,300)','[300,400)','[400,2000)']
blab=['36–99','100–149','150–174','175–199','200–249','250–299','300–399','400+']
o=[]
genes=[('MBNL2','Falls in the striatum from 150 CAGs',0.5,1.1,[0.6,0.8,1.0]),('HTT','Falls in the cortex past 250 CAGs',0.7,1.2,[0.8,1.0,1.2]),('EZH2','PRC2; rises in the cortex',0.5,2.1,[0.5,1.0,1.5,2.0]),('CREBBP','CBP; rises in the cortex',0.9,1.4,[1.0,1.2,1.4]),('TCERG1','Rises in the striatum only',0.8,1.4,[0.8,1.0,1.2,1.4])]
PW=259;PG=24;PH=200;Y0=60
for i,(g,title,lo,hi,ticks) in enumerate(genes):
    x0=i*(PW+PG)
    X=lambda j: x0+60+j*(PW-72)/7
    Y=lambda v: Y0+PH-(v-lo)/(hi-lo)*PH
    t(x0,20,g,17,w=600); t(x0,40,title,12.5,fill=Q)
    for v in ticks:
        o.append(f'<line x1="{x0+44}" y1="{Y(v):.1f}" x2="{x0+PW}" y2="{Y(v):.1f}" stroke="{GRID}"/>')
        t(x0+38,Y(v)+4,f'{v:.1f}',11.5,'end',Q)
    o.append(f'<line x1="{x0+44}" y1="{Y(1):.1f}" x2="{x0+PW}" y2="{Y(1):.1f}" stroke="{Q}" stroke-width="1.5"/>')
    for reg,c,dx in (('striatum',STR,-4),('cortex',CTX,4)):
        pts=[r for r in D['cand'] if r['GENE']==g and r['region']==reg]
        pts.sort(key=lambda r: bins.index(r['bin']))
        path=' '.join(('M' if j==0 else 'L')+f'{X(bins.index(r["bin"]))+dx:.1f} {Y(float(r["mean"])):.1f}' for j,r in enumerate(pts))
        o.append(f'<path d="{path}" fill="none" stroke="{c}" stroke-width="2"/>')
        for r in pts:
            x=X(bins.index(r['bin']))+dx
            l=max(lo,float(r['lo']));h=min(hi,float(r['hi']))
            o.append(f'<line x1="{x:.1f}" y1="{Y(l):.1f}" x2="{x:.1f}" y2="{Y(h):.1f}" stroke="{c}" stroke-width="1.5" opacity="0.6"/>')
            o.append(f'<circle cx="{x:.1f}" cy="{Y(float(r["mean"])):.1f}" r="4" fill="{c}"/>')
    for j,b in enumerate(blab):
        o.append(f'<text x="{X(j):.1f}" y="{Y0+PH+16}" text-anchor="end" font-size="11" fill="{Q}" transform="rotate(-40 {X(j):.1f} {Y0+PH+16})">{b}</text>')
B=''.join(o)
open(os.path.join(HERE, 'A.svg'), 'w').write(A)
open(os.path.join(HERE, 'B.svg'), 'w').write(B)
print('wrote A.svg and B.svg beside this script')
