"""The two drawn expansion panels of the first board.

The curves are drawn from a stated growth rule rather than fitted to data,
which the board caption says.

This writes an SVG fragment, not a whole board. The fragment was pasted into
the hand-written board beside it, which is the source of record the page is
assembled from; the script is kept so the drawing can be regenerated.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import math, random

INK='#1d1d1b'; Q='#5c5a55'; GRID='#e4dfd4'; BLUE='#2f6db5'; ORANGE='#d9822b'; ACC='#4f7a4a'
PH=300; Y0=56; PW=600
K0=0.10; SIG=0.8
def A_of(L0):
    u0=(L0-35)/10.0
    return (10/0.7)*(u0**-0.7 - 11.5**-0.7)
T50=lambda L0: A_of(L0)/K0
def share(L0,t):
    if L0<=35 or t<=0: return 0.0
    z=math.log(t/T50(L0))/SIG
    return 0.5*(1+math.erf(z/math.sqrt(2)))
onset=lambda cag: 649.0/(cag-30)

def traj(L0,k,years=72,dt=0.25):
    L=L0; out=[(0,L)]; t=0
    while t<years:
        L += k*((L-35)/10.0)**1.7*dt; t += dt
        if L>260: out.append((t,260)); break
        out.append((t,L))
    return out

def build(x0,L0,seed):
    X=lambda t: x0+70+t/72*(PW-90)
    Y=lambda L: Y0+PH-(min(L,250)-35)/215*PH
    Ys=lambda s: Y0+PH-s*PH
    o=[f'<text x="{x0}" y="22" font-size="17" font-weight="600" fill="{INK}">Born with {L0} repeats</text>',
       f'<text x="{x0}" y="42" font-size="13" fill="{Q}">Twenty-four neurons from one person, each expanding at its own speed</text>']
    o.append(f'<rect x="{x0+70}" y="{Y(90):.1f}" width="{PW-90}" height="{Y(35)-Y(90):.1f}" fill="#efe8d8" opacity="0.55"/>')
    for L in (50,100,150,200,250):
        o.append(f'<line x1="{x0+70}" y1="{Y(L):.1f}" x2="{x0+PW}" y2="{Y(L):.1f}" stroke="{GRID}"/>')
        o.append(f'<text x="{x0+62}" y="{Y(L)+4:.1f}" text-anchor="end" font-size="12" fill="{Q}">{L}</text>')
    crossed=0; rnd=random.Random(seed)
    for i in range(24):
        k=K0*math.exp(rnd.gauss(0,SIG)); p=traj(L0,k)
        hot=max(l for _,l in p)>=150
        if hot: crossed+=1
        d=' '.join(('M' if j==0 else 'L')+f'{X(t):.1f} {Y(l):.1f}' for j,(t,l) in enumerate(p))
        o.append(f'<path d="{d}" fill="none" stroke="{ORANGE if hot else BLUE}" stroke-width="{2 if hot else 1.4}" opacity="{0.95 if hot else 0.5}"/>')
    o.append(f'<line x1="{x0+70}" y1="{Y(150):.1f}" x2="{x0+PW}" y2="{Y(150):.1f}" stroke="{INK}" stroke-width="2" stroke-dasharray="4 3"/>')
    o.append(f'<text x="{x0+PW}" y="{Y(150)-8:.1f}" text-anchor="end" font-size="12" font-weight="600" fill="{INK}">150: this cell begins to change</text>')
    # accumulation curve on a right-hand scale
    pts=[(t, share(L0,t)) for t in [j*0.5 for j in range(1,145)]]
    d=' '.join(('M' if j==0 else 'L')+f'{X(t):.1f} {Ys(s):.1f}' for j,(t,s) in enumerate(pts))
    o.append(f'<path d="{d}" fill="none" stroke="{ACC}" stroke-width="3" stroke-dasharray="7 4"/>')
    for s in (0.5,1.0):
        o.append(f'<text x="{x0+PW+8}" y="{Ys(s)+4:.1f}" font-size="11.5" fill="{ACC}">{int(s*100)}%</text>')
    o.append(f'<text x="{x0+PW+8}" y="{Ys(0)+4:.1f}" font-size="11.5" fill="{ACC}">0%</text>')
    sh=share(L0,onset(L0)); ax=X(onset(L0))
    o.append(f'<line x1="{ax:.1f}" y1="{Y0-6}" x2="{ax:.1f}" y2="{Y(35):.1f}" stroke="{INK}" stroke-width="1.5" stroke-dasharray="2 3"/>')
    o.append(f'<text x="{ax:.1f}" y="{Y(35)+20:.1f}" text-anchor="middle" font-size="12" font-weight="600" fill="{INK}">{onset(L0):.0f}</text>')
    o.append(f'<line x1="{x0+70}" y1="{Y(35):.1f}" x2="{x0+PW}" y2="{Y(35):.1f}" stroke="{Q}" stroke-width="1.5"/>')
    for t in (0,20,40,60):
        o.append(f'<text x="{X(t):.1f}" y="{Y(35)+20:.1f}" text-anchor="middle" font-size="12" fill="{Q}">{t}</text>')
    o.append(f'<text x="{x0+70+(PW-90)/2:.1f}" y="{Y(35)+40:.1f}" text-anchor="middle" font-size="13" fill="{INK}">Age in years</text>')
    o.append(f'<circle cx="{X(0):.1f}" cy="{Y(L0):.1f}" r="4" fill="{INK}"/>')
    o.append(f'<text x="{X(0)+8:.1f}" y="{Y(L0)-10:.1f}" font-size="12" font-weight="600" fill="{INK}">starts at {L0}</text>')
    o.append(f'<text x="{x0}" y="{Y0+PH+68}" font-size="13" font-weight="600" fill="{INK}">{crossed} of these 24 cells pass 150 within a lifetime. Diagnosis is expected at {onset(L0):.0f}, with {sh*100:.0f}% of cells past 150 by then</text>')
    return ''.join(o), crossed

a,ca=build(0,42,17); b,cb=build(720,50,17)
head=(f'<text x="34" y="{Y0+PH//2}" transform="rotate(-90 34 {Y0+PH//2})" text-anchor="middle" font-size="13" fill="{INK}">Repeat length in that one cell</text>'
      f'<text x="754" y="{Y0+PH//2}" transform="rotate(-90 754 {Y0+PH//2})" text-anchor="middle" font-size="13" fill="{INK}">Repeat length in that one cell</text>')
legend=(f'<line x1="0" y1="436" x2="26" y2="436" stroke="{ACC}" stroke-width="3" stroke-dasharray="7 4"/>'
        f'<text x="34" y="440" font-size="13" fill="{INK}">Dashed green line, right-hand scale: the share of that person’s neurons that have passed 150 by that age</text>')
foot=(f'<text x="0" y="466" font-size="13" fill="{Q}">Shaded band: the slow phase, up to about 90 repeats. Each line is one neuron. The shapes and the spread follow the two-stage dynamic the two papers describe; the paths are drawn, not measured.</text>'
      f'<text x="0" y="486" font-size="13" fill="{Q}">Expected age at diagnosis is 649/(CAG−30) years, from the CAG-Age Product as Enroll-HD standardizes it (CAP = 100 at the expected age of diagnosis; Warner et al. 2020).</text>')
open(os.path.join(HERE, 'panel.svg'), 'w').write(head+a+b+legend+foot)

# ---- third chart: accumulation only, six starting lengths
CW=1392; CH=330; cY0=66; cPH=210
X=lambda t: 100+t/80*(CW-200)
Y=lambda s: cY0+cPH-s*cPH
o=[f'<text x="0" y="22" font-size="17" font-weight="600" fill="{INK}">Below about 36 repeats nothing happens; above it, every added repeat moves the whole curve left</text>',
   f'<text x="0" y="42" font-size="13" fill="{Q}">Share of one person’s neurons that have passed 150, by age, for six starting lengths</text>']
for s in (0,0.25,0.5,0.75,1.0):
    o.append(f'<line x1="100" y1="{Y(s):.1f}" x2="{CW-100}" y2="{Y(s):.1f}" stroke="{GRID}"/>')
    o.append(f'<text x="92" y="{Y(s)+4:.1f}" text-anchor="end" font-size="12" fill="{Q}">{int(s*100)}%</text>')
COL={40:'#2f6db5',50:'#d9822b',60:'#b8473a'}
o.append(f'<line x1="100" y1="{Y(0):.1f}" x2="{CW-100}" y2="{Y(0):.1f}" stroke="#7b5ea7" stroke-width="3"/>')
o.append(f'<text x="170" y="{Y(0)-12:.1f}" font-size="12.5" font-weight="600" fill="#7b5ea7">starts at 10, 20 or 30: stable, no cell reaches 150</text>')
for L0 in (40,50,60):
    pts=[(t,share(L0,t)) for t in [j*0.5 for j in range(1,161)]]
    d=' '.join(('M' if j==0 else 'L')+f'{X(t):.1f} {Y(s):.1f}' for j,(t,s) in enumerate(pts))
    o.append(f'<path d="{d}" fill="none" stroke="{COL[L0]}" stroke-width="2.5"/>')
    oa=onset(L0); sh=share(L0,oa)
    o.append(f'<circle cx="{X(oa):.1f}" cy="{Y(sh):.1f}" r="5" fill="{COL[L0]}"/>')
    o.append(f'<line x1="{X(oa):.1f}" y1="{Y(sh):.1f}" x2="{X(oa):.1f}" y2="{Y(0):.1f}" stroke="{COL[L0]}" stroke-width="1" stroke-dasharray="2 3"/>')
    o.append(f'<text x="{X(oa):.1f}" y="{Y(sh)-12:.1f}" text-anchor="middle" font-size="11.5" font-weight="600" fill="{COL[L0]}">diagnosis at {oa:.0f}</text>')
    yl=share(L0,79)
    o.append(f'<text x="{CW-96}" y="{Y(yl)+4:.1f}" font-size="12.5" font-weight="600" fill="{COL[L0]}">starts at {L0}</text>')
for t in (0,20,40,60,80):
    o.append(f'<text x="{X(t):.1f}" y="{Y(0)+20:.1f}" text-anchor="middle" font-size="12" fill="{Q}">{t}</text>')
o.append(f'<text x="{(100+CW-100)/2:.1f}" y="{Y(0)+42:.1f}" text-anchor="middle" font-size="13" fill="{INK}">Age in years</text>')
o.append(f'<text x="34" y="{cY0+cPH//2}" transform="rotate(-90 34 {cY0+cPH//2})" text-anchor="middle" font-size="13" fill="{INK}">Share of neurons past 150</text>')
o.append(f'<text x="0" y="{cY0+cPH+68}" font-size="13" fill="{Q}">Dots mark the expected age at diagnosis, 649/(CAG−30) years. The curves follow from the drawn expansion speeds and are not fitted to onset data; that they put diagnosis at a similar share in all three cases is a</text>')
o.append(f'<text x="0" y="{cY0+cPH+88}" font-size="13" fill="{Q}">consequence of the drawn model, not evidence for it. Lengths of 36 to 39 are shown here as stable, though in that range the disease sometimes appears late in life.</text>')
open(os.path.join(HERE, 'panel2.svg'), 'w').write(''.join(o))
print(ca,cb,[round(share(L,onset(L))*100) for L in (40,42,50,60)])
