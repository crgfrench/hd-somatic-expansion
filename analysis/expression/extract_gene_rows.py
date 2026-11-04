"""Pull a chosen set of gene rows out of a dense R expression matrix.

The released matrices are dense base-R doubles. The cortex one is
39,560 x 78,855, which is 25 GB if read whole, so rds_stream walks it a block
of columns at a time and each block is reduced immediately to the rows asked
for. What comes out is a CSV of counts per nucleus for those genes, plus that
nucleus's total, which is small enough to work with anywhere.

The column barcodes are read from the matrix's trailing dimnames attribute.
They are needed for the striatum, whose released summary carries fewer rows
than the matrix has columns, so nuclei there cannot be matched by position.
"""
import argparse, numpy as np, pandas as pd
from rds_stream import Reader, read_attributes

ap = argparse.ArgumentParser()
ap.add_argument('--matrix', required=True, help='raw.dge.rds for one region')
ap.add_argument('--genes', required=True, help='the matrix row names, one per line')
ap.add_argument('--rows', required=True, help='the genes to keep, one per line')
ap.add_argument('--out', required=True)
ap.add_argument('--block', type=int, default=256)
a = ap.parse_args()

genes = [l.rstrip('\n') for l in open(a.genes)]
wanted = [l.strip() for l in open(a.rows) if l.strip()]
index = {g: i for i, g in enumerate(genes)}
missing = [w for w in wanted if w not in index]
if missing:
    print('not in this matrix, skipped: %s' % ', '.join(missing))
keep_names = [w for w in wanted if w in index]
rows = np.array([index[w] for w in keep_names])
nrow = len(genes)

r = Reader(a.matrix)
r.header()
t, has_attr = r.flags()
assert t == 14, 'expected a double vector, got SEXP type %d' % t
ncol = r.length() // nrow
print('%d genes x %d nuclei, keeping %d rows' % (nrow, ncol, len(rows)), flush=True)

keep = np.zeros((len(rows), ncol))
total = np.zeros(ncol)
j = 0
while j < ncol:
    k = min(a.block, ncol - j)
    blk = np.frombuffer(r.f.read(nrow * 8 * k), dtype='>f8').reshape(k, nrow).T
    keep[:, j:j + k] = blk[rows, :]
    total[j:j + k] = blk.sum(axis=0)
    j += k

attrs = read_attributes(r) if has_attr else {}
dn = attrs.get('dimnames')
cells = list(dn[1]) if dn and len(dn) > 1 else ['col%d' % i for i in range(ncol)]
out = pd.DataFrame({'CELL_BARCODE': cells, 'umis': total})
for i, g in enumerate(keep_names):
    out[g] = keep[i, :].astype(int)
out.to_csv(a.out, index=False)
print('wrote %s, %d rows' % (a.out, len(out)), flush=True)
