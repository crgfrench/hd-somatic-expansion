import numpy as np, time
from rds_stream import reduce_matrix
f = "HD_cortex/data/deepdive_dge_data/caudate/SPN/raw.dge.rds"
nrow = 39350
import gzip, struct
# one group for everything, so num[:,0] is the row sums
t0=time.time()
g = np.zeros(9426, dtype=int)
num, colsum, ncol, genes = reduce_matrix(f, nrow, g, 1)
print("streamed %d x %d in %.0f s" % (nrow, ncol, time.time()-t0))
np.save("spn_rowsums_stream.npy", num[:,0]); np.save("spn_colsums_stream.npy", colsum)
print("rowsum head:", num[:5,0], "total:", num[:,0].sum())
print("colsum head:", colsum[:5])
print("genes parsed:", (genes[:3] if genes else None))
