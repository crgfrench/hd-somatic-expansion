"""First pass over a matrix: gene names, cell barcodes and per-column totals."""
import numpy as np, pandas as pd, sys, os, time, pickle
from rds_stream import reduce_matrix, peek_shape
path, nrow, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t0 = time.time()
R = reduce_matrix(path, nrow=nrow)
pickle.dump(dict(genes=R["genes"], cells=R["cells"], colsum=R["colsum"], ncol=R["ncol"]),
            open(out, "wb"))
print("%s: %d x %d in %.0f s; genes %s cells %s"
      % (out, nrow, R["ncol"], time.time()-t0,
         len(R["genes"]) if R["genes"] else None,
         len(R["cells"]) if R["cells"] else None))
