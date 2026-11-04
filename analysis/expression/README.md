# Repair-gene expression against expansion achieved

The expression matrices in the archive (`raw.dge.rds`) are **dense** R matrices, not
sparse. The cortex one is 39,560 genes by 78,855 nuclei, which is 1.4 GB gzipped on
disk and 25 GB once R loads it. Only a high-memory machine can hold it.

`rds_stream.py` avoids that. R writes a matrix one column at a time, so the file can be
read as a stream and each nucleus added into its group as it arrives, with peak memory
of one column plus the result. The cortex matrix reduces in 53 seconds in under a
gigabyte. `validate_reader.py` checks it against R: row and column sums match exactly,
and the column sums equal the transcript counts recorded in the released per-nucleus
table.

`scan_matrix.py` makes a first pass for gene names, cell barcodes and per-column totals.
Nuclei are matched to the summary table **by barcode**, never by position: the striatal
summary has 6,359 rows against the matrix's 9,426 columns, and the two gene lists differ
between regions (39,560 genes in cortex, 39,350 in striatum).

`repair_genes_vs_expansion.py` does the analysis. Each donor's expansion is the 90th
percentile repeat length in its neurons less the length inherited. Each gene's per-donor
level is standardized: nuclei are binned by how far their own repeat has moved from the
inherited length, a rate is computed inside each bin, and the bins are recombined with
one set of weights shared by every donor, so each donor's figure describes the same
mixture of cells. Without this, the share of a donor's nuclei still near its inherited
length correlates with expansion at −0.98 and contaminates everything.

All 17,516 expressed genes go through the identical treatment, which is what makes the
result interpretable: FAN1 at rho −0.92 ranks first of 17,516 in cortex, and none of the
1,935 genes within a factor of 1.3 of its abundance comes close. Holding inherited
length constant it is −0.82, p = 0.004. The striatum shows nothing at −0.03, on six
donors of whom five inherited 42 repeats.

Cause and consequence cannot be separated here. FAN1 protein has been reported lower in
HD cortex than in controls, so a decline driven by disease stage would produce the same
correlation.

## Running it

Python with numpy, pandas and scipy. R is needed only to convert the two `summary.rds`
files to CSV. Extract `donor_metadata.txt` and the two `deepdive_dge_data` directories
from the Zenodo archive into the working directory, then:

```
python3 scan_matrix.py <path>/cortex/glutamatergic/raw.dge.rds 39560 cortex_scan.pkl
python3 scan_matrix.py <path>/caudate/SPN/raw.dge.rds          39350 spn_scan.pkl
python3 repair_genes_vs_expansion.py
```
