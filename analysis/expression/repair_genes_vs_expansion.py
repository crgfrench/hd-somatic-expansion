import numpy as np, pandas as pd, pickle, time
from scipy.stats import rankdata, spearmanr
from rds_stream import reduce_matrix

meta = pd.read_csv("HD_cortex/data/donor_metadata.txt", sep="\t")
meta = meta[meta.STATUS == "pwHD"][["DONOR","AGE","HD_CAG_LENGTH"]]
BRK, MIN = [-np.inf, 2, 5, 10, 20, 40, np.inf], 20
PANEL = ["FAN1","MSH2","MSH3","MSH6","MLH1","MLH3","PMS1","PMS2","LIG1",
         "OGG1","NEIL1","POLB","XRCC1","PARP1","TDP1","HTT"]

def run(path, summary_csv, scan_pkl, label):
    s = pickle.load(open(scan_pkl, "rb"))
    genes = np.array(s["genes"]); cells = list(s["cells"])
    cpos = {b: i for i, b in enumerate(cells)}
    sd = pd.read_csv(summary_csv, low_memory=False)
    sd["col"] = sd.CELL_BARCODE.map(cpos)                 # matched by barcode, never by position
    assert sd.col.notna().all(), "%d summary rows missing from the matrix" % sd.col.isna().sum()
    sd["col"] = sd.col.astype(int)
    sd = sd.merge(meta, on="DONOR", how="left")
    sd["gain"] = sd.CAG_LENGTH - sd.HD_CAG_LENGTH
    ok = sd.gain.notna()
    d = sd[ok].groupby("DONOR").agg(n=("gain","size"), age=("AGE","first"),
          inherited=("HD_CAG_LENGTH","first"), expansion=("CAG_LENGTH", lambda v: np.quantile(v,.9)))
    d["expansion"] -= d["inherited"]
    don = d[d.n >= 200].sort_index()

    sd["bin"] = pd.cut(sd.gain, BRK, right=False)
    use = sd[ok & sd.DONOR.isin(don.index)].copy()
    cnt = use.pivot_table(index="bin", columns="DONOR", values="col", aggfunc="size",
                          fill_value=0, observed=False)
    usable = [b for b in cnt.index if cnt.loc[b].min() >= MIN]
    use = use[use.bin.isin(usable)]
    print("\n== %s: %d donors, %d nuclei of %d columns, bins kept: %s"
          % (label, len(don), len(use), len(cells), ", ".join(str(b) for b in usable)), flush=True)
    print(use.pivot_table(index="DONOR", columns="bin", values="col", aggfunc="size",
                          fill_value=0, observed=False)[usable].to_string(), flush=True)

    pairs = list(zip(use.DONOR, use.bin.astype(str)))
    levels = sorted(set(pairs)); lix = {p: i for i, p in enumerate(levels)}
    group = np.full(len(cells), -1, dtype=np.int64)
    group[use.col.values] = [lix[p] for p in pairs]

    t0 = time.time()
    R = reduce_matrix(path, group=group, ngroup=len(levels), nrow=len(genes))
    print("streamed in %.0f s" % (time.time()-t0), flush=True)
    den = np.zeros(len(levels)); np.add.at(den, group[use.col.values], R["colsum"][use.col.values])
    wN = use.bin.astype(str).value_counts()
    w = np.array([wN[str(b)] for b in usable], float); w /= w.sum()
    donors = list(don.index)
    std = np.zeros((len(genes), len(donors)))
    for bi, b in enumerate(usable):
        cols = [lix[(dn, str(b))] for dn in donors]
        std += (R["num"][:, cols] / den[cols]) * w[bi]
    tpm = R["num"].sum(axis=1) / R["num"].sum() * 1e6
    return don, donors, std, tpm, genes

def rho_all(don, donors, std, tpm, genes):
    keep = (tpm >= 1) & ((std > 0).sum(axis=1) >= len(donors) - 1)
    m = std[keep]
    rk = np.apply_along_axis(rankdata, 1, m)
    y = rankdata(don.expansion.values).astype(float); y -= y.mean(); y /= np.sqrt((y**2).sum())
    rx = rk - rk.mean(axis=1, keepdims=True); rx /= np.sqrt((rx**2).sum(axis=1, keepdims=True))
    return pd.DataFrame({"gene": genes[keep], "tpm": tpm[keep], "rho": rx @ y})

def report(A, label):
    print("\n---- %s: %d genes standardised ----" % (label, len(A)), flush=True)
    print(A[A.gene.isin(PANEL)].sort_values("rho")
           .to_string(index=False, float_format=lambda v: "%8.3f" % v), flush=True)
    print("\n%-7s %8s %9s %12s %16s" % ("gene","rho","TPM","genome-wide","same-abundance"))
    for g in ["FAN1","MSH2","NEIL1","PARP1","PMS2","MSH3"]:
        if g not in set(A.gene): continue
        o = float(A.loc[A.gene==g,"rho"].iloc[0]); t = float(A.loc[A.gene==g,"tpm"].iloc[0])
        pool = A[(A.tpm > t/1.3) & (A.tpm < t*1.3) & (A.gene != g)]
        print("%-7s %+8.2f %9.1f %11.2f%% %11.2f%%  (n=%d)" %
              (g, o, t, 100*np.mean(np.abs(A.rho) >= abs(o)),
               100*np.mean(np.abs(pool.rho) >= abs(o)), len(pool)), flush=True)

D = "HD_cortex/data/deepdive_dge_data"
res = {}
for key, path, csv, pkl in (
    ("cortex",   D+"/cortex/glutamatergic/raw.dge.rds", "cortex_glutamatergic_summary.csv", "cortex_scan.pkl"),
    ("striatum", D+"/caudate/SPN/raw.dge.rds",          "caudate_SPN_summary.csv",          "spn_scan.pkl")):
    don, donors, std, tpm, genes = run(path, csv, pkl, key)
    print(don.to_string(), flush=True)
    A = rho_all(don, donors, std, tpm, genes); A.to_csv("rho_%s_allgenes.csv" % key, index=False)
    gi = {g: i for i, g in enumerate(genes)}
    pd.DataFrame({g: std[gi[g]] for g in PANEL if g in gi}, index=donors).join(don) \
      .to_csv("%s_standardised.csv" % key)
    report(A, key.upper()); res[key] = A; del std

b = res["cortex"].merge(res["striatum"], on="gene", suffixes=("_cortex","_striatum"))
b[["gene","rho_cortex","rho_striatum"]].to_csv("rho_both_regions.csv", index=False)
print("\n---- both regions ----")
print(b[b.gene.isin(PANEL)][["gene","rho_cortex","rho_striatum"]]
      .sort_values("rho_cortex").to_string(index=False, float_format=lambda v: "%8.3f" % v))
print("\nacross %d genes in both, the two sets of rho correlate at %.2f"
      % (len(b), spearmanr(b.rho_cortex, b.rho_striatum).statistic))
print("\nDONE", flush=True)
