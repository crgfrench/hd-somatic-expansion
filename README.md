# Somatic CAG expansion in Huntington's disease: a re-analysis of released data

Code and figures from a re-analysis of the data released with Kashin et al.,
*Cell* (2026) (cortex) and Handsaker et al., *Cell* 188:623 (2025)
(striatum). Everything here was computed from the public archives; no new
experimental data is involved.

Craig W. French, Princeton, New Jersey · cf3357@princeton.edu ·
ORCID 0009-0001-8847-3018

## The figures

`figures/Repeat_Disease_Thresholds.pdf` opens in the browser from here, one page
per board, each page sized to the board it holds. `figures/Repeat_Disease_Thresholds.html` is the same content as a
single self-contained page, which a browser renders once downloaded. Seven
boards:

1. Repeat lengths across 47 expansion disorders, on a logarithmic scale, with
   the Huntington's disease row marked and the 150-CAG cell-level threshold drawn
   against them.
2. What other repeat disorders suggest about where a threshold can come from.
3. The distribution of the shortest disease-causing length across those 47.
4. The re-analysis of the two papers' released data, including which phase D
   genes sit together on the chromosome.
5. Expansion timing: what the donor records allow to be measured, and what
   slowing the expansion would buy.
6. FAN1 and the pace of expansion: the MutSβ–FAN1 competition, FAN1 against
   expansion across the 11 cortex donors, where that correlation sits among all
   17,516 genes, and the two regions compared.
7. When the silencing fails: four tests that score each nucleus by its own
   repeat length rather than by its donor.

## What was run

### Per-gene onset fits — `notebooks/01_phaseC_gene_onsets.ipynb`

For each phase C gene, a breakpoint Poisson regression of per-nucleus counts
against CAG length, with an offset for that donor-and-subtype reference rate.
Ninety-five percent intervals come from the profile deviance with an
overdispersion scale.

Striatal start points cluster between 100 and 170 CAGs. Cortical ones spread
from the shortest lengths measured, and 71 rising genes have no detectable start
point — they change steadily from 40 CAGs upward. 117 genes are shared between
the regions, their onsets correlated at Spearman 0.23.

### PRC2 — `analysis/prc2_enrichment.py`

Fisher exact tests of the phase C and phase D gene sets (p < 10⁻¹⁰ in the
released per-gene statistics) against the H3K27me3-marked and
knockout-induced lists of von Schimmelmann et al., *Nat Neurosci* 19:1321
(2016). Phase D genes are largely the set a neuron switches on when PRC2 is
removed: 18.4% of cortical phase D genes against 0.47% of the background.
Removing every phase D gene from both sides leaves a cortical phase C
enrichment resting on three genes.

### Where the phase genes sit — `analysis/phase_gene_clustering.py`

Whether the phase C and phase D gene sets are bunched along the chromosome,
which a domain-level failure of polycomb silencing would produce. The answer is
specific rather than general: 11 of the 23 genes named HOXA to HOXD with a
paralog number in the cortex statistics
are phase D genes against 0.22 expected, all four clusters with their embedded
lncRNAs, and 4 of 17 in the striatum against 0.54. The only other multi-gene run
is MEG3 with MEG8 at the imprinted DLK1-DIO3 locus. Set the clusters aside and
the remaining phase D genes are no more bunched than chance.

The script reports three tests and the order matters. The spatial test against a
null drawn from every tested gene is biased and is printed only so the bias is
visible: the phase sets run near 190 transcripts per million against 17 for all
tested genes, while clustered paralog families are quiet. Matching the null on
abundance changes cortical phase C from p = 0.66 to p = 0.002. It cannot be
applied to the HOX result, because all 11 of those genes sit below the 1 TPM
floor of the abundance table, so the gene-family count is the evidence there.

A GRCh38 gene table with columns `symbol`, `chr` and `start` is needed, and an
abundance table measured in the same region as the statistics it matches.

### Repair genes and expansion — `analysis/expression/`

Whether a donor's own DNA-repair gene levels track how far that donor's repeat
travelled. Across the 11 cortex donors FAN1 gives Spearman rho = −0.92, ranking first of
17,516 genes put through identical standardization; holding inherited length constant it
is −0.82, p = 0.004. The striatum shows nothing, on six donors with almost no variation
in inherited length. `analysis/expression/README.md` has the method, the controls and
the limits, including the reverse-causation problem a post-mortem snapshot cannot settle.

### When the silencing fails — `analysis/silencing_tests.py`

Four tests on the per-nucleus data, each scoring a nucleus by the repeat it
carries rather than by the donor it came from.

**A.** All four HOX clusters are above their own baseline across 140–199 CAGs.
Each cluster is taken relative to its own rate below 100 copies, since HOXD's
baseline is four times the others and an absolute scale reads that difference as
a difference in onset. Pooled over 140–199 the factors are 3.7, 2.0, 5.3 and 6.5
for HOXA to HOXD, on 27, 9, 25 and 121 transcripts, each above 1 at p = 1.6e-08,
0.04, 4.1e-11 and 1.1e-55; beyond 260 all four run 50 to 130 times baseline.
Four polycomb domains lifting inside one 60-copy window, which contains the 150
both papers identify, is what a failure of the shared silencing would produce. A
common onset length is past what these counts can settle: the four factors differ
from one another at p = 0.0001, and the two bands where HOXB records nothing
expected 1.5 transcripts to begin with.

Panel A plots fold change for that reason, and `figures/src/hox_cluster_bands.json`
holds the per-cluster transcript counts by band that it draws.

**B.** At a fixed repeat length, a nucleus's HOX level still tracks how far it
has gone in phase C: past 230 CAGs, rho = +0.32 with the rising genes and −0.25
with the identity genes, holding the released control gene set, read depth and
repeat length constant, on 891 nuclei. The control set is what separates this
from a drift in transcriptome quality — across the same bands it moves +0.13
while the identity genes move −0.34. Phase C and phase D therefore advance together rather than
as two clocks both reading off repeat length, though a correlation of this kind
leaves open whether one process drives both.

**C.** The MutSβ-to-FAN1 balance does not predict expansion between cell types
of one person, though it appears to between people. Measured only in nuclei
whose repeat has moved under 10 copies, and compared across subtypes inside a
donor, the mean correlation is −0.20 over six donors, opposite to the sign
the mechanism predicts; in the one donor whose subtypes span a wide range of
expansion it is −0.28 at p = 0.0005, with subtype size explaining neither
variable. This was meant to answer the sample size behind the FAN1 result and
does not. Two features of the design limit what the null result rules out.
Restricting the measurement to barely expanded nuclei samples a different
fraction of each subtype, since in a fast-expanding subtype those nuclei are its
stragglers. And cell types differ in chromatin state, transcription rate and
much else, so the comparison credits MutSβ and FAN1 with whatever those other
differences do.

**D.** The two regions derepress the same genes at different repeat lengths.
Across the 42 HOX genes present in both expression matrices — the 39 named
HOXA to HOXD with a paralog number, plus HOXA10-AS, HOXA11-AS and HOXC13-AS —
levels in nuclei carrying 260 copies or more correlate at rho = +0.86. But at 140–199 CAGs the cortex sits at 5.2 times its
own baseline (p = 3 × 10⁻⁶⁸) while the striatum sits at 0.46 — and that striatal
band holds 18 million transcripts, enough that a cortex-sized rise would have
shown about 75 HOX counts where 7 were seen. Both reach the same place past 260
CAGs, 102 times baseline against 61.

Gene rows are pulled from the matrices by
`analysis/expression/extract_gene_rows.py`, which also recovers the column
barcodes, needed because the striatal summary has 6,359 rows against 9,426
matrix columns.

### Subtype thresholds — `notebooks/02_subtype_thresholds.ipynb`, `analysis/subtype_thresholds.py`

Whether a neuron subtype's threshold tracks its expression of candidate
regulators. The result is null across 28 subtypes and the design is
underpowered; it is included for completeness.

### Expansion rate — `notebooks/03`, `analysis/expansion_*.py`, `analysis/age_versus_expansion.py`

An attempt to estimate the exponent *p* in d(repeat)/dt ∝ ((repeat − 35)/10)^*p*
from the per-nucleus length distributions, using age at death as elapsed time.

The exponent is not identifiable in this cohort, and the reason is the
substantive result. Across the 54 HD donors in `donor_metadata.txt`, age at
death and inherited CAG length correlate at Spearman rho = −0.79,
p = 2 × 10⁻¹². Age at death is an outcome of the expansion rather than an
independent exposure time, so the donors with the most elapsed years are the
ones whose repeats grow slowest. Within the 11 donors with per-nucleus
lengths, age predicts the expansion achieved at rho = −0.11 (p = 0.74) while
inherited length predicts it at rho = +0.85 (p = 0.001). Measuring the rate
would need donors who died of unrelated causes across a span of ages, or
repeated samples from living people.

`expansion_exponent_likelihood.py` fits the full length distribution and pins at
the edge of the grid; `expansion_exponent_from_age.py` identifies the exponent
instead by requiring the implied rate to be age-independent, and
`expansion_diagnostics.py` shows that the striatal fit is worse with age than
without it. Both are kept so the failure can be inspected.

### Slowing and delay — `analysis/selection_and_delay.py`

In any growth law, rate and time appear only as a product, so cutting the rate by
a fraction *f* multiplies every remaining interval by 1/(1 − *f*). The anchor is
observed clinical onset: 649/(CAG − 30) years, the Enroll-HD calibration of the
CAG-Age Product. A person born with 42 repeats has an expected onset at 54, and
treatment begun at 30 needs 29% slowing to buy ten years. Begun at 50 the same
ten years needs 71%. This holds whatever the exponent is.

## Reproducing the page

The seven boards live in `figures/src/`. Four are written by hand and are the
source of record for themselves; three are generated. From a clean checkout:

```
cd figures/src
python3 board_timing.py          # writes Timing.dc.html
python3 board_fan1.py            # writes FAN1.dc.html
python3 board_silencing.py       # writes Silencing.dc.html
python3 panel_hox.py             # rewrites one panel inside Reanalysis.dc.html
python3 build_standalone_page.py # assembles all seven into figures/…html
python3 build_pdf.py             # renders that to PDF, one board per page
```

Each generator reads a JSON file beside it that an analysis script wrote, so the
drawings follow from the data rather than being redrawn by hand: `board_timing.py`
from `donor_table.json`, `board_fan1.py` from `repair_genes.json` as written by
`analysis/expression/repair_genes_vs_expansion.py`, `board_silencing.py` from
`silencing_data.json` as written by `analysis/silencing_tests.py` together with
`hox_cluster_bands.json`, and
`panel_hox.py` from `hox_clusters.json` as written by
`analysis/phase_gene_clustering.py`, rewriting the region between the
`hox:start` and `hox:end` markers in the re-analysis board. That panel orders the
genes within each cluster by paralog number, which is the order they occupy on
the chromosome, rather than by the annotation start: the gene table keeps one row
per symbol, and for HOXC6 that row belongs to an upstream isoform whose start
falls between HOXC10 and HOXC9. Assembling the
page needs Python and Node; the PDF step also needs Playwright with a browser
(`npm install playwright && npx playwright install chromium`). Running it this
way reproduces `figures/Repeat_Disease_Thresholds.html` byte for byte.

Two further scripts, `board_reanalysis.py` and `board_expansion_panels.py`,
write SVG fragments rather than whole boards. Those fragments were pasted into
the two hand-written boards they belong to, and the scripts are kept so the
drawings can be regenerated. `board_reanalysis.py` needs a per-gene onset
summary as JSON, which notebook 01 produces and which is not committed here.

## Running it

Python 3.11 with `numpy`, `scipy` and `pandas`. R with `Matrix` for the
notebooks, which were run on Colab against the Zenodo archive.

```
python3 analysis/prc2_enrichment.py /path/to/input/files
```

With no argument the scripts read from `data/`. The input files each script
needs are listed in `data/PUT_INPUT_FILES_HERE.txt`. Two of them,
`per_nucleus_repeat_lengths.csv.gz` and `phaseC_gene_onsets.csv`, are written
by notebook 01 rather than downloaded.


## Sources

- Kashin, Lee et al., *Cell* (2026), doi:10.1016/j.cell.2026.10.001. Zenodo record 18202664.
- Handsaker et al., *Cell* 188:623 (2025).
- von Schimmelmann et al., *Nat Neurosci* 19:1321 (2016), supplementary tables 1 and 3.
- Depienne and Mandel, *Am J Hum Genet* 108:764 (2021), for the disease table.
