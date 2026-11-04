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
single self-contained page, which a browser renders once downloaded. Six
boards:

1. Repeat lengths across 47 expansion disorders, on a logarithmic scale, with
   the Huntington's disease row marked and the 150-CAG cell-level threshold drawn
   against them.
2. What other repeat disorders suggest about where a threshold can come from.
3. The distribution of the shortest disease-causing length across those 47.
4. The re-analysis of the two papers' released data.
5. Expansion timing: what the donor records allow to be measured, and what
   slowing the expansion would buy.
6. FAN1 and the pace of expansion: the MutSβ–FAN1 competition, FAN1 against
   expansion across the 11 cortex donors, where that correlation sits among all
   17,516 genes, and the two regions compared.

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

### Repair genes and expansion — `analysis/expression/`

Whether a donor's own DNA-repair gene levels track how far that donor's repeat
travelled. Across the 11 cortex donors FAN1 gives Spearman rho = −0.92, ranking first of
17,516 genes put through identical standardization; holding inherited length constant it
is −0.82, p = 0.004. The striatum shows nothing, on six donors with almost no variation
in inherited length. `analysis/expression/README.md` has the method, the controls and
the limits, including the reverse-causation problem a post-mortem snapshot cannot settle.

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
ones whose repeats grow slowest. Within the 11 donors carrying per-nucleus
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
a fraction *f* multiplies every remaining interval by 1/(1 − *f*). Anchored on
observed clinical onset — 649/(CAG − 30) years, the Enroll-HD calibration of the
CAG-Age Product — a person born with 42 repeats has an expected onset at 54, and
treatment begun at 30 needs 29% slowing to buy ten years. Begun at 50 the same
ten years needs 71%. This holds whatever the exponent is.

## Running it

Python 3.11 with `numpy`, `scipy` and `pandas`. R with `Matrix` for the
notebooks, which were run on Colab against the Zenodo archive.

```
python3 analysis/prc2_enrichment.py /path/to/input/files
```

With no argument the scripts read from `data/`. The input files each script
needs are listed in `data/PUT_INPUT_FILES_HERE.txt`. Two of them —
`per_nucleus_repeat_lengths.csv.gz` and `phaseC_gene_onsets.csv` — are written
by notebook 01 rather than downloaded.

## Rebuilding the figures

`figures/src/board_timing.py` writes the expansion-timing board and
`figures/src/board_fan1.py` the FAN1 board, the latter from `repair_genes.json`
as written by `analysis/expression/repair_genes_vs_expansion.py`.
`build_standalone_page.py` assembles the six boards into the single HTML page,
expanding the disease table through `expand_template.js` (Node). `build_pdf.py`
then renders that page to PDF, one board per page at the board's own height,
through `build_pdf.js` (Playwright) and pypdf.

## Sources

- Kashin, Lee et al., *Cell* (2026), doi:10.1016/j.cell.2026.10.001. Zenodo record 18202664.
- Handsaker et al., *Cell* 188:623 (2025).
- von Schimmelmann et al., *Nat Neurosci* 19:1321 (2016), supplementary tables 1 and 3.
- Depienne and Mandel, *Am J Hum Genet* 108:764 (2021), for the disease table.
