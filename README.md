# An interferon-stimulated gene programme marks palbociclib resistance in ER-positive breast cancer models

[![DOI](https://zenodo.org/badge/DOI/PASTE_DOI_HERE.svg)](https://doi.org/PASTE_DOI_HERE)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Analysis code, result tables and figures for the manuscript

> Thakar SB, Kumbhar NM, Sonawane KD. *An interferon-stimulated gene programme marks
> palbociclib resistance in ER-positive breast cancer models: nested machine-learning
> validation and structure-based repurposing against STAT1.* (under review)

Everything in `results/` and `figures/` was produced by the code in `code/` from public data.
No value in the manuscript was entered by hand.

---

## What this study found

Three independent palbociclib-resistant ER+ model systems converge on 378 reproducible
expression changes dominated by interferon-stimulated genes, with STAT1 as the most central
hub protein. A 17-gene signature built from them reaches an apparent cross-validated AUC of
**1.00** — and **0.65** once feature selection is nested inside each training fold. The
derived risk score is prognostic in TCGA, fails to replicate in METABRIC, and does not predict
response to ribociclib in the FELINE trial.

The negative validation results are reported in full. The gap between apparent and nested
performance is the point.

| Result | Value |
|---|---|
| Robust genes across three models | 378 (264 up, 114 down) |
| Leading Hallmark set (up) | interferon gamma response, 30/146, *P* = 2.5 × 10⁻²¹ |
| Network hubs | 14, all interferon-stimulated; STAT1 degree 44 |
| Apparent AUC | 1.00 |
| **Nested leave-one-dataset-out AUC** | **0.65 (0.39–0.90)** |
| TCGA risk score | HR 1.53 per SD (1.19–1.96), *P* = 7.9 × 10⁻⁴ |
| METABRIC | HR 0.92 per SD (0.85–1.00) — does not replicate |
| FELINE, baseline score vs ribociclib response | AUC 0.63, *P* = 0.42 — does not predict |
| Redocking validation (positive control) | 1.24 Å RMSD |

---

## Repository layout

```
code/      analysis scripts (see ANALYSIS.md for the run order)
results/   every result table, including Supplementary_Tables.xlsx
figures/   Figures 1-7 and S1-S3, as PNG and PDF
docking/   prepared STAT1 receptor and the docked poses of the leading compounds
qpcr/      validated RT-qPCR assay design, protocol and ddCq analysis script
```

`ANALYSIS.md` documents the data sources, the run order and the statistical methods.

## Quick start

```bash
git clone https://github.com/sambhajithakar/palbociclib-resistance-isg.git
cd palbociclib-resistance-isg
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

The scripts expect the raw downloads in a `data/raw/` directory; `ANALYSIS.md` lists every
accession and the exact endpoint each file came from. The result tables in `results/` let you
reproduce every number and figure in the paper without re-downloading anything.

## Data sources

All public. GEO: GSE130437, GSE222367, GSE229235, GSE93204, GSE158724.
TCGA-BRCA via UCSC Xena. METABRIC via cBioPortal. STRING v12.0. RCSB PDB 1BF5. ChEMBL.
The underlying datasets remain under their own terms of use; see the note in `LICENSE`.

## Citing

Please cite both the paper and this archive. `CITATION.cff` carries machine-readable metadata;
GitHub renders a "Cite this repository" button from it.

## Licence

MIT for the code (see `LICENSE`). Result tables and figures are released under
CC-BY-4.0. Redistributed third-party data remain under their original licences.
