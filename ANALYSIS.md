# Palbociclib resistance in ER+ breast cancer — analysis code and results

Code and result tables for the manuscript *"An interferon-stimulated gene programme marks
palbociclib resistance in ER-positive breast cancer models: nested machine-learning validation
and structure-based repurposing against STAT1."*

Everything in `results/` and `figures/` was produced by the code in `code/` from the public
datasets listed below. No value in the manuscript was entered by hand.

## Data sources

| Source | Identifier | What was taken | How |
|---|---|---|---|
| GEO | GSE130437 | NCBI-generated raw counts (GRCh38.p13); 6 ER+ MCF7 samples | `ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts` |
| GEO | GSE222367 | NCBI-generated raw counts; 31 palbociclib-related samples | same |
| GEO | GSE229235 | authors' raw count matrix (supplementary file) | GEO FTP |
| GEO | GSE93204 | series matrix + GPL6480 annotation (NeoPalAna) | GEO FTP |
| UCSC Xena | TCGA-BRCA | HiSeqV2 expression, clinical matrix, curated survival | `tcga-xena-hub.s3.amazonaws.com` |
| cBioPortal | brca_metabric | mRNA for the candidate genes, patient and sample clinical data | REST API |
| STRING | v12.0 | interaction network of the 378 robust genes (score ≥ 0.4) | REST API |
| Enrichr | — | MSigDB Hallmark 2020, GO BP 2023, KEGG 2021, Reactome 2022 | `maayanlab.cloud/Enrichr/geneSetLibrary` |
| RCSB PDB | 1BF5 | STAT1 dimer bound to DNA | `files.rcsb.org` |
| ChEMBL | — | 3,311 approved small molecules (max_phase = 4) | REST API |
| GEO | GSE158724 | FELINE trial: 10x cancer-cell counts, cell-type labels, sample metadata | GEO FTP |
| NCBI | RefSeq / BLAST | transcript records and `refseq_rna` specificity search for qPCR primer design | E-utilities, BLAST URL API |

Sample assignments were read from the GEO **sample records**, not guessed from titles. Excluded:
ER-negative MDA-MB-231 samples (GSE130437), abemaciclib-resistant samples (GSE222367) and
triple-negative xenografts (GSE229235).

## Run order

| Step | Script | Produces |
|---|---|---|
| 1 | `s1_deg.py` | per-dataset differential expression, robust rank aggregation, centred matrix |
| 2 | `s2_enrich_network.py` | GO/KEGG/Hallmark over-representation, preranked GSEA, STRING network, hubs |
| 3 | `s3_ml.py` | LASSO, SVM-RFE, RF-Boruta, gradient boosting, core signature, apparent AUCs |
| 4 | `s7_nested.py` | **fully nested leave-one-dataset-out validation** |
| 5 | `s4_survival.py` | TCGA and METABRIC cohorts, LASSO-Cox score, KM, Cox, time-AUC, nomogram, DCA |
| 6 | `s5_neopalana.py` | paired on-treatment scores in the NeoPalAna trial |
| 6b | `s8_feline.py` | **external validation in ribociclib-treated FELINE patients** |
| 7 | `d1_receptor.py` | STAT1 SH2 receptor, biological dimer, grid centre |
| 8 | `d2_ligands.py` | ligand library preparation (RDKit + Meeko) |
| 9 | `d5_validation.py` | redocking validation of the grid |
| 10 | `d3_screen.py` | stage-1 virtual screen |
| 11 | `d4_stage2.py` | stage-2 triplicate re-docking and contact analysis |
| 12 | `f0…f6_*.py` | Figures 1–7 |

`stats_utils.py` and `surv_utils.py` hold the statistical methods implemented directly:
empirical Bayes moderated *t* (Smyth 2004), robust rank aggregation (Kolde 2012),
hypergeometric over-representation, preranked GSEA, Cox regression (Breslow partial likelihood),
LASSO-Cox by proximal gradient with cross-validated deviance, Kaplan–Meier, log-rank,
Harrell's C, IPCW time-dependent AUC and survival decision-curve analysis. These were checked
against simulated data with known effects before use.

## Environment

```
python 3.11
numpy scipy pandas scikit-learn matplotlib networkx
rdkit 2026.03.6   vina 1.2.7   meeko 0.8.0   gemmi 0.7.5   openmm 8.6.1   pdbfixer 1.12.0
```

Seed 2026 throughout. The screen was run on 2 CPU cores; on more cores, raise `EXH` in
`d3_screen.py` from 4 to 8 and relax the library filter in the same file.

## Key results

| Result | Value |
|---|---|
| Robust genes across three models | 378 (264 up, 114 down) |
| Leading Hallmark set (up) | interferon gamma response, 30/146, *P* = 2.5 × 10⁻²¹ |
| Depleted sets | E2F targets (NES −2.39), G2-M checkpoint (−1.94) |
| Network hubs | 14, all interferon-stimulated; STAT1 degree 44 |
| Core signature | 17 genes (≥ 2 of 3 algorithms) |
| Apparent AUC | 1.00 |
| **Nested leave-one-dataset-out AUC** | **0.65 (0.39–0.90)** |
| TCGA risk score | HR 1.53 per SD (1.19–1.96), *P* = 7.9 × 10⁻⁴ |
| METABRIC | HR 0.92 per SD (0.85–1.00) — does not replicate |
| Added value over clinical model | none (C 0.5986 vs 0.5984, *P* = 0.79) |
| NeoPalAna, E2F targets C1D1 → C1D15 | −0.16, *P* = 0.0039 |
| NeoPalAna, resistance score | unchanged, *P* = 0.45 |
| **FELINE, baseline resistance score vs response (ribociclib)** | **AUC 0.63, *P* = 0.42 — does not predict** |
| FELINE, interferon hub / risk score at baseline | AUC 0.44 / 0.53 |
| FELINE, E2F targets day 0 → day 14 | letrozole −0.47 (*P* = 0.014); + ribociclib −0.10 (*P* = 0.095) |
| Redocking validation | 1.24 Å (16 Å box centred on the phosphate subsite) |

## qPCR validation assay (`qpcr/`)

A complete, specificity-verified RT-qPCR design for the eight core genes and four candidate
reference genes: primer thermodynamics from SantaLucia nearest-neighbour parameters,
exon-junction placement from RefSeq records, NCBI BLAST against `refseq_rna` and an in-silico
PCR test that found **no off-target gene for any of the twelve pairs**. Includes the protocol,
a 384-well plate layout and an efficiency-corrected ΔΔCq analysis script with MIQE quality
gates. **This is a design, not data — no amplification was performed.**

## What is still missing

Wet-lab work, which cannot be done computationally: running the qPCR assay above in a
palbociclib-resistant line, and viability or biophysical testing of any docking hit. The
docking results are hypothesis-generating only — Vina scores are not binding affinities, and
the SH2 phosphotyrosine groove is a hard target for small molecules.
