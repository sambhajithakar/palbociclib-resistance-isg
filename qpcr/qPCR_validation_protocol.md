# RT–qPCR validation protocol for the interferon-stimulated core signature

**Status: this is a validated experimental *design*, not experimental data.**
The primers below were designed and checked computationally (nearest-neighbour
thermodynamics, exon-junction placement, NCBI BLAST against `refseq_rna`, and an
in-silico PCR test). No amplification has been performed. The protocol is written
so that it can be executed directly in the laboratory and the resulting data
analysed with the accompanying script.

---

## 1. Purpose

To confirm, at the transcript level in an independent palbociclib-resistant
derivative, the direction of change reported for the core genes of the
interferon-stimulated programme identified by the integrative analysis.

**Primary hypothesis.** *STAT1*, *ISG15*, *LGALS3BP* and *PARP12* are expressed at
higher levels in palbociclib-resistant ER+ cells than in the isogenic
palbociclib-sensitive parental line.

**Secondary.** *CASP4*, *SH3TC1*, *LAMB1* and *MARCKS* change in the same
direction as in the discovery meta-analysis.

## 2. Cell models

| Line | Description | Selection |
|---|---|---|
| MCF7-par | parental ER+, palbociclib-sensitive | — |
| MCF7-PR | palbociclib-resistant derivative | stepwise palbociclib, ≥6 months, maintained in 1 µM |
| T47D-par | parental ER+, palbociclib-sensitive | — |
| T47D-PR | palbociclib-resistant derivative | as above |

Resistance must be confirmed by a dose–response assay (palbociclib 0–10 µM, 6 days)
and the resistance index reported as IC50(PR)/IC50(par). Withdraw palbociclib for
72 h before RNA extraction so that the measurement reflects the resistant state
rather than acute drug exposure. Three independent biological replicates
(separate passages, separate extractions) per line.

## 3. RNA and cDNA

1. Extract total RNA at 70–80 % confluence (column-based kit with on-column DNase I).
2. Assess integrity: A260/A280 1.9–2.1, A260/A230 > 2.0; RIN ≥ 8 where a
   bioanalyser is available.
3. Reverse-transcribe 1 µg total RNA in 20 µL with random hexamers plus oligo-dT.
4. For every RNA sample prepare a **no-reverse-transcriptase (−RT)** control. This is
   mandatory for *ISG15*, whose amplicon cannot span an exon junction (see §4).
5. Dilute cDNA 1:10 in nuclease-free water; use 2 µL per 10 µL reaction.

## 4. Primers

All sequences 5′→3′. Tm calculated with the SantaLucia (1998) unified
nearest-neighbour parameters at 0.2 µM primer, 50 mM Na⁺, 3 mM Mg²⁺, 0.8 mM dNTP,
with the Owczarzy (2004) salt correction.

| Gene | Role | RefSeq | Forward | Rev | Amplicon (bp) | Tm F/R (°C) | Junction | Off-target genes |
|---|---|---|---|---|---|---|---|---|
| PARP12 | target | NM_022750.4 | GCAAACCTGCAATACCAAG | CAGACCTTCTGATACTCTTCC | 129 | 60.5 / 60.5 | primer spans | 0 |
| CASP4 | target | NM_033306.3 | TGCTTTCTGCTCTTCAACG | CTGTACCTTCCGAAATACTTCC | 139 | 61.7 / 61.6 | primer spans | 0 |
| SH3TC1 | target | NM_001410712.1 | TTTCATCAGTGGGCTCTTAG | CTCGGAAGGTCATCTCTTC | 151 | 60.4 / 60.4 | primer spans | 0 |
| LAMB1 | target | NM_002291.3 | GTGGTCACTACATTTGCTC | GGACGGAATGTCTTGAAAG | 140 | 59.3 / 59.3 | primer spans | 0 |
| STAT1 | target | NM_001437284.1 | GTATTACTCCAGGCCAAAGG | GGGTGAACTTCAGACACAG | 105 | 61.4 / 61.5 | primer spans | 0 |
| ISG15 | target | NM_005101.4 | GCAGATCACCCAGAAGATC | AGAGGTTCGTCGCATTTG | 156 | 61.1 / 61.1 | **none — use −RT** | 0 |
| LGALS3BP | target | NM_005567.4 | TGGTCTGCACCAATGAAAC | CTGGCTGTCAAAGATCTGG | 89 | 61.7 / 61.6 | primer spans | 0 |
| MARCKS | target | NM_002356.7 | TGCCCAGTTCTCCAAGAC | GCCGTTTACCTTCACGTG | 124 | 62.9 / 62.4 | amplicon spans | 0 |
| ACTB | reference | NM_001101.5 | CAAGATCATTGCTCCTCCTG | CATACTCCTGCTTGCTGATC | 107 | 61.7 / 61.7 | primer spans | 0 |
| GAPDH | reference | NM_001357943.2 | ATGACAACAGCCTCAAGATC | CTGTGGTCATGAGTCCTTC | 114 | 61.0 / 61.0 | primer spans | 0 |
| RPLP0 | reference | NM_053275.4 | GAAATCCTGAGTGATGTGC | GGAGATGTTGAGCATGTTC | 87 | 59.4 / 59.4 | primer spans | 0 |
| TBP | reference | NM_001172085.2 | AAGACCATTGCACTTCGTG | CTGGACTGTTCTTCACTCTTG | 152 | 61.7 / 61.7 | primer spans | 0 |

*ISG15* is a two-exon gene whose coding sequence lies almost entirely within
exon 2; no primer pair meeting the thermodynamic constraints can span the single
junction. The −RT control therefore carries the burden of excluding genomic
amplification for this gene.

### 4.1 In-silico specificity test

The 24 oligonucleotides were submitted to NCBI BLAST (`refseq_rna`, *Homo sapiens*,
`blastn-short`, word size 7). Hits with ≥ 95 % identity over the full primer
length were treated as annealing-competent. A pair was scored as producing an
off-target product only when the forward primer annealed to the plus strand and
the reverse primer to the minus strand of the same transcript within 50–2000 nt.
Sixty-seven transcripts satisfied this condition across the twelve pairs; NCBI
`esummary` annotation confirmed that **every one is an isoform of the intended
gene**. No pair amplifies a transcript of any other gene. The isoform counts
(e.g. 12 *STAT1*, 25 *SH3TC1* transcripts) reflect broad isoform coverage, which is
desirable for total-gene quantification.

## 5. Reaction and cycling

10 µL reactions in 384-well format (scale to 20 µL for 96-well):

| Component | Volume | Final |
|---|---|---|
| 2× SYBR Green master mix | 5.0 µL | 1× |
| Forward primer (10 µM) | 0.2 µL | 200 nM |
| Reverse primer (10 µM) | 0.2 µL | 200 nM |
| cDNA (1:10) | 2.0 µL | — |
| Nuclease-free water | 2.6 µL | — |

Cycling: 95 °C 2 min; 40 × (95 °C 15 s, 60 °C 30 s with acquisition);
melt curve 65 → 95 °C in 0.5 °C increments.

Three technical replicates per biological replicate. Include on every plate:
no-template control (NTC), −RT control, and an inter-run calibrator (a fixed
aliquot of pooled parental cDNA) so that plates can be scaled to one another.

## 6. Mandatory quality controls (MIQE)

1. **Standard curve** per assay from a 5-point, 10-fold serial dilution of pooled
   cDNA, in triplicate. Accept only assays with efficiency 90–110 %
   (slope −3.58 to −3.10) and R² ≥ 0.99. Report E and R² for every assay.
2. **Melt curve**: a single peak per assay; the NTC must show no peak at the
   product Tm, or a Cq at least 5 cycles later than the least abundant sample.
3. **−RT control**: Cq ≥ 35 or ≥ 5 cycles above the matched +RT reaction.
4. **Reference-gene stability**: rank *ACTB*, *GAPDH*, *RPLP0* and *TBP* by geNorm M
   and NormFinder stability across all samples; normalise to the geometric mean of
   the two most stable. Do **not** assume *ACTB* or *GAPDH* is stable — CDK4/6
   inhibition slows proliferation and can alter both.
5. Technical replicates with SD(Cq) > 0.5 are re-run.
6. Cq > 35 is reported as "not detected" and excluded from quantification.

## 7. Analysis

Relative quantification by the ΔΔCq method with efficiency correction
(Pfaffl 2001), parental line as calibrator, geometric mean of the two most stable
reference genes as normaliser. Statistics on log2 fold-change across the three
biological replicates: two-sided one-sample *t*-test against 0 within each cell
model, and a paired test across models; Benjamini–Hochberg correction across the
eight target genes. The run script `ddct_analysis.py` in this folder performs all
of this from a plate-export CSV.

Expected direction, from the discovery meta-analysis: all eight targets up in the
resistant line. The validation is considered successful if at least the four
primary genes are significantly up (FDR < 0.05) in both cell models with a
log2 fold-change concordant in sign with the meta-analysis estimate.

## 8. Reporting

Report per MIQE: cell provenance and passage, RNA quality metrics, RT enzyme and
priming, primer sequences with RefSeq accession and amplicon length, master mix
and instrument, efficiency and R² per assay, Cq of NTC and −RT, reference-gene
selection method, number of biological and technical replicates, and the
normalisation and statistical model.

## 9. Files in this folder

| File | Content |
|---|---|
| `qpcr_primers_final.csv` | the twelve validated pairs with thermodynamics and specificity |
| `qpcr_primers_all_candidates.csv` | up to three alternative pairs per gene |
| `qpcr_insilico_pcr.csv` | in-silico PCR result per pair |
| `blast_primers_tabular.txt` | raw NCBI BLAST output (RID BB04ERFB014) |
| `plate_layout_template.csv` | 384-well layout for one biological replicate |
| `ddct_analysis.py` | efficiency-corrected ΔΔCq analysis with QC gates |
| `design_primers.py` | the primer-design code |
| `insilico_pcr.py` | the BLAST parsing / in-silico PCR code |

## References

SantaLucia J (1998) *PNAS* 95:1460–1465.
Owczarzy R *et al.* (2004) *Biochemistry* 43:3537–3554.
von Ahsen N *et al.* (2001) *Clin Chem* 47:1956–1961.
Pfaffl MW (2001) *Nucleic Acids Res* 29:e45.
Vandesompele J *et al.* (2002) *Genome Biol* 3:research0034.
Andersen CL *et al.* (2004) *Cancer Res* 64:5245–5250.
Bustin SA *et al.* (2009) *Clin Chem* 55:611–622 (MIQE).
