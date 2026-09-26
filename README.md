# Pharmacogenomic Mapping of an NSCLC Quiescence Signature

Analysis code and derived results for the Discover Oncology manuscript mapping the
688-gene quiescence signature (Foglietta et al., *Int J Mol Sci* 2022;23(17):9869,
PMC9456317; 676 genes recovered) against DepMap 26Q1, GDSC2, PRISM 20Q2, and
TCGA-LUAD/LUSC, plus a pre-registered proliferation-stratified follow-up.

## Layout

- `manuscript/draft.md` — full manuscript
- `figures/` — Figures 1–8 (PNG + PDF)
- `scripts/08_proliferation_stratified.py` — quartile stratification + dependency/drug reanalysis
- `scripts/09_figure7_update.py` — regenerates Figure 7 evidence table
- `data/signature/signature_676.csv` — recovered signature (Table S1)
- `data/depmap/` — enrichment CSVs, result JSONs, CCLE TPM JSONs, DepMap Model metadata
- `data/gdsc/PortalCompounds.csv`, `data/prism/*cell_line_info.csv`,
  `data/prism/*treatment_info.csv` — compound/cell-line metadata
- `data/tcga/*.json`, `data/depmap/*.json` — specificity, survival, control summaries

## Large inputs (not committed; fetch to reproduce)

- `data/depmap/CRISPRGeneEffect_26Q1.csv` (440 MB) —
  `https://huggingface.co/datasets/ChanghaoKan/crispr-depmap/resolve/main/CRISPRGeneEffect_26Q1.csv`
- `data/gdsc/GDSC2_AUC_Matrix.csv` (3.9 MB) — same mirror, `GDSC2_AUC_Matrix.csv`
- `data/prism/prism_20q2_secondary_dose_response.csv` (290 MB) —
  `https://ndownloader.figshare.com/files/36794595`
- TCGA STAR counts — GDC API `POST /data` with file UUIDs
  (`cases.project.project_id` TCGA-LUAD/TCGA-LUSC, STAR-Counts workflow)
- CCLE TPM — cBioPortal API, study `ccle_broad_2025`
- MDPI supplement zip — `https://mdpi-res.com/d_attachment/ijms/ijms-23-09869/article_deploy/ijms-23-09869-s001.zip`

## Environment

Python 3.14, pandas 2.3.3, numpy 2.4.6, scipy 1.17.1, lifelines 0.30.3,
matplotlib 3.10.8, statsmodels 0.14.6.

## Key results

| Analysis | P-value | Call |
|---|---|---|
| DepMap NSCLC enrichment (permutation) | 0.8821 | Null |
| GDSC2 drug sensitivity | 0.5450 | Null |
| PRISM drug sensitivity | <0.0001 | Nominal (promiscuity caveat) |
| TCGA specificity / survival | 1.0 / 0.92, 0.47 | Null |
| Proliferation-stratified dependency (24 vs 24) | 0.5315 | Null |
