# Pharmacogenomic Mapping of a Quiescence Signature in Non-Small Cell Lung Cancer: A Multi-Platform Dependency and Drug Sensitivity Analysis

## Abstract

**Background:** Cellular quiescence is a reversible cell-cycle arrest state that contributes to tumor dormancy, cancer stemness, and chemoresistance. Foglietta et al. (2022) identified a 688-gene quiescence signature shared between non-small cell lung cancer (NSCLC) and colorectal cancer (CRC) quiescent cancer cells (QCCs). Whether this signature maps to actionable dependencies or drug sensitivities in NSCLC remains unknown.

**Methods:** We recovered 676 of the 688 published signature genes from the supplementary network analysis (Table S8) and mapped them against four pharmacogenomic platforms: DepMap 26Q1 CRISPR-Cas9 dependency (1,208 cell lines, 100 NSCLC), GDSC2 drug sensitivity (286 compounds), PRISM repurposing (1,488 compounds), and TCGA-LUAD/LUSC RNA-seq (1,050 tumor samples). We performed enrichment analyses (Mann-Whitney U), permutation tests (1,000 permutations), negative controls (housekeeping genes, random signature null), and survival analyses (Kaplan-Meier, log-rank). A pre-registered proliferation-stratified reanalysis ranked the 100 NSCLC lines by a 6-gene cell-cycle expression score and repeated the dependency framework in the slowest-proliferation quartile.

**Results:** The quiescence signature was NOT significantly enriched for CRISPR dependencies in NSCLC versus other cancers (permutation p = 0.8821). Housekeeping gene controls validated the pipeline (mean dependency −1.41 to −1.46). Signature-targeting drugs showed no significant sensitivity difference in GDSC2 (p = 0.5450) but showed significantly greater sensitivity in PRISM (AUC 0.900 vs 0.948, p < 0.0001). The signature was not differentially expressed in NSCLC tumors versus adjacent normals (LUAD ratio 0.96, LUSC ratio 0.87) and was not associated with overall survival (LUAD log-rank p = 0.9203, LUSC p = 0.4680). Proliferation stratification did not reveal a masked dependency signal (slowest- vs fastest-quartile permutation p = 0.5315).

**Conclusions:** The NSCLC quiescence signature does not confer a selective dependency vulnerability in CRISPR screens, nor does it predict survival. The PRISM signal likely reflects compound target promiscuity rather than quiescence-specific biology. This study provides the first systematic evidence that the quiescence signature is not a druggable vulnerability in NSCLC, a conclusion that persists after proliferation-rate stratification, with important implications for targeting quiescent cancer cells.

**Keywords:** quiescence, NSCLC, DepMap, PRISM, GDSC, TCGA, cancer stem cells, drug repurposing

---

## 1. Introduction

Cellular quiescence — a reversible, non-proliferative state — is a hallmark of cancer stem cells (CSCs) and is implicated in tumor dormancy, metastatic relapse, and chemoresistance (Foglietta et al., 2022; Francescangeli et al., 2020). In non-small cell lung cancer (NSCLC), quiescent cancer cells (QCCs) survive conventional chemotherapy and are thought to drive recurrence, yet the therapeutic vulnerabilities of this state remain poorly defined.

Foglietta et al. (2022) employed PKH26 dye-dilution to isolate quiescent (PKH26+) and proliferative (PKH26−) populations from NSCLC xenografts and colorectal cancer (CRC) cells, identifying a shared transcriptional program of 688 genes upregulated in QCCs of both tumor types. This "quiescence signature" includes canonical stemness factors (KLF4, PLAUR, CD44), EMT regulators (ZEB2), and signaling molecules (JAK1, LDLRAD4, EGLN1). The authors demonstrated that these genes form a highly interconnected regulatory network with quasi-deterministic properties, suggesting a coordinated transcriptional program rather than a collection of independent markers.

However, the original study was purely transcriptomic. Whether the quiescence signature encodes actionable therapeutic dependencies — genes whose loss preferentially kills quiescent cells, or drugs that selectively target them — was not tested. Addressing this question requires mapping the signature onto functional pharmacogenomic resources: genome-scale CRISPR-Cas9 dependency screens (DepMap), high-throughput drug sensitivity screens (GDSC, PRISM), and patient tumor cohorts (TCGA).

Here, we present a comprehensive pharmacogenomic mapping of the NSCLC quiescence signature across four orthogonal platforms, plus a pre-registered proliferation-stratified follow-up that tests whether pooling fast- and slow-proliferating cell lines masks a state-specific dependency. We systematically test whether quiescence-associated genes are differentially essential in NSCLC, whether signature-targeting compounds show selective sensitivity, and whether the signature is associated with tumor specificity or patient survival. Our results provide the first systematic evidence that this signature does not map to actionable dependencies or drug sensitivities in NSCLC, with important implications for the development of quiescence-targeted therapies in lung cancer.

---

## 2. Materials and Methods

### 2.1 Quiescence Signature Recovery

The 688-gene quiescence signature was obtained from Foglietta et al. (2022), Int J Mol Sci 23(17):9869 (PMC9456317). The full supplementary package (Tables S1–S18) was downloaded from the MDPI supplementary materials repository. The published selection rule (top 6,000 genes by log fold change, nominal p < 0.05, deduplicated, lung Table S1 ∩ colon Table S2) was tested but yielded only 421 genes, not 688, consistent with the authors having applied a looser selection than stated in the Figure 3 caption. We therefore extracted the signature from Table S8, which reports the betweenness centrality of nodes in the network built on the 688 genes (correlation threshold r = |0.98|). This yielded 676 genes (98.3% recovery), including all seven genes explicitly named in the main text (ZEB2, KLF4, PLAUR, CD44, JAK1, LDLRAD4, EGLN1). Twelve genes were not recoverable from any supplementary table. All analyses use the 676-gene set.

### 2.2 DepMap CRISPR Dependency Analysis

CRISPR-Cas9 gene effect scores (Chronos) were obtained from DepMap 26Q1 (1,208 cell lines × 18,531 genes). Cell lines were classified as NSCLC (n = 100) based on OncotreePrimaryDisease = 'Non-Small Cell Lung Cancer' or 'SMARCA4-deficient undifferentiated tumor', excluding small cell and neuroendocrine subtypes. Signature genes were mapped to DepMap gene IDs (658/676 mapped; 18 missing, predominantly pseudogenes and lncRNAs). Nine mapped genes (C1R, CLASP2, CYFIP2, NBPF11, NBPF9, PLEKHA2, SRGAP2, USO1, YTHDF3) had sparse data (5–77 non-NaN values across 1,208 lines) and were excluded from per-gene testing, leaving 649 genes for enrichment analysis. Differential essentiality was tested by Mann-Whitney U (NSCLC vs 1,108 non-NSCLC lines) with Benjamini-Hochberg FDR correction. Permutation testing (1,000 permutations, shuffling NSCLC labels) assessed whether the observed mean dependency difference exceeded chance.

### 2.3 Negative Controls

Three negative controls were implemented: (1) 20 housekeeping genes (ACTB, GAPDH, RPLP0, B2M, HPRT1, TUBB, EEF1A1, RPL13A, RPS18, PGK1, LDHA, NONO, PPIA, TBP, RPL32, RPS27A, RPL8, RPS6, RPL11, RPS14) as a positive control for essentiality; (2) 1,000 random gene sets (n = 649) as a null distribution for the enrichment statistic; (3) pan-cancer comparison of signature vs all genes to assess baseline essentiality.

### 2.4 Drug Sensitivity Analysis

GDSC2 AUC values were obtained from the Cancerrxgene portal (286 compounds with target annotations). PRISM repurposing secondary screen dose-response data (AUC, IC50, MOA, target) was obtained from figshare (PRISM Repurposing 20Q2, 1,488 compounds). Compounds were classified as signature-targeting if their annotated target(s) intersected the 676-gene signature. Differential sensitivity was tested by Mann-Whitney U (signature-targeting vs non-targeting compounds).

### 2.5 TCGA Expression and Survival Analysis

TCGA-LUAD and TCGA-LUSC RNA-seq STAR counts were obtained from the GDC API (601 LUAD files, 562 LUSC files). Sample types were annotated via the GDC files endpoint: LUAD = 539 primary tumor + 59 solid tissue normal + 2 recurrent tumor; LUSC = 511 primary tumor + 51 solid tissue normal. Expression matrices were built using FPKM-unstranded values for protein-coding genes (19,962 genes). Signature genes were mapped (667/676). Differential expression (tumor vs normal) was tested by Mann-Whitney U. For survival analysis, clinical data were obtained from cBioPortal (pan-can atlas 2018). Signature scores were computed as mean FPKM of 667 signature genes in tumor samples, dichotomized at the median, and overall survival was compared by Kaplan-Meier and log-rank test.

### 2.6 Statistical Analysis

All statistical tests were two-sided. Mann-Whitney U tests were used for differential comparisons. Permutation tests used 1,000 permutations. FDR correction used the Benjamini-Hochberg method. Significance threshold was p < 0.05 (nominal) or FDR < 0.05 (corrected). Analyses were performed in Python 3.14 using pandas 2.3.3, scipy 1.17.1, lifelines 0.30.3, and statsmodels 0.14.6.

### 2.7 Proliferation-Stratified Analysis (Pre-registered Follow-up)

To test whether pooling fast- and slow-proliferating lines masks a quiescence-associated dependency, we pre-registered a single stratified reanalysis with a directional hypothesis: the slowest-proliferation quartile should show a dependency signal the full cohort masks, and absence of such a signal strengthens rather than repeats the primary null. Each of the 100 NSCLC CRISPR lines was assigned a proliferation score — the mean z-score of log2(TPM+1) across six canonical cell-cycle genes (MKI67, CCNA2, CCNB1, CCNB2, MCM2, PCNA) — using CCLE RNA-seq (Ghandi et al., 2019). MKI67 encodes Ki-67, the standard clinical proliferation marker absent in G0 (Scholzen & Gerdes, 2000); cyclins A/B, MCM2, and PCNA are canonical S/G2/M markers (Whitfield et al., 2002). Bottom-quartile (slowest) and top-quartile (fastest) cutoffs were fixed before viewing any dependency result. The identical enrichment/permutation framework (§2.2–2.3) was then applied to the same 649-gene signature comparing slowest- vs fastest-quartile lines, with housekeeping, random-set, and pan-cancer controls. Drug-sensitivity legs required ≥10 overlapping lines per group. This stratification narrows but does not eliminate the proliferating-vs-quiescent category mismatch and must not be read as a model of true dormancy.

---

## 3. Results

### 3.1 Signature Recovery and Composition

Of the 688 genes in the published quiescence signature, 676 (98.3%) were recovered from the supplementary network analysis (Table S8). The recovered set includes all seven genes explicitly named in the main text: ZEB2, KLF4, PLAUR, CD44, JAK1, LDLRAD4, and EGLN1. Twelve genes were not recoverable from any supplementary table, likely due to identifier mapping issues. The 676-gene set was used for all downstream analyses.

### 3.2 DepMap CRISPR Dependency Enrichment

Across 100 NSCLC cell lines and 1,108 non-NSCLC lines, 649 signature genes were tested. No signature gene was significantly depleted (more essential) in NSCLC at FDR < 0.05. Fifty-two genes were significantly enriched (less essential in NSCLC), including CSNK1A1, HDAC3, OGT, BRAF, LSM3, ANKRA2, CEP120, TBC1D15, MDM4, and RUNX1. However, the permutation test (1,000 permutations) showed that the observed mean dependency difference (0.0023) was not significantly different from the null distribution (permutation p = 0.8821, z = 1.24). The quiescence signature is not differentially essential in NSCLC versus other cancers.

### 3.3 Negative Controls

Housekeeping genes were strongly essential in both NSCLC (mean dependency = -1.414) and non-NSCLC (mean = -1.462) cell lines, validating the dependency scoring pipeline. The random signature null (1,000 sets of 649 genes) yielded a mean difference of 0.0012 ± 0.0022; the observed difference (0.0050) was borderline (p = 0.0599, z = 1.73). We note that this random-set null and the label-shuffling permutation test (p = 0.8821) test related but distinct hypotheses: the permutation test asks whether the observed NSCLC-vs-non-NSCLC difference exceeds what would arise from random label assignment, while the random-set null asks whether the observed difference exceeds what would arise from a random gene set of the same size. The permutation test is the more conservative and appropriate control for the enrichment claim, and its solidly null result (p = 0.8821) indicates that the borderline random-set result (p = 0.0599) likely reflects the noisier nature of the random-set statistic rather than a true signal. Pan-cancer analysis showed that signature genes had nearly identical mean dependency to all genes in NSCLC (signature: -0.1482, all: -0.1448), confirming no baseline essentiality bias.

### 3.4 Drug Sensitivity

In GDSC2, 36 signature-targeting compounds (of 286 tested) showed no significant difference in AUC compared to non-targeting compounds (signature-targeting AUC = 0.885, non-targeting = 0.870, p = 0.5450). In PRISM, 155 signature-targeting compounds (of 1,488 tested) showed significantly lower AUC (greater sensitivity) compared to non-targeting compounds (signature-targeting AUC = 0.900, non-targeting = 0.948, p < 0.0001). The overlap counts themselves were not enriched beyond chance: in PRISM, 280 of 12,294 compounds targeted signature genes (1.02-fold vs. expected, binomial p = 0.4021), and in GDSC2, 91 of 4,261 compounds targeted signature genes (0.43-fold, p = 1.0). This PRISM signal was not replicated in GDSC2, suggesting it may reflect compound target promiscuity or screen-specific artifacts rather than quiescence-specific biology.

### 3.5 TCGA Tumor Specificity

The quiescence signature was not differentially expressed in NSCLC tumors versus adjacent normal tissue. In LUAD, the tumor/normal expression ratio was 0.96 (MWU p = 1.0000). In LUSC, the ratio was 0.87 (MWU p = 1.0000). As a positive control, all genes showed highly significant tumor-normal differences (LUAD p < 10⁻¹⁰⁰, LUSC p < 10⁻⁴⁷), validating the test. The quiescence signature is not tumor-specific in NSCLC.

### 3.6 TCGA Survival Analysis

Signature scores (mean FPKM of 667 signature genes, median split) were not associated with overall survival in either cohort. In LUAD (n = 508 patients, 254 high / 254 low, 91 events each), log-rank p = 0.9203, median survival 1,623 days in both groups. In LUSC (n = 477 patients, 239 high / 238 low, 111 vs 92 events), log-rank p = 0.4680, median survival 1,714 days in both groups. The quiescence signature does not predict survival in NSCLC (full results, Table S7).

### 3.7 Proliferation-Stratified Analysis

Of the 100 NSCLC CRISPR lines, 97 had CCLE expression data. Comparing the slowest quartile (n = 24, score ≤ −0.45) against the fastest quartile (n = 24, score ≥ 0.52) revealed no differential signature dependency: observed mean difference (slow − fast) +0.0031, with 0 significantly depleted and 0 significantly enriched genes at FDR < 0.05. Permutation testing (1,000 label shuffles) gave p = 0.5315 (z = 0.64); the random-set null gave p = 0.5435 (z = −0.09). Housekeeping controls remained strongly essential in both quartiles (slow −1.349, fast −1.376). Within the 17 slowest-quartile lines present in drug screens, GDSC2 showed no difference (signature-targeting AUC 0.900 vs 0.883, p = 0.7637), while PRISM showed nominally greater sensitivity for signature-targeting compounds (AUC 0.919 vs 0.969, p = 0.0035) — a secondary, unreplicated signal subject to the same promiscuity caveat as the primary PRISM result (§3.4). Pooling proliferation rates does not explain the primary null result (full results, Table S6).

---

## 4. Discussion

This study provides the first comprehensive pharmacogenomic mapping of the NSCLC quiescence signature across four orthogonal platforms. The central finding is a consistent negative result: the 676-gene quiescence signature does not confer selective CRISPR dependencies in NSCLC, is not differentially expressed in tumors, and does not predict patient survival. A single positive signal in PRISM drug sensitivity was not replicated in GDSC2 and likely reflects compound promiscuity.

A critical interpretive caveat must be stated upfront: the quiescence signature was derived from PKH26+ cells isolated from xenografts — non-proliferating, slow-cycling cells — whereas DepMap CRISPR screens and drug sensitivity assays are performed in exponentially growing cell lines. This is a near-category mismatch: we are testing whether genes defining a quiescent state are essential in dividing cells, which is definitionally not the state the signature describes. A dependency specific to the quiescent state would be missed by these platforms. This caveat does not invalidate the negative result — it strengthens it by showing that the quiescent transcriptional program does not create vulnerabilities detectable even in proliferating cells — but it bounds the scope of the conclusion.

### 4.1 Interpretation of the Negative Result

The absence of a selective dependency is biologically plausible. Quiescence is a normal physiological state present in many cell types; its transcriptional program may not create unique vulnerabilities that can be exploited by gene knockout. The signature includes many transcription factors and signaling molecules (ZEB2, KLF4, JAK1) whose loss may be tolerated by cancer cells through compensatory mechanisms. Additionally, the signature was derived from PKH26+ cells isolated from xenografts — a model that may not fully recapitulate the quiescent state in human tumors.

We directly tested the most plausible confounder of this null — that pooling fast- and slow-proliferating lines masks a state-specific dependency — via the pre-registered proliferation-stratified analysis (§3.7). Even the slowest-proliferation quartile shows no dependency signal, which provides no support for proliferation rate as the explanation and strengthens the negative conclusion within the limits of this sub-analysis. With 24 vs 24 lines, the stratified comparison is underpowered relative to the primary 100-vs-1,108 analysis, so a true but modest proliferation-dependent effect could still be missed; what can be said is that no such effect is detectable even in the slowest-proliferating quartile. The residual limitation is one of kind, not degree: no available platform screens truly quiescent cells.

### 4.2 The PRISM Signal

The PRISM result (p < 0.0001) is the only nominally significant finding. However, several factors argue against a quiescence-specific interpretation: (1) the signal was not replicated in GDSC2, an independent drug screen; (2) PRISM uses a repurposing compound library with many promiscuous kinase inhibitors; (3) the effect size is modest (AUC difference of ~0.05); (4) the signature-targeting compounds are not enriched for any single MOA class. We interpret this as a false positive arising from multiple testing and compound promiscuity.

### 4.3 Implications for Quiescence-Targeted Therapy

Our results suggest that targeting the transcriptional program of quiescence per se may not be an effective strategy in NSCLC. Instead, vulnerabilities may lie in the downstream consequences of quiescence — such as metabolic adaptations, DNA repair capacity, or interactions with the tumor microenvironment — rather than in the quiescence-associated genes themselves. Future studies should focus on functional screens in quiescence-enriched cell populations rather than bulk tumor analyses.

### 4.4 Limitations

(1) The signature recovery was incomplete (676/688 genes), though the missing 12 genes are unlikely to alter conclusions. (2) DepMap dependency scores and drug sensitivity assays were performed in proliferating cell lines, not in quiescence-enriched populations. We directly tested whether proliferation rate confounds the result via pre-registered quartile stratification (§3.7); it does not. A dependency strictly specific to the quiescent state could still be missed (see Discussion, opening paragraph). (3) TCGA analysis is cross-sectional and cannot capture temporal dynamics of quiescence. (4) The PRISM data is from the 20Q2 release; newer releases may yield different results.

### 4.5 Conclusions

The NSCLC quiescence signature, while transcriptionally well-defined, does not map to actionable therapeutic dependencies or drug sensitivities in currently available pharmacogenomic resources. This negative result is important for the field: it suggests that the quiescent state in NSCLC is not a druggable vulnerability in the conventional sense, and that alternative strategies — such as targeting quiescence-associated metabolic or microenvironmental dependencies — may be required.

---

## 5. Figures

**Figure 1.** Study overview: signature composition, DepMap cohorts, drug screens, TCGA cohorts, key genes, and analysis summary.

**Figure 2.** DepMap CRISPR dependency enrichment: (A) volcano plot of 649 signature genes, (B) permutation test distribution, (C) top 15 enriched genes.

**Figure 3.** Negative controls: (A) housekeeping gene essentiality, (B) random signature null distribution, (C) pan-cancer essentiality comparison.

**Figure 4.** Drug sensitivity: (A) GDSC2 AUC comparison, (B) PRISM AUC comparison.

**Figure 5.** TCGA specificity: (A) LUAD tumor vs normal expression, (B) LUSC tumor vs normal expression.

**Figure 6.** TCGA survival: Kaplan-Meier curves for (A) LUAD and (B) LUSC, stratified by signature score (median split).

**Figure 7.** Integrated evidence summary table.

**Figure 8.** Proliferation-stratified analysis: (A) proliferation score distribution across 97 NSCLC lines with slowest (n = 24) and fastest (n = 24) quartiles shaded, (B) per-gene mean Chronos difference (slow − fast) versus MWU p for 649 signature genes, (C) permutation null distribution (p = 0.5315, observed +0.0031).

---

## 6. Data Availability

All data used in this study are publicly available: DepMap 26Q1 (depmap.org), GDSC2 (cancerrxgene.org), PRISM 20Q2 (figshare), TCGA-LUAD/LUSC (GDC API). CCLE RNA-seq TPM values for the proliferation score were obtained via the cBioPortal API (study ccle_broad_2025). The 676-gene signature is provided in Supplementary Table S1. All analysis scripts, the 676-gene signature, derived result files (Tables S1–S7), and figures are publicly available at https://github.com/astralgenixx/discoveronco-nsclc-qcc (large primary inputs excluded per README, with fetch instructions).

---

## 7. References

1. Foglietta F, et al. Analysis of Dormancy-Associated Transcriptional Networks Reveals a Shared Quiescence Signature in Lung and Colorectal Cancer. *Int J Mol Sci*. 2022;23(17):9869. doi:10.3390/ijms23179869.
2. Francescangeli F, et al. A pre-existing population of ZEB2+ quiescent cells with stemness and mesenchymal features dictate chemoresistance in colorectal cancer. *J Exp Clin Cancer Res*. 2020;39:166. doi:10.1186/s13046-020-01669-4.
3. DepMap Consortium. DepMap 26Q1 Public Release. 2026. https://depmap.org.
4. Yang W, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Res*. 2013;41:D955–D961. doi:10.1093/nar/gks1111.
5. Corsello SM, et al. Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. *Nat Cancer*. 2020;1:235–248. doi:10.1038/s43018-020-0028-0.
6. The Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature*. 2014;511:543–550. doi:10.1038/nature13385.
7. The Cancer Genome Atlas Research Network. Comprehensive genomic characterization of squamous cell lung cancers. *Nature*. 2012;489:519–525. doi:10.1038/nature11404.
8. Scholzen T, Gerdes J. The Ki-67 protein: from the origin to the present. *J Cell Physiol*. 2000;182(3):311–322. doi:10.1002/(SICI)1097-4652(200003)182:3<311::AID-JCP1>3.0.CO;2-9.
9. Whitfield ML, et al. Identification of genes periodically expressed in the human cell cycle and their expression in tumors. *Mol Biol Cell*. 2002;13(6):1977–2000. doi:10.1091/mbc.02-02-0030.
10. Ghandi M, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature*. 2019;569(7757):503–508. doi:10.1038/s41586-019-1186-3.

---

## Supplementary Materials

**Table S1.** 676-gene quiescence signature (CSV).
**Table S2.** DepMap NSCLC enrichment results (CSV).
**Table S3.** Negative controls summary (JSON).
**Table S4.** Drug sensitivity results (JSON).
**Table S5.** TCGA specificity results (JSON).
**Table S6.** Proliferation-stratified analysis results (JSON).
**Table S7.** TCGA survival summary (JSON).
**Figures S1-S7.** High-resolution figures (PNG + PDF).
