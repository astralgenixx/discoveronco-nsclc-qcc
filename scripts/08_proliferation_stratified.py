"""Proliferation-stratified reanalysis (§3.7 addendum) — PRE-REGISTERED DESIGN.

Hypothesis (directional, falsifiable): if the null CRISPR result in the full
NSCLC cohort (§3.2) is driven by pooling fast- and slow-proliferating lines,
the slowest-proliferation quartile should show a detectable dependency signal
(low-prolif lines more dependent on signature genes) that the pooled analysis
masks. A null here closes the door on the 'wrong cell-state' explanation.

Pre-specified choices (fixed BEFORE looking at any dependency result):
  PROLIF_GENES = [MKI67, CCNA2, CCNB1, CCNB2, MCM2, PCNA] (Entrez 4288, 890,
    891, 9133, 4171, 5111) — Ki-67 is the clinical proliferation marker
    (Scholzen & Gerdes, J Cell Physiol 2000); cyclins/MCM2/PCNA are canonical
    periodically expressed cell-cycle genes (Whitfield et al., Mol Biol Cell
    2002). Expression: CCLE RNA-seq TPM via cBioPortal ccle_broad_2025.
  Score = mean of per-gene z-scores of log2(TPM+1), z-scored across the 100
    NSCLC lines. Rank; bottom quartile (n=25) = low-prolif, top (n=25) = high.
  Same 649-gene signature set and same framework as §2.2–2.3 (MWU + FDR,
    housekeeping check, 1,000-permutation label shuffle, 1,000 random-set null).
  Drug leg only if >=10 overlapping lines per group in GDSC2/PRISM; else
    abandon explicitly.

This does NOT model true quiescence — it stratifies proliferating lines by a
transcriptional proliferation proxy. Stated as such in Methods/Discussion.
"""

import json
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/sabujsamaddar/P_for_PhD/discoveronco/results/tcga"
print(
    "versions:",
    __import__("pandas").__version__,
    __import__("numpy").__version__,
    __import__("scipy").__version__,
    matplotlib.__version__,
)

# ---------- 1. Proliferation scores ----------
ENTREZ2SYM = {
    4288: "MKI67",
    890: "CCNA2",
    891: "CCNB1",
    9133: "CCNB2",
    4171: "MCM2",
    5111: "PCNA",
}
tpms = {}
for ez, sym in ENTREZ2SYM.items():
    recs = json.load(open(f"{BASE}/data/depmap/ccle_tpm_{ez}.json"))
    tpms[sym] = {r["sampleId"]: r["value"] for r in recs}
tpm_df = pd.DataFrame(tpms)  # CCLE sampleId x 6 genes
print("CCLE TPM samples:", tpm_df.shape)

model = pd.read_csv(f"{BASE}/data/depmap/Model_26Q1.csv")
NSCLC_DISEASES = {
    "Non-Small Cell Lung Cancer",
    "SMARCA4-deficient undifferentiated tumor",
}
crispr = pd.read_csv(f"{BASE}/data/depmap/CRISPRGeneEffect_26Q1.csv", index_col=0)
crispr.columns = [c.split(" (")[0].strip().upper() for c in crispr.columns]
# Identical cohort to §3.2: strict-NSCLC models WITH CRISPR data (n=100)
nsclc = model[
    model["OncotreePrimaryDisease"].isin(NSCLC_DISEASES)
    & model["ModelID"].isin(crispr.index)
].copy()
print("NSCLC CRISPR cohort:", len(nsclc))

# Map DepMap NSCLC lines -> CCLE TPM via CCLEName (= cBioPortal sampleId)
nsclc = nsclc[nsclc["CCLEName"].isin(tpm_df.index)].copy()
print("NSCLC lines with TPM:", len(nsclc))
sub = np.log2(tpm_df.loc[nsclc["CCLEName"], list(ENTREZ2SYM.values())] + 1)
z = (sub - sub.mean()) / sub.std(ddof=0)
score = z.mean(axis=1)
nsclc = nsclc.set_index("CCLEName")
nsclc["prolif_score"] = score
nsclc = nsclc.sort_values("prolif_score").reset_index()
n = len(nsclc)
q = n // 4
low = nsclc.iloc[:q]  # slowest quartile
high = nsclc.iloc[-q:]  # fastest quartile
print(
    f"low quartile n={len(low)} score<= {low['prolif_score'].max():.2f}; "
    f"high quartile n={len(high)} score>= {high['prolif_score'].min():.2f}"
)
print("LOW:", ", ".join(low["CellLineName"].tolist()))
print("HIGH:", ", ".join(high["CellLineName"].tolist()))
low_ids = low["ModelID"].tolist()
high_ids = high["ModelID"].tolist()

# ---------- 2. Dependency data (same 649 genes) ----------
enr = pd.read_csv(f"{BASE}/data/depmap/depmap_nsclc_enrichment.csv")
genes649 = enr["gene"].tolist() if "gene" in enr.columns else enr.iloc[:, 0].tolist()
print("signature genes reused:", len(genes649))
HOUSE = [
    "ACTB",
    "GAPDH",
    "RPLP0",
    "B2M",
    "HPRT1",
    "TUBB",
    "EEF1A1",
    "RPL13A",
    "RPS18",
    "PGK1",
    "LDHA",
    "NONO",
    "PPIA",
    "TBP",
    "RPL32",
    "RPS27A",
    "RPL8",
    "RPS6",
    "RPL11",
    "RPS14",
]

use_genes = [g for g in genes649 if g in crispr.columns]
use_house = [g for g in HOUSE if g in crispr.columns]
print(f"mapped: {len(use_genes)}/649 signature, {len(use_house)}/20 housekeeping")
L = crispr.loc[[i for i in low_ids if i in crispr.index], use_genes]
H = crispr.loc[[i for i in high_ids if i in crispr.index], use_genes]
print("low/high rows in CRISPR:", L.shape[0], H.shape[0])
Lh = crispr.loc[L.index, use_house]
Hh = crispr.loc[H.index, use_house]

# ---------- 3. Same framework: MWU per gene + FDR ----------
from statsmodels.stats.multitest import multipletests

pvals, diffs = [], []
for g in use_genes:
    a = L[g].dropna().values
    b = H[g].dropna().values
    if len(a) < 3 or len(b) < 3:
        pvals.append(np.nan)
        diffs.append(np.nan)
        continue
    _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    pvals.append(p)
    diffs.append(float(np.nanmean(a) - np.nanmean(b)))
pvals = np.array(pvals)
mask = ~np.isnan(pvals)
fdr = np.full_like(pvals, np.nan)
fdr[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
sig_dep = [
    g
    for g, f, d in zip(use_genes, fdr, diffs)
    if not np.isnan(f) and f < 0.05 and d < 0
]  # more essential in low
sig_enr = [
    g
    for g, f, d in zip(use_genes, fdr, diffs)
    if not np.isnan(f) and f < 0.05 and d > 0
]
obs_diff = float(np.nanmean(L.values) - np.nanmean(H.values))
print(
    f"obs mean diff (low-high): {obs_diff:.4f}; "
    f"depleted(FDR,low<high): {len(sig_dep)}; enriched: {len(sig_enr)}"
)

# ---------- 4. Negative controls ----------
# (a) housekeeping
hk_low = float(np.nanmean(Lh.values))
hk_high = float(np.nanmean(Hh.values))
# (b) permutation: shuffle low/high labels 1000x
rng = np.random.default_rng(7)
both = pd.concat([L, H])
labs = np.array([0] * len(L) + [1] * len(H))
nulls = np.empty(1000)
for i in range(1000):
    perm = rng.permutation(labs)
    nulls[i] = np.nanmean(both.values[perm == 0]) - np.nanmean(both.values[perm == 1])
perm_p = float((np.sum(np.abs(nulls) >= abs(obs_diff)) + 1) / (1000 + 1))
perm_z = float((obs_diff - nulls.mean()) / nulls.std(ddof=1))
# (c) random gene sets 1000 x 649
universe = [c for c in crispr.columns if c not in set(use_genes) | set(use_house)]
rand_diffs = np.empty(1000)
for i in range(1000):
    gs = rng.choice(universe, size=len(use_genes), replace=False)
    rand_diffs[i] = np.nanmean(crispr.loc[L.index, gs].values) - np.nanmean(
        crispr.loc[H.index, gs].values
    )
rand_p = float((np.sum(np.abs(rand_diffs) >= abs(obs_diff)) + 1) / 1001)
rand_z = float((obs_diff - rand_diffs.mean()) / rand_diffs.std(ddof=1))
# (d) baseline: signature vs all genes within low-prolif subset
sig_low = float(np.nanmean(L.values))
all_low = float(np.nanmean(crispr.loc[L.index].values))
print(f"housekeeping low={hk_low:.3f} high={hk_high:.3f}")
print(f"permutation p={perm_p:.4f} z={perm_z:.2f}")
print(f"random-set p={rand_p:.4f} z={rand_z:.2f}")
print(f"low-subset signature={sig_low:.4f} all-genes={all_low:.4f}")

res = {
    "n_nsclc_tpm": int(n),
    "n_low": int(len(L)),
    "n_high": int(len(H)),
    "prolif_genes": list(ENTREZ2SYM.values()),
    "low_ids": low_ids,
    "high_ids": high_ids,
    "low_names": low["CellLineName"].tolist(),
    "high_names": high["CellLineName"].tolist(),
    "n_genes_tested": int(len(use_genes)),
    "obs_mean_diff_low_minus_high": obs_diff,
    "n_depleted_FDR": len(sig_dep),
    "depleted_genes": sig_dep,
    "n_enriched_FDR": len(sig_enr),
    "enriched_genes": sig_enr[:15],
    "housekeeping": {"low": hk_low, "high": hk_high, "n_genes": len(use_house)},
    "permutation": {
        "n": 1000,
        "p_two_sided": perm_p,
        "z": perm_z,
        "null_mean": float(nulls.mean()),
        "null_std": float(nulls.std(ddof=1)),
    },
    "random_set_null": {
        "n": 1000,
        "p_two_sided": rand_p,
        "z": rand_z,
        "null_mean": float(rand_diffs.mean()),
        "null_std": float(rand_diffs.std(ddof=1)),
    },
    "baseline_low_subset": {"signature": sig_low, "all_genes": all_low},
}
json.dump(res, open(f"{BASE}/data/depmap/proliferation_stratified.json", "w"), indent=2)
pd.DataFrame(
    {"gene": use_genes, "mean_diff_low_minus_high": diffs, "mwu_p": pvals, "fdr": fdr}
).to_csv(f"{BASE}/data/depmap/proliferation_stratified_enrichment.csv", index=False)

# ---------- 5. Drug leg: overlap check (min n>=10/group) ----------
sig676 = set(
    pd.read_csv(f"{BASE}/data/signature/signature_676.csv")["Symbol"]
    .str.strip()
    .str.upper()
)


def targets_hit(spec):
    if not isinstance(spec, str):
        return False
    return any(
        t.strip().upper() in sig676
        for t in spec.replace(",", ";").split(";")
        if t.strip()
    )


drug = {}
# GDSC2: per-compound mean AUC across low-quartile lines; MWU sig vs non
gdsc = pd.read_csv(f"{BASE}/data/gdsc/GDSC2_AUC_Matrix.csv", index_col=0)
g_low = [i for i in low_ids if i in gdsc.index]
portal = pd.read_csv(f"{BASE}/data/gdsc/PortalCompounds.csv")
pmap = dict(zip(portal["CompoundID"], portal["GeneSymbolOfTargets"]))
compounds = [c for c in gdsc.columns if c in pmap]
if len(g_low) >= 10:
    m = gdsc.loc[g_low, compounds].mean(axis=0)
    s = [c for c in compounds if targets_hit(pmap[c])]
    o = [c for c in compounds if not targets_hit(pmap[c])]
    _, p = stats.mannwhitneyu(m[s], m[o], alternative="two-sided")
    drug["gdsc2_low_only"] = {
        "n_lines": len(g_low),
        "n_sig": len(s),
        "n_non": len(o),
        "sig_auc": float(m[s].mean()),
        "non_auc": float(m[o].mean()),
        "mwu_p": float(p),
    }
else:
    drug["gdsc2_low_only"] = {"abandoned": True, "n_lines": len(g_low)}
print("GDSC2 low-only:", drug["gdsc2_low_only"])
# PRISM: per-compound (name) mean AUC across low-quartile lines
pr = pd.read_csv(
    f"{BASE}/data/prism/prism_20q2_secondary_dose_response.csv",
    usecols=["depmap_id", "name", "target", "auc"],
)
pr = pr[pr["depmap_id"].isin(low_ids)]
if pr["depmap_id"].nunique() >= 10:
    means, cls = {}, {}
    for nm, grp in pr.groupby("name"):
        means[nm] = grp["auc"].mean()
        cls[nm] = targets_hit(grp["target"].iloc[0])
    ms = pd.Series(means)
    s = [k for k in means if cls[k]]
    o = [k for k in means if not cls[k]]
    _, p = stats.mannwhitneyu(ms[s], ms[o], alternative="two-sided")
    drug["prism_low_only"] = {
        "n_lines": int(pr["depmap_id"].nunique()),
        "n_sig": len(s),
        "n_non": len(o),
        "sig_auc": float(ms[s].mean()),
        "non_auc": float(ms[o].mean()),
        "mwu_p": float(p),
    }
else:
    drug["prism_low_only"] = {
        "abandoned": True,
        "n_lines": int(pr["depmap_id"].nunique()),
    }
print("PRISM low-only:", drug["prism_low_only"])
res["drug_low_only"] = drug
json.dump(res, open(f"{BASE}/data/depmap/proliferation_stratified.json", "w"), indent=2)

# ---------- 6. Figure 8 ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
a = ax[0]
a.hist(nsclc["prolif_score"], bins=20, color="0.75", edgecolor="white")
a.axvspan(
    nsclc["prolif_score"].min(),
    low["prolif_score"].max(),
    color="steelblue",
    alpha=0.25,
    label=f"low quartile (n={len(low)})",
)
a.axvspan(
    high["prolif_score"].min(),
    nsclc["prolif_score"].max(),
    color="firebrick",
    alpha=0.2,
    label=f"high quartile (n={len(high)})",
)
a.set_xlabel("Proliferation score (mean z of log2 TPM+1)")
a.set_ylabel("NSCLC cell lines")
a.legend(fontsize=8)
a.set_title("A. Proliferation stratification")
b = ax[1]
b.scatter(
    diffs, -np.log10(np.where(np.isnan(pvals), 1, pvals)), s=8, color="0.4", alpha=0.6
)
b.axhline(-np.log10(0.05), ls="--", lw=1, color="grey")
b.set_xlabel("Mean Chronos diff (low − high)")
b.set_ylabel("-log10 MWU p")
b.set_title(f"B. Low vs high (n={len(use_genes)} genes)")
c = ax[2]
c.hist(nulls, bins=30, color="0.75", edgecolor="white", label="null")
c.axvline(obs_diff, color="steelblue", lw=2, label=f"observed {obs_diff:.4f}")
c.set_xlabel("Mean diff under label shuffle")
c.set_ylabel("Permutations")
c.legend(fontsize=8)
c.set_title(f"C. Permutation p={perm_p:.4f}")
fig.suptitle("Proliferation-stratified dependency analysis (pre-registered)")
fig.tight_layout()
fig.savefig(f"{BASE}/figures/figure8_proliferation_stratified.png", dpi=200)
fig.savefig(f"{BASE}/figures/figure8_proliferation_stratified.pdf")
print("Figure 8 saved.")
print("DONE")
