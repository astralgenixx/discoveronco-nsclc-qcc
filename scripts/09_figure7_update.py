"""Regenerate Figure 7 (Integrated Evidence Summary) with the proliferation-stratified row.

Original Figure 7 generator (07_generate_figures.py) is superseded; this script
reproduces the same table style and appends the §3.7 result for parity.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROWS = [
    ("DepMap NSCLC enrichment", "0 depleted / 52 enriched", "0.8821", "NO"),
    ("Housekeeping control", "Mean dependency -1.41 (NSCLC)", "\u2014", "PASS"),
    ("Random signature null", "Observed 0.0050", "0.0599", "BORDERLINE"),
    ("GDSC2 drug sensitivity", "AUC 0.885 vs 0.870", "0.5450", "NO"),
    ("PRISM drug sensitivity", "AUC 0.900 vs 0.948", "0.0000", "YES"),
    ("TCGA LUAD specificity", "Ratio 0.96", "1.0000", "NO"),
    ("TCGA LUSC specificity", "Ratio 0.87", "1.0000", "NO"),
    ("TCGA LUAD survival", "Median 1623 vs 1623 days", "0.9203", "NO"),
    ("TCGA LUSC survival", "Median 1714 vs 1714 days", "0.4680", "NO"),
    (
        "Proliferation-stratified dependency",
        "0 depleted / 0 enriched (slow vs fast)",
        "0.5315",
        "NO",
    ),
]

STATUS_COLOR = {
    "NO": "#e74c3c",
    "PASS": "#3498db",
    "BORDERLINE": "#f39c12",
    "YES": "#2ecc71",
}

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("off")
ax.set_title("Integrated Evidence Summary", fontsize=16, fontweight="bold", pad=20)

table_data = [["Analysis", "Result", "P-value", "Significant?"]] + [
    list(r) for r in ROWS
]
table = ax.table(cellText=table_data, loc="center", colWidths=[0.30, 0.32, 0.16, 0.22])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.9)

for (r, c), cell in table.get_celld().items():
    if r == 0:
        cell.set_facecolor("#f2f2f2")
        cell.get_text().set_fontweight("bold")
    elif c == 3 and r > 0:
        cell.set_facecolor(STATUS_COLOR[ROWS[r - 1][3]])

plt.tight_layout()
plt.savefig("figures/figure7_evidence_summary.png", dpi=150, bbox_inches="tight")
plt.savefig("figures/figure7_evidence_summary.pdf", bbox_inches="tight")
print("Figure 7 regenerated with", len(ROWS), "rows.")
