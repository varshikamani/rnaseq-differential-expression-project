"""
RNA-seq differential expression: breast cancer vs healthy-donor platelets (GEO: GSE68086)

Run:  python rnaseq_analysis.py path/to/GSE68086_TEP_data_matrix.txt
Needs: pandas numpy matplotlib seaborn pydeseq2
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

data_path = sys.argv[1] if len(sys.argv) > 1 else "GSE68086_TEP_data_matrix.txt"
out = "."  # run from the repository root

# 1. Load counts and keep the two groups of interest
counts = pd.read_csv(data_path, sep="\t", index_col=0)
breast_cols = [c for c in counts.columns if "Breast" in c]
healthy_cols = [c for c in counts.columns if "HD" in c]
subset = counts[breast_cols + healthy_cols]
print("Breast:", len(breast_cols), "Healthy:", len(healthy_cols), "Genes:", subset.shape[0])

# 2. Filter low-count genes (>= 10 total counts across samples)
keep = subset.sum(axis=1) >= 10
filtered = subset[keep]
print("Genes after filtering:", filtered.shape[0])

# 3. Metadata (order must match the count matrix columns)
metadata = pd.DataFrame(
    {"group": ["breast"] * len(breast_cols) + ["healthy"] * len(healthy_cols)},
    index=filtered.columns,
)
assert list(filtered.columns) == list(metadata.index)

# 4. DESeq2 model (design = ~ group), healthy donors as reference
dds = DeseqDataSet(counts=filtered.T, metadata=metadata, design="~group", ref_level=("group", "healthy"), quiet=True)
dds.deseq2()

# 5. Variance-stabilizing transformation for PCA / heatmaps
dds.vst(use_design=False)
vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)

# 6. PCA
centered = vst - vst.mean()
u, s, vt = np.linalg.svd(centered.values, full_matrices=False)
var = s**2 / (s**2).sum() * 100
pcs = pd.DataFrame(u[:, :2] * s[:2], index=vst.index, columns=["PC1", "PC2"])
pcs["group"] = metadata["group"]
plt.figure(figsize=(6, 5))
sns.scatterplot(data=pcs, x="PC1", y="PC2", hue="group", palette={"breast": "#d95f02", "healthy": "#1b9e77"}, s=60)
plt.xlabel(f"PC1 ({var[0]:.0f}% variance)")
plt.ylabel(f"PC2 ({var[1]:.0f}% variance)")
plt.title("PCA of VST-transformed counts")
plt.tight_layout()
plt.savefig(f"{out}/plots/pca_plot.png", dpi=200)
plt.close()

# 7. Sample-to-sample distance heatmap
from scipy.spatial.distance import pdist, squareform
dist = pd.DataFrame(squareform(pdist(vst.values)), index=vst.index, columns=vst.index)
lut = {"breast": "#d95f02", "healthy": "#1b9e77"}
colors = metadata["group"].map(lut)
g = sns.clustermap(dist, cmap="viridis_r", row_colors=colors, col_colors=colors, xticklabels=False, yticklabels=False, figsize=(7, 7))
g.fig.suptitle("Sample-to-sample distances", y=1.02)
g.savefig(f"{out}/plots/sample_distance_heatmap.png", dpi=200)
plt.close("all")

# 8. Differential expression test: breast vs healthy
stat = DeseqStats(dds, contrast=["group", "breast", "healthy"], quiet=True)
stat.summary()
res = stat.results_df.copy()
res = res.sort_values("padj")
res.to_csv(f"{out}/results/deseq2_all_genes.csv")
sig = res[(res["padj"] < 0.05)]
sig.to_csv(f"{out}/results/significant_genes_padj_0.05.csv")
print("Significant genes (padj < 0.05):", len(sig),
      "| up in breast:", int((sig["log2FoldChange"] > 0).sum()),
      "| down in breast:", int((sig["log2FoldChange"] < 0).sum()))

# 9. Volcano plot
v = res.dropna(subset=["padj"]).copy()
v["neglog10"] = -np.log10(v["padj"].clip(lower=1e-300))
v["status"] = "not significant"
v.loc[(v.padj < 0.05) & (v.log2FoldChange > 1), "status"] = "up in breast (padj<0.05, |LFC|>1)"
v.loc[(v.padj < 0.05) & (v.log2FoldChange < -1), "status"] = "down in breast (padj<0.05, |LFC|>1)"
plt.figure(figsize=(7, 5.5))
sns.scatterplot(data=v, x="log2FoldChange", y="neglog10", hue="status",
                palette={"not significant": "#bbbbbb", "up in breast (padj<0.05, |LFC|>1)": "#d95f02", "down in breast (padj<0.05, |LFC|>1)": "#1f78b4"},
                s=8, linewidth=0)
plt.axhline(-np.log10(0.05), ls="--", c="k", lw=0.6)
plt.axvline(1, ls="--", c="k", lw=0.6); plt.axvline(-1, ls="--", c="k", lw=0.6)
plt.xlabel("log2 fold change (breast vs healthy)")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot")
plt.legend(fontsize=7, loc="upper right")
plt.tight_layout()
plt.savefig(f"{out}/plots/volcano_plot.png", dpi=200)
plt.close()

# 10. Heatmap of the top 30 genes (by adjusted p-value)
top = res.dropna(subset=["padj"]).head(30).index
mat = vst[top]
mat = (mat - mat.mean()) / mat.std()
order = metadata.sort_values("group").index
g = sns.clustermap(mat.loc[order].T, cmap="vlag", col_cluster=False, row_cluster=True,
                   col_colors=metadata.loc[order, "group"].map(lut), xticklabels=False, figsize=(9, 8))
g.fig.suptitle("Top 30 differentially expressed genes (z-scored VST)", y=1.02)
g.savefig(f"{out}/plots/top30_genes_heatmap.png", dpi=200)
plt.close("all")
print("Done.")
