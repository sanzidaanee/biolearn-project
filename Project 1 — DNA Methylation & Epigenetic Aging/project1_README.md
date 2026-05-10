# DNA Methylation Analysis for Epigenetic Aging
### Mapping Human Aging by DNA Methylation Visualizations Using GEO Datasets

This project applies Horvath's epigenetic clock to a large public DNA methylation dataset (GSE40279) to visualize and quantify how methylation patterns change with age. Three complementary visualizations — UMAP, Delta Beta distribution, and Lollipop plot — are used to characterize age-related epigenetic changes at both the global and site-specific level.

---

## Background

Aging is characterized by systematic changes at CpG sites across the genome — a process captured by **epigenetic clocks**. Horvath's Clock, one of the most widely validated of these models, uses ~353 CpG sites trained across multiple tissues to estimate **DNA methylation age (DNAmAge)**. Because it is tissue-agnostic and well-benchmarked, it is ideal for whole blood datasets like GSE40279 and serves as a reference point for comparing newer clocks such as PhenoAge and GrimAge.

---

## Dataset

| Field | Details |
|---|---|
| **GEO Accession** | [GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) |
| **Sample type** | Whole blood |
| **Data type** | Genome-wide DNA methylation (β-values) |
| **Age range** | ~19–101 years |
| **Loaded via** | Biolearn `DataLibrary` |

---

## Analysis Pipeline

### 1. Load Data and Apply Horvath's Clock

```python
from biolearn.data_library import DataLibrary
from biolearn.model_gallery import ModelGallery

# Load GSE40279 from GEO via Biolearn
geo = DataLibrary().get("GSE40279")
gse40279 = geo.load()

# Extract beta values and metadata
beta = gse40279.dnam       # CpGs x Samples
meta = gse40279.metadata   # Samples x metadata

# Apply Horvath v1 epigenetic clock
mg = ModelGallery()
clock = mg.get("Horvathv1")
pred = clock.predict(gse40279)
```

Horvath's clock imputes any missing CpG values and returns a predicted `DNAmAge` for each sample. These predictions are added to the dataset for downstream visualization.

---

### 2. UMAP Visualization of Methylation Profiles Colored by DNAmAge

```python
import scanpy as sc
import anndata as ad

# Build AnnData object (samples x CpGs)
beta_T = beta.T
adata = ad.AnnData(X=beta_T.values, obs=meta.reset_index(drop=True),
                   var=pd.DataFrame(index=beta_T.columns))
adata.obs["DNAmAge"] = pred["Predicted"].values

# Scale, PCA, UMAP
sc.pp.scale(adata)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata)

sc.pl.umap(adata, color=["DNAmAge"], cmap="plasma", size=30,
           title="UMAP: DNA methylation profiles colored by DNAmAge")
```

**What the plot shows:**  
Each point is a sample embedded in 2D using UMAP. Points are coloured by predicted DNAmAge (purple = ~20 years, yellow = ~90 years).

**Key findings:**
- Two distinct clusters emerge, reflecting global shifts in methylation state between age groups
- Within the right-side cluster, DNAmAge shows a clear gradient from younger (purple/blue, bottom) to older (yellow, top)
- The alignment between UMAP structure and DNAmAge confirms that Horvath's clock captures meaningful biological variation, not noise
- Younger samples (20s–30s) form epigenetically distinct states from older samples (60s–90s)

---

### 3. Delta Beta (Δβ) Distribution Plot

```python
# Split samples into young (<40) and old (>70)
young_idx = adata.obs["age"] < 40
old_idx   = adata.obs["age"] > 70

beta_df   = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
mean_young = beta_df.loc[young_idx].mean(axis=0)
mean_old   = beta_df.loc[old_idx].mean(axis=0)
delta_beta = mean_old - mean_young

# Histogram of Δβ
sns.histplot(delta_beta, bins=50, kde=True)
plt.axvline(0, color="red", linestyle="--")
plt.xlabel("Δβ (Old – Young)")
plt.ylabel("Number of CpGs")
```

**What the plot shows:**  
The distribution of per-CpG methylation differences (Δβ = mean_old − mean_young) across all genome-wide sites. Positive Δβ = hypermethylation with age; negative Δβ = hypomethylation with age.

**Key findings:**
- The distribution is roughly bell-shaped and centred near Δβ ≈ 0 — most CpG sites remain stable across the lifespan
- Tails in both directions indicate a subset of sites with strong, systematic age-related changes
- The distribution is approximately symmetric, suggesting a balance between sites gaining and losing methylation
- CpGs in the tails are the primary contributors to epigenetic clocks — they capture aging signal rather than noise
- Hypermethylated sites (Δβ > 0) are typically enriched in CpG islands and promoters, associated with transcriptional silencing with age
- Hypomethylated sites (Δβ < 0) are typically found in intergenic or repetitive regions, linked to genomic instability in older individuals

---

### 4. Lollipop Plot of Top CpGs by Δβ

```python
# Top 15 CpGs by absolute Δβ
delta_df = pd.DataFrame({"CpG": beta_df.columns, "DeltaBeta": delta_beta})
delta_df["abs_db"] = delta_df["DeltaBeta"].abs()
delta_df = delta_df.sort_values("abs_db", ascending=False)
top = delta_df.head(15)

plt.figure(figsize=(12, 6))
plt.vlines(x=range(15), ymin=0, ymax=top["DeltaBeta"], color='skyblue', linewidth=2)
plt.scatter(range(15), top["DeltaBeta"], color='blue', s=50)
plt.xticks(range(15), top["CpG"], rotation=90, fontsize=8)
plt.xlabel("CpG sites (Top Δβ)")
plt.ylabel("Δβ (Old – Young)")
```

**What the plot shows:**  
Each lollipop represents a single CpG site — the stem height shows the magnitude of Δβ, and the dot marks the exact value. Sites are ranked by absolute Δβ.

**Key findings:**
- `cg16867657` (ELOVL2 gene) shows the strongest hypermethylation with age — one of the most reproducible aging markers across datasets
- `cg10501210`, `cg19283806` show strong hypomethylation — consistent with the Δβ histogram
- Both hyper- and hypomethylated CpGs are represented among the top sites
- Many of these top CpGs are directly part of Horvath's 353-site clock or show behaviours consistent with known clock CpGs
- Hypermethylated CpGs in this set occur in promoter/CpG island regions (linked to gene silencing); hypomethylated CpGs occur in enhancers or repetitive elements (linked to genomic instability)

---

## Key Conclusions

- Horvath's clock successfully captures chronological aging in GSE40279 whole blood data, visible as a structured gradient in the UMAP embedding — distinct epigenetic states correspond to distinct age groups
- While the vast majority of CpG sites are stable across the lifespan, a subset undergoes strong, directional methylation changes with age — these are the sites that epigenetic clocks exploit
- The lollipop plot identifies the most age-sensitive individual CpG sites, including `cg16867657` (ELOVL2), a well-established aging biomarker replicated across multiple independent cohorts

---

## Requirements

```bash
pip install biolearn scanpy anndata pandas numpy matplotlib seaborn umap-learn
```

| Library | Purpose |
|---|---|
| `biolearn` | Load GEO datasets, apply Horvath's clock |
| `scanpy` / `anndata` | UMAP pipeline (PCA, neighbors, embedding) |
| `pandas` / `numpy` | Data manipulation, Δβ calculation |
| `matplotlib` / `seaborn` | Visualization (UMAP, histogram, lollipop) |

---

## Files

```
.
+-- code.ipynb     # Full analysis notebook (data loading → UMAP → Δβ → lollipop)
+-- report.md      # Written analysis with plots and biological interpretation
+-- README.md      # This file
```

---

## References

1. Trapp, A., Kerepesi, C., & Gladyshev, V. N. (2021). Profiling epigenetic age in single cells. *Nature Aging*, 1(12), 1189–1201.
2. Garagnani, P., et al. (2012). Methylation of ELOVL2 gene as a new epigenetic marker of age. *Aging Cell*, 11(6), 1132–1134.

---

*Part of the Harvard Aging Initiative × Biomarkers of Aging Consortium Fall 2025 collaboration. See the [main repository README](../README.md) for the full project overview.*
