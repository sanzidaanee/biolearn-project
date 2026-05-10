
# Mapping Human Aging by DNA Methylation Visualizations Using GEO datasets

## Introduction

Aging can be characterized by several canonical hallmarks, including epigenetic alterations at CG dinucleotides (CpG sites) [1]. Machine learning models based on DNA methylation (DNAm) levels, termed ‘epigenetic clocks’ have revolutionized the aging field [1]. These clocks are the estimators of chronological age, and predict various measures of biological aging and disease risk with remarkable precision indicating the strong conservation of epigenetic aging patterns across species. 

Horvath’s Clock is one of the widely validated epigenetic clocks, trained across multiple tissues. It uses ~353 CpG sites, many of which are conserved across tissues, making it ideal for a general dataset like GSE40279. It provides a benchmark for newer clocks (PhenoAge, GrimAge, etc.) that can be compared against Horvath’s original model.

## Exploring Age-Related Epigenetic Changes by Analyzing GEO Datasets

The GEO dataset (GSE40279), which contains genome-wide DNA methylation (β-values) from whole blood samples across different ages. By using the Biolearn platform, Horvath DNA methylation clock (Horvathv1), an epigenetic clock model based on methylation levels at a defined set of CpG sites is applied. This model predicted DNA methylation age (DNAmAge) for each sample by imputing missing values and computing an epigenetic age estimate. These predicted ages are added into the dataset and visualized using UMAP (Uniform Manifold Approximation and Projection for Dimension Reduction), which allows one to notice how methylation-based biological aging relates to the overall structure. 

To further explore site-specific age-related changes,a Δβ plot is generated , which shows the distribution of methylation differences (Δβ) between younger and older individuals. This highlights whether CpG sites tend to gain or lose methylation with age. Finally, a lollipop plot  is created of the top CpGs ranked by Δβ, which provides a more detailed, gene-level perspective of the sites most strongly associated with aging.


The results of this analysis can be used to quantify biological age. The DNAm clocks of Horvath estimate the "biological age" of a tissue based on its methylation signatures. Additionally, we can use the Horvath clock to determine whether methylation patterns in GEO datasets correspond to biological aging trajectories to assess aging-related epigenetic changes. Also, the predicted and actual ages are compared to identify samples or groups with accelerated or reduced epigenetic aging.


## Results

### Visualization Strategies for Age-Associated Methylation: UMAP, Δβ, and Lollipop Plots

To analyze DNA methylation patterns in aging, we performed UMAP, Lollipop, and Δβ Plots by using GEO datasets. 

### UMAP Visualization of DNA Methylation Profiles Colored by Age

#### The plot shows

- Data: The plot is a 2D UMAP embedding of genome-wide DNA methylation profiles from GSE40279.

- Coloring: Each point is a sample, colored by predicted DNA methylation age (DNAmAge) from Horvath’s clock. The color scale runs from ~20 (purple/blue, younger) to ~90 (yellow, older).

- Axes: UMAP1 and UMAP2 are non-linear dimensions summarizing high-dimensional methylation variation. Clusters can reflect shared epigenetic patterns.



<img width="598" height="426" alt="1" src="https://github.com/user-attachments/assets/50375c67-481b-4785-9bc8-0bce130c4fc2" />




### Main observations

- It shows two clear clusters, one in left (UMAP1 negative side) and another one on the right (UMAP1 positive side). Within the right-side cluster, DNAmAge colors show a progression from purple/blue (younger) at the bottom to yellow (older) at the top.


- The separation and gradient suggest that methylation profiles systematically change with age.

### Biological Insights

- Horvath’s clock successfully reflects chronological age in GSE40279, visible as a smooth gradient in the UMAP embedding.

- Younger samples (20s–30s) cluster differently from older samples (60s–90s). This suggests that global methylation changes accumulate with age and create distinct “epigenetic states.”

- The fact that UMAP structure aligns with DNAmAge supports that the clock captures meaningful biological variation, not just noise.


## Delta Beta (Δβ) Plot Highlighting Epigenetic Changes with Aging


### The plot shows

- X-axis (Δβ = Old − Young): The difference in DNA methylation (β-values) at each CpG site between older individuals and younger individuals.

   - Positive Δβ: Higher methylation in the older group (hypermethylation with age).

  - Negative Δβ: Lower methylation in the older group (hypomethylation with age).


- Y-axis: The number of CpGs showing that Δβ value.


- Distribution: Roughly bell-shaped, centered near 0, with most CpGs showing only small changes between old and young groups.


<img width="854" height="553" alt="2" src="https://github.com/user-attachments/assets/470cc497-ccd4-427d-98f4-00835ec80326" />




### Main observations

- Most CpG sites do not show large methylation differences between young and old individuals in this graph.

- Some CpG sites show strong age-related methylation changes — both hypermethylation (positive tail) and hypomethylation (negative tail).

- The distribution appears in this graph relatively symmetric, suggesting a balance between sites gaining and losing methylation with age.


### Biological interpretation

- The majority of CpGs remain stable across the lifespan, which is why the peak is so high near Δβ ≈ 0.

- The CpGs in the tails are the ones that systematically change with age. These sites are likely contributors to epigenetic clocks like Horvath’s clock, since they capture aging signals rather than noise.

- Hypermethylated sites (Δβ > 0) which are often enriched in CpG islands and promoters, leading to repression of gene expression with age while in hypomethylated sites (Δβ < 0) in intergenic or repetitive regions, leading to increased genomic instability in older individuals.


## Mapping Age-Related DNA Methylation Changes via Lollipop Plot

### Plot shows
- X-axis: Individual CpG sites (Illumina probe IDs, e.g., cg16867657, cg10501210).

- Y-axis (Δβ = Old − Young): The methylation difference at each CpG between old and young individuals.

  - Positive values: CpGs gain methylation (hypermethylation) with age.
  - Negative values: CpGs lose methylation (hypomethylation) with age.

- Lollipop stems: Each stem represents the magnitude of methylation difference for that CpG site.


<img width="1135" height="563" alt="3" src="https://github.com/user-attachments/assets/c0abf80a-a6e0-4e5d-bbaf-44a6f6a3ae86" />



### Main observations

- Large Δβ values: Some CpGs (e.g., cg16867657, cg22454769) show strong hypermethylation in older individuals. Others (e.g., cg10501210, cg19283806) show strong hypomethylation.

- Top CpGs represent the most age-sensitive sites in GSE40279, where methylation differences are strongest between age groups.
 Both hypermethylation and hypomethylation are seen, consistent with the Δβ distribution histogram shown earlier.

### Biological interpretation
- CpG sites like cg16867657 is a well-known aging marker in the ELOVL2 gene [2] and are among the most reproducible sites for predicting age across datasets.
- Hypermethylated CpGs often occur in promoters/ CpG islands, linked to transcriptional silencing with age.
- Hypomethylated CpGs often occur in enhancers or repetitive elements, linked to loss of genomic stability.
- Horvath’s model selects CpGs that show consistent, predictable methylation shifts with age. Many CpGs in this plot (like cg16867657) are either directly part of Horvath’s 353 clock CpGs or behave similarly.


## Conclusion

- Horvath's epigenetic clock estimates that methylation profiles separate into distinct clusters and gradients corresponding to age. A gradient of DNA methylation patterns is observed between younger and older individuals, demonstrating how epigenetic patterns evolve with age.

- Δβ distribution shows that while most CpG sites remain stable between young and old groups, a subset undergoes significant hyper- or hypomethylation with age. These age-sensitive CpGs are responsible for Horvath's clock and explain how methylation profiles reflect biological aging.

- The lollipop plot identifies the CpG sites with the largest methylation differences between young and old individuals. These epigenetic signals of aging include hypermethylated CpGs (like cg16867657) and hypomethylated CpGs.



## References

1. Trapp, A., Kerepesi, C., & Gladyshev, V. N. (2021). Profiling epigenetic age in single cells. Nature Aging, 1(12), 1189-1201.

2. Garagnani, P., Bacalini, M. G., Pirazzini, C., Gori, D., Giuliani, C., Mari, D., ... & Franceschi, C. (2012). Methylation of ELOVL 2 gene as a new epigenetic marker of age. Aging cell, 11(6), 1132-1134.



