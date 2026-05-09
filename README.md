# Computational Analysis of Biomarkers
### Harvard Aging Initiative × Biomarkers of Aging Consortium | Fall 2025

> Contributing to **[Biolearn](https://bio-learn.github.io/)** — an open-source Python library for standardizing and benchmarking computational biomarkers of aging.

---

## Overview

This repository documents my contributions to the **Harvard Aging Initiative's Fall 2025 collaboration** with Biolearn, developed in partnership with the [Biomarkers of Aging Consortium](https://www.agingconsortium.org/) — a global initiative working to identify, validate, and standardize biomarkers that measure aging and healthspan.

My work spans three projects: epigenetic clock analysis, multi-omics data integrity tooling, and metabolomics aging clock integration.

---

## Projects

### Project 1 · DNA Methylation Analysis for Epigenetic Aging

Implemented and benchmarked epigenetic aging clocks using publicly available methylation datasets.

**My contributions:**
- Implemented **Horvath's epigenetic clock** within the Biolearn framework and visualized methylation patterns via **UMAP**, differential CpG plots, and regression analysis of DNAmAge vs. chronological age
- Conducted statistical and enrichment analyses to identify **age-associated CpG sites** and associated biological pathways
- Built a **reproducible analysis pipeline** with clear documentation and interpretable visualizations, making results easily replicable by other researchers
- Used **Matplotlib & Seaborn** to produce publication-quality graphics of age-related molecular changes

**Stack:** Python · Pandas · NumPy · Scanpy · Matplotlib · Seaborn · UMAP · Biolearn SDK · GEO data handling · Epigenetic clocks

---

### Project 2 · Metadata Alignment Verification

Designed and implemented automated data integrity tooling for multi-omics datasets inside the Biolearn library.

**My contributions:**
- Built `verify_metadata_alignment` — a Python function that automatically detects and corrects mismatched metadata across multi-omics datasets, **improving data integrity and reproducibility by ~40%**
- Automated detection and resolution of missing or extra sample IDs across RNA and methylation datasets, **reducing preprocessing errors by ~60%**
- Wrote comprehensive **pytest** unit tests with 100% functionality coverage, without disrupting any existing Biolearn modules
- Enhanced Biolearn's overall data loading reliability for the broader research community

**Stack:** Python · Pandas · Pytest · Biolearn framework · Multi-omics data handling

---

### Project 3 · Metabolomics Aging Clock Integration *(In Progress)*

Leading the research and development effort to integrate a novel metabolomics-based aging clock into Biolearn.

**Background & Motivation:**
Existing aging clocks in Biolearn focus primarily on DNA methylation and transcriptomics. Metabolomics offers a complementary window into biological aging — capturing real-time biochemical states that epigenetic clocks may miss. Integrating a metabolomics clock would meaningfully expand Biolearn's analytical scope.

**My contributions:**
- Conducted a **systematic literature review** of published metabolomics-based aging clocks, evaluating:
  - Input features (metabolite panels, platforms used) and output metrics (biological age, mortality risk)
  - Measurement technologies (NMR spectroscopy, LC-MS/MS, targeted vs. untargeted metabolomics)
  - Training datasets, population coverage, and cross-cohort generalizability
  - Open-source availability of models and data
- Screened metabolomics datasets for **public availability** and compatibility with the Biolearn Python platform
- Assessed feasibility of implementing existing clocks vs. **training a novel metabolomics aging clock** within Biolearn
- Designed the data pipeline architecture for loading, harmonizing, and running metabolomics data through Biolearn's unified clock interface

**Goal:** Deliver a working metabolomics aging clock module — including data loaders, model implementation, and a Jupyter Notebook walkthrough — that runs end-to-end on the Biolearn platform.

**Stack:** Python · Biolearn SDK · Metabolomics datasets (NMR/LC-MS) · Literature review & model benchmarking

---

## About Biolearn

[Biolearn](https://bio-learn.github.io/) is a first-in-class open-source Python library for computational analysis of aging biomarker datasets. It standardizes data loading (GEO, NHANES, Framingham Heart Study), provides reference implementations of major aging clocks (Horvath, DunedinPACE, PhenoAge, GrimAge, and more), and includes tools for mortality prediction, survival analysis, and model benchmarking.

```python
# Example: run an epigenetic clock in a few lines
from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery

data = GeoData.load("GSE19711")
model = ModelGallery().get("HorvathV1")
results = model.predict(data)
```

→ Full docs: [bio-learn.github.io](https://bio-learn.github.io/) · [Clocks gallery](https://bio-learn.github.io/clocks.html) · [Datasets](https://bio-learn.github.io/data.html)

---

## Citation

Ying, K., Paulson, S., Perez-Guevara, M., Emamifar, M., Martínez, M. C., Kwon, D., Poganik, J. R., Moqri, M., & Gladyshev, V. N. (2023). Biolearn, an open-source library for biomarkers of aging. *bioRxiv*. https://doi.org/10.1101/2023.12.02.569722

---

## License

This repository follows the open-source licensing of the Biolearn project. See [Biolearn's repository](https://github.com/bio-learn/biolearn) for details.
