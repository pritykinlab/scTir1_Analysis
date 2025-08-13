# 🧬 Treg Foxp3 Degradation — High-Throughput Genomics Analysis

This repository contains a collection of Jupyter notebooks and associated datasets for the comprehensive analysis of **regulatory T cells (Tregs)** undergoing **acute Foxp3 degradation**, across diverse biological contexts for the following publication: (link TBD). Leveraging both **single-cell** and **bulk RNA-sequencing** technologies, these analyses uncover dynamic transcriptional changes following the targeted loss of *Foxp3*


---

## 📁 Repository Structure

The repository is organized into **single-cell**, **bulk transcriptomics**, and **shared utility code**, each addressing a distinct role in the analysis workflow.

### 🔬 Single-Cell RNA-seq Analyses

These notebooks explore the effects on transcription in diverse Treg contexts following Foxp3 degradation:

- **Foxp3 Degradation in Adult Secondary Lymphoid Organs (SLO)**  
  - `scrna-filter-d0.ipynb`  (Initial filtering of D0 cells)
  - `scrna-postfilter.ipynb` (Complete analysis)

- **Foxp3 Degradation Across Multiple Organs**  
  - `scrna-organ.ipynb`

- **Foxp3 Degradation in the Tumor Microenvironment**  
  - `scrna-tumor.ipynb`

- **Foxp3 Degradation in Adoptive Transfer Settings**  
  - `scrna-transfer.ipynb`

---

### 🧪 Bulk RNA-seq Analyses

Bulk RNA-seq analyses include both newly generated and reanalyzed datasets:

- `bulk_rna_data/`  
  Contains gene expression matrices, differential expression analyses, and metadata relevant to Foxp3 degradation experiments.

---

### 🧰 Reusable Functions and Utilities

- `code/`  
  A collection of shared functions and modules used across the notebooks, including:
  - Data loading and preprocessing utilities  
  - Plotting functions  
  - Statistical tests and helper methods

These utilities ensure modularity and reproducibility across the analysis workflows.

---

## 🛠️ Requirements

To run the notebooks, ensure you have a Python environment with the following key packages:

- `scanpy`
- `scanpy`
- `pandas`, `numpy`, `matplotlib`, `seaborn`
- `scikit-learn`

We recommend setting up a virtual environment using `venv` or `conda` to manage dependencies.

---

## 📣 Questions & Feedback

For questions, feedback, or other concerns, feel free to reach out via e-mail or by posting a GitHub issue.

---
