# 🧬 Treg Foxp3 Degradation — High-Throughput Genomics Analysis

This repository contains a collection of Jupyter notebooks and associated datasets for the comprehensive analysis of **regulatory T cells (Tregs)** undergoing **acute Foxp3 degradation**, across diverse biological contexts. Leveraging both **single-cell** and **bulk RNA-sequencing** technologies, these analyses uncover dynamic transcriptional changes following the targeted loss of *Foxp3*, a master regulator of Treg identity and function.

---

## 📁 Repository Structure

The repository is organized into **single-cell** and **bulk transcriptomics** components, each addressing a distinct experimental context.

### 🔬 Single-Cell RNA-seq Analyses

These notebooks explore the cellular heterogeneity and transcriptional dynamics of Tregs following Foxp3 degradation:

- **Foxp3 Degradation in Adult Secondary Lymphoid Organs (SLO)**  
  - `scrna-filter-d0.ipynb`  
  - `scrna-postfilter.ipynb`

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

## 💡 Purpose

The goal of this project is to dissect the downstream consequences of **acute Foxp3 loss** in Tregs, with a focus on:

- Loss of Treg identity  
- Emergence of effector-like transcriptional programs  
- Tissue- and context-specific regulatory rewiring  
- Implications for tumor immunity and autoimmunity

---

## 🛠️ Requirements

To run the notebooks, ensure you have a Python environment with the following key packages:

- `scanpy`
- `anndata`
- `pandas`, `numpy`, `matplotlib`, `seaborn`
- `scikit-learn`

We recommend setting up a virtual environment using `venv` or `conda` to manage dependencies.

---

## 📣 Citation & Acknowledgments

If you find this repository useful in your research, please cite the forthcoming publication (preprint link TBD) and acknowledge the authors.

For questions, feedback, or collaboration opportunities, feel free to reach out.

---
