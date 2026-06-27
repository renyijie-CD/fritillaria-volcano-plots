# Volcano plot reproducibility files

This repository provides processed volcano-plot data tables and the Python script used for volcano-plot generation and feature-screening visualization in the revised manuscript.

## Contents

- `volcano_plot.py`  
  Python script used to generate volcano plots.

- `volcano_plot_data_tables/`  
  Processed volcano-plot data tables for the four pairwise comparisons:
  - `volcano_stats_FW-B_vs_FW-F.xlsx`
  - `volcano_stats_F-B_vs_F-F.xlsx`
  - `volcano_stats_FW-B_vs_F-B.xlsx`
  - `volcano_stats_FW-F_vs_F-F.xlsx`

- `preview_plots_png/`  
  Preview PNG figures generated using `volcano_plot.py`.

- `run_examples.txt`  
  Example command lines for reproducing the volcano plots.

## Data description

The processed volcano-plot data tables contain feature-level statistics used for volcano-plot generation, including feature ID, m/z, retention time, log2FC, fold change, and p value.

These tables are intended to reproduce the volcano plots and feature-screening visualization from processed data. They are not raw LC-MS vendor files and are not full sample-by-feature peak-intensity matrices.

## Screening threshold

The volcano plots use the following threshold:

- `|log2FC| >= 1`
- `p < 0.05`

Features meeting these criteria are shown using a continuous color scale according to log2FC. Non-selected features are shown in grey.

## Notes on compound names

Compound names in the volcano-plot data tables are software-generated preliminary annotations. Final representative compounds reported in Tables S8-S11 of the Supplementary Material were manually curated using accurate mass, retention time, MS/MS evidence, reference standards when available, VIP values, peak quality, and the in-house FCB-related compound library.

## Software requirements

The script requires Python 3 and the following packages:

```bash
pip install numpy pandas matplotlib openpyxl
```

On Windows, if multiple Python versions are installed, use:

```bash
python -m pip install numpy pandas matplotlib openpyxl
```

or

```bash
py -m pip install numpy pandas matplotlib openpyxl
```
