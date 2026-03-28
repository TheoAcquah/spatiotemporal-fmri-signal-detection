# spatiotemporal-fmri-signal-detection
Unlike traditional GLM approaches, this framework captures:  Temporal evolution of brain activity Distributed activation patterns Improved detection under autocorrelation


# Spatiotemporal fMRI Signal Detection Framework

This repository provides code and sample outputs for the paper:

"Advancing Temporal Dynamics in Spatial Random Field Theory: A Framework for fMRI Signal Detection"

## 🔍 Overview
We propose a spatiotemporal signal detection framework based on:
- Time-adaptive Gaussian random fields
- FAR(1) temporal modeling
- Global test statistics (Xmax, Ymax)

## 📊 Dataset
Data used in this study are publicly available from:
OpenNeuro ds000114 (Gorgolewski et al., 2013)

Link: https://openneuro.org/datasets/ds000114

## 📁 Repository Structure
- `code/` – analysis pipelines (GLM, Xmax, Ymax)
- `data/` – sample processed data
- `outputs/` – figures and tables from the paper
- `docs/` – workflow diagram

## ▶️ How to Run
```bash
pip install -r requirements.txt
python code/xmax/run_xmax.py




spatiotemporal-fmri-detection/
│
├── README.md
├── requirements.txt
├── LICENSE
│
├── data/
│   ├── raw/                # (optional: DO NOT upload full dataset)
│   ├── processed/
│   │   └── sample_subject/ # small demo subset
│   └── description.md
│
├── code/
│   ├── preprocessing/
│   ├── glm/
│   ├── xmax/
│   ├── ymax/
│   └── visualization/
│
├── outputs/
│   ├── figures/
│   └── tables/
│
└── docs/
    └── workflow.png
