# statistics_sandbox

Sandbox for exploratory statistics, meta-analysis workflows, and teaching notebooks across multiple tools (Jupyter, Quarto/R, Stata).

## Project Inventory

| Path | Type | Primary entry points |
|---|---|---|
| `/` (repo root) | Jupyter + Quarto | `scheme simulation.ipynb`, `visualizing poor mans p value.ipynb`, `Representative Survey Data.qmd` |
| `ASA Baseline Use MA/` | Stata meta-analysis | `ASA BL MA.do`, `ASA BL MA.xls` |
| `Conway Thorax supplement and code/` | R Markdown analysis | `TcCO2 meta-analysis.Rmd`, `data.Rdata` |
| `Heterotopic Ossification/` | Stata meta-analysis | `HO MA.do`, `HO MA.xls` |
| `ICU Glucose/` | Stata meta-analysis | `GluMA.do`, `O2MA.do` |
| `Medicaid Provider Data/` | Quarto + Python EDA | `medicaid-provider-spending-eda.qmd`, `medicaid-provider-spending.parquet` |
| `NLP Workshop/` | Jupyter workshop notebook | `NLP Workshop.ipynb` |
| `NRH Vent Proj/` | Quarto project | `code/NRH Vent.qmd`, `code/sci_poster.qmd` |
| `O2 Target Meta-Analysis/` | Quarto/R project | `O2 Target Meta-Analysis.qmd`, `O2 Target Meta-Analysis.Rproj` |
| `Obesity NIV Reintubation/` | Stata meta-analysis | `ObesNIV MA.do`, `ObesNIV MA.xls` |
| `PFT sim and figs/` | Stata simulation | `PFT sim.do` |
| `Pepe Diagnostic Testing Book/` | Stata diagnostics toolkit | `Code/*.ado`, `Datasets/*.dta` |
| `RCT Duplicate/` | Stata + R reproducibility analysis | `RCT Duplicate Analysis.do`, `DUPLICATE-advances_statistics_code_reprod.R` |
| `Software Carpentry Workshop/` | Quarto training materials | `_lessons_live/*.qmd`, `Quarto Docs/Test Doc.qmd` |
| `Steroids PNa/` | Stata meta-analysis | `Steroids PNa.do`, `Steroids PNa MA.xls` |
| `TcCO2/` | Stata analysis | `TcCO2.do`, `TcCo2 data table.xlsx` |
| `dication_workflow/` | Jupyter dictation pipeline | `dictation_notebook.ipynb` |
| `diagnosis_llms/` | **New transferred exploratory notebooks** | `notebooks/display_reasoning.ipynb`, `notebooks/multi_class_cont.ipynb`, `notebooks/test_notebook.ipynb` |

## New Notebook Transfer (2026-03-02)

Transferred bundle was unpacked to:

- `diagnosis_llms/notebooks/`
- `diagnosis_llms/src/diagnosis_llms/`
- `diagnosis_llms/env/`
- `diagnosis_llms/README_TRANSFER.md`

`diagnosis_llms/pyproject.toml` was added as a local project marker so notebook path bootstrap logic resolves correctly.

## Execution Guide

### 1) Jupyter notebooks (general)

From repo root:

```bash
jupyter lab
```

Then open any `.ipynb` notebook. For projects with custom dependencies, use a dedicated virtual environment first.

### 2) `diagnosis_llms` transferred notebooks

From repo root:

```bash
cd "diagnosis_llms"
python3 -m venv .venv
source .venv/bin/activate
pip install -r env/requirements.notebook-transfer.txt
python -m ipykernel install --user --name notebook-transfer --display-name "Python (notebook-transfer)"
jupyter lab notebooks
```

Notes:

- `test_notebook.ipynb` is a stabilized default-execution version for this repo (heavy/external cells are skipped).
- Full exploratory version is preserved as `diagnosis_llms/notebooks/test_notebook_full_original.ipynb`.
- Running the full-original notebook may require API keys (`OPENAI_API_KEY`, `ANTHROPIC_API_KEY`), shell tools, network downloads, and substantial compute.

### 3) Quarto projects (`.qmd`)

Typical workflow:

```bash
quarto render "O2 Target Meta-Analysis/O2 Target Meta-Analysis.qmd"
quarto render "Medicaid Provider Data/medicaid-provider-spending-eda.qmd"
```

Prerequisites: Quarto CLI, plus required R/Python dependencies per project.

### 4) Stata projects (`.do`)

Open the corresponding `.do` file in Stata and run it from the project directory so relative paths resolve.

Examples:

- `ASA Baseline Use MA/ASA BL MA.do`
- `TcCO2/TcCO2.do`
- `RCT Duplicate/RCT Duplicate Analysis.do`

### 5) R Markdown scripts (`.Rmd`)

Example:

```bash
Rscript -e "rmarkdown::render('Conway Thorax supplement and code/TcCO2 meta-analysis.Rmd')"
```

Install needed R packages before rendering.

## Notes

- This repository is multi-project and does not currently define one global environment for every folder.
- Some project folders include generated outputs and local environments (for example `Medicaid Provider Data/.venv`).
- Paths contain spaces; quote paths in shell commands.

DETAIN-IVF simulation notebook launch (Binder):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/reblocke/statistics_sandbox/HEAD?urlpath=%2Fdoc%2Ftree%2Fscheme+simulation.ipynb)
