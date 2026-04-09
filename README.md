# Resubmission Analysis Workflow

This repository contains the notebooks, scripts, configuration files, and exported analysis outputs used to retrain the revision models and regenerate the manuscript figures, statistics tables, and source-data tables.

The workflow is organized around two entry styles:

- interactive notebook reruns for step-by-step analysis
- shell and Python scripts for non-interactive batch execution

## Repository Layout

```text
resubmission/
├── README.md
├── environment-parnet_clean.yml
├── environment-tfmodisco.yml
├── requirements.txt
├── configs/
├── data/
├── notebooks/
├── results/
└── scripts/
```

## Environments

Two environments are used throughout the workflow:

- `parnet_clean` for metadata assembly, model training, UMAP export, CAM export, figure generation, and most notebooks
- `tfmodisco` only for motif discovery and motif report generation

### Main Analysis Environment

Create the main environment:

```bash
conda env create -f environment-parnet_clean.yml
conda activate resubmission-parnet-clean
```

Install the Python packages exported from the original `parnet_clean` environment:

```bash
pip install -r requirements.txt
```

Install `ir_toolkit` from GitHub:

```bash
pip install git+https://github.com/melonheader/ir_toolkit.git
```

Notes:

- `requirements.txt` is a package snapshot from the environment used for the resubmission analyses.
- A machine-local dependency, `rinalmo @ file:///...`, is not included in `requirements.txt` because it is not portable.
- GPU-specific packages are kept in `requirements.txt` for fidelity with the original runs.

### tf-modisco Environment

Create the dedicated motif-analysis environment separately:

```bash
conda env create -f environment-tfmodisco.yml
conda activate resubmission-tfmodisco
```

This environment provides the `modisco-lite` Python package together with the MEME suite tools used by `tomtom`.

## Running The Notebooks

Launch Jupyter from the main analysis environment:

```bash
conda activate resubmission-parnet-clean
jupyter lab
```

Use the `resubmission-parnet-clean` kernel for all analysis notebooks. The tf-modisco step is executed through the shell wrapper and uses the dedicated `tfmodisco` environment internally.

Notebook guide:

- `notebooks/1.reassemble_metadata.ipynb`
  Rebuilds `data/metadata_selected.csv` from the assembled source tables.
- `notebooks/2.hl_revision_workflow.ipynb`
  Main workflow notebook for run validation, retraining, UMAP export, CAM export, tf-modisco execution, and per-run plotting.
- `notebooks/3.remake_fig3.ipynb`
  Regenerates Figure 3 panels and associated source tables.
- `notebooks/4.remake_sfig4.ipynb`
  Regenerates Supplementary Figure 4 panels and associated source tables.
- `notebooks/5.remake_sfig5.ipynb`
  Regenerates Supplementary Figure 5 analyses and tables.
- `notebooks/6.detained_introns_analysis.ipynb`
  Runs the detained-intron follow-up analyses.
- `notebooks/7.intron_motif_enrichment.ipynb`
  Runs motif-enrichment analyses for intron subsets.

## Main Retraining Workflow

The half-life revision workflow is defined by:

- `configs/hl_revision_runs.json`
- `scripts/hl_revision_pipeline.py`
- `scripts/run_hl_revision_workflow.sh`
- `scripts/run_all_hl_revision_workflow.sh`
- `scripts/run_modisco_report.sh`
- `notebooks/2.hl_revision_workflow.ipynb`

### Inspect Configured Runs

```bash
python \
  resubmission/scripts/hl_revision_pipeline.py list-runs
```

### Train One Configured Run

```bash
python \
  resubmission/scripts/hl_revision_pipeline.py train \
  --run hl_revised_50percgap
```

Train a single fold:

```bash
python \
  resubmission/scripts/hl_revision_pipeline.py train \
  --run hl_revised_50percgap \
  --folds 1
```

### Export UMAP Embeddings

```bash
python \
  resubmission/scripts/hl_revision_pipeline.py export-umap \
  --run hl_revised_50percgap
```

### Export CAM / tf-modisco Inputs

```bash
python \
  resubmission/scripts/hl_revision_pipeline.py export-modisco-inputs \
  --run hl_revised_50percgap \
  --cam-modes final_logit_linearized,branch_signed
```

### Run tf-modisco

```bash
resubmission/scripts/run_modisco_report.sh \
  --run hl_revised_50percgap \
  --cam-mode final_logit_linearized \
  --motif-db /path/to/pwms_all_motifs_ids.meme
```

Run the branch-signed interpretation in parallel naming space:

```bash
resubmission/scripts/run_modisco_report.sh \
  --run hl_revised_50percgap \
  --cam-mode branch_signed \
  --motif-db /path/to/pwms_all_motifs_ids.meme
```

### Export Per-Run Plots

```bash
python \
  resubmission/scripts/export_hl_revision_plots.py \
  --run hl_revised_50percgap \
  --cam-mode final_logit_linearized
```

## Shell Wrappers

For a single run, the main wrapper is:

```bash
resubmission/scripts/run_hl_revision_workflow.sh \
  --run hl_revised_50percgap \
  --train \
  --umap \
  --modisco-inputs \
  --modisco-report \
  --plots \
  --cam-modes final_logit_linearized,branch_signed \
  --motif-db /path/to/pwms_all_motifs_ids.meme
```

For all configured runs:

```bash
resubmission/scripts/run_all_hl_revision_workflow.sh \
  --train \
  --cam-modes final_logit_linearized,branch_signed
```

With no stage flags, `run_all_hl_revision_workflow.sh` regenerates downstream outputs by default:

- UMAP embeddings
- tf-modisco inputs
- tf-modisco reports
- plot exports

## Regenerating Figures, Statistics, And Source Tables

The figure-remaking notebooks and plot-export scripts write outputs into:

- `results/plots/` for manuscript figure panels, source CSVs, and statistics tables
- `results/models/*/plot_exports/` for per-run exported plotting tables and figures

To consolidate CSV statistics tables into a single workbook:

```bash
python resubmission/scripts/merge_csv_to_xlsx.py \
  --input-dir resubmission/results/plots \
  --output-xlsx resubmission/results/plots/statistics_tables_unified.xlsx
```

## Raw Data Processing Helpers

The repository also includes the shell scripts used to download and process the PRJNA608890 raw sequencing data under:

- `scripts/raw_data/PRJNA608890/`

Included helpers:

- `scripts/raw_data/PRJNA608890/download_sra.sh`
  Downloads SRA accessions and converts them to compressed FASTQ files.
- `scripts/raw_data/PRJNA608890/trim_reads.sh`
  Runs read trimming for the WGBS sample set.
- `scripts/raw_data/PRJNA608890/process_WGBS.sh`
  Builds the filtered Bismark genome index, aligns WGBS reads, deduplicates alignments, and extracts methylation calls.
- `scripts/raw_data/PRJNA608890/process_ssDRIP.sh`
  Trims, aligns, and deduplicates the ssDRIP-seq reads.
- `scripts/raw_data/PRJNA608890/process_ssDRIP_macs3.sh`
  Calls strand-specific ssDRIP peaks with `macs3`.

These scripts were copied from the original raw-data processing workspace and are kept as the exact command-line helpers used for download and preprocessing.

## Recommended Order

For a full rerun, the usual order is:

1. Run `notebooks/1.reassemble_metadata.ipynb` to rebuild `metadata_selected.csv`.
2. Run `notebooks/2.hl_revision_workflow.ipynb` or the shell wrappers to retrain models and export UMAP, CAM, tf-modisco, and per-run plotting outputs.
3. Run `notebooks/3.remake_fig3.ipynb`, `notebooks/4.remake_sfig4.ipynb`, and `notebooks/5.remake_sfig5.ipynb` to regenerate the manuscript figures and source tables.
4. Run `scripts/merge_csv_to_xlsx.py` to collect the exported statistics CSVs into one workbook if needed.

## Key Outputs

- `data/metadata_selected.csv`
- `results/plots/*.svg`
- `results/plots/*_source.csv`
- `results/plots/*stats*.csv`
- `results/models/*/plot_exports/*`
