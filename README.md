# Supporting repository for the publication: The periphery of nuclear speckles defines a spatially and temporally regulated compartment of long-lived intron-retained RNAs that resolves during mitosis


This repository contains the Python notebooks, R scripts, R Markdown notebooks, shell wrappers, configuration files, and exported analysis outputs used to regenerate manuscript figures and train models to predict IR and IR HL.

## Repository Layout

```text
.
├── README.md
├── environment-parnet_clean.yml
├── environment-tfmodisco.yml
├── environment-r.yml
├── requirements.txt
├── configs/
├── data/
├── notebooks/
├── R/
├── results/
└── scripts/
```

## Environments

Three environments are used in this repository:

- `parnet-clean` for metadata assembly, model training, UMAP export, CAM export, figure generation, and Python notebooks
- `tfmodisco` for motif discovery and motif report generation
- `r` for the R analyses under `R/`

### Main Analysis Environment

Create the main environment:

```bash
conda env create -f environment-parnet_clean.yml
conda activate parnet-clean
```

Install the Python packages exported from the original `parnet-clean` environment:

```bash
pip install -r requirements.txt
```

Install `ir_toolkit` from GitHub:

```bash
pip install git+https://github.com/melonheader/ir_toolkit.git
```

### tf-modisco Environment

Create the dedicated motif-analysis environment separately:

```bash
conda env create -f environment-tfmodisco.yml
conda activate tfmodisco
```

This environment provides `modisco-lite` together with the MEME suite tools used by `tomtom`.

### R Analysis Environment

Create the R environment for the analyses in `R/`:

```bash
conda env create -f environment-r.yml
conda activate r
```

## Running The Python Notebooks

Launch Jupyter from the main analysis environment:

```bash
conda activate parnet-clean
jupyter lab
```

Use the `parnet-clean` kernel for the notebooks in `notebooks/`. The tf-modisco step is executed through the shell wrapper and uses the dedicated `tfmodisco` environment internally.

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

## Running The R Workflows

The R code is organized by figure:

- `R/Figure_1/`
  Figure 1 length, GC content, and nuclear enrichment analysis.
- `R/Figure_2/`
  Figure 2 stability analysis notebook.
- `R/Figure_4/`
  Figure 4 TSA-seq overlap analysis.
- `R/Figure_4_distance/RANDOM_SPOTS/`
  Random-spot null-distribution analyses for iPS and HUVEC distance-to-speckle measurements.

These workflows assume they are launched from their local `scripts/` directory. Example commands:

```bash
conda activate r

cd R/Figure_1/scripts
Rscript Figure_1_length_gc_nuclear_enrichment.R

cd ../../Figure_2/scripts
Rscript -e "rmarkdown::render('Figure_2.Rmd')"

cd ../../Figure_4/scripts
Rscript Figure_4_TSAseq_IR_stability.R

cd ../../Figure_4_distance/RANDOM_SPOTS/scripts
Rscript FINAL_iPS_IR_RNAS_vs_random_spots_KS.R
Rscript FINAL_HUVEC_IR_RNAS_vs_random_spots_KS.R
```

## Model Retraining Workflow

The workflow is defined by:

- `configs/hl_revision_runs.json`
- `scripts/hl_revision_pipeline.py`
- `scripts/run_hl_revision_workflow.sh`
- `scripts/run_all_hl_revision_workflow.sh`
- `scripts/run_modisco_report.sh`
- `notebooks/2.hl_revision_workflow.ipynb`

### Inspect Configured Runs

```bash
python scripts/hl_revision_pipeline.py list-runs
```

### Train One Configured Run

```bash
python scripts/hl_revision_pipeline.py train --run hl_revised_50percgap
```

## Shell Wrappers

For a single run, the main wrapper is:

```bash
scripts/run_hl_revision_workflow.sh \
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
scripts/run_all_hl_revision_workflow.sh \
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
python scripts/merge_csv_to_xlsx.py \
  --input-dir results/plots \
  --output-xlsx results/plots/statistics_tables_unified.xlsx
```

## Raw Data Processing Helpers

The repository also includes the shell scripts used to download and process the PRJNA608890 raw sequencing data under `scripts/raw_data/PRJNA608890/`.

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
