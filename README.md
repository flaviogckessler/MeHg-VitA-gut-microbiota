# MeHg–Vitamin A Gut Microbiota

This repository supports the reproducibility of a multigenerational gut microbiota study investigating methylmercury (MeHg) exposure, vitamin A supplementation, and their combined effects. It brings together curated metadata, processed 16S rRNA sequencing products, analysis results, and figure-generation workflows in a portable project structure.

## Study overview

The experimental groups are:

- `CN`: control;
- `VitA`: vitamin A;
- `MeHg`: methylmercury;
- `VitA+MeHg`: combined vitamin A and methylmercury.

The microbiota data were generated in two sequencing batches:

| Batch | Biological material | Time point | Platform |
|---|---|---|---|
| `batch01_iontorrent_dams-f0_f1-pnd21` | F0 dams and F1 offspring | F1 at postnatal day 21 | Ion Torrent |
| `batch02_illumina-miseq_f1-f2_pnd365` | F1 and F2 offspring | Postnatal day 365 | Illumina MiSeq |

Sex, generation, treatment, sequencing batch, analysis inclusion, and data-deposition information are standardized in [`data/metadata/sample_manifest_all_batches.csv`](data/metadata/sample_manifest_all_batches.csv).

## Repository status

This repository is under active reproducibility curation.

- Figure 3 panels A-K and Figure 4 panels A-D have portable, validated workflows.
- The PICRUSt2 terminal workflow, exact batch-specific inputs, analysis-ready outputs, and weighted-NSTI quality metrics have been harmonized and validated against the supplied pipeline products; Figure 6 generation remains pending.
- Additional main and supplementary figure workflows will be migrated and validated progressively.
- Legacy scripts already used as provenance for portable workflows are isolated under [`scripts/legacy_scripts/`](scripts/legacy_scripts/). Unused legacy scripts remain directly under `scripts/` until their corresponding analyses are curated.

## Repository structure

| Path | Contents |
|---|---|
| [`config/`](config/) | Reserved for project-wide configuration files. |
| [`data/`](data/) | Data organization and provenance documentation. |
| [`data/metadata/`](data/metadata/) | Canonical sample manifests, preserved source metadata, SRA information, and validation summaries. |
| [`data/processed/`](data/processed/) | Processed ASV tables, taxonomy, representative sequences, BIOM files, relative abundances, and read tracking. |
| [`figures/`](figures/) | Main-figure workflows, validation scripts, and lightweight generated PNG outputs. |
| [`results/`](results/) | Alpha-diversity and taxon-abundance tables, functional-prediction exports, and legacy summary tables. |
| [`scripts/`](scripts/) | Analysis and curation scripts, including unused legacy scripts awaiting migration. |
| [`scripts/legacy_scripts/`](scripts/legacy_scripts/) | Unchanged legacy references already mapped to portable figure and PICRUSt2 workflows. |
| [`supplementary_figures/`](supplementary_figures/) | Reserved for supplementary-figure workflows and outputs. |
| [`docs/`](docs/) | Reserved for extended project documentation. |

Canonical data and metadata are stored once under `data/`. Figure scripts should read these files directly rather than duplicate them within figure directories.

## Quick start

The validated Figure 3 workflows require R and the `tidyverse`, `rstatix`, and `ggpubr` packages. The Figure 4 workflows use base R only. Start R from the repository root and run:

```r
source("figures/Figure_3/validate_figure_3_abc.R")
source("figures/Figure_3/generate_figure_3_abc.R")
generate_figure_3_abc()
source("figures/Figure_3/validate_figure_3_dg.R")
source("figures/Figure_3/generate_figure_3_dg.R")
generate_figure_3_dg()
source("figures/Figure_3/validate_figure_3_hk.R")
source("figures/Figure_3/generate_figure_3_hk.R")
generate_figure_3_hk()
source("figures/Figure_4/validate_figure_4.R")
source("figures/Figure_4/generate_figure_4_a.R")
generate_figure_4_a()
source("figures/Figure_4/generate_figure_4_b.R")
generate_figure_4_b()
source("figures/Figure_4/generate_figure_4_c.R")
generate_figure_4_c()
source("figures/Figure_4/generate_figure_4_d.R")
generate_figure_4_d()
```

The validation scripts check sample inclusion, abundance values, statistical methods, legacy checkpoints where retained, plot construction, and output generation. Generated panels are written to [`figures/Figure_3/output/`](figures/Figure_3/output/) and [`figures/Figure_4/output/`](figures/Figure_4/output/).

See [`figures/README.md`](figures/README.md) for figure sources and generation details.

The reconstructed Linux-terminal PICRUSt2 workflow is documented in [`scripts/04_functional_prediction/picrust2_workflow.md`](scripts/04_functional_prediction/picrust2_workflow.md). A full rerun is optional and time-consuming; the retained EC, KO, and MetaCyc matrices were validated line-by-line against the supplied original pipeline outputs.

## Data organization and analysis rules

- Use `data/metadata/sample_manifest_all_batches.csv` as the canonical sample-level metadata table.
- Join sample-level data using `sample_id`.
- Honor the `include_primary_analysis` field when selecting samples.
- Samples `10F` and `22F` are preserved for provenance but excluded from the primary analysis because of low sequencing quality.
- The unexpected batch 02 count-table column `t3` must not be analyzed unless valid metadata are supplied.
- Do not duplicate canonical data or metadata inside figure directories.
- Generate new figure outputs as lightweight, readable PNG files; do not generate new TIFF files.

Further details are available in the [`data`](data/README.md), [`metadata`](data/metadata/README.md), [`processed data`](data/processed/README.md), [`figures`](figures/README.md), and [`results`](results/README.md) documentation.

## Raw sequence availability

Raw FASTQ files are not stored in this repository. They have been deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663):

- Dams F0 and F1-PND21: `SAMN61406324`–`SAMN61406377`;
- F1-PND365 and F2-PND365: `SAMN61407672`–`SAMN61407719`.

The records are currently on hold and are intended for public release with the manuscript. The repository does not infer individual sample-to-BioSample assignments where these were not supplied.

## Figure-output policy

Generated figures are stored as compact PNG files suitable for inspection and reproducibility checks. Manually supplied publication/reference artifacts may be retained in their original format, but figure-generation scripts must not create new TIFF files. See [`figures/README.md`](figures/README.md).

## Citation and license

Citation information will be added when the associated manuscript is available. A repository license has not yet been assigned; one should be selected before public release.
