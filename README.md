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

- Figure 3 panels A-K have portable, validated workflows.
- Additional main and supplementary figure workflows will be migrated and validated progressively.
- Scripts retained directly under `scripts/` are legacy source scripts. Some still contain historical paths, exploratory code, or obsolete comments and should not be considered portable until they are migrated into the corresponding analysis or figure directory.

## Repository structure

| Path | Contents |
|---|---|
| [`config/`](config/) | Reserved for project-wide configuration files. |
| [`data/`](data/) | Data organization and provenance documentation. |
| [`data/metadata/`](data/metadata/) | Canonical sample manifests, preserved source metadata, SRA information, and validation summaries. |
| [`data/processed/`](data/processed/) | Processed ASV tables, taxonomy, representative sequences, BIOM files, relative abundances, and read tracking. |
| [`figures/`](figures/) | Main-figure workflows, validation scripts, and lightweight generated PNG outputs. |
| [`results/`](results/) | Alpha-diversity and taxon-abundance tables, functional-prediction exports, and legacy summary tables. |
| [`scripts/`](scripts/) | Analysis scripts, including legacy scripts awaiting reproducibility migration. |
| [`supplementary_figures/`](supplementary_figures/) | Reserved for supplementary-figure workflows and outputs. |
| [`docs/`](docs/) | Reserved for extended project documentation. |

Canonical data and metadata are stored once under `data/`. Figure scripts should read these files directly rather than duplicate them within figure directories.

## Quick start

The validated Figure 3 workflows require R and the `tidyverse`, `rstatix`, and `ggpubr` packages. Start R from the repository root and run:

```r
source("figures/Figure_3/validate_figure_3_abc.R")
source("figures/Figure_3/generate_figure_3_abc.R")
source("figures/Figure_3/validate_figure_3_dg.R")
source("figures/Figure_3/generate_figure_3_dg.R")
source("figures/Figure_3/validate_figure_3_hk.R")
source("figures/Figure_3/generate_figure_3_hk.R")
```

The validation scripts check sample inclusion, abundance values, statistical methods, plot construction, and output generation. The generation scripts write readable PNG panels to [`figures/Figure_3/output/`](figures/Figure_3/output/).

See [`figures/README.md`](figures/README.md) for figure sources and generation details.

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
