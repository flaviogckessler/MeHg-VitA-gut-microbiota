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

- Figure 3 panels A-K, Figure 4 panels A-D, Figure 5 panel D, Figure 6 panels A-C, Figure S1 panels A-H, all 75 taxon-abundance panels in Figures S2-S5, panels E of Figures 7 and 8, and all ten KO z-score heatmaps in Figures S8-S16 have portable, tested workflows. The 19 KO-level source tables shared by Figures 7, 8, and S8-S16 also have a portable, regression-tested generator. Their KEGG pathway panels remain retained manual exports prepared with the generated z-score tables as the coloring reference; supplied final composites are retained for visual validation.
- Figure S3 generated panels correct the erroneous treatment ordering in the historical PND21 exports. Their authoritative order is `CN`, `VitA`, `MeHg`, `VitA+MeHg` with blue, green, red, and black fills; the historical composite is retained only for provenance.
- Figure S6 is a retained manual, qualitative synthesis made by visually reviewing Figure 5 and Figures S2-S5. Its interpretive construction and reproducibility boundary are documented in the [`supplementary-figure guide`](supplementary_figures/README.md#figure-s6); it is not presented as the output of a statistical or scripted workflow.
- The PICRUSt2 terminal workflow, exact batch-specific inputs, analysis-ready outputs, and weighted-NSTI quality metrics have been harmonized and validated against the supplied pipeline products. Figure 6 uses only the canonical KO predictions and canonical sample manifest.
- Additional main and supplementary figure workflows will be migrated and validated progressively.
- Legacy scripts already used as provenance for portable workflows are isolated under [`scripts/legacy_scripts/`](scripts/legacy_scripts/). Unused legacy scripts remain directly under `scripts/` until their corresponding analyses are curated.

## Repository structure

| Path | Contents |
|---|---|
| [`config/`](config/) | Reserved for project-wide configuration files. |
| [`data/`](data/) | Data organization and provenance documentation. |
| [`data/metadata/`](data/metadata/) | Canonical sample manifests, preserved source metadata, SRA information, and validation summaries. |
| [`data/processed/`](data/processed/) | Processed ASV tables, taxonomy, representative sequences, BIOM files, relative abundances, and read tracking. |
| [`data/reference/`](data/reference/) | Small, versioned reference objects required by portable analysis workflows. |
| [`figures/`](figures/) | Main-figure workflows, validation scripts, and lightweight generated PNG outputs. |
| [`results/`](results/) | Alpha-diversity and taxon-abundance tables, functional-prediction exports, and legacy summary tables. |
| [`scripts/`](scripts/) | Analysis and curation scripts, including unused legacy scripts awaiting migration. |
| [`scripts/legacy_scripts/`](scripts/legacy_scripts/) | Unchanged legacy references already mapped to portable figure and PICRUSt2 workflows. |
| [`supplementary_figures/`](supplementary_figures/) | Supplementary-figure workflows, validation scripts, retained composites, and generated panel outputs. |
| [`docs/`](docs/) | Reserved for extended project documentation. |

Canonical data and metadata are stored once under `data/`. Figure scripts should read these files directly rather than duplicate them within figure directories.

## Quick start

The validated Figure 3 and Figure S1 alpha-diversity workflows and the Figures S2-S5 taxon-abundance workflows require R and the `tidyverse`, `rstatix`, and `ggpubr` packages. The Figure 4 workflows use base R only. Start R from the repository root and run:

```r
source("figures/Figure_3/validate_figure_3_abc.R")
source("figures/Figure_3/generate_figure_3_abc.R")
generate_figure_3_abc()
source("figures/Figure_3/validate_figure_3_dg.R")
source("figures/Figure_3/generate_figure_3_dg.R")
generate_figure_3_dg()
source("supplementary_figures/Figure_S1/validate_figure_s1_ah.R")
source("supplementary_figures/Figure_S1/generate_figure_s1_ah.R")
generate_figure_s1_ah()
source("supplementary_figures/Figure_S2/validate_figure_s2_abc.R")
source("supplementary_figures/Figure_S2/generate_figure_s2_abc.R")
generate_figure_s2_abc()
source("supplementary_figures/Figure_S3/validate_figure_s3_aq.R")
source("supplementary_figures/Figure_S3/generate_figure_s3_aq.R")
generate_figure_s3_aq()
source("supplementary_figures/Figure_S4/validate_figure_s4_a_an.R")
source("supplementary_figures/Figure_S4/generate_figure_s4_a_an.R")
generate_figure_s4_a_an()
source("supplementary_figures/Figure_S5/validate_figure_s5_ao.R")
source("supplementary_figures/Figure_S5/generate_figure_s5_ao.R")
generate_figure_s5_ao()
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
source("figures/Figure_6/generate_figure_6_a.R")
generate_figure_6_a()
source("figures/Figure_6/generate_figure_6_b.R")
generate_figure_6_b()
source("figures/Figure_6/generate_figure_6_c.R")
generate_figure_6_c()
source("scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R")
generate_kegg_pathway_ko_z_scores()
source("scripts/04_functional_prediction/validate_kegg_pathway_ko_z_scores.R")
validate_kegg_pathway_ko_z_scores()
```

The Figure 6 generators additionally require `ggpicrust2`, `ALDEx2`, `ggplot2`, `patchwork`, and `scales`; its seeded ALDEx2 step can take several minutes per panel. The KO z-score source-data workflow uses base R plus the versioned KEGG reference stored under `data/reference/`; it does not require an installed `ggpicrust2` package. Generated figure outputs are written to each figure's `output/` directory; Figures S3-S5 further separate individual panels into `output/females/` and `output/male/`. Generated source tables stay with their batch under `results/functional_prediction/`.

See the [main-figure guide](figures/README.md) and
[supplementary-figure guide](supplementary_figures/README.md) for figure sources
and generation details. The supplied Figures S1-S5 final composites
were assembled manually in Adobe Photoshop and are retained as visual
references rather than R-generated outputs.

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

Main figures use `figures/Figure_<number>/output/`; supplementary figures with
portable panel workflows use `supplementary_figures/Figure_S<number>/output/`.
Supplied assembled composites use `Figure_<id>_FinalFigure.png`, while portable
panel outputs use `Figure<id><panel>_generated.png`. Descriptive incoming
filenames should be normalized to this convention once the figure number and
panel mapping are confirmed. Scientific cohort conflicts must be documented
and resolved rather than hidden by renaming an artifact.

## Citation and license

Citation information will be added when the associated manuscript is available. A repository license has not yet been assigned; one should be selected before public release.
