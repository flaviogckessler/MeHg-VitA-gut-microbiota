# MeHg–Vitamin A Gut Microbiota

This repository supports the reproducibility of a multigenerational gut microbiota study investigating methylmercury (MeHg) exposure, vitamin A supplementation, and their combined effects. It brings together curated metadata, processed 16S rRNA sequencing products, analysis results, and figure-generation workflows in a portable project structure. Analyses are implemented in R [@r-core], with method-specific software cited below and in the figure guides.

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

## Analysis workflow

[![Bioinformatics analysis workflow](figures/Workflow/output/Analysis_Workflow_generated.svg)](figures/Workflow/output/Analysis_Workflow_generated.png)

The workflow begins with the manually defined 2 × 2 exposure design and F0–F2 cohorts. Two external 16S rRNA gene amplicon-sequencing batches are processed with DADA2, assigned SILVA 138.1 taxonomy, and combined with curated metadata to create analysis-ready ASV counts, taxonomy, and representative sequences.

From the analysis-ready data, the repository separates the analyses into two parallel branches:

- **Taxonomic ecology:** community composition, alpha diversity, differential-taxon testing, and the *Dorea*–mercury association produce Figures 3–4, 5D, and S1–S5. Figure S6 is a manual qualitative derivative of these taxonomic results.
- **Predicted functional profiles:** PICRUSt2 predicted KO/EC/MetaCyc profiles support parallel KO PCA, differential testing of predicted KEGG pathway profiles, KO heatmap and KEGG-map interpretation, and targeted predicted mercury-related ortholog analyses. These outputs map to Figures 6–8, S7–S16, and S18. Behavioral Figure S17 is outside this workflow.

Card labels and border patterns distinguish scripted, reconstructed, hybrid, and manual steps. The scientific DAG, wording, evidence status, execution status, visual mapping, and figure-output mappings are available in the machine-readable [`workflow_spec.json`](figures/Workflow/workflow_spec.json); an accessible text description is provided in [`Analysis_Workflow_alt_text.md`](figures/Workflow/output/Analysis_Workflow_alt_text.md). The rendering entry point and rasterizer are retained as implementation references, while the standalone validator and shared validation-helper source files are intentionally retained only in local working copies for now and excluded from version control. Consequently, end-to-end regeneration currently requires those local companion files.

## Repository status

This repository is under active reproducibility curation.

- The two batch-specific FASTQ workflows have been consolidated into a portable, linear [`DADA2 preprocessing script`](scripts/01_sequence_processing/preprocess_16s.R). Its configuration, versioned SILVA 138.1 Zenodo references, archived-data contracts, and batch-specific output serializers have been validated without raw reads. A full rerun and automatic comparison with the archived products remain pending public SRA release and an accession-to-sample mapping.
- Figure 3 panels A-K, Figure 4 panels A-D, Figure 5 panel D, Figure 6 panels A-C, Figure S1 panels A-H, all 75 taxon-abundance panels in Figures S2-S5, panels E of Figures 7 and 8, and all ten KO z-score heatmaps in Figures S8-S16 have portable, tested workflows. The 19 KO-level source tables shared by Figures 7, 8, and S8-S16 also have a portable, regression-tested generator. Their KEGG pathway panels remain retained manual exports prepared with the generated z-score tables as the coloring reference; supplied final composites are retained for visual validation.
- Figure S3 generated panels and final composite use the same authoritative treatment order and colors as the other abundance panels: `CN` (`#98F5FF`), `VitA` (`#B3EE3A`), `MeHg` (`#CD2626`), and `VitA+MeHg` (`#4D4D4D`).
- Figure S6 is a manual, qualitative synthesis made by visually reviewing the current Figure 5 and Figures S2-S5 composites. Its interpretive construction and reproducibility boundary are documented in the [`supplementary-figure guide`](supplementary_figures/README.md#figure-s6); it is not presented as the output of a statistical or scripted workflow.
- The PICRUSt2 terminal workflow, exact batch-specific inputs, analysis-ready outputs, and weighted-NSTI quality metrics have been harmonized and validated against the supplied pipeline products. Figure 6 uses only the canonical KO predictions and canonical sample manifest.
- Historical scripts cited by the provenance documentation are isolated under [`scripts/legacy_scripts/`](scripts/legacy_scripts/) and preserved unchanged. Unreferenced exploratory scripts and placeholder-only directories were removed during the repository-wide review.

## Repository structure

| Path | Contents |
|---|---|
| [`data/`](data/) | Data organization and provenance documentation. |
| [`data/metadata/`](data/metadata/) | Canonical sample manifests, preserved source metadata, SRA information, and validation summaries. |
| [`data/processed/`](data/processed/) | Processed ASV tables, taxonomy, representative sequences, BIOM files, relative abundances, and read tracking. |
| [`data/reference/`](data/reference/) | Small, versioned reference objects required by portable analysis workflows. |
| [`data/external/`](data/external/) | Download and checksum documentation plus the ignored local location for the SILVA 138.1 DADA2 reference FASTAs. |
| [`figures/`](figures/) | Main-figure workflows and lightweight generated outputs. |
| [`figures/Workflow/`](figures/Workflow/) | Machine-readable analysis-workflow specification plus the canonical SVG, matching PNG, accessibility text, schema, and rendering sources. |
| [`results/`](results/) | Alpha-diversity and taxon-abundance tables, functional-prediction exports, and legacy summary tables. |
| [`scripts/`](scripts/) | Maintained preprocessing and functional-analysis workflows plus documented historical provenance. |
| [`scripts/legacy_scripts/`](scripts/legacy_scripts/) | Unchanged legacy references already mapped to portable figure and PICRUSt2 workflows. |
| [`supplementary_figures/`](supplementary_figures/) | Supplementary-figure workflows, retained composites, and generated panel outputs. |
| [`bibliography.bib`](bibliography.bib) | Shared bibliography for software, databases, and analysis methods cited throughout the documentation. |

Canonical data and metadata are stored once under `data/`. Figure scripts should read these files directly rather than duplicate them within figure directories.

## Quick start

Sequence preprocessing is documented in the [`16S preprocessing guide`](scripts/01_sequence_processing/README.md). Raw source sequencing data are deposited in the SRA and are not distributed here; after downloading the DOI-pinned SILVA references and after SRA release, each batch is processed by setting `DADA2_BATCH_ID` and `RAW_FASTQ_DIR` and running `Rscript scripts/01_sequence_processing/preprocess_16s.R`. With the retained local validator companion present, validation-only mode checks the script, canonical manifests, archived processed-data contracts, reference configuration, and any downloaded reference checksums without FASTQ data. A full run automatically compares the generated batch-specific products with `data/processed/` and fails if their scientific contents differ.

The validated Figure 3 workflows use `tidyverse`, `scales`, `rstatix`, `ggpubr`, and `rlang`; Figures S1-S5 use `dplyr`, `ggplot2`, `rstatix`, and `ggpubr` [@tidyverse; @dplyr; @ggplot2; @scales; @rstatix; @ggpubr; @rlang]. The Figure 4 workflows use only packages distributed with R. Start R from the repository root and run:

```r
source("figures/Figure_3/generate_figure_3_abc.R")
generate_figure_3_abc()
source("figures/Figure_3/generate_figure_3_dg.R")
generate_figure_3_dg()
source("supplementary_figures/Figure_S1/generate_figure_s1_ah.R")
generate_figure_s1_ah()
source("supplementary_figures/Figure_S2/generate_figure_s2_abc.R")
generate_figure_s2_abc()
source("supplementary_figures/Figure_S3/generate_figure_s3_aq.R")
generate_figure_s3_aq()
source("supplementary_figures/Figure_S4/generate_figure_s4_a_an.R")
generate_figure_s4_a_an()
source("supplementary_figures/Figure_S5/generate_figure_s5_ao.R")
generate_figure_s5_ao()
source("figures/Figure_3/generate_figure_3_hk.R")
generate_figure_3_hk()
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
```

Standalone validator files used during repository curation are intentionally local-only and therefore are not included in the public quick-start commands.

The reconstructed historical Figure 6 environment uses `ggpicrust2` 2.0.0 and `ALDEx2` 1.28.1, with R 4.2.x (likely R 4.2.1) retained as contextual evidence; see the [main-figure guide](figures/README.md) for the provenance and its session-snapshot limitation. The generators additionally require `ggplot2`, `patchwork`, and `scales` [@ggpicrust2; @aldex2; @ggplot2; @patchwork; @scales], and the seeded ALDEx2 step can take several minutes per panel. The KO z-score source-data workflow uses base R plus the versioned KEGG reference stored under `data/reference/` [@kegg]; it does not require an installed `ggpicrust2` package. Generated figure outputs are written to each figure's `output/` directory; Figures S3-S5 further separate individual panels into `output/females/` and `output/male/`. Generated source tables stay with their batch under `results/functional_prediction/`.

See the [supplementary-figure guide](supplementary_figures/README.md) for its
figure sources and generation details. The supplied Figures S1-S5 final composites
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

Raw FASTQ files are not stored in this repository. They have been deposited in the NCBI Sequence Read Archive [@sra] under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663):

- Dams F0 and F1-PND21: `SAMN61406324`–`SAMN61406377`;
- F1-PND365 and F2-PND365: `SAMN61407672`–`SAMN61407719`.

The records are currently on hold and are intended for public release with the manuscript. The repository does not infer individual sample-to-BioSample assignments where these were not supplied.

The deposited material is the source sequencing data (raw reads), not the analysis source code. The versioned preprocessing source code is [`scripts/01_sequence_processing/preprocess_16s.R`](scripts/01_sequence_processing/preprocess_16s.R). A full raw-read rerun also requires the accession-to-`sample_id` mapping that was not included in the supplied deposition metadata.

## Figure-output policy

Generated figures are stored as compact PNG files suitable for inspection and reproducibility checks. The repository-wide analysis workflow is additionally committed as an accessible canonical SVG with a visually matching PNG. Manually supplied publication/reference artifacts may be retained in their original format, but figure-generation scripts must not create new TIFF files. See [`figures/README.md`](figures/README.md).

Main figures use `figures/Figure_<number>/output/`; supplementary figures with
portable panel workflows use `supplementary_figures/Figure_S<number>/output/`.
Supplied assembled composites use `Figure_<id>_FinalFigure.png`, while portable
panel outputs use `Figure<id><panel>_generated.png`. Descriptive incoming
filenames should be normalized to this convention once the figure number and
panel mapping are confirmed. Scientific cohort conflicts must be documented
and resolved rather than hidden by renaming an artifact.

## License

The repository is distributed under the [`LICENSE`](LICENSE) file.

## References

- <a id="ref-aldex2"></a>Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: Characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome, 2*, 15. [https://doi.org/10.1186/2049-2618-2-15](https://doi.org/10.1186/2049-2618-2-15)
- <a id="ref-rlang"></a>Henry, L., & Wickham, H. (2026). *rlang: Functions for base types and core R and Tidyverse features* (Version 1.3.0) [R package]. [https://doi.org/10.32614/CRAN.package.rlang](https://doi.org/10.32614/CRAN.package.rlang)
- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-ggpubr"></a>Kassambara, A., & Erdey, L. (2026). *ggpubr: ggplot2 based publication ready plots* [R package]. [https://doi.org/10.32614/CRAN.package.ggpubr](https://doi.org/10.32614/CRAN.package.ggpubr)
- <a id="ref-rstatix"></a>Kassambara, A. (2026). *rstatix: Pipe-friendly framework for basic statistical tests* [R package]. [https://doi.org/10.32614/CRAN.package.rstatix](https://doi.org/10.32614/CRAN.package.rstatix)
- <a id="ref-sra"></a>Leinonen, R., Sugawara, H., Shumway, M., & International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. *Nucleic Acids Research, 39*(Database issue), D19–D21. [https://doi.org/10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)
- <a id="ref-patchwork"></a>Pedersen, T. L. (2026). *patchwork: The composer of plots* [R package]. [https://doi.org/10.32614/CRAN.package.patchwork](https://doi.org/10.32614/CRAN.package.patchwork)
- <a id="ref-scales"></a>Posit Software, PBC. (2026). *scales: Scale functions for visualization* [R package]. [https://doi.org/10.32614/CRAN.package.scales](https://doi.org/10.32614/CRAN.package.scales)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
- <a id="ref-ggpicrust2"></a>Yang, C., Mai, J., Cao, X., Burberry, A., Cominelli, F., & Zhang, L. (2023). ggpicrust2: An R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics, 39*(8), btad470. [https://doi.org/10.1093/bioinformatics/btad470](https://doi.org/10.1093/bioinformatics/btad470)
- <a id="ref-dplyr"></a>Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2026). *dplyr: A grammar of data manipulation* (Version 1.2.1) [R package]. [https://doi.org/10.32614/CRAN.package.dplyr](https://doi.org/10.32614/CRAN.package.dplyr)
- <a id="ref-ggplot2"></a>Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd ed.). Springer. [https://doi.org/10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
- <a id="ref-tidyverse"></a>Wickham, H., Averick, M., Bryan, J., et al. (2019). Welcome to the Tidyverse. *Journal of Open Source Software, 4*(43), 1686. [https://doi.org/10.21105/joss.01686](https://doi.org/10.21105/joss.01686)
