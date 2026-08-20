# Supplementary figures

Standalone validator source files used during repository curation are retained locally but intentionally excluded from version control for now. The commands below therefore document public generators only.

## Software, R packages, and sample-size convention

`N` denotes biological samples (animals), not taxa or KOs. Microbiome counts are reported after applying `include_primary_analysis`. Packages distributed with R are cited collectively as R [@r-core].

| Figure(s) | Software or required R packages | Analysis N |
|---|---|---|
| S1 | `dplyr`, `ggplot2`, `ggpubr`, and `rstatix`, plus base-R `stats` and `utils` [@dplyr; @ggplot2; @ggpubr; @rstatix; @r-core] | Panels A/E: 18; B/F: 34; C/G: 24; D/H: 24. |
| S2 | `dplyr`, `ggplot2`, `ggpubr`, and `rstatix`, plus base-R `stats`, `tools`, and `utils` [@dplyr; @ggplot2; @ggpubr; @rstatix; @r-core] | 18 dams. |
| S3 | Same as S2 | 34 F1-PND21 animals (18 female, 16 male). |
| S4 | Same as S2 | 24 F1-PND365 animals (12/sex). |
| S5 | Same as S2 | 24 F2-PND365 animals (12/sex). |
| S6 | Adobe Photoshop/manual synthesis; no R packages [@adobe-photoshop] | No new observations; inherits S2-S5 cohort N values (18, 34, 24, and 24). |
| S7 | `ggplot2` plus base-R `stats` and `utils` [@ggplot2; @r-core] | Panels A-D: 18, 34, 24, and 24. |
| S8-S14 | Base-R `graphics`, `grDevices`, `stats`, and `utils`; KEGG for manual pathway views [@r-core; @kegg] | 34 F1-PND21 animals per figure. |
| S15 | Base-R `graphics`, `grDevices`, `stats`, and `utils`; KEGG for manual pathway views [@r-core; @kegg] | Panel I: 34 F1-PND21; panel J: 24 F1-PND365. |
| S16 | Base-R `graphics`, `grDevices`, `stats`, and `utils`; KEGG for manual pathway views [@r-core; @kegg] | 24 F2-PND365 animals. |
| S17 | GraphPad Prism [@graphpad-prism]; no R packages | Endpoint-specific N is transcribed below from non-empty Prism values. |
| S18 | `ggplot2` plus base-R `grDevices`, `grid`, and `stats` [@ggplot2; @r-core] | Cohort columns: 18, 34, 24, and 24. |

## Figure S1

Figure S1 reports Observed richness (panels A-D) and InvSimpson diversity (panels E-H) for the same four cohorts used in Figure 3. Panels A/E use 18 F0 dams (`CN` 5, `VitA` 4, `MeHg` 5, `VitA+MeHg` 4); panels B/F use 34 F1-PND21 offspring (`CN` 10, `VitA` 7, `MeHg` 10, `VitA+MeHg` 7); panels C/G and D/H each use 24 animals, with three animals per sex-by-treatment cell.

The official generator is [`generate_figure_s1_ah.R`](Figure_S1/generate_figure_s1_ah.R). During curation, a local-only validator checked the canonical alpha-diversity inputs, cohort membership, statistical results, output dimensions, and supplied composite. The generator requires `dplyr`, `ggplot2`, `ggpubr`, and `rstatix`, together with base-R `stats` and `utils`.

```bash
Rscript supplementary_figures/Figure_S1/generate_figure_s1_ah.R
```

## Figure S2

Figure S2 contains three individual Dams (F0) relative-abundance boxplots. Panels are assigned from left to right in the supplied composite. The portable workflow reads the canonical batch-01 [`taxon-relative-abundance table`](../results/taxon_relative_abundance/batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv) and joins its sample columns to the batch-01 [`sample manifest`](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv). It honors `include_primary_analysis` and selects the 18 F0 dams: 5 `CN`, 4 `VitA`, 5 `MeHg`, and 4 `VitA+MeHg` samples. No Figure S2-specific source-data copy is required because the historical `AbundanciasTaxons_Todos.csv` values are already preserved under the repository-compliant canonical filename.

| Panel | Canonical taxon row | Display label | Significant FDR-adjusted Dunn comparisons |
|------------------|------------------|------------------|------------------|
| A | `g__[Eubacterium] ventriosum group` | `g__.Eubacterium..ventriosum.group` | `CN`-`MeHg`; `VitA`-`MeHg`; `MeHg`-`VitA+MeHg` |
| B | `g__[Ruminococcus] gauvreauii group` | `g__.Ruminococcus..gauvreauii.group` | `CN`-`VitA+MeHg`; `MeHg`-`VitA+MeHg` |
| C | `g__Oribacterium` | `g__Oribacterium` | `MeHg`-`VitA+MeHg` |

All panels use a four-group Kruskal-Wallis test followed by all six pairwise Dunn tests with Benjamini-Hochberg FDR adjustment. Only adjusted comparisons at or below 0.05 are annotated. This deliberately nonparametric workflow reproduces the tests and labels in the historical source panels; it does not choose a method from a normality test.

[`generate_figure_s2_abc.R`](Figure_S2/generate_figure_s2_abc.R) uses the shared [`figure_s2_s5_abundance_common.R`](figure_s2_s5_abundance_common.R) workflow and generates the three individual 600 x 500 PNGs [`FigureS2A_generated.png`](Figure_S2/output/FigureS2A_generated.png), [`FigureS2B_generated.png`](Figure_S2/output/FigureS2B_generated.png), and [`FigureS2C_generated.png`](Figure_S2/output/FigureS2C_generated.png). It does not assemble the publication layout. The supplied [`Figure_S2_FinalFigure.png`](Figure_S2/output/Figure_S2_FinalFigure.png) was assembled manually in Adobe Photoshop and is retained only as a visual reference. The three original 600 x 500 TIFF panel exports were inspected from the archived `Alberto/images/Abundancias/NONparametric/Dams_F0_NONparametric/` directory; the portable workflow creates PNG files only.

Generate all three panels from the repository root:

```bash
Rscript supplementary_figures/Figure_S2/generate_figure_s2_abc.R
```

During curation, a local-only validator checked source schemas, taxon identity, F0 membership, treatment counts, finite relative abundances, the exact Kruskal-Wallis results, all six Dunn comparisons per panel, the expected FDR-adjusted significant pairs, PNG names and dimensions, generated panels, and the supplied composite.

## Figures S3-S5

Figures S3-S5 extend the same portable taxon-relative-abundance workflow to offspring cohorts. Together they generate 72 individual 600 x 500 PNG panels; they do not assemble the supplied Photoshop composites.

| Figure | Cohort | N (female/male) | Female panels | Male panels | Numerical source | Metadata source | Generator |
|---|---|---:|---|---|---|---|---|
| S3 | F1-PND21 | 18/16 | [`A-J`](Figure_S3/output/females/) | [`K-Q`](Figure_S3/output/male/) | [`batch01` abundance](../results/taxon_relative_abundance/batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv) | [`batch01` manifest](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv) | [`generator`](Figure_S3/generate_figure_s3_aq.R) |
| S4 | F1-PND365 | 12/12 | [`A-AC`](Figure_S4/output/females/) | [`AD-AN`](Figure_S4/output/male/) | [`batch02` abundance](../results/taxon_relative_abundance/batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv) | [`batch02` manifest](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv) | [`generator`](Figure_S4/generate_figure_s4_a_an.R) |
| S5 | F2-PND365 | 12/12 | [`A-H`](Figure_S5/output/females/) | [`I-O`](Figure_S5/output/male/) | [`batch02` abundance](../results/taxon_relative_abundance/batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv) | [`batch02` manifest](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv) | [`generator`](Figure_S5/generate_figure_s5_ao.R) |

For each figure, female generated panels are stored under `output/females/` and male generated panels under `output/male/`. The manually assembled final composite remains directly under `output/`. The shared generator derives the destination from each panel's manifest sex; local-only curation checks verified the sex-specific paths and absence of generated panels at the output root.

The canonical abundance files are exact copies of the historical `AbundanciasTaxons_Todos.csv` and `AbundanciasTaxons_1year.csv` inputs, so no new figure-specific data copies or renaming were needed. Together with Figure S2, these generators use the shared [`figure_s2_s5_abundance_common.R`](figure_s2_s5_abundance_common.R) renderer, join abundance columns to the relevant manifest, honor `include_primary_analysis`, and select generation and sex explicitly. Each taxon is tested with a four-group Kruskal-Wallis test followed by all six Dunn tests with Benjamini-Hochberg FDR adjustment; only adjusted comparisons at or below 0.05 are drawn. Every workflow uses the common treatment order `CN`, `VitA`, `MeHg`, `VitA+MeHg` and the shared colors `#98F5FF`, `#B3EE3A`, `#CD2626`, `#4D4D4D`, respectively. Figure S3 generated panels and [`Figure_S3_FinalFigure.png`](Figure_S3/output/Figure_S3_FinalFigure.png) use this same order and palette.

Generate the panels from the repository root:

```bash
Rscript supplementary_figures/Figure_S3/generate_figure_s3_aq.R
Rscript supplementary_figures/Figure_S4/generate_figure_s4_a_an.R
Rscript supplementary_figures/Figure_S5/generate_figure_s5_ao.R
```

Local-only curation checks independently reproduced the Kruskal-Wallis statistics and checked the exact omnibus p-values, Dunn comparison signatures, adjusted significance counts, source checksums and schemas, cohort membership and treatment counts, common treatment order and colors, output names, sex-specific directories and dimensions, generated PNGs, and final-composite dimensions. The authoritative composites are [`Figure S3`](Figure_S3/output/Figure_S3_FinalFigure.png), [`Figure S4`](Figure_S4/output/Figure_S4_FinalFigure.png), and [`Figure S5`](Figure_S5/output/Figure_S5_FinalFigure.png).

## Figure S6

Figure S6 is a manually prepared qualitative synthesis, not the output of a computational figure-generation workflow. The author made it after side-by-side visual inspection of the current final [`Figure 5`](../figures/Figure_5/output/Figure_5_FinalFigure.png) and [`Figures S2`](Figure_S2/output/Figure_S2_FinalFigure.png), [`S3`](Figure_S3/output/Figure_S3_FinalFigure.png), [`S4`](Figure_S4/output/Figure_S4_FinalFigure.png), and [`S5`](Figure_S5/output/Figure_S5_FinalFigure.png).

To construct the summary, the author reviewed the treatment patterns and displayed comparisons in those figures within each available cohort and sex, then manually placed the taxa judged to be associated with `VitA`, `MeHg`, or `VitA+MeHg` into the corresponding box. Green, red, and black borders identify those three association categories, respectively. Dams (F0) are shown without a sex subdivision; F1-PND21, F1-PND365, and F2-PND365 are divided into females and males. A taxon may occur in more than one box when the visual interpretation supported more than one association, and an empty box means that no taxon was manually assigned to that cohort/sex/category combination.

This classification is descriptive and interpretive. Figure S6 does not add a formal statistical test, pooled analysis, quantitative selection threshold, or automated rule beyond the analyses already displayed in its source figures. Consequently, it should be reproduced by repeating the documented visual review and manual transcription rather than by running an R script. The authoritative output is the 1836 x 1281 raster [`Figure_S6_FinalFigure.png`](Figure_S6/Figure_S6_FinalFigure.png); no editable layout file or Figure S6 generator is available.

Figure S6 adds no observations of its own. Its evidence base inherits the included cohort sizes shown in Figures S2-S5: N=18 F0 dams, N=34 F1-PND21, N=24 F1-PND365, and N=24 F2-PND365.

## Figure S7

Figure S7 is a four-panel PCA of the untransformed PICRUSt-predicted KEGG Ortholog (KO) scores. The official portable workflow is [`generate_figure_s7_ad.R`](Figure_S7/generate_figure_s7_ad.R). It generates the individual plots used in the supplied, manually assembled [`Figure_S7_FinalFigure.png`](Figure_S7/output/Figure_S7_FinalFigure.png).

| Panel | Cohort | Numerical source | Metadata source | Generated PNG |
|---------------|---------------|---------------|---------------|---------------|
| A | Dams (F0) | [`batch01` KO predictions](../results/functional_prediction/batch01_iontorrent_dams-f0_f1-pnd21/kegg_ortholog_abundance.tsv) | [`batch01` manifest](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv) | [`FigureS7A_generated.png`](Figure_S7/output/FigureS7A_generated.png) |
| B | F1-PND21 | [`batch01` KO predictions](../results/functional_prediction/batch01_iontorrent_dams-f0_f1-pnd21/kegg_ortholog_abundance.tsv) | [`batch01` manifest](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv) | [`FigureS7B_generated.png`](Figure_S7/output/FigureS7B_generated.png) |
| C | F1-PND365 | [`batch02` KO predictions](../results/functional_prediction/batch02_illumina-miseq_f1-f2_pnd365/kegg_ortholog_abundance.tsv) | [`batch02` manifest](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv) | [`FigureS7C_generated.png`](Figure_S7/output/FigureS7C_generated.png) |
| D | F2-PND365 | [`batch02` KO predictions](../results/functional_prediction/batch02_illumina-miseq_f1-f2_pnd365/kegg_ortholog_abundance.tsv) | [`batch02` manifest](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv) | [`FigureS7D_generated.png`](Figure_S7/output/FigureS7D_generated.png) |

The generator joins matrix columns to the cohort manifests, honors `include_primary_analysis`, and therefore excludes low-quality batch01 samples `10F` and `22F` and the undocumented batch02 column `t3`. It applies base-R `prcomp()` to samples by all KOs with centering and no variance scaling, matching the legacy analysis. Treatment is encoded by the original ggplot2 discrete colors and 68% normal ellipses; sex is encoded by point shape for the offspring cohorts. The panel headings in the supplied composite were added during manual assembly and are not repeated inside the generated plot areas.

Run the generator from the repository root:

```bash
Rscript supplementary_figures/Figure_S7/generate_figure_s7_ad.R
```

The same script validates the canonical source schemas and dimensions, cohort membership (18, 34, 24, and 24 samples), finite non-negative KO scores, treatment ellipse covariance, deterministic PCA orientation, the four reference variance pairs, each 1800 x 900 generated PNG, and the supplied 2082 x 1179 composite. No new source-data file is needed: the canonical, descriptively named PICRUSt2 matrices and sample manifests are used directly.

## Figures S8-S16

Figures S8-S16 follow the same reproducibility boundary as main Figures 7 and 8. The KEGG pathway views were prepared manually, while their KO z-score heatmaps now have portable base-R renderers. The supplied final composites are retained as visual references under each figure's `output/` directory.

Figures S8-S14 each use the 34 included F1-PND21 samples (`CN` 10, `VitA` 7, `MeHg` 10, `VitA+MeHg` 7). Figure S15 uses those 34 samples for panel I and 24 F1-PND365 samples for panel J (six per treatment, three per sex-by-treatment cell). Figure S16 uses 24 F2-PND365 samples with the same balanced allocation.

| Figure | Pathway | KEGG map | Manuscript cohort | Heatmap panel(s) | Generator | Supplied composite |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| S8 | Non-homologous end-joining | [`map03450`](https://www.kegg.jp/pathway/map03450) | F1-PND21 | E | [`generate_figure_s8_e.R`](Figure_S8/generate_figure_s8_e.R) | [`Figure_S8_FinalFigure.png`](Figure_S8/output/Figure_S8_FinalFigure.png) |
| S9 | Oxidative phosphorylation | [`map00190`](https://www.kegg.jp/pathway/map00190) | F1-PND21 | E | [`generate_figure_s9_e.R`](Figure_S9/generate_figure_s9_e.R) | [`Figure_S9_FinalFigure.png`](Figure_S9/output/Figure_S9_FinalFigure.png) |
| S10 | Riboflavin metabolism | [`map00740`](https://www.kegg.jp/pathway/map00740) | F1-PND21 | E | [`generate_figure_s10_e.R`](Figure_S10/generate_figure_s10_e.R) | [`Figure_S10_FinalFigure.png`](Figure_S10/output/Figure_S10_FinalFigure.png) |
| S11 | Nicotinate and nicotinamide metabolism | [`map00760`](https://www.kegg.jp/pathway/map00760) | F1-PND21 | E | [`generate_figure_s11_e.R`](Figure_S11/generate_figure_s11_e.R) | [`Figure_S11_FinalFigure.png`](Figure_S11/output/Figure_S11_FinalFigure.png) |
| S12 | Taurine and hypotaurine metabolism | [`map00430`](https://www.kegg.jp/pathway/map00430) | F1-PND21 | E | [`generate_figure_s12_e.R`](Figure_S12/generate_figure_s12_e.R) | [`Figure_S12_FinalFigure.png`](Figure_S12/output/Figure_S12_FinalFigure.png) |
| S13 | Toluene degradation | [`map00623`](https://www.kegg.jp/pathway/map00623) | F1-PND21 | E | [`generate_figure_s13_e.R`](Figure_S13/generate_figure_s13_e.R) | [`Figure_S13_FinalFigure.png`](Figure_S13/output/Figure_S13_FinalFigure.png) |
| S14 | Caprolactam degradation | [`map00930`](https://www.kegg.jp/pathway/map00930) | F1-PND21 | E | [`generate_figure_s14_e.R`](Figure_S14/generate_figure_s14_e.R) | [`Figure_S14_FinalFigure.png`](Figure_S14/output/Figure_S14_FinalFigure.png) |
| S15 | Proteasome | [`map03050`](https://www.kegg.jp/pathway/map03050) | F1-PND21 and F1-PND365 | I and J | [`generate_figure_s15_ij.R`](Figure_S15/generate_figure_s15_ij.R) | [`Figure_S15_FinalFigure.png`](Figure_S15/output/Figure_S15_FinalFigure.png) |
| S16 | Penicillin and cephalosporin biosynthesis | [`map00311`](https://www.kegg.jp/pathway/map00311) | F2-PND365 | E | [`generate_figure_s16_e.R`](Figure_S16/generate_figure_s16_e.R) | [`Figure_S16_FinalFigure.png`](Figure_S16/output/Figure_S16_FinalFigure.png) |

### Sources and rendering

-   [`figure_s8_s16_common.R`](figure_s8_s16_common.R) contains the shared source checks, reference KO orders, complete-linkage row clustering, deterministic dendrogram orientation, blue-white-red coloring, and PNG renderer. The thin per-figure scripts above are independently executable.
-   The clean F1-PND21 inputs are in [`batch01_iontorrent_dams-f0_f1-pnd21/kegg_pathway_ko_z_scores/`](../results/functional_prediction/batch01_iontorrent_dams-f0_f1-pnd21/kegg_pathway_ko_z_scores/). Figure S15J and the supplied Figure S16 source are in [`batch02_illumina-miseq_f1-f2_pnd365/kegg_pathway_ko_z_scores/`](../results/functional_prediction/batch02_illumina-miseq_f1-f2_pnd365/kegg_pathway_ko_z_scores/).
-   The source tables are generated from canonical PICRUSt2 KO predictions and the canonical manifest by [`generate_kegg_pathway_ko_z_scores.R`](../scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R). Local-only regression checks compared all 19 clean tables with the preserved legacy exports during curation.
-   The legacy workflow applied `na.omit()` before rendering. The new renderer makes this behavior explicit: three undefined constant-profile rows are omitted from S9 and seven from S13, and the exact KO IDs are validated.
-   Panels A-D in S8-S14 and S16 were prepared manually at `https://www.kegg.jp/pathway/mapxxxxx`, one export for each treatment in the order `CN`, `VitA`, `MeHg`, and `VitA+MeHg`. Their pathway boxes use the matching KO z-score column and the shared blue-white-red scale. Figure S15 analogously contains four PND21 pathway views (A-D) and four PND365 views (E-H). These KEGG panels and the assembled composites are retained reference outputs, not products of the R heatmap scripts.

### Output and filename convention

Each figure uses `supplementary_figures/Figure_S#/output/`, matching the organization of the main-figure workflows. Supplied composites are named `Figure_S#_FinalFigure.png`; generated heatmaps are named `FigureS#<panel>_generated.png`. Figure S15 therefore contains two generated files, `FigureS15I_generated.png` and `FigureS15J_generated.png`.

### Run

Run any individual generator from the repository root, for example:

```bash
Rscript supplementary_figures/Figure_S8/generate_figure_s8_e.R
Rscript supplementary_figures/Figure_S15/generate_figure_s15_ij.R
```

Generate all ten heatmap panels in one R session:

``` r
source("supplementary_figures/figure_s8_s16_common.R")
generate_all_supplementary_heatmaps()
```

During curation, a local-only end-to-end validator regenerated every panel in a temporary directory and checked the per-figure scripts, source schemas, total and finite KO sets, explicit undefined-row omissions, row-wise z-score invariants, reference row orders, dendrogram compatibility, generated dimensions, generated panels, and supplied-composite dimensions. The verified coverage is ten heatmaps containing 175 finite KO rows and ten explicitly omitted undefined rows.

## Figure S17

-   Source: [`Figure_S17_SourceData.pzfx`](Figure_S17_SourceData.pzfx).
-   Generation: the full figure was generated in GraphPad Prism and exported as [`Figure_S17.png`](Figure_S17.png).
-   N values below are transcribed from non-empty Prism values and are shown as male/female:

| Cohort and endpoint | `CN` | `VitA` | `MeHg` | `VitA+MeHg` |
|---|---|---|---|---|
| F1 PND180 open-field total distance, center time, and freezing time | 8/8 | 8/7 | 7/8 | 8/8 |
| F1 PND180 plus-maze open-arm time | 7/7 | 7/7 | 7/7 | 7/7 |
| F2 PND28 open-field total distance and freezing time | 4/4 | 4/4 | 4/4 | 4/4 |
| F2 PND28 open-field center time | 4/4 | 4/4 | 4/4 | 3/4 |
| F2 PND32 plus-maze open-arm time | 4/4 | 4/4 | 4/4 | 4/4 |
| F2 PND180 plus-maze open-arm time | 4/4 | 8/7 | 8/7 | 8/8 |
| F2 PND180 open-field total distance and freezing time | 4/4 | 7/8 | 7/8 | 8/8 |
| F2 PND180 open-field center time | 4/4 | 7/8 | 7/8 | 7/8 |

## Figure S18

Figure S18 is the repository-generated 4 x 4 comparison of four mercuric-related KEGG Orthologs across F0 dams, F1-PND21, F1-PND365, and F2-PND365. KEGG Orthologs are rows and cohort datasets are columns; no panel letters are assigned. The official portable workflow is [`generate_figure_s18.R`](Figure_S18/generate_figure_s18.R).

| KO row | Annotation |
|------------------------------------|------------------------------------|
| `K00520` (`merA`) | mercuric reductase [EC:1.16.1.1] |
| `K08363` (`merT`) | mercuric ion transport protein |
| `K08364` (`merP`) | periplasmic mercuric ion binding protein |
| `K08365` (`merR`) | MerR family transcriptional regulator, mercuric resistance operon regulatory protein |

| Column | Cohort    | Batch                                 | Included samples |
|--------------|--------------|-------------------------------|-------------:|
| 1      | F0-Dams   | `batch01_iontorrent_dams-f0_f1-pnd21` |               18 |
| 2      | F1-PND21  | `batch01_iontorrent_dams-f0_f1-pnd21` |               34 |
| 3      | F1-PND365 | `batch02_illumina-miseq_f1-f2_pnd365` |               24 |
| 4      | F2-PND365 | `batch02_illumina-miseq_f1-f2_pnd365` |               24 |

The numerical sources are the canonical PICRUSt2 [`kegg_ortholog_abundance.tsv`](../results/functional_prediction/batch01_iontorrent_dams-f0_f1-pnd21/kegg_ortholog_abundance.tsv) matrix for batch 01 and the corresponding [`kegg_ortholog_abundance.tsv`](../results/functional_prediction/batch02_illumina-miseq_f1-f2_pnd365/kegg_ortholog_abundance.tsv) matrix for batch 02. Cohorts are selected through the batch-specific [`sample_manifest.csv`](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv) and [`sample_manifest.csv`](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv) files. The generator honors `include_primary_analysis`, so low-quality samples `10F` and `22F` are excluded. KO selection and annotation labels come from the preserved manual KEGG REST result [`kegg_api_find_ko_mercuric_orthologs.tsv`](../data/reference/kegg/kegg_api_find_ko_mercuric_orthologs.tsv).

Within each dataset-by-KO cell, the boxplots compare the four exposure groups after pooling sexes; sex remains visible through point shape. The script fits a one-way treatment ANOVA and performs all six Tukey HSD pairwise exposure comparisons per cell. Only adjusted comparisons below 0.05 receive brackets and asterisks (`*` p \< 0.05, `**` p \< 0.01, `***` p \< 0.001, and `****` p \< 0.0001). The current data contain no Tukey-adjusted comparison below 0.05 among the 96 evaluated comparisons.

The official 2700 x 2400 PNG is [`FigureS18_generated.png`](Figure_S18/output/FigureS18_generated.png). Generate it from the repository root with:

```bash
Rscript supplementary_figures/Figure_S18/generate_figure_s18.R
```

The script validates source schemas, matrix dimensions, cohort membership, sample counts, required sex-by-treatment replication, KO presence in both matrices, finite non-negative abundances, the expected 96 post-hoc comparisons, and output dimensions. It reports the 16 omnibus ANOVA values and any adjusted significant comparisons when run directly.

## References

- <a id="ref-adobe-photoshop"></a>Adobe Inc. (2026). *Adobe Photoshop* [Computer software]. [https://www.adobe.com/products/photoshop.html](https://www.adobe.com/products/photoshop.html)
- <a id="ref-graphpad-prism"></a>GraphPad Software. (2026). *GraphPad Prism* [Computer software; repository source files record Prism 8.0.1]. [https://www.graphpad.com/](https://www.graphpad.com/)
- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-ggpubr"></a>Kassambara, A., & Erdey, L. (2026). *ggpubr: ggplot2 based publication ready plots* [R package]. [https://doi.org/10.32614/CRAN.package.ggpubr](https://doi.org/10.32614/CRAN.package.ggpubr)
- <a id="ref-rstatix"></a>Kassambara, A. (2026). *rstatix: Pipe-friendly framework for basic statistical tests* [R package]. [https://doi.org/10.32614/CRAN.package.rstatix](https://doi.org/10.32614/CRAN.package.rstatix)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
- <a id="ref-dplyr"></a>Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2026). *dplyr: A grammar of data manipulation* (Version 1.2.1) [R package]. [https://doi.org/10.32614/CRAN.package.dplyr](https://doi.org/10.32614/CRAN.package.dplyr)
- <a id="ref-ggplot2"></a>Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd ed.). Springer. [https://doi.org/10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
