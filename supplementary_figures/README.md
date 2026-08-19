# Supplementary figures

## Figure S1

Figure S1 complements the Shannon-diversity results in main Figure 3 panels
D-G by showing Observed richness and InvSimpson diversity for the same four
canonical cohorts. The portable workflow reads the batch-specific
[`alpha-diversity tables`](../results/alpha_diversity/) and joins them to the
corresponding canonical [`sample manifests`](../data/metadata/). It honors
`include_primary_analysis`, excluding low-quality samples `10F` and `22F`, and
does not analyze undocumented batch 02 sample `t3`.

| Panel | Metric | Cohort | Samples | Statistical workflow |
|---|---|---|---:|---|
| A | Observed | Dams (F0) | 18 | Within-sex one-way ANOVA and Tukey HSD |
| B | Observed | F1-PND21 | 34 | Between-sex Wilcoxon test; within-sex Dunn tests with FDR adjustment |
| C | Observed | F1-PND365 | 24 | Between-sex Wilcoxon test; within-sex Dunn tests with FDR adjustment |
| D | Observed | F2-PND365 | 24 | Between-sex Welch t-test; within-sex one-way ANOVA and Tukey HSD |
| E | InvSimpson | Dams (F0) | 18 | Within-sex one-way ANOVA and Tukey HSD |
| F | InvSimpson | F1-PND21 | 34 | Between-sex Welch t-test; within-sex one-way ANOVA and Tukey HSD |
| G | InvSimpson | F1-PND365 | 24 | Between-sex Welch t-test; within-sex one-way ANOVA and Tukey HSD |
| H | InvSimpson | F2-PND365 | 24 | Between-sex Welch t-test; within-sex one-way ANOVA and Tukey HSD |

[`generate_figure_s1_ah.R`](Figure_S1/generate_figure_s1_ah.R) implements the
same treatment order, colors, sex facets, test labels, and PNG dimensions used
for Figure 3 panels D-G. It generates
`FigureS1A_generated.png`-`FigureS1H_generated.png` in
[`Figure_S1/output/`](Figure_S1/output/). The matching
[`validate_figure_s1_ah.R`](Figure_S1/validate_figure_s1_ah.R) checks the
canonical sources, manifest join, exact cohort and exclusion sets, panel
mapping, statistical-test coverage, PNG names and dimensions, committed panel
outputs, and supplied composite.

The supplied [`Figure_S1_FinalFigure.png`](Figure_S1/output/Figure_S1_FinalFigure.png)
was assembled manually in Adobe Photoshop from the individual panels and is
retained as a visual reference; the R workflow does not recreate that
composite. The reference labels F1-PND365 as `n = 25`, reflecting the original
inclusion of `t3`. The portable panels deliberately report `n = 24` because
`t3` has no canonical metadata row; this follows the repository-wide cohort
rule and makes the difference auditable.

Generate and validate the panels from the repository root:

``` powershell
Rscript supplementary_figures/Figure_S1/generate_figure_s1_ah.R
Rscript supplementary_figures/Figure_S1/validate_figure_s1_ah.R
```

## Figure S2

Figure S2 contains three individual Dams (F0) relative-abundance boxplots.
Panels are assigned from left to right in the supplied composite. The portable
workflow reads the canonical batch-01
[`taxon-relative-abundance table`](../results/taxon_relative_abundance/batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv)
and joins its sample columns to the batch-01
[`sample manifest`](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv).
It honors `include_primary_analysis` and selects the 18 F0 dams: 5 `CN`, 4
`VitA`, 5 `MeHg`, and 4 `VitA+MeHg` samples. No Figure S2-specific source-data
copy is required because the historical `AbundanciasTaxons_Todos.csv` values
are already preserved under the repository-compliant canonical filename.

| Panel | Canonical taxon row | Display label | Significant FDR-adjusted Dunn comparisons |
|---|---|---|---|
| A | `g__[Eubacterium] ventriosum group` | `g__.Eubacterium..ventriosum.group` | `CN`-`MeHg`; `VitA`-`MeHg`; `MeHg`-`VitA+MeHg` |
| B | `g__[Ruminococcus] gauvreauii group` | `g__.Ruminococcus..gauvreauii.group` | `CN`-`VitA+MeHg`; `MeHg`-`VitA+MeHg` |
| C | `g__Oribacterium` | `g__Oribacterium` | `MeHg`-`VitA+MeHg` |

All panels use a four-group Kruskal-Wallis test followed by all six pairwise
Dunn tests with Benjamini-Hochberg FDR adjustment. Only adjusted comparisons
at or below 0.05 are annotated. This deliberately nonparametric workflow
reproduces the tests and labels in the historical source panels; it does not
choose a method from a normality test.

[`generate_figure_s2_abc.R`](Figure_S2/generate_figure_s2_abc.R) generates the
three individual 600 x 500 PNGs
[`FigureS2A_generated.png`](Figure_S2/output/FigureS2A_generated.png),
[`FigureS2B_generated.png`](Figure_S2/output/FigureS2B_generated.png), and
[`FigureS2C_generated.png`](Figure_S2/output/FigureS2C_generated.png). It does
not assemble the publication layout. The supplied
[`Figure_S2_FinalFigure.png`](Figure_S2/output/Figure_S2_FinalFigure.png) was
assembled manually in Adobe Photoshop and is retained only as a visual
reference. The three original 600 x 500 TIFF panel exports were inspected from
the archived `Alberto/images/Abundancias/NONparametric/Dams_F0_NONparametric/`
directory; the portable workflow creates PNG files only.

Generate and validate all three panels from the repository root:

``` powershell
Rscript supplementary_figures/Figure_S2/generate_figure_s2_abc.R
Rscript supplementary_figures/Figure_S2/validate_figure_s2_abc.R
```

[`validate_figure_s2_abc.R`](Figure_S2/validate_figure_s2_abc.R) checks source
schemas, taxon identity, F0 membership, treatment counts, finite relative
abundances, the exact Kruskal-Wallis results, all six Dunn comparisons per
panel, the expected FDR-adjusted significant pairs, PNG names and 600 x 500
dimensions, committed generated panels, and the supplied composite.

## Figures S8-S16

Figures S8-S16 follow the same reproducibility boundary as main Figures 7 and 8. The KEGG pathway views were prepared manually, while their KO z-score heatmaps now have portable base-R renderers. The supplied final composites are retained as visual references under each figure's `output/` directory.

| Figure | Pathway | KEGG map | Manuscript cohort | Heatmap panel(s) | Generator | Supplied composite |
|----|----|----|----|----|----|----|
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
-   The source tables are generated from canonical PICRUSt2 KO predictions and the canonical manifest by [`generate_kegg_pathway_ko_z_scores.R`](../scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R). [`validate_kegg_pathway_ko_z_scores.R`](../scripts/04_functional_prediction/validate_kegg_pathway_ko_z_scores.R) checks all 19 clean tables against the preserved legacy exports.
-   The legacy workflow applied `na.omit()` before rendering. The new renderer makes this behavior explicit: three undefined constant-profile rows are omitted from S9 and seven from S13, and the exact KO IDs are validated.
-   Panels A-D in S8-S14 and S16 were prepared manually at `https://www.kegg.jp/pathway/mapxxxxx`, one export for each treatment in the order `CN`, `VitA`, `MeHg`, and `VitA+MeHg`. Their pathway boxes use the matching KO z-score column and the shared blue-white-red scale. Figure S15 analogously contains four PND21 pathway views (A-D) and four PND365 views (E-H). These KEGG panels and the assembled composites are retained reference outputs, not products of the R heatmap scripts.

### Output and filename convention

Each figure uses `supplementary_figures/Figure_S#/output/`, matching the organization of the main-figure workflows. Supplied composites are named `Figure_S#_FinalFigure.png`; generated heatmaps are named `FigureS#<panel>_generated.png`. Figure S15 therefore contains two generated files, `FigureS15I_generated.png` and `FigureS15J_generated.png`.

### Run and validate

Run any individual generator from the repository root, for example:

``` powershell
Rscript supplementary_figures/Figure_S8/generate_figure_s8_e.R
Rscript supplementary_figures/Figure_S15/generate_figure_s15_ij.R
```

Generate all ten heatmap panels in one R session:

``` r
source("supplementary_figures/figure_s8_s16_common.R")
generate_all_supplementary_heatmaps()
```

Run the complete rendering validation:

``` powershell
Rscript supplementary_figures/validate_figure_s8_s16.R
```

The validator regenerates every panel in a temporary directory and checks all per-figure scripts, source schemas, total and finite KO sets, explicit undefined-row omissions, row-wise z-score invariants, reference row orders, dendrogram compatibility, generated dimensions, committed generated panels, and supplied-composite dimensions. Current coverage is ten heatmaps containing 175 finite KO rows and ten explicitly omitted undefined rows.

## Figure S17

-   Source: [`Figure_S17_SourceData.pzfx`](Figure_S17_SourceData.pzfx).
-   Generation: the full figure was generated in GraphPad Prism and exported as [`Figure_S17.png`](Figure_S17.png).

## Figure S18

Figure S18 is the repository-generated 4 x 4 comparison of four
mercuric-related KEGG Orthologs across F0 dams, F1-PND21, F1-PND365, and
F2-PND365. KEGG Orthologs are rows and cohort datasets are columns; no panel
letters are assigned. The
official portable workflow is
[`generate_figure_s18.R`](Figure_S18/generate_figure_s18.R).

| KO row | Annotation |
|---|---|
| `K00520` (`merA`) | mercuric reductase [EC:1.16.1.1] |
| `K08363` (`merT`) | mercuric ion transport protein |
| `K08364` (`merP`) | periplasmic mercuric ion binding protein |
| `K08365` (`merR`) | MerR family transcriptional regulator, mercuric resistance operon regulatory protein |

| Column | Cohort | Batch | Included samples |
|---|---|---|---:|
| 1 | F0-Dams | `batch01_iontorrent_dams-f0_f1-pnd21` | 18 |
| 2 | F1-PND21 | `batch01_iontorrent_dams-f0_f1-pnd21` | 34 |
| 3 | F1-PND365 | `batch02_illumina-miseq_f1-f2_pnd365` | 24 |
| 4 | F2-PND365 | `batch02_illumina-miseq_f1-f2_pnd365` | 24 |

The numerical sources are the canonical PICRUSt2
[`kegg_ortholog_abundance.tsv`](../results/functional_prediction/batch01_iontorrent_dams-f0_f1-pnd21/kegg_ortholog_abundance.tsv)
matrix for batch 01 and the corresponding
[`kegg_ortholog_abundance.tsv`](../results/functional_prediction/batch02_illumina-miseq_f1-f2_pnd365/kegg_ortholog_abundance.tsv)
matrix for batch 02. Cohorts are selected through the batch-specific
[`sample_manifest.csv`](../data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv)
and
[`sample_manifest.csv`](../data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv)
files. The generator honors `include_primary_analysis`, so low-quality samples
`10F` and `22F` are excluded. KO selection and annotation labels come from the
preserved manual KEGG REST result
[`kegg_api_find_ko_mercuric_orthologs.tsv`](../data/reference/kegg/kegg_api_find_ko_mercuric_orthologs.tsv).

Within each dataset-by-KO cell, the boxplots compare the four exposure groups
after pooling sexes; sex remains visible through point shape. The script fits a
one-way treatment ANOVA and performs all six Tukey HSD pairwise exposure
comparisons per cell. Only adjusted comparisons below 0.05 receive brackets and
asterisks (`*` p < 0.05, `**` p < 0.01, `***` p < 0.001, and `****` p <
0.0001). The current data contain no Tukey-adjusted comparison below 0.05 among
the 96 evaluated comparisons.

The official 2700 x 2400 PNG is
[`FigureS18_generated.png`](Figure_S18/output/FigureS18_generated.png).
Generate it from the repository root with:

``` powershell
Rscript supplementary_figures/Figure_S18/generate_figure_s18.R
```

The script validates source schemas, matrix dimensions, cohort membership,
sample counts, required sex-by-treatment replication, KO presence in both
matrices, finite non-negative abundances, the expected 96 post-hoc comparisons,
and output dimensions. It reports the 16 omnibus ANOVA values and any adjusted
significant comparisons when run directly.
