# Analysis workflow accessibility text

## Alt text

Workflow for a 2 × 2 vitamin A and methylmercury study. Sequential sample preparation, 16S rRNA gene sequencing, DADA2 processing, SILVA taxonomy, and data curation lead to two parallel lanes: taxonomic ecology and PICRUSt2 predicted functional profiles. PICRUSt2 feeds four independent analyses, including KO PCA, predicted KEGG pathway testing, scripted KO heatmaps with manual map coloring and composites, and predicted mercury-related orthologs. Border patterns identify scripted, restructured, hybrid, and manual steps.

## Caption

Sequential 16S rRNA gene processing produces analysis-ready ASV data, which feed parallel taxonomic analyses and PICRUSt2 predicted functional-profile analyses.

### Sample exclusions

- `batch01_iontorrent_dams_f0_f1_pnd21`: 10F, 22F — Low sequencing quality (include_primary_analysis must be true).
- `batch02_illumina_miseq_f1_f2_pnd365`: t3 — Undocumented count-table column without valid metadata (Exclude from downstream analyses unless valid metadata are supplied).

### Scripted/manual boundaries

- **Study design and cohorts:** The study design and cohort definitions were predetermined by the research project and supplied by the author as authoritative contextual information; they are not restructured analysis outputs.
- **Sample preparation and 16S rRNA gene amplicon sequencing:** Sample preparation and sequencing were performed manually outside the maintained analysis-code workflow.
- **Analysis-ready data:** Canonical metadata and documented sample-inclusion decisions are curated inputs.
- **PICRUSt2 prediction:** The terminal workflow is restructured; retained predicted profiles are the authoritative numerical products.
- **KEGG pathway interrogation:** KO z-score tables and heatmaps are scripted; KEGG pathway views are manually colored from matching treatment z-score columns and final composites are manually assembled.
- **Manual qualitative analysis:** The author visually reviewed source figures and manually assigned taxa to association boxes; Figure S6 adds no formal statistical test, pooled analysis, quantitative selection threshold, or automated rule.

### Figure notes

- Exact trimming parameters and ASV counts are documented in the sequence-processing workflow rather than displayed in the figure.
- Taxonomic manuscript mapping: Figures 3–4, 5D, S1–S5; manual taxonomic derivative: Figure S6.
- Functional manuscript mapping: Figures 6–8, S7–S16, S18.
- Figure S17 contains behavioral phenotyping and is excluded from this bioinformatics workflow.
