# PICRUSt2 functional predictions

These directories contain the analysis-ready, unstratified PICRUSt2 predictions used for downstream figure generation. Files are separated by sequencing batch and named by their contents.

| File | PICRUSt2 pipeline source | Contents |
|---|---|---|
| `enzyme_commission_abundance.tsv` | `EC_metagenome_out/pred_metagenome_unstrat.tsv.gz` | Enzyme Commission abundance by sample |
| `kegg_ortholog_abundance.tsv` | `KO_metagenome_out/pred_metagenome_unstrat.tsv.gz` | KEGG Ortholog abundance by sample |
| `metacyc_pathway_abundance.tsv` | `pathways_out/path_abun_unstrat.tsv.gz` | MetaCyc pathway abundance by sample |
| `weighted_nsti.tsv` | `KO_metagenome_out/weighted_nsti.tsv.gz` | Per-sample weighted NSTI quality metric |
| `kegg_ortholog_abundance_legacy_ko2.tsv` | legacy export; exact transformation undocumented | Preserved `KO2` checkpoint, not a primary PICRUSt2 pipeline export |

The EC, KO, and MetaCyc files previously curated from `Data_metadata/picrust/` were validated against the newly supplied original PICRUSt2 pipeline directories. All six files match the compressed pipeline outputs line-for-line after the EC and KO first-column header was renamed from `function` to `ID`. The MetaCyc files required no content change. The weighted NSTI files are direct decompressions of the pipeline outputs.

## Batch-level validation

| Batch | Samples in PICRUSt2 output | EC rows | KO rows | MetaCyc rows | Downstream exclusions |
|---|---:|---:|---:|---:|---|
| `batch01_iontorrent_dams-f0_f1-pnd21` | 54 | 2,242 | 6,710 | 431 | `10F`, `22F` via `include_primary_analysis` |
| `batch02_illumina-miseq_f1-f2_pnd365` | 49 | 2,106 | 6,026 | 407 | undocumented `t3` |

The batch 02 PICRUSt2 output contains `t3`, although the canonical metadata describe only 48 samples. The column is preserved for provenance but must not be analyzed without valid metadata. Downstream scripts should join the matrix columns to [`data/metadata/sample_manifest_all_batches.csv`](../../data/metadata/sample_manifest_all_batches.csv) and retain only documented primary-analysis samples.

The exact input tables and representative sequences are stored with descriptive `picrust2_input_*` names under [`data/processed/`](../../data/processed/). See the portable [`picrust2_workflow.md`](../../scripts/04_functional_prediction/picrust2_workflow.md) and its unchanged [`legacy command log`](../../scripts/legacy_scripts/Script_picrust_mehg.txt) for command and provenance details.
