# PICRUSt2 functional-prediction workflow

This is the portable Linux-terminal reconstruction of the commands preserved in [`Script_picrust_mehg.txt`](../legacy_scripts/Script_picrust_mehg.txt). Run it from the repository root. It uses the exact ASV tables and representative-sequence files that generated the retained PICRUSt2 outputs; it does not substitute the broader canonical ASV tables.

## Historical scope and limitations

| Batch | Material | Retained PICRUSt2 input | Samples in output | Input ASVs | ASVs in `seqtab_norm` |
|---|---|---:|---:|---:|---:|
| `batch01_iontorrent_dams-f0_f1-pnd21` | F0 dams and F1-PND21 offspring | original `mestrado` run | 54 | 1,684 | 1,683 |
| `batch02_illumina-miseq_f1-f2_pnd365` | F1-PND365 and F2-PND365 offspring | final filtered `outfiles` run | 49, including undocumented `t3` | 2,739 | 2,679 |

The batch 02 directory also contained an earlier 11,704-ASV attempt that reported 3,413 poorly aligned sequences. That attempt did not generate the retained aggregate output and is not part of this workflow. The final filtered run used `exprabund_1year_outfiles.txt` and `exprasv_outfiles.fna`; its supplied normalized output contains 2,679 ASV rows, 60 fewer than the input. This corrects the overwritten legacy note that reported 58 exclusions. The supplied batch 02 output contains 49 sample columns because it includes `t3`, despite the legacy note stating 48 samples.

The exact historical PICRUSt2 and BIOM package versions were not recorded. Preserved EPA-NG logs show that PICRUSt2 ran from a Python 3.9 Conda environment, but that does not identify the package version. The archived outputs are therefore the authoritative results for figure generation. Record the complete environment specification if the pipeline is rerun.

## Prerequisites and input verification

Activate a Linux Conda environment that provides `biom` and `picrust2_pipeline.py`, then verify the archived inputs before running the time-consuming steps:

```bash
set -euo pipefail
export LC_ALL=C

repo_root="$(git rev-parse --show-toplevel)"
work_root="${TMPDIR:-/tmp}/mehg_vita_picrust2"
mkdir -p "$work_root"

conda activate picrust2
command -v biom
command -v picrust2_pipeline.py
conda list --explicit > "$work_root/conda_environment_explicit.txt"

printf '%s  %s\n' \
  '4afd72de23737678e6638d3820122b16a067f549902ed03693ffbb9a33ce5de2' \
  "$repo_root/data/processed/batch01_iontorrent_dams-f0_f1-pnd21/picrust2_input_asv_counts.tsv" \
  | sha256sum --check -
printf '%s  %s\n' \
  'd7c665fab8df749858bd65318bf162306d1c4a4426fb2f2bac566ed8dffe02df' \
  "$repo_root/data/processed/batch01_iontorrent_dams-f0_f1-pnd21/picrust2_input_representative_sequences.fasta" \
  | sha256sum --check -
printf '%s  %s\n' \
  '8a5c7d3182276186c2772c0e793e1d78a89644e2efb23e37ff99857e4e9052e6' \
  "$repo_root/data/processed/batch02_illumina-miseq_f1-f2_pnd365/picrust2_input_asv_counts.tsv" \
  | sha256sum --check -
printf '%s  %s\n' \
  'd4be00bb3a0225e386b0cb85a6801ae569d6ce38817fcb5c3ba9cf7ecd8d909b' \
  "$repo_root/data/processed/batch02_illumina-miseq_f1-f2_pnd365/picrust2_input_representative_sequences.fasta" \
  | sha256sum --check -
```

## Reconstructed pipeline

The function below keeps all expensive pipeline products under `work_root`, refuses to reuse a pre-existing batch work directory, and creates harmonized, analysis-ready exports without overwriting the retained repository results.

```bash
run_picrust2_batch() {
  local batch_id="$1"
  local asv_table="$2"
  local sequences="$3"
  local batch_work="$work_root/$batch_id"
  local feature_table="$batch_work/feature_table.biom"
  local pipeline_out="$batch_work/picrust2_pipeline"
  local harmonized_out="$batch_work/harmonized"

  test -s "$asv_table"
  test -s "$sequences"
  if [[ -e "$batch_work" ]]; then
    printf 'Refusing to overwrite existing batch work directory: %s\n' "$batch_work" >&2
    return 1
  fi

  mkdir -p "$batch_work" "$harmonized_out"

  biom convert \
    --input-fp "$asv_table" \
    --output-fp "$feature_table" \
    --to-hdf5 \
    --table-type 'OTU table'

  picrust2_pipeline.py \
    --study_fasta "$sequences" \
    --input "$feature_table" \
    --output "$pipeline_out" \
    --processes 1 \
    --verbose

  gzip -cd "$pipeline_out/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    | awk 'BEGIN { FS=OFS="\t" } NR==1 { $1="ID" } { print }' \
    > "$harmonized_out/enzyme_commission_abundance.tsv"
  gzip -cd "$pipeline_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
    | awk 'BEGIN { FS=OFS="\t" } NR==1 { $1="ID" } { print }' \
    > "$harmonized_out/kegg_ortholog_abundance.tsv"
  gzip -cd "$pipeline_out/pathways_out/path_abun_unstrat.tsv.gz" \
    > "$harmonized_out/metacyc_pathway_abundance.tsv"
  gzip -cd "$pipeline_out/KO_metagenome_out/weighted_nsti.tsv.gz" \
    > "$harmonized_out/weighted_nsti.tsv"
}

batch01='batch01_iontorrent_dams-f0_f1-pnd21'
run_picrust2_batch \
  "$batch01" \
  "$repo_root/data/processed/$batch01/picrust2_input_asv_counts.tsv" \
  "$repo_root/data/processed/$batch01/picrust2_input_representative_sequences.fasta"

batch02='batch02_illumina-miseq_f1-f2_pnd365'
run_picrust2_batch \
  "$batch02" \
  "$repo_root/data/processed/$batch02/picrust2_input_asv_counts.tsv" \
  "$repo_root/data/processed/$batch02/picrust2_input_representative_sequences.fasta"
```

The long options above are equivalent to the legacy short options `-s`, `-i`, `-o`, and `-p`. A full rerun is expected to take hours; the recorded elapsed times were 6,581.13 seconds for batch 01 and 7,987.52 seconds for the final batch 02 run.

## Compare a rerun with the retained results

PICRUSt2 and reference-database version changes can alter predictions. Compare a completed rerun with the archived result matrices before using it in a figure:

```bash
for batch_id in "$batch01" "$batch02"; do
  for result_file in \
    enzyme_commission_abundance.tsv \
    kegg_ortholog_abundance.tsv \
    metacyc_pathway_abundance.tsv \
    weighted_nsti.tsv
  do
    cmp \
      "$work_root/$batch_id/harmonized/$result_file" \
      "$repo_root/results/functional_prediction/$batch_id/$result_file"
  done
done
```

No full rerun was performed during repository curation. Instead, all six retained EC, KO, and MetaCyc matrices were compared line-by-line with the supplied compressed pipeline outputs and matched exactly after changing only the EC/KO first-column header from `function` to `ID`. Both retained `weighted_nsti.tsv` files are direct decompressions of the corresponding pipeline outputs.

## Downstream sample selection

PICRUSt2 was run on all supplied columns. Figure-generation scripts must join samples to [`data/metadata/sample_manifest_all_batches.csv`](../../data/metadata/sample_manifest_all_batches.csv), honor `include_primary_analysis`, exclude batch 01 samples `10F` and `22F`, and exclude batch 02 column `t3` unless valid metadata are supplied.
