[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$sourceRoot = (Resolve-Path (Join-Path $repoRoot '..')).Path
$sourceData = Join-Path $sourceRoot 'Data_metadata'
$sourceImages = Join-Path $sourceRoot 'Alberto\images'

$batch01 = 'batch01_iontorrent_dams-f0_f1-pnd21'
$batch02 = 'batch02_illumina-miseq_f1-f2_pnd365'
$sraBioProject = 'PRJNA1489663'
$sraReleaseStatus = 'on hold; public upon publication of the manuscript'

$metadataRoot = Join-Path $repoRoot 'data\metadata'
$processedRoot = Join-Path $repoRoot 'data\processed'
$resultsRoot = Join-Path $repoRoot 'results'

$directories = @(
    (Join-Path $metadataRoot $batch01),
    (Join-Path $metadataRoot $batch02),
    (Join-Path $processedRoot $batch01),
    (Join-Path $processedRoot $batch02),
    (Join-Path $resultsRoot "functional_prediction\$batch01"),
    (Join-Path $resultsRoot "functional_prediction\$batch02"),
    (Join-Path $resultsRoot 'alpha_diversity'),
    (Join-Path $resultsRoot 'taxon_relative_abundance'),
    (Join-Path $resultsRoot 'legacy_summary_tables')
)

foreach ($directory in $directories) {
    New-Item -ItemType Directory -Force -Path $directory | Out-Null
}

function Copy-DataFile {
    param(
        [Parameter(Mandatory = $true)][string]$Source,
        [Parameter(Mandatory = $true)][string]$Destination
    )

    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        throw "Required source file not found: $Source"
    }

    Copy-Item -LiteralPath $Source -Destination $Destination -Force
}

# Preserve source metadata exactly, then create standardized analysis manifests below.
Copy-DataFile -Source (Join-Path $sourceData 'Metadata_mehg.csv') -Destination (Join-Path $metadataRoot "$batch01\metadata_original_semicolon.csv")
Copy-DataFile -Source (Join-Path $sourceData 'DamsF0_F1-PND21\sample-metadata.txt') -Destination (Join-Path $metadataRoot "$batch01\sample_metadata_legacy.tsv")
Copy-DataFile -Source (Join-Path $sourceData 'metadata.csv') -Destination (Join-Path $metadataRoot "$batch02\metadata_original.csv")

$processedCopies = @(
    @{ Source = 'DamsF0_F1-PND21\seqtab.txt'; Destination = "$batch01\asv_counts.tsv" },
    @{ Source = 'DamsF0_F1-PND21\tax.txt'; Destination = "$batch01\asv_taxonomy.tsv" },
    @{ Source = 'DamsF0_F1-PND21\tax_m.txt'; Destination = "$batch01\asv_taxonomy_renamed_ids.tsv" },
    @{ Source = 'DamsF0_F1-PND21\seqtab-lineage.txt'; Destination = "$batch01\asv_counts_with_taxonomy.tsv" },
    @{ Source = 'DamsF0_F1-PND21\rep-seqs.fna'; Destination = "$batch01\representative_sequences.fasta" },
    @{ Source = 'DamsF0_F1-PND21\otu_biom.biom'; Destination = "$batch01\feature_table.biom" },
    @{ Source = 'DamsF0_F1-PND21\reads_per_sample.txt'; Destination = "$batch01\read_tracking.csv" },
    @{ Source = 'F1_F2_PND365\seqtab.txt'; Destination = "$batch02\asv_counts.tsv" },
    @{ Source = 'F1_F2_PND365\relative_abundance.txt'; Destination = "$batch02\asv_relative_abundance_with_taxonomy.tsv" },
    @{ Source = 'F1_F2_PND365\tax.txt'; Destination = "$batch02\asv_taxonomy.tsv" },
    @{ Source = 'F1_F2_PND365\rep-seqs.fna'; Destination = "$batch02\representative_sequences.fasta" },
    @{ Source = 'F1_F2_PND365\otu_biom.biom'; Destination = "$batch02\feature_table.biom" }
)

foreach ($copy in $processedCopies) {
    Copy-DataFile -Source (Join-Path $sourceData $copy.Source) -Destination (Join-Path $processedRoot $copy.Destination)
}

$functionalCopies = @(
    @{ Source = 'picrust\Results_mestrado\EC_PND21_Dams.tsv'; Destination = "$batch01\enzyme_commission_abundance.tsv" },
    @{ Source = 'picrust\Results_mestrado\KO_PND21_Dams.tsv'; Destination = "$batch01\kegg_ortholog_abundance.tsv" },
    @{ Source = 'picrust\Results_mestrado\KO2_PND21_Dams.tsv'; Destination = "$batch01\kegg_ortholog_abundance_legacy_ko2.tsv" },
    @{ Source = 'picrust\Results_mestrado\Path_MetaCyc_PND21_Dams.tsv'; Destination = "$batch01\metacyc_pathway_abundance.tsv" },
    @{ Source = 'picrust\Results_1year\EC_1year_F1_F2.tsv'; Destination = "$batch02\enzyme_commission_abundance.tsv" },
    @{ Source = 'picrust\Results_1year\KO_1year_F1_F2.tsv'; Destination = "$batch02\kegg_ortholog_abundance.tsv" },
    @{ Source = 'picrust\Results_1year\KO2_1year_F1_F2.tsv'; Destination = "$batch02\kegg_ortholog_abundance_legacy_ko2.tsv" },
    @{ Source = 'picrust\Results_1year\Path_MetaCyc_1year_F1_F2.tsv'; Destination = "$batch02\metacyc_pathway_abundance.tsv" }
)

foreach ($copy in $functionalCopies) {
    Copy-DataFile -Source (Join-Path $sourceData $copy.Source) -Destination (Join-Path $resultsRoot "functional_prediction\$($copy.Destination)")
}

$summaryCopies = @(
    @{ Source = 'ZScore_Todos.csv'; Destination = 'all_batches_z_scores.csv' }
)

foreach ($copy in $summaryCopies) {
    Copy-DataFile -Source (Join-Path $sourceData $copy.Source) -Destination (Join-Path $resultsRoot "legacy_summary_tables\$($copy.Destination)")
}

$alphaDiversityCopies = @(
    @{ Source = 'Alpha_filtered_Metadata_F0F1.csv'; Destination = 'batch01_iontorrent_dams-f0_f1-pnd21_alpha_diversity_with_metadata.csv' },
    @{ Source = 'alpha_diversity.csv'; Destination = 'batch02_illumina-miseq_f1-f2_pnd365_alpha_diversity_with_metadata.csv' }
)

foreach ($copy in $alphaDiversityCopies) {
    Copy-DataFile -Source (Join-Path $sourceData $copy.Source) -Destination (Join-Path $resultsRoot "alpha_diversity\$($copy.Destination)")
}

$taxonAbundanceCopies = @(
    @{ Source = 'AbundanciasTaxons_Todos.csv'; Destination = 'batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv' },
    @{ Source = 'AbundanciasTaxons_1year.csv'; Destination = 'batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv' }
)

foreach ($copy in $taxonAbundanceCopies) {
    Copy-DataFile -Source (Join-Path $sourceImages $copy.Source) -Destination (Join-Path $resultsRoot "taxon_relative_abundance\$($copy.Destination)")
}

$expectedBatch01 = @(
    'M1F', 'M2F', 'M3F', 'M4F', 'M5F', 'M12F', 'M13F', 'M14F', 'M15F',
    'M16F', 'M6F', 'M7F', 'M8F', 'M9F', 'M18F', 'M19F', 'M21F', 'M22F',
    '25F', '26F', '27F', '28F', '30F', '37F', '38F', '39F', '40F', '41F',
    '31F', '33F', '34F', '36F', '42F', '44F', '46F', '47F', '1F', '2F', '3F',
    '5F', '6F', '13F', '14F', '15F', '16F', '17F', '10F', '11F', '7F', '8F',
    '19F', '20F', '22F', '23F'
)

$expectedBatch02 = @(
    't1', 't2', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13',
    't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24',
    't25', 't26', 't27', 't28', 't29', 't30', 't31', 't32', 't33', 't34', 't35',
    't36', 't37', 't38', 't39', 't40', 't41', 't42', 't43', 't44', 't45', 't46',
    't47', 't48', 't49'
)

$batch01SourceMetadata = Import-Csv -LiteralPath (Join-Path $sourceData 'Metadata_mehg.csv') -Delimiter ';'
$batch02SourceMetadata = Import-Csv -LiteralPath (Join-Path $sourceData 'metadata.csv')

$batch01Manifest = foreach ($sampleId in $expectedBatch01) {
    $row = $batch01SourceMetadata | Where-Object sampleID -eq $sampleId
    if (@($row).Count -ne 1) {
        throw "Batch 01 sample '$sampleId' has $(@($row).Count) metadata matches; expected exactly one."
    }

    $include = $sampleId -notin @('10F', '22F')
    [pscustomobject][ordered]@{
        sample_id = $sampleId
        batch_id = $batch01
        sequencing_platform = 'Ion Torrent'
        generation = $row.Geracao
        life_stage = if ($row.Geracao -eq 'F0') { 'dam' } else { 'offspring' }
        collection_day_pnd = if ($row.Geracao -eq 'F1') { 21 } else { '' }
        treatment = $row.Grupo
        sex = $row.Sexo
        include_primary_analysis = $include
        exclusion_reason = if ($include) { '' } else { 'low sequencing quality' }
        raw_data_repository = 'NCBI Sequence Read Archive (SRA)'
        sra_bioproject = $sraBioProject
        biosample_accession_range = 'SAMN61406324-SAMN61406377'
        raw_data_publicly_available = $false
        raw_data_release_status = $sraReleaseStatus
        source_metadata_file = 'Metadata_mehg.csv'
        legacy_observation = ''
    }
}

$batch02Manifest = foreach ($sampleId in $expectedBatch02) {
    $row = $batch02SourceMetadata | Where-Object sampleID -eq $sampleId
    if (@($row).Count -ne 1) {
        throw "Batch 02 sample '$sampleId' has $(@($row).Count) metadata matches; expected exactly one."
    }

    [pscustomobject][ordered]@{
        sample_id = $sampleId
        batch_id = $batch02
        sequencing_platform = 'Illumina MiSeq'
        generation = $row.generation
        life_stage = 'offspring'
        collection_day_pnd = 365
        treatment = $row.treatment
        sex = $row.sex
        include_primary_analysis = $true
        exclusion_reason = ''
        raw_data_repository = 'NCBI Sequence Read Archive (SRA)'
        sra_bioproject = $sraBioProject
        biosample_accession_range = 'SAMN61407672-SAMN61407719'
        raw_data_publicly_available = $false
        raw_data_release_status = $sraReleaseStatus
        source_metadata_file = 'metadata.csv'
        legacy_observation = $row.obs
    }
}

$batch01Manifest |
    Export-Csv -LiteralPath (Join-Path $metadataRoot "$batch01\sample_manifest.csv") -NoTypeInformation -Encoding UTF8
$batch02Manifest |
    Export-Csv -LiteralPath (Join-Path $metadataRoot "$batch02\sample_manifest.csv") -NoTypeInformation -Encoding UTF8
@($batch01Manifest) + @($batch02Manifest) |
    Export-Csv -LiteralPath (Join-Path $metadataRoot 'sample_manifest_all_batches.csv') -NoTypeInformation -Encoding UTF8

$sraAvailability = @(
    [pscustomobject][ordered]@{
        batch_id = $batch01
        sequencing_platform = 'Ion Torrent'
        biological_material = 'Dams_F0 and F1-PND21'
        sra_bioproject = $sraBioProject
        biosample_accession_start = 'SAMN61406324'
        biosample_accession_end = 'SAMN61406377'
        biosample_count = 54
        release_status = $sraReleaseStatus
    },
    [pscustomobject][ordered]@{
        batch_id = $batch02
        sequencing_platform = 'Illumina MiSeq'
        biological_material = 'F1-PND365 and F2-PND365'
        sra_bioproject = $sraBioProject
        biosample_accession_start = 'SAMN61407672'
        biosample_accession_end = 'SAMN61407719'
        biosample_count = 48
        release_status = $sraReleaseStatus
    }
)

$sraAvailability |
    Export-Csv -LiteralPath (Join-Path $metadataRoot 'sra_data_availability.csv') -NoTypeInformation -Encoding UTF8

function Get-TabHeaderSampleIds {
    param(
        [Parameter(Mandatory = $true)][string]$Path,
        [Parameter(Mandatory = $true)][int]$FirstSampleColumn
    )

    $header = Get-Content -LiteralPath $Path -TotalCount 1
    $fields = $header.Split("`t")
    return @($fields[$FirstSampleColumn..($fields.Count - 1)])
}

$batch01ProcessedIds = Get-TabHeaderSampleIds -Path (Join-Path $processedRoot "$batch01\asv_counts.tsv") -FirstSampleColumn 1
$batch02ProcessedIds = Get-TabHeaderSampleIds -Path (Join-Path $processedRoot "$batch02\asv_counts.tsv") -FirstSampleColumn 1

$missingBatch01 = @($expectedBatch01 | Where-Object { $_ -notin $batch01ProcessedIds })
$extraBatch01 = @($batch01ProcessedIds | Where-Object { $_ -notin $expectedBatch01 })
$missingBatch02 = @($expectedBatch02 | Where-Object { $_ -notin $batch02ProcessedIds })
$extraBatch02 = @($batch02ProcessedIds | Where-Object { $_ -notin $expectedBatch02 })

$validationRows = @(
    [pscustomobject][ordered]@{
        batch_id = $batch01
        expected_samples = $expectedBatch01.Count
        metadata_samples = $batch01SourceMetadata.Count
        processed_table_samples = $batch01ProcessedIds.Count
        primary_analysis_samples = @($batch01Manifest | Where-Object include_primary_analysis -eq $true).Count
        excluded_samples = '10F;22F'
        missing_from_processed_table = $missingBatch01 -join ';'
        unexpected_in_processed_table = $extraBatch01 -join ';'
        sra_bioproject = $sraBioProject
        biosample_accession_range = 'SAMN61406324-SAMN61406377'
        sra_release_status = $sraReleaseStatus
    },
    [pscustomobject][ordered]@{
        batch_id = $batch02
        expected_samples = $expectedBatch02.Count
        metadata_samples = $batch02SourceMetadata.Count
        processed_table_samples = $batch02ProcessedIds.Count
        primary_analysis_samples = $batch02Manifest.Count
        excluded_samples = ''
        missing_from_processed_table = $missingBatch02 -join ';'
        unexpected_in_processed_table = $extraBatch02 -join ';'
        sra_bioproject = $sraBioProject
        biosample_accession_range = 'SAMN61407672-SAMN61407719'
        sra_release_status = $sraReleaseStatus
    }
)

$validationRows |
    Export-Csv -LiteralPath (Join-Path $metadataRoot 'data_validation_summary.csv') -NoTypeInformation -Encoding UTF8

Write-Output "Organized data under: $repoRoot"
Write-Output "Batch 01: $($expectedBatch01.Count) samples; 2 flagged for exclusion (10F, 22F)."
Write-Output "Batch 02: $($expectedBatch02.Count) metadata samples; processed table extras: $($extraBatch02 -join ', ')."
Write-Output "Raw sequencing data: NCBI SRA BioProject $sraBioProject ($sraReleaseStatus)."
