# External reference data

Large third-party databases belong under `data/external/` and are ignored by Git. Only this provenance and download guide is versioned.

## SILVA 138.1 for DADA2

The 16S preprocessing workflow uses the open dataset [Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2](https://doi.org/10.5281/zenodo.4587955), version 1, published by Michael R. McLaren and Benjamin J. Callahan. The workflow requires:

| Local path below `data/external/silva/` | Size reported by Zenodo | MD5 |
|---|---:|---|
| `silva_nr99_v138.1_train_set.fa.gz` | 137.3 MB | `6b41db7139834c71171f8ce5b5918fc6` |
| `silva_species_assignment_v138.1.fa.gz` | 78.7 MB | `f21c2d97c79ff07c17949a9622371a4c` |

Download URLs:

- `https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1`
- `https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1`
- optional license copy: `https://zenodo.org/records/4587955/files/SILVA_LICENSE.txt?download=1`

The complete R download and checksum commands are in the [`16S preprocessing guide`](../../scripts/01_sequence_processing/README.md#reference-download). Do not decompress the FASTA files: DADA2 reads the `.fa.gz` files directly. If the files are kept elsewhere, set `SILVA_TAXONOMY_FASTA` and `SILVA_SPECIES_FASTA` to their paths.

Zenodo reports the SILVA 138.1 material under CC BY 4.0. Retain the license file with any redistributed local copy, and cite the Zenodo DOI, SILVA database, and DADA2 method in publications as requested by the record.

The preprocessing script refuses a local reference whose MD5 differs from the values published by Zenodo. It uses `dada2::assignTaxonomy()` with the taxonomy training FASTA and `dada2::addSpecies()` with the species-assignment FASTA, matching the historical scripts.

The Zenodo record states that this v1 reference has known label problems affecting 3 of 895 families and 59 of 3,936 genera. The affected labels are listed in the linked [`bad-taxa.csv`](https://github.com/mikemc/dada2-reference-databases/blob/main/silva-138.1/v1/bad-taxa.csv). Version 138.1 is retained because it is the study reference; do not silently substitute a newer SILVA release. Record any deliberate reference update as a new analysis version.

The reference database, formatted dataset, and classifier are cited as SILVA, the Zenodo 138.1 DADA2 dataset, and DADA2 [@silva; @silva-dada2-138-1; @dada2].

## References

- <a id="ref-dada2"></a>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- <a id="ref-silva-dada2-138-1"></a>McLaren, M. R., & Callahan, B. J. (2021). *Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2* (Version 1) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.4587955](https://doi.org/10.5281/zenodo.4587955)
- <a id="ref-silva"></a>Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. *Nucleic Acids Research, 41*(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
