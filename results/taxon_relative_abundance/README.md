# Taxon relative abundance

- `batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv` preserves the values from the original source file named `AbundanciasTaxons_Todos.csv`.
- `batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv` preserves the values from the original source file named `AbundanciasTaxons_1year.csv`.

Each table contains taxonomic descriptors in rows and sample-level relative abundances in columns. These files are the canonical repository locations used at runtime; figure workflows select the required taxa and join sample identifiers to the canonical manifests in [`data/metadata/`](../../data/metadata/).
