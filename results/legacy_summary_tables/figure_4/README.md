# Figure 4 legacy validation checkpoints

These unchanged `Signif_*.csv` tables preserve the significant-taxon lists and treatment means used by the legacy Figure 4 heatmap workflow. They are validation checkpoints, not primary analysis inputs. The portable generators recompute Kruskal-Wallis p-values, treatment means, and z-scores from the canonical abundance tables and sample manifests.

| Checkpoint | Figure 4 subset |
|---|---|
| [`Signif_dams.csv`](Signif_dams.csv) | Panel A, F0 dams |
| [`Signif_offs_F1Female.csv`](Signif_offs_F1Female.csv) | Panel B, F1 PND21 females |
| [`Signif_offs_F1Male.csv`](Signif_offs_F1Male.csv) | Panel B, F1 PND21 males |
| [`Signif_1year_F1Female.csv`](Signif_1year_F1Female.csv) | Panel C, F1 PND365 females |
| [`Signif_1year_F1Male.csv`](Signif_1year_F1Male.csv) | Panel C, F1 PND365 males |
| [`Signif_1year_F2Female.csv`](Signif_1year_F2Female.csv) | Panel D, F2 PND365 females |
| [`Signif_1year_F2Male.csv`](Signif_1year_F2Male.csv) | Panel D, F2 PND365 males |

The checkpoints are consumed by [`validate_figure_4.R`](../../../figures/Figure_4/validate_figure_4.R). Canonical runtime inputs are documented in the main [`figures/README.md`](../../../figures/README.md).
