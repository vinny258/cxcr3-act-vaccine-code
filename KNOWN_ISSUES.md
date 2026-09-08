# Known issues and reproducibility notes

Checked 7 September 2026 by re-running each script from a clean copy and comparing the output
against the archived figure. This records honestly what does and does not reproduce.

## Reproduces exactly

| Panel | Result |
|---|---|
| Fig. 3C | 100% pixel-identical. The DE tables are byte-identical to the archived ones |
| fig. S4F | 100% pixel-identical, pathway list identical |
| fig. S8A | Reproduces via the two-step route below |
| Fig. 4H | Day 2 and Day 3 volcanoes 100% identical; Day 1 98.8%, combined 95.8% |

## Reproduces with cosmetic differences

**Fig. 4I, Fig. 4J** (97.0% and 92.4% identical). The scripts originally called
`windowsFonts(Arial = windowsFont("Arial"))`, which exists only on Windows. That is now guarded
so the scripts run anywhere, using Arial when available and the default family otherwise. On a
machine without Arial the glyph widths differ slightly. Data, layout and labels are unchanged.

## Does not reproduce exactly

**fig. S8C** (95.7% identical). Two of the twenty pathways shown differ, both at the rank-10
cutoff: the regenerated plot gains `coagulation` and `generation of second messenger molecules`
and loses `il2 stat5 signaling` and `extracellular matrix organization`. The adjusted P value
scale also shifts about tenfold, which indicates a different number of gene sets was tested.

The cause is the GSEA implementation, not the gene sets. fig. S4F uses `fgsea()` directly with
`msigdbr`'s current API and reproduces bit-for-bit; fig. S8C uses `clusterProfiler::GSEA()` and
the deprecated `msigdbr(category=, subcategory=)` arguments. **Pin clusterProfiler and
enrichplot** and update that call to `collection=` / `subcollection=`.

## Needs a large input that is not in this repository

| Panel | Needs |
|---|---|
| Fig. 3C | `shared data/wt_D0_D1.rds`, the annotated lymph node Seurat object (1.6 GB) |
| fig. S8A step 1 | The three GSE307143 matrix, features and barcode files |
| fig. S8B | `AllTumors_seurat.rds` |
| Fig. 3D, Fig. 3H, fig. S5A-C, fig. S10B-C | Saved `.rds` objects from upstream steps |

## fig. S8A: why it is split in two

`FigS8A_01_build_annotate.R` clusters and annotates; `FigS8A_02_plot_violin.R` draws the panel
from the saved object.

The annotation maps **cluster numbers** to cell type names. Cluster numbering is not stable
across Seurat versions: re-running under Seurat 5.5.0 produced 16 clusters where the original
produced 17, so the label for cluster 16 (Neutrophils) was never applied and that population
silently disappeared from the figure. Step 1 now stops with an explicit error if the map does
not apply completely. Step 2 is deterministic and reproduces the published panel from the
archived annotated object.

Step 2 also pins three things that were originally done by hand in Illustrator: the row order,
the relabelling of "Proliferating cells" to "Proliferating T/NK cells", and the restriction to
Cxcl10 and Cxcl9 (the original call included `Il21r`).

## Package version drift

The methods state R 4.3.3 to 4.5.3, Seurat 5.1.0 to 5.4.0, DESeq2 1.50.2, fgsea 1.36.2. The
checking run used R 4.6.0, Seurat 5.5.0, DESeq2 1.52.0, fgsea 1.38.0. Most panels are
insensitive to this. The exceptions are fig. S8C above, and any script that depends on cluster
numbering.

## Output filenames

Four scripts write a filename that differs from the one filed in the lab archive, because the
file was renamed when archived:

| Script writes | Archived as |
|---|---|
| `Fig3C_pseudobulk_volcano_d1_vs_d0.pdf` | `Fig3C_pseudobulk.pdf` |
| `cxcl10_boxplot.pdf` | `Fig4I_cxcl10_boxplot.pdf` |
| `FigS8B_GSE121861_B16F10_violin.pdf` | `FigS8B_GSE121861_violin.pdf` |
| `GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome_dotplot.pdf` | `FigS8C_GSEA_ChAd_d1_vs_PBS_d1.pdf` |

The contents are the same. `helpers/script_panel_map.json` gives the authoritative mapping.
