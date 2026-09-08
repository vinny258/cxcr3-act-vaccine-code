# Analysis code

Code for the figures of:

**CXCR3 on transferred T cells enables vaccine-boosted adoptive T-cell therapy in solid tumors**
Pereira Almeida et al.

## What is here

20 scripts covering 30 figure panels. Each was taken from the most
recent version in the lab's per-panel archive or working script tree, so this is the code as it
stands, not as first written.

| Script | Figure panels |
|---|---|
| `R/fig3/Fig3C_pseudobulk_DE.R` | Fig3C |
| `R/fig3/Fig3D_percluster_DE_plots.R` | Fig3D |
| `R/fig4/Fig4I_cxcl10_boxplot.R` | Fig4I, FigS9E, FigS9F, FigS9G, FigS9H, FigS9I, FigS9J, FigS9K, FigS9L |
| `R/fig4/Fig4J_cxcl9_boxplot.R` | Fig4J |
| `R/figS10/FigS10B_UMAP.R` | FigS10B |
| `R/figS10/FigS10C_UMAP.R` | FigS10C |
| `R/figS10/FigS10D_transferred_clustering.R` | FigS10D |
| `R/figS10/FigS10E_transferred_clustering.R` | FigS10E |
| `R/figS10/FigS10F_transferred_clustering.R` | FigS10F |
| `R/figS4/FigS4A_volcano_ChAd.R` | FigS4A |
| `R/figS4/FigS4F_pseudobulk_GSEA.R` | FigS4F |
| `R/figS5/UMAP_Exp333_Tet_cluster_abs_stats.R` | FigS5C |
| `R/figS5/UMAP_Exp333_Tet_cluster_stats.R` | FigS5B |
| `R/figS5/UMAP_Exp333_Tet_clusters.R` | Fig3H, Fig3I, FigS5A |
| `R/figS8/B16F10_GSE307143_Cxcl10.R` | FigS8A |
| `R/figS8/FigS8A_01_build_annotate.R` | FigS8A |
| `R/figS8/FigS8A_02_plot_violin.R` | FigS8A |
| `R/figS8/FigS8B_GSE121861_B16F10_violin.R` | FigS8B |
| `R/figS8/FigS8C_GSEA_ChAd_d1_vs_PBS_d1.R` | FigS8C |
| `R/figS9/biomarker_analysis.R` | Fig4H, FigS8D, FigS8E, FigS9E, FigS9F, FigS9G, FigS9H, FigS9I, FigS9J, FigS9K, FigS9L |

Some scripts feed several panels: one analysis produces several outputs that were placed in
different figures. That mapping is also in `helpers/script_panel_map.json`.

## What is not here

- **File-management utilities.** Around 26 Python scripts that copy files into the per-panel
  archive and update its index. They are lab housekeeping, not analysis, and produce no result
  in the paper.
- **Large inputs.** Seurat objects, FASTQs and count matrices. See "Getting the data".
- **Superseded versions.** Earlier drafts of scripts, kept in the lab archive under
  `_superseded_*` folders.

## Getting the data

| Data | Where |
|---|---|
| Bulk RNA-seq, B16F10-OVA tumours | NCBI GEO, accession `[GSE PENDING]` |
| scRNA-seq, vaccine-draining lymph node | NCBI GEO, accession `[GSE PENDING]` |
| Public reanalysis: syngeneic tumour atlas | GEO GSE307143 |
| Public reanalysis: mouse syngeneic tumour models | GEO GSE121861 |
| Public reanalysis: human ChAdOx1 blood | `[accession pending, see paper]` |

Tabulated values behind every figure panel are in data S1 of the supplementary materials.

## Running it

Each script resolves its own directory and writes outputs beside itself. Override with the
`FIG_OUTDIR` environment variable.

    Rscript R/figS8/FigS8A_02_plot_violin.R

Some scripts need a large input that is not in this repository. They stop with a clear message
saying what is missing.

**Read `KNOWN_ISSUES.md` before you conclude that something is broken.** Several panels do not
reproduce byte-for-byte under current package versions, and the reasons are documented.

## Environment

`helpers/ENVIRONMENT.md` and `helpers/package_versions.csv` record the R version and the exact
package versions used to run and check this code, alongside the versions stated in the paper's
methods. Where they differ, see `KNOWN_ISSUES.md`.

## Layout

    R/fig3/ … R/figS10/     analysis scripts, one directory per figure
    helpers/                environment record, script-to-panel map, source-data export
