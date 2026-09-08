# Bundled intermediate objects

Small R objects needed to redraw figure panels without re-running the upstream analysis.
They are included because re-deriving them is not reliable: cluster numbering is not stable
across Seurat versions, and a rebuild of fig. S8A silently lost a cell population. See
`../KNOWN_ISSUES.md`.

| File | Size | Used by |
|---|---|---|
| `UMAP_Exp333_Tet_combined_df.rds` | 428 KB | Fig. 3H, Fig. 3I, fig. S5A, fig. S5B, fig. S5C |
| `FigS10B_umap_data.rds` | 2.3 MB | fig. S10B |
| `FigS10C_umap_data.rds` | 2.3 MB | fig. S10C |
| `de_merged.rds`, `dot_df.rds`, `pct_long.rds` | 4 KB each | Fig. 3D |

Larger objects are not bundled here:

| Object | Size | Where to get it |
|---|---|---|
| `wt_D0_D1.rds` | 1.6 GB | Deposited with the lymph node scRNA-seq in GEO, accession `[GSE PENDING]` |
| `B16F10_GSE307143_annotated.rds` | 329 MB | Deposited in this Zenodo record |
| `AllTumors_seurat.rds` | 164 MB | Deposited in this Zenodo record |

Point `FIG_DATADIR` at the folder holding whichever of these a script needs.
