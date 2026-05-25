# Visium HD Analysis System — Agent Guide

## Project Overview

A single-file R Shiny application for end-to-end analysis of 10x Genomics Visium HD spatial transcriptomics data. Built on Seurat v5 with support for multi-sample integration, ROI-based analysis, cell type annotation, differential expression, enrichment, trajectory inference, and cell-cell communication.

**Location:** `/home/winston/Documents/visiumHD_ready/`

**Files:**
| File | Lines | Purpose |
|---|---|---|
| `visium_HD_analysis.R` | ~1555 | Complete Shiny app (UI + Server in one file) |
| `USER_MANUAL.md` | ~436 | Traditional Chinese user manual (14 sections) |

## Architecture

Monolithic Shiny app using `shinydashboard`. No separate `ui.R`/`server.R` — everything lives in `visium_HD_analysis.R`.

**Structure within visium_HD_analysis.R:**
- Lines 1-91: Package loading, constants, utility functions
- Lines 92-405: UI definition (dashboardPage with 9 tabs)
- Lines 407-1555: Server logic (reactive values, event handlers, render functions)

**State Management:**
All app state is held in a single `reactiveValues()` object (`rv`):
```
rv$datasets        # Named list of dataset objects (seurat, organism, bin_size, source)
rv$active_dataset  # Name string of currently active dataset
rv$organism        # "human" or "mouse"
rv$roi_a_cells     # Cell IDs for ROI-A
rv$roi_b_cells     # Cell IDs for ROI-B
rv$de_results      # Data frame of DE results
rv$enrich_results  # clusterProfiler enrichment result object
rv$ref_seurat      # scRNA-seq reference Seurat object (for RCTD)
rv$slingshot_res   # Slingshot trajectory result (SDS object)
rv$cellchat_obj    # CellChat object
```

Helper `active_obj()` reactive returns the Seurat object for the active dataset.

## 9 Analysis Tabs

### Tab 1: Data Loading (`tab_load`)
- **Space Ranger load:** Browse or type path → auto-detects `spatial/`, `binned_outputs/`, `.h5` files
- **RDS load:** Upload pre-saved Seurat objects
- **RCTD reference:** Upload scRNA-seq reference RDS for deconvolution
- **Multi-dataset table:** Set active, remove datasets
- **Bin sizes:** 2µm, 8µm (default), 16µm, or custom
- **Key functions:** `auto_detect_dir()`, `Load10X_Spatial()`

### Tab 2: QC & Preprocessing (`tab_qc`)
- **QC filtering:** min nCount (default 100), min nFeature (default 50), max percent.mt (default 25%)
- **Normalization:** LogNormalize (default) or SCTransform
- **Sketch workflow:** Auto-enabled for >100K cells; LeverageScore subsampling to 50K cells
- **Pipeline:** Normalize → FindVariableFeatures → ScaleData → PCA → FindNeighbors → FindClusters → UMAP
- **Multi-dataset integration:** Merge → Normalize → PCA → IntegrateLayers (Harmony/CCA/RPCA) → Cluster → UMAP
- **QC plots:** Violin plot (nCount, nFeature, percent.mt) + scatter plots

### Tab 3: ROI Selection & Dim Reduction (`tab_roi`)
- **Spatial lasso:** Plotly interactive plot with lasso selection tool
- **ROI-A / ROI-B:** Save lasso selections as named cell sets
- **Dim reduction plots:** UMAP, t-SNE, PCA — colorable by any metadata column
- **Performance:** Auto-downsamples to 50K points for spatial plot when >50K cells
- **Coordinate detection:** Tries `GetTissueCoordinates()`, falls back to UMAP embeddings

### Tab 4: Spatial Projection & Cell Type Annotation (`tab_spatial`)
- **SpatialDimPlot:** Project clusters/cell types onto tissue image
- **Gene expression maps:** SpatialFeaturePlot for individual genes
- **Three annotation methods:**
  1. **Manual:** Select cluster → assign cell type name
  2. **Marker-based:** Input marker genes → AddModuleScore → threshold → assign cell type
  3. **RCTD:** spacexr reference-based deconvolution → writes `first_type` and `cell_type` columns

### Tab 5: Differential Expression (`tab_de`)
- **Four comparison modes:**
  1. `roi` — All ROI-A cells vs all ROI-B cells
  2. `roi_celltype` — Same cell type in ROI-A vs ROI-B
  3. `group` — Arbitrary metadata column + two groups
  4. `all` — FindAllMarkers for every cluster
- **Parameters:** min log2FC (default 0.25), heatmap top N (default 30)
- **Output:** DT table (p_val, avg_log2FC, pct.1, pct.2, p_val_adj) + heatmap (DoHeatmap, downsampled to 500 cells)

### Tab 6: Functional Enrichment (`tab_enrich`)
- **Types:** GO BP, GO MF, GO CC, KEGG, GSEA (GO BP)
- **Gene ID conversion:** SYMBOL → ENTREZID via org.Hs.eg.db / org.Mm.eg.db
- **Output:** Dotplot (dotplot from enrichplot) + result table
- **p-value cutoff:** Default 0.05, adjustable

### Tab 7: Trajectory Inference (`tab_traj`)
- **Engine:** slingshot
- **Inputs:** Cluster column, optional start cluster, PCA dimensions
- **Advanced mode:** Subset to specific ROI + cell type before running slingshot
- **Output:** Dual-panel — UMAP scatter + PCA with slingshot curves overlaid (black paths)

### Tab 8: Cell-Cell Communication (`tab_cellchat`)
- **Engine:** CellChat
- **Inputs:** Cell type column, interaction DB (Secreted Signaling / ECM-Receptor / Cell-Cell Contact), scope (all / ROI-A / ROI-B)
- **Pipeline:** createCellChat → subsetDB → identifyOverExpressedGenes → identifyOverExpressedInteractions → computeCommunProb → filterCommunication → computeCommunProbPathway → aggregateNet
- **Visualizations:** Circle plot, bubble plot, heatmap
- **Output table:** source, target, ligand, receptor, pathway, probability

### Tab 9: Report Export (`tab_report`)
- **Methods section:** Auto-generated English text suitable for paper submission
- **Includes:** All parameters, package versions, citations (Seurat, clusterProfiler, slingshot, CellChat)
- **Download:** Timestamped .txt file
- **Session info:** R version + all package versions

## Sidebar: Global Visualization Controls
- Point size (0.1-5, default 0.8)
- Opacity (0.1-1, default 0.9)
- Color palette: 高對比 (custom 40-color), Set1, Dark2, Viridis, Plasma
- Resolution selector: Appears when multiple Spatial assays exist
- Refresh button: Re-renders all plots with new settings

## Key Constants & Utilities

```r
PALETTE_CAT      # 40-color categorical palette (indigo → rose)
ORG_DB           # Species-specific config: orgDb, mt pattern, rb pattern
  human: org.Hs.eg.db, "^MT-", "^RP[SL]"
  mouse: org.Mm.eg.db, "^mt-", "^Rp[sl]"
get_palette(n)   # Returns n colors from selected palette
safe_run(expr)   # tryCatch wrapper with showNotification
auto_detect_dir() # Resolves Space Ranger path → correct outs/ directory
```

## Package Dependencies

**CRAN:** shiny, shinydashboard, shinyFiles, shinyjs, plotly, ggplot2, dplyr, tidyr, tibble, RColorBrewer, viridis, scales, pheatmap, DT, patchwork

**Bioconductor:** Seurat, SeuratObject, hdf5r, clusterProfiler, enrichplot, DOSE, org.Hs.eg.db, org.Mm.eg.db, AnnotationDbi, slingshot, SingleCellExperiment

**GitHub:** CellChat (jinworks/CellChat), spacexr (dmcable/spacexr, optional for RCTD)

**Other:** harmony

## System Requirements

| Component | Minimum | Recommended |
|---|---|---|
| R | 4.2.0 | 4.3+ |
| RAM | 16 GB | 64+ GB (Visium HD datasets often >10GB) |
| CPU | 4 cores | 8+ cores |
| Storage | SSD | — |

## Launch

```r
shiny::runApp("visium_HD_analysis.R")
# or with custom port:
shiny::runApp("visium_HD_analysis.R", port = 3838, host = "0.0.0.0")
```

## Recommended Analysis Workflow

**Single sample:**
Tab 1 (Load) → Tab 2 (QC + Normalize + Cluster) → Tab 3 (ROI Selection) → Tab 4 (Cell Annotation) → Tab 5 (DE) → Tab 6 (Enrichment) → Tab 7 (Trajectory) → Tab 8 (CellChat) → Tab 9 (Report)

**Multi-sample:**
Load + QC each sample separately → Tab 2 (Merge + Integrate with Harmony) → Continue same pipeline

## Known Patterns & Pitfalls

1. **Sketch workflow** is critical for datasets >100K cells — uses Seurat v5 ProjectData to project clusters/UMAP back to full dataset
2. **Spatial coordinate detection** has multiple fallbacks: GetTissueCoordinates → UMAP embeddings. Column name detection regex: `^(x|imagerow|imagecol)`
3. **RCTD requires spacexr** — not installed by default, must be manually added via devtools
4. **CellChat requires cell_type column** to exist before running — user must complete Tab 4 annotation first
5. **Multi-dataset integration** creates a new dataset named `integrated_Nds` and sets it as active
6. **Downsampling for visualizations:** Spatial plots cap at 50K points, heatmaps subsample to 500 cells
7. **org.Hs.eg.db / org.Mm.eg.db** organism selection affects QC gene patterns AND enrichment analysis
8. **Bin size auto-extraction:** If user browses to a `square_NNnum/` directory, bin size is auto-detected from the folder name
9. **Max request size** is set to 2GB (`shiny.maxRequestSize = 2048 * 1024^2`)
10. **All plots respond** to sidebar visualization controls via `input$refresh_plot` trigger
