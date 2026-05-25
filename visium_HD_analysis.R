# =============================================================================
# Visium HD 分析系統 — 完整獨立 Shiny App
# 支援：檔案瀏覽/載入、QC、降維分群、ROI 圈選、差異分析、富集、軌跡
# =============================================================================

# --- 套件載入 -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyFiles)
  library(shinyjs)
  library(plotly)
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
  library(hdf5r)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  library(pheatmap)
  library(DT)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(slingshot)
  library(SingleCellExperiment)
  library(CellChat)
  library(harmony)
})

options(shiny.maxRequestSize = 2048 * 1024^2)

# --- 常數與工具函數 -----------------------------------------------------------
PALETTE_CAT <- c(
  "#6366f1","#8b5cf6","#a855f7","#d946ef","#ec4899",
  "#f43f5e","#ef4444","#f97316","#f59e0b","#eab308",
  "#84cc16","#22c55e","#10b981","#14b8a6","#06b6d4",
  "#0ea5e9","#3b82f6","#2563eb","#4f46e5","#7c3aed",
  "#9333ea","#c026d3","#db2777","#e11d48","#dc2626",
  "#ea580c","#d97706","#ca8a04","#65a30d","#16a34a",
  "#059669","#0d9488","#0891b2","#0284c7","#1d4ed8",
  "#4338ca","#6d28d9","#7e22ce","#a21caf","#be185d"
)

ORG_DB <- list(
  human = list(orgDb = "org.Hs.eg.db", mt = "^MT-", rb = "^RP[SL]"),
  mouse = list(orgDb = "org.Mm.eg.db", mt = "^mt-", rb = "^Rp[sl]")
)

get_palette <- function(n, style = "categorical") {
  if (style == "categorical") {
    if (n <= length(PALETTE_CAT)) PALETTE_CAT[1:n]
    else colorRampPalette(PALETTE_CAT)(n)
  } else viridis::viridis(n)
}

safe_run <- function(expr, msg, session = getDefaultReactiveDomain()) {
  tryCatch(expr, error = function(e) {
    showNotification(paste0("Error in ", msg, ": ", e$message),
                     type = "error", duration = 10)
    NULL
  })
}

# --- Auto-detect directory helper ---
auto_detect_dir <- function(data_dir) {
  if (!dir.exists(file.path(data_dir, "spatial")) &&
      !dir.exists(file.path(data_dir, "binned_outputs"))) {
    if (dir.exists(file.path(data_dir, "outs"))) {
      return(file.path(data_dir, "outs"))
    } else if (basename(data_dir) == "binned_outputs") {
      return(dirname(data_dir))
    } else if (grepl("^square_[0-9]+um$", basename(data_dir))) {
      if (basename(dirname(data_dir)) == "binned_outputs") {
        return(dirname(dirname(data_dir)))
      }
    }
  }
  data_dir
}

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Visium HD 分析系統", titleWidth = 300),

  # --- Sidebar ---------------------------------------------------------------
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "main_menu",
      menuItem("1. 載入資料", tabName = "tab_load", icon = icon("folder-open")),
      menuItem("2. QC 與前處理", tabName = "tab_qc", icon = icon("filter")),
      menuItem("3. ROI 圈選與降維", tabName = "tab_roi", icon = icon("crop")),
      menuItem("4. 空間投影與註解", tabName = "tab_spatial", icon = icon("map-marker-alt")),
      menuItem("5. 差異分析與熱圖", tabName = "tab_de", icon = icon("fire")),
      menuItem("6. 富集分析", tabName = "tab_enrich", icon = icon("flask")),
      menuItem("7. 軌跡分析", tabName = "tab_traj", icon = icon("route")),
      menuItem("8. 細胞通訊", tabName = "tab_cellchat", icon = icon("link")),
      menuItem("9. 分析報告", tabName = "tab_report", icon = icon("file-alt"))
    ),
    hr(),
    h4("🎨 視覺化控制", style = "padding-left:15px; color: white;"),
    sliderInput("pt_size", "點大小:", min = 0.1, max = 5, value = 0.8, step = 0.1),
    sliderInput("pt_alpha", "不透明度:", min = 0.1, max = 1, value = 0.9, step = 0.1),
    selectInput("color_pal", "調色盤:",
                choices = c("高對比" = "cat", "Set1" = "Set1", "Dark2" = "Dark2",
                            "Viridis" = "viridis", "Plasma" = "plasma")),
    uiOutput("resolution_selector"),
    actionButton("refresh_plot", "🔄 重繪圖形", icon = icon("sync"),
                 style = "margin:10px 15px; background:#f39c12; color:white; border:none;")
  ),

  # --- Body ------------------------------------------------------------------
  dashboardBody(
    useShinyjs(),
    tabItems(

      # == Tab 1: Data Loading ==
      tabItem(tabName = "tab_load",
        h2("資料載入與管理"),
        fluidRow(
          box(width = 6, title = tagList(icon("folder-open"), " Space Ranger 載入"),
              status = "primary", solidHeader = TRUE,
              textInput("ds_name", "Dataset 名稱", placeholder = "e.g., skin_8um"),
              radioButtons("organism", "物種", choices = c("Human" = "human", "Mouse" = "mouse"),
                           selected = "human", inline = TRUE),
              radioButtons("bin_size", "Bin Size",
                           choices = c("2µm" = "002", "8µm" = "008", "16µm" = "016", "Custom" = "custom"),
                           selected = "008", inline = TRUE),
              conditionalPanel("input.bin_size == 'custom'",
                textInput("custom_bin", "自訂 Bin Code", placeholder = "e.g., 032")),
              tags$label("Space Ranger 輸出目錄"),
              tags$div(style = "display:flex; gap:8px; margin-bottom:8px;",
                shinyDirButton("dir_choose", "Browse...",
                               title = "選擇 Space Ranger 輸出目錄",
                               icon = icon("folder-open"), class = "btn-info"),
                tags$span(style = "flex:1; font-size:12px; padding:6px 10px; background:#f0f0f0; border-radius:6px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;",
                  textOutput("sel_dir_text", inline = TRUE))
              ),
              textInput("data_dir", "或手動輸入路徑", placeholder = "/path/to/spaceranger/outs"),
              uiOutput("dir_preview"),
              radioButtons("data_fmt", "資料格式",
                           choices = c("HDF5 (.h5)" = "h5", "Matrix (.mtx)" = "mtx"),
                           selected = "h5", inline = TRUE),
              hr(),
              actionButton("btn_load", "載入 Dataset", icon = icon("upload"),
                           class = "btn-primary btn-block")
          ),
          box(width = 6, title = tagList(icon("file-import"), " RDS 載入"),
              status = "info", solidHeader = TRUE,
              fileInput("rds_file", "上傳 RDS 檔案", accept = c(".rds", ".RDS")),
              textInput("rds_name", "Dataset 名稱", placeholder = "e.g., my_analysis"),
              radioButtons("rds_org", "物種", choices = c("Human" = "human", "Mouse" = "mouse"),
                           selected = "human", inline = TRUE),
              hr(),
              actionButton("btn_load_rds", "載入 RDS", icon = icon("file-upload"),
                           class = "btn-primary btn-block"),
              hr(),
              h4("scRNA-seq Reference (RCTD 用)"),
              fileInput("ref_rds_file", "上傳 Reference RDS", accept = c(".rds", ".RDS")),
              textInput("ref_label_col", "Cell Type 欄位名稱", value = "cell_type",
                        placeholder = "e.g., subclass_label")
          )
        ),
        fluidRow(
          box(width = 12, title = tagList(icon("database"), " 已載入 Datasets"),
              DT::dataTableOutput("ds_table"),
              tags$div(style = "display:flex; gap:10px; margin-top:8px;",
                actionButton("btn_set_active", "設為作用中", icon = icon("check"), class = "btn-info"),
                actionButton("btn_remove_ds", "移除", icon = icon("trash-alt"), class = "btn-danger"))
          )
        ),
        fluidRow(
          box(width = 12, title = "作用中 Dataset 摘要", collapsible = TRUE,
              verbatimTextOutput("ds_summary"))
        )
      ),

      # == Tab 2: QC ==
      tabItem(tabName = "tab_qc",
        h2("品質控制與前處理"),
        fluidRow(
          box(width = 4, title = "QC 過濾參數", status = "warning", solidHeader = TRUE,
              sliderInput("qc_ncount", "最低 nCount:", min = 0, max = 5000, value = 100, step = 50),
              sliderInput("qc_nfeat", "最低 nFeature:", min = 0, max = 2000, value = 50, step = 10),
              sliderInput("qc_mt", "最高 percent.mt (%):", min = 0, max = 100, value = 25, step = 1),
              actionButton("btn_filter", "套用 QC 過濾", icon = icon("filter"), class = "btn-warning btn-block"),
              hr(),
              h4("前處理 Pipeline"),
              radioButtons("norm_method", "Normalization 方式:",
                           choices = c("LogNormalize" = "log", "SCTransform" = "sct"),
                           selected = "log", inline = TRUE),
              checkboxInput("use_sketch", "大資料集使用 Sketch Workflow", value = TRUE),
              numericInput("sketch_n", "Sketch cells:", value = 50000, min = 10000, step = 10000),
              numericInput("cluster_res", "Clustering Resolution:", value = 0.8, min = 0.1, max = 5, step = 0.1),
              numericInput("n_pcs", "PCA dims:", value = 30, min = 5, max = 100, step = 5),
              actionButton("btn_process", "執行 Normalize + Cluster + UMAP",
                           icon = icon("play"), class = "btn-danger btn-block"),
              hr(),
              h4("多資料集整合 (Batch Correction)"),
              helpText("需先載入 ≥2 個 datasets 並各自完成 QC 過濾。"),
              radioButtons("integ_method", "Integration 方法:",
                           choices = c("Harmony" = "harmony", "CCA" = "cca", "RPCA" = "rpca"),
                           selected = "harmony", inline = TRUE),
              actionButton("btn_integrate", "合併 + 整合所有 Datasets",
                           icon = icon("object-group"), class = "btn-success btn-block")
          ),
          box(width = 8, title = "QC 視覺化", status = "info",
              plotOutput("qc_vln", height = "300px"),
              plotOutput("qc_scatter", height = "300px")
          )
        )
      ),

      # == Tab 3: ROI & Dim Reduction ==
      tabItem(tabName = "tab_roi",
        h2("ROI 圈選與降維分析"),
        fluidRow(
          box(width = 4, title = "ROI 操作", status = "success", solidHeader = TRUE,
              helpText("使用 Plotly 的 Lasso 工具在空間圖上圈選 ROI。"),
              actionButton("btn_set_roi_a", "儲存為 ROI-A", icon = icon("save"), class = "btn-success btn-block"),
              actionButton("btn_set_roi_b", "儲存為 ROI-B", icon = icon("save"), class = "btn-info btn-block"),
              hr(),
              textOutput("roi_a_info"),
              textOutput("roi_b_info"),
              hr(),
              selectInput("dim_color", "降維圖上色依據:",
                          choices = c("seurat_clusters"), selected = "seurat_clusters"),
              actionButton("btn_update_dim_choices", "更新欄位選單", icon = icon("sync"))
          ),
          box(width = 8, title = "空間座標 (Lasso 圈選)", status = "primary",
              plotlyOutput("spatial_lasso", height = "500px"))
        ),
        fluidRow(
          tabBox(width = 12, title = "降維圖",
            tabPanel("UMAP", plotOutput("plot_umap", height = "550px")),
            tabPanel("t-SNE", plotOutput("plot_tsne", height = "550px")),
            tabPanel("PCA", plotOutput("plot_pca", height = "550px"))
          )
        )
      ),

      # == Tab 4: Spatial & Annotation ==
      tabItem(tabName = "tab_spatial",
        h2("空間投影與細胞型態註解"),
        fluidRow(
          box(width = 4, title = "註解控制", status = "primary", solidHeader = TRUE,
              selectInput("annot_col", "Cell Type 欄位:", choices = NULL),
              actionButton("btn_refresh_annot", "更新欄位", icon = icon("sync")),
              hr(),
              h4("RCTD Deconvolution"),
              helpText("需先在 Tab 1 上傳 scRNA-seq Reference RDS。"),
              actionButton("btn_run_rctd", "執行 RCTD", icon = icon("dna"), class = "btn-primary btn-block"),
              hr(),
              h4("手動標註"),
              selectInput("manual_cluster", "選擇 Cluster:", choices = NULL),
              textInput("manual_label", "標註名稱", placeholder = "e.g., Mast cell"),
              actionButton("btn_manual_annot", "套用標註", icon = icon("tag"), class = "btn-success"),
              hr(),
              h4("以 Markers 定義細胞群"),
              textAreaInput("marker_genes", "輸入 Marker Genes (逗號分隔):",
                            placeholder = "e.g., TPSAB1, CPA3, KIT", rows = 3),
              textInput("marker_cell_name", "細胞群名稱:", placeholder = "e.g., Mast cell"),
              sliderInput("marker_threshold", "Module Score 閾值:", min = -1, max = 2, value = 0, step = 0.1),
              actionButton("btn_marker_annot", "依 Markers 定義細胞", icon = icon("bullseye"), class = "btn-info btn-block")
          ),
          box(width = 8, title = "空間投影", status = "info",
              plotOutput("plot_spatial_proj", height = "600px"),
              hr(),
              textInput("sp_gene", "查看基因空間表現:", placeholder = "e.g., KRT15"),
              plotOutput("plot_spatial_gene", height = "400px"))
        )
      ),

      # == Tab 5: DE & Heatmap ==
      tabItem(tabName = "tab_de",
        h2("差異基因表達分析"),
        fluidRow(
          box(width = 4, title = "DE 參數", status = "danger", solidHeader = TRUE,
              selectInput("de_mode", "比較模式:",
                          choices = c("ROI-A vs ROI-B" = "roi",
                                      "跨 ROI 同一細胞群比較" = "roi_celltype",
                                      "指定群組比較" = "group",
                                      "Find All Markers" = "all")),
              conditionalPanel("input.de_mode == 'roi_celltype'",
                selectInput("de_celltype", "選擇細胞群:", choices = NULL)),
              conditionalPanel("input.de_mode == 'group'",
                selectInput("de_ident_col", "Identity 欄位:", choices = NULL),
                selectInput("de_group1", "Group 1:", choices = NULL),
                selectInput("de_group2", "Group 2:", choices = NULL)
              ),
              sliderInput("de_logfc", "min log2FC:", min = 0, max = 2, value = 0.25, step = 0.05),
              sliderInput("de_top_n", "Heatmap Top N:", min = 5, max = 100, value = 30, step = 5),
              actionButton("btn_run_de", "執行差異分析", icon = icon("play"), class = "btn-danger btn-block")
          ),
          box(width = 8, title = "差異分析結果",
              DT::dataTableOutput("de_table"),
              hr(),
              plotOutput("de_heatmap", height = "600px"))
        )
      ),

      # == Tab 6: Enrichment ==
      tabItem(tabName = "tab_enrich",
        h2("功能富集分析"),
        fluidRow(
          box(width = 4, title = "富集參數", status = "primary", solidHeader = TRUE,
              selectInput("enrich_type", "分析類型:",
                          choices = c("GO (BP)" = "go_bp", "GO (MF)" = "go_mf", "GO (CC)" = "go_cc",
                                      "KEGG" = "kegg", "GSEA (GO)" = "gsea_go")),
              sliderInput("enrich_pval", "p-value cutoff:", min = 0.001, max = 0.1, value = 0.05, step = 0.005),
              numericInput("enrich_show", "顯示前 N 條:", value = 20, min = 5, max = 50),
              actionButton("btn_run_enrich", "執行富集分析", icon = icon("flask"), class = "btn-primary btn-block")
          ),
          box(width = 8, title = "富集結果",
              plotOutput("enrich_plot", height = "500px"),
              hr(),
              DT::dataTableOutput("enrich_table"))
        )
      ),

      # == Tab 7: Trajectory ==
      tabItem(tabName = "tab_traj",
        h2("軌跡推斷分析 (Slingshot)"),
        fluidRow(
          box(width = 4, title = "軌跡參數", status = "info", solidHeader = TRUE,
              selectInput("traj_cluster_col", "Cluster 欄位:", choices = NULL),
              selectInput("traj_start", "起始 Cluster (optional):", choices = c("Auto" = "")),
              numericInput("traj_dims", "使用 PCA dims:", value = 10, min = 2, max = 50),
              checkboxInput("traj_roi_subset", "僅分析特定 ROI 的特定細胞群", value = FALSE),
              conditionalPanel("input.traj_roi_subset",
                selectInput("traj_roi_choice", "選擇 ROI:", choices = c("ROI-A" = "a", "ROI-B" = "b", "兩者合併" = "both")),
                selectInput("traj_celltype", "篩選細胞群:", choices = NULL)
              ),
              actionButton("btn_run_traj", "執行軌跡分析", icon = icon("route"), class = "btn-info btn-block")
          ),
          box(width = 8, title = "軌跡投影",
              plotOutput("traj_plot", height = "600px"))
        )
      ),

      # == Tab 8: CellChat ==
      tabItem(tabName = "tab_cellchat",
        h2("細胞通訊分析 (CellChat)"),
        fluidRow(
          box(width = 4, title = "CellChat 參數", status = "warning", solidHeader = TRUE,
              selectInput("cc_celltype_col", "Cell Type 欄位:", choices = NULL),
              selectInput("cc_db", "Interaction Database:",
                          choices = c("Secreted Signaling" = "Secreted Signaling",
                                      "ECM-Receptor" = "ECM-Receptor",
                                      "Cell-Cell Contact" = "Cell-Cell Contact")),
              radioButtons("cc_scope", "分析範圍:",
                           choices = c("全部 cells" = "all", "僅 ROI-A" = "roi_a",
                                       "僅 ROI-B" = "roi_b"),
                           selected = "all", inline = TRUE),
              actionButton("btn_run_cellchat", "執行 CellChat", icon = icon("link"), class = "btn-warning btn-block"),
              hr(),
              selectInput("cc_plot_type", "視覺化類型:",
                          choices = c("Circle Plot" = "circle", "Bubble Plot" = "bubble",
                                      "Heatmap" = "heatmap")),
              numericInput("cc_top_n", "顯示前 N 條 pathway:", value = 20, min = 5, max = 50)
          ),
          box(width = 8, title = "CellChat 結果",
              plotOutput("cc_plot", height = "600px"),
              hr(),
              DT::dataTableOutput("cc_table"))
        )
      ),

      # == Tab 9: Report Export ==
      tabItem(tabName = "tab_report",
        h2("分析報告與版本資訊"),
        fluidRow(
          box(width = 12, title = tagList(icon("file-alt"), " Methods Section (投稿用)"),
              status = "primary", solidHeader = TRUE,
              helpText("點擊「產出分析報告」按鈕生成完整的 Methods 段落和軟體版本資訊，可直接用於論文投稿。"),
              actionButton("btn_gen_report", "產出分析報告", icon = icon("file-alt"),
                           class = "btn-primary"),
              downloadButton("btn_dl_methods", "下載 Methods 文字",
                             class = "btn-success", style = "margin-left:10px;"),
              hr(),
              h4("Methods (English)"),
              verbatimTextOutput("methods_text"),
              hr(),
              h4("分析參數摘要"),
              DT::dataTableOutput("params_table"),
              hr(),
              h4("軟體版本資訊 (sessionInfo)"),
              verbatimTextOutput("session_info")
          )
        )
      )
    ) # end tabItems
  ) # end dashboardBody
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # --- Shared Reactive Values -----------------------------------------------
  rv <- reactiveValues(
    datasets       = list(),
    active_dataset = NULL,
    organism       = "human",
    roi_a_cells    = NULL,
    roi_b_cells    = NULL,
    de_results     = NULL,
    enrich_results = NULL,
    ref_seurat     = NULL,
    slingshot_res  = NULL,
    cellchat_obj   = NULL
  )

  # Helper: get active Seurat object
  active_obj <- reactive({
    req(rv$active_dataset)
    rv$datasets[[rv$active_dataset]]$seurat
  })

  # --- Resolution Selector --------------------------------------------------
  output$resolution_selector <- renderUI({
    sobj <- active_obj()
    if (is.null(sobj)) return(NULL)
    assays <- names(sobj@assays)
    spatial_assays <- assays[grepl("^Spatial", assays)]
    if (length(spatial_assays) <= 1) return(NULL)
    selectInput("active_assay", "解析度 (Assay):", choices = spatial_assays,
                selected = DefaultAssay(sobj))
  })

  observeEvent(input$active_assay, {
    sobj <- active_obj()
    req(sobj, input$active_assay)
    if (input$active_assay %in% names(sobj@assays)) {
      DefaultAssay(sobj) <- input$active_assay
      rv$datasets[[rv$active_dataset]]$seurat <- sobj
    }
  })

  # --- shinyFiles: Directory Browser ----------------------------------------
  volumes <- c(Home = path.expand("~"), getVolumes()())
  shinyDirChoose(input, "dir_choose", roots = volumes, session = session,
                 restrictions = system.file(package = "base"))

  observeEvent(input$dir_choose, {
    if (!is.integer(input$dir_choose)) {
      chosen <- parseDirPath(volumes, input$dir_choose)
      if (length(chosen) > 0 && nchar(chosen) > 0) {
        updateTextInput(session, "data_dir", value = as.character(chosen))
      }
    }
  })

  output$sel_dir_text <- renderText({
    d <- trimws(input$data_dir)
    if (is.null(d) || d == "") "尚未選擇目錄" else d
  })

  # --- Directory Structure Preview ------------------------------------------
  output$dir_preview <- renderUI({
    d <- trimws(input$data_dir)
    if (is.null(d) || d == "" || !dir.exists(d)) return(NULL)
    d <- auto_detect_dir(d)
    has_sp  <- dir.exists(file.path(d, "spatial"))
    has_bin <- dir.exists(file.path(d, "binned_outputs"))
    h5s     <- list.files(d, pattern = "filtered.*\\.h5$", recursive = TRUE)
    bin_dirs <- character(0)
    if (has_bin) {
      bin_dirs <- list.dirs(file.path(d, "binned_outputs"), recursive = FALSE, full.names = FALSE)
      bin_dirs <- bin_dirs[grepl("^square_", bin_dirs)]
    }
    icon_ok <- icon("check-circle", style = "color:#22c55e; margin-right:6px;")
    icon_no <- icon("times-circle", style = "color:#ef4444; margin-right:6px;")
    items <- tagList(
      tags$div(style = "font-size:12px; margin-bottom:4px;",
               if (has_sp) icon_ok else icon_no,
               tags$span(paste0("spatial/ ", if (has_sp) "✓" else "✗"))),
      tags$div(style = "font-size:12px; margin-bottom:4px;",
               if (length(h5s) > 0) icon_ok else icon_no,
               tags$span(if (length(h5s) > 0) paste0("H5: ", paste(h5s, collapse = ", ")) else "No .h5 files")),
      tags$div(style = "font-size:12px; margin-bottom:4px;",
               if (has_bin) icon_ok else icon_no,
               tags$span(if (has_bin) "binned_outputs/ (Visium HD)" else "No binned_outputs/ (standard Visium)"))
    )
    if (length(bin_dirs) > 0) {
      labels <- gsub("square_0*(\\d+)um", "\\1 µm", bin_dirs)
      items <- tagList(items,
        tags$div(style = "font-size:12px; margin-top:6px; font-weight:600; color:#8b5cf6;",
                 "Available Bin Sizes: ", paste(labels, collapse = ", ")))
    }
    tags$div(style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px; padding:10px; margin-top:8px;",
             tags$div(style = "font-weight:600; margin-bottom:6px;", icon("search"), " 目錄結構"), items)
  })

  # --- Load Space Ranger Data -----------------------------------------------
  observeEvent(input$btn_load, {
    req(input$data_dir)
    ds_name <- trimws(input$ds_name)
    if (ds_name == "") ds_name <- paste0("dataset_", length(rv$datasets) + 1)
    bin_size <- input$bin_size
    if (bin_size == "custom") {
      bin_size <- trimws(input$custom_bin)
      if (bin_size == "") { showNotification("請輸入 bin size", type = "error"); return() }
    }
    data_dir <- auto_detect_dir(trimws(input$data_dir))
    if (!dir.exists(data_dir)) {
      showNotification(paste("目錄不存在:", data_dir), type = "error"); return()
    }
    # Extract bin from directory name if user pointed to a specific bin folder
    if (grepl("^square_[0-9]+um$", basename(trimws(input$data_dir)))) {
      extracted <- sub("^square_0*([0-9]+)um$", "\\1", basename(trimws(input$data_dir)))
      bin_size <- sprintf("%03d", as.numeric(extracted))
    }

    withProgress(message = "載入中...", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Reading Space Ranger output...")
        bin_dir <- file.path(data_dir, "binned_outputs", paste0("square_", bin_size, "um"))
        if (dir.exists(bin_dir)) {
          num_bin <- suppressWarnings(as.numeric(bin_size))
          bin_arg <- if (!is.na(num_bin)) num_bin else bin_size
          sobj <- Load10X_Spatial(data.dir = data_dir, bin.size = c(bin_arg),
                                 assay = "Spatial", slice = ds_name)
          def_a <- DefaultAssay(sobj)
          if (def_a != "Spatial") {
            args <- list(sobj); args[[def_a]] <- "Spatial"; names(args)[1] <- ""
            sobj <- do.call(RenameAssays, args); DefaultAssay(sobj) <- "Spatial"
          }
        } else {
          h5f <- list.files(data_dir, pattern = "filtered.*\\.h5$", full.names = TRUE)
          if (length(h5f) > 0 && input$data_fmt == "h5") {
            sobj <- Load10X_Spatial(data.dir = data_dir, filename = basename(h5f[1]),
                                   assay = "Spatial", slice = ds_name)
          } else {
            sobj <- Load10X_Spatial(data.dir = data_dir, assay = "Spatial", slice = ds_name)
          }
        }
        incProgress(0.5, detail = "Computing QC metrics...")
        org <- ORG_DB[[input$organism]]
        sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = org$mt)
        sobj[["percent.rb"]] <- PercentageFeatureSet(sobj, pattern = org$rb)
        sobj[["nCount_Spatial"]]   <- colSums(GetAssayData(sobj, layer = "counts"))
        sobj[["nFeature_Spatial"]] <- colSums(GetAssayData(sobj, layer = "counts") > 0)
        incProgress(0.3, detail = "Finalizing...")
        rv$datasets[[ds_name]] <- list(seurat = sobj, organism = input$organism,
                                       bin_size = bin_size, source = "spaceranger")
        rv$active_dataset <- ds_name
        rv$organism <- input$organism
        showNotification(paste0("'", ds_name, "' 載入成功! ", ncol(sobj), " cells, ", nrow(sobj), " genes"), type = "message")
        incProgress(0.1, detail = "Done!")
      }, error = function(e) {
        showNotification(paste("載入錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # --- Load RDS -------------------------------------------------------------
  observeEvent(input$btn_load_rds, {
    req(input$rds_file)
    ds_name <- trimws(input$rds_name)
    if (ds_name == "") ds_name <- tools::file_path_sans_ext(input$rds_file$name)
    withProgress(message = "載入 RDS...", value = 0, {
      tryCatch({
        incProgress(0.3)
        sobj <- readRDS(input$rds_file$datapath)
        if (!inherits(sobj, "Seurat")) stop("非 Seurat 物件")
        org <- ORG_DB[[input$rds_org]]
        if (!"percent.mt" %in% colnames(sobj@meta.data))
          sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = org$mt)
        rv$datasets[[ds_name]] <- list(seurat = sobj, organism = input$rds_org,
                                       bin_size = "unknown", source = "rds")
        rv$active_dataset <- ds_name
        rv$organism <- input$rds_org
        showNotification(paste0("RDS '", ds_name, "' 載入成功!"), type = "message")
      }, error = function(e) {
        showNotification(paste("RDS 錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # --- Load Reference for RCTD ---------------------------------------------
  observeEvent(input$ref_rds_file, {
    req(input$ref_rds_file)
    tryCatch({
      rv$ref_seurat <- readRDS(input$ref_rds_file$datapath)
      showNotification("Reference 載入成功!", type = "message")
    }, error = function(e) {
      showNotification(paste("Reference 錯誤:", e$message), type = "error")
    })
  })

  # --- Dataset Table --------------------------------------------------------
  output$ds_table <- DT::renderDataTable({
    ds_info <- lapply(names(rv$datasets), function(nm) {
      ds <- rv$datasets[[nm]]
      data.frame(Name = nm, Cells = ncol(ds$seurat), Genes = nrow(ds$seurat),
                 Organism = ds$organism, BinSize = ds$bin_size, Source = ds$source,
                 Active = ifelse(nm == rv$active_dataset, "✓", ""), stringsAsFactors = FALSE)
    })
    if (length(ds_info) == 0) return(data.frame(Name = character(), Cells = integer()))
    do.call(rbind, ds_info)
  }, selection = "single", options = list(pageLength = 5, dom = "t"), rownames = FALSE)

  observeEvent(input$btn_set_active, {
    sel <- input$ds_table_rows_selected
    if (is.null(sel)) { showNotification("請先選擇 dataset", type = "warning"); return() }
    nms <- names(rv$datasets)
    rv$active_dataset <- nms[sel]
    rv$organism <- rv$datasets[[nms[sel]]]$organism
  })

  observeEvent(input$btn_remove_ds, {
    sel <- input$ds_table_rows_selected
    if (is.null(sel)) return()
    removed <- names(rv$datasets)[sel]
    rv$datasets[[removed]] <- NULL
    if (rv$active_dataset == removed) {
      remaining <- names(rv$datasets)
      rv$active_dataset <- if (length(remaining) > 0) remaining[1] else NULL
    }
  })

  output$ds_summary <- renderPrint({
    req(rv$active_dataset, rv$datasets[[rv$active_dataset]])
    sobj <- rv$datasets[[rv$active_dataset]]$seurat
    cat("Active:", rv$active_dataset, "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n")
    print(sobj)
  })

  # --- QC Plots -------------------------------------------------------------
  output$qc_vln <- renderPlot({
    sobj <- active_obj(); req(sobj)
    feats <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), colnames(sobj@meta.data))
    if (length(feats) == 0) return(NULL)
    VlnPlot(sobj, features = feats, pt.size = 0, ncol = length(feats)) & theme_minimal()
  })

  output$qc_scatter <- renderPlot({
    sobj <- active_obj(); req(sobj)
    if (!"percent.mt" %in% colnames(sobj@meta.data)) return(NULL)
    p1 <- FeatureScatter(sobj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
    p2 <- FeatureScatter(sobj, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
    p1 + p2
  })

  # --- QC Filter ------------------------------------------------------------
  observeEvent(input$btn_filter, {
    sobj <- active_obj(); req(sobj)
    n_before <- ncol(sobj)
    sobj <- subset(sobj,
      nCount_Spatial   > input$qc_ncount &
      nFeature_Spatial > input$qc_nfeat &
      percent.mt       < input$qc_mt)
    rv$datasets[[rv$active_dataset]]$seurat <- sobj
    showNotification(paste0("QC 過濾完成: ", n_before, " → ", ncol(sobj), " cells"), type = "message")
  })

  # --- Processing Pipeline --------------------------------------------------
  observeEvent(input$btn_process, {
    sobj <- active_obj(); req(sobj)
    n_cells <- ncol(sobj)
    use_sketch <- input$use_sketch && n_cells > 100000

    withProgress(message = "前處理中...", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Normalizing...")
        if (input$norm_method == "sct") {
          sobj <- SCTransform(sobj, verbose = FALSE)
        } else {
          sobj <- NormalizeData(sobj)
          sobj <- FindVariableFeatures(sobj)
          sobj <- ScaleData(sobj)
        }

        if (use_sketch) {
          incProgress(0.2, detail = paste0("Sketching ", input$sketch_n, " cells..."))
          sobj <- SketchData(sobj, ncells = min(input$sketch_n, n_cells),
                             method = "LeverageScore", sketched.assay = "sketch")
          DefaultAssay(sobj) <- "sketch"
          sobj <- FindVariableFeatures(sobj)
          sobj <- ScaleData(sobj)
          sobj <- RunPCA(sobj, reduction.name = "pca.sketch", npcs = input$n_pcs)
          incProgress(0.2, detail = "Clustering...")
          sobj <- FindNeighbors(sobj, reduction = "pca.sketch", dims = 1:input$n_pcs)
          sobj <- FindClusters(sobj, resolution = input$cluster_res)
          sobj <- RunUMAP(sobj, reduction = "pca.sketch", reduction.name = "umap.sketch",
                          return.model = TRUE, dims = 1:input$n_pcs)
          incProgress(0.2, detail = "Projecting...")
          orig_assay <- grep("^Spatial", names(sobj@assays), value = TRUE)[1]
          sobj <- ProjectData(sobj, assay = orig_assay,
                              full.reduction = "full.pca.sketch",
                              sketched.assay = "sketch",
                              sketched.reduction = "pca.sketch",
                              umap.model = "umap.sketch", dims = 1:input$n_pcs,
                              refdata = list(seurat_clusters = "seurat_clusters"))
          DefaultAssay(sobj) <- orig_assay
        } else {
          incProgress(0.2, detail = "PCA...")
          sobj <- RunPCA(sobj, npcs = input$n_pcs)
          incProgress(0.2, detail = "Clustering...")
          sobj <- FindNeighbors(sobj, dims = 1:input$n_pcs)
          sobj <- FindClusters(sobj, resolution = input$cluster_res)
          sobj <- RunUMAP(sobj, dims = 1:input$n_pcs)
        }

        incProgress(0.1, detail = "Done!")
        rv$datasets[[rv$active_dataset]]$seurat <- sobj
        showNotification("前處理完成!", type = "message")
      }, error = function(e) {
        showNotification(paste("前處理錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # --- Multi-Dataset Integration --------------------------------------------
  observeEvent(input$btn_integrate, {
    if (length(rv$datasets) < 2) {
      showNotification("需至少 2 個 datasets 才能整合", type = "warning"); return()
    }
    withProgress(message = "整合中...", value = 0, {
      tryCatch({
        incProgress(0.1, detail = "Merging datasets...")
        obj_list <- lapply(names(rv$datasets), function(nm) {
          s <- rv$datasets[[nm]]$seurat
          s$sample_id <- nm
          s
        })
        merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(rv$datasets))

        incProgress(0.2, detail = "Normalizing...")
        if (input$norm_method == "sct") {
          merged <- SCTransform(merged, verbose = FALSE)
        } else {
          merged <- NormalizeData(merged)
          merged <- FindVariableFeatures(merged)
          merged <- ScaleData(merged)
        }

        incProgress(0.2, detail = "Running PCA...")
        merged <- RunPCA(merged, npcs = input$n_pcs)

        incProgress(0.3, detail = paste0("Integrating (", input$integ_method, ")..."))
        if (input$integ_method == "harmony") {
          merged <- IntegrateLayers(merged, method = HarmonyIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "integrated",
                                   verbose = FALSE)
        } else if (input$integ_method == "cca") {
          merged <- IntegrateLayers(merged, method = CCAIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "integrated",
                                   verbose = FALSE)
        } else {
          merged <- IntegrateLayers(merged, method = RPCAIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "integrated",
                                   verbose = FALSE)
        }
        merged <- FindNeighbors(merged, reduction = "integrated", dims = 1:input$n_pcs)
        merged <- FindClusters(merged, resolution = input$cluster_res)
        merged <- RunUMAP(merged, reduction = "integrated", dims = 1:input$n_pcs)

        incProgress(0.1, detail = "Done!")
        ds_name <- paste0("integrated_", length(rv$datasets), "ds")
        rv$datasets[[ds_name]] <- list(seurat = merged, organism = rv$organism,
                                       bin_size = "mixed", source = "integrated")
        rv$active_dataset <- ds_name
        showNotification(paste0("整合完成! '", ds_name, "' — ", ncol(merged), " cells"), type = "message")
      }, error = function(e) {
        showNotification(paste("整合錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # --- ROI Lasso Selection --------------------------------------------------
  output$spatial_lasso <- renderPlotly({
    sobj <- active_obj(); req(sobj)
    coords <- tryCatch(GetTissueCoordinates(sobj), error = function(e) NULL)
    if (is.null(coords) || nrow(coords) == 0) {
      # Fallback: use embeddings if no spatial coords
      red <- intersect(c("umap", "umap.sketch", "full.umap.sketch"), names(sobj@reductions))
      if (length(red) == 0) return(plotly_empty() %>% layout(title = "請先執行前處理"))
      emb <- Embeddings(sobj, reduction = red[1])[, 1:2]
      coords <- data.frame(x = emb[, 1], y = emb[, 2], cell = rownames(emb))
    } else {
      coords$cell <- rownames(coords)
      cnames <- colnames(coords)
      x_col <- cnames[grep("^(x|imagerow|imagecol)", cnames, ignore.case = TRUE)][1]
      y_col <- cnames[grep("^(y|imagerow|imagecol)", cnames, ignore.case = TRUE)]
      y_col <- y_col[y_col != x_col][1]
      if (is.na(x_col) || is.na(y_col)) { x_col <- cnames[1]; y_col <- cnames[2] }
      coords$x <- coords[[x_col]]; coords$y <- coords[[y_col]]
    }
    # Add cluster info if available
    if ("seurat_clusters" %in% colnames(sobj@meta.data)) {
      coords$cluster <- sobj$seurat_clusters[coords$cell]
    } else {
      coords$cluster <- "1"
    }
    # Downsample for performance
    if (nrow(coords) > 50000) {
      idx <- sample(nrow(coords), 50000)
      coords <- coords[idx, ]
    }
    plot_ly(coords, x = ~x, y = ~y, key = ~cell, color = ~cluster,
            colors = get_palette(length(unique(coords$cluster))),
            type = "scattergl", mode = "markers",
            marker = list(size = 3, opacity = 0.7),
            source = "lasso_src") %>%
      layout(dragmode = "lasso",
             yaxis = list(scaleanchor = "x", scaleratio = 1, autorange = "reversed"),
             title = "使用 Lasso 工具圈選 ROI")
  })

  observeEvent(input$btn_set_roi_a, {
    sel <- event_data("plotly_selected", source = "lasso_src")
    if (is.null(sel) || nrow(sel) == 0) {
      showNotification("請先在空間圖上圈選區域", type = "warning"); return()
    }
    rv$roi_a_cells <- sel$key
    showNotification(paste0("ROI-A 設定完成: ", length(rv$roi_a_cells), " cells"), type = "message")
  })

  observeEvent(input$btn_set_roi_b, {
    sel <- event_data("plotly_selected", source = "lasso_src")
    if (is.null(sel) || nrow(sel) == 0) {
      showNotification("請先在空間圖上圈選區域", type = "warning"); return()
    }
    rv$roi_b_cells <- sel$key
    showNotification(paste0("ROI-B 設定完成: ", length(rv$roi_b_cells), " cells"), type = "message")
  })

  output$roi_a_info <- renderText({
    if (is.null(rv$roi_a_cells)) "ROI-A: 尚未設定"
    else paste0("ROI-A: ", length(rv$roi_a_cells), " cells")
  })
  output$roi_b_info <- renderText({
    if (is.null(rv$roi_b_cells)) "ROI-B: 尚未設定"
    else paste0("ROI-B: ", length(rv$roi_b_cells), " cells")
  })

  # --- Update dropdown choices from metadata --------------------------------
  observeEvent(input$btn_update_dim_choices, {
    sobj <- active_obj(); req(sobj)
    cols <- colnames(sobj@meta.data)
    updateSelectInput(session, "dim_color", choices = cols,
                      selected = if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1])
  })

  # --- Dim Reduction Plots --------------------------------------------------
  make_dim_plot <- function(reduction) {
    input$refresh_plot
    sobj <- active_obj(); req(sobj)
    avail <- names(sobj@reductions)
    # Map to available reduction
    red_map <- list(
      umap = c("umap", "umap.sketch", "full.umap.sketch"),
      tsne = c("tsne"),
      pca  = c("pca", "pca.sketch", "full.pca.sketch")
    )
    candidates <- red_map[[reduction]]
    found <- intersect(candidates, avail)
    if (length(found) == 0) return(NULL)
    color_by <- input$dim_color
    if (is.null(color_by) || !color_by %in% colnames(sobj@meta.data)) color_by <- "seurat_clusters"
    if (!color_by %in% colnames(sobj@meta.data)) return(NULL)

    isolate({
      n_groups <- length(unique(sobj@meta.data[[color_by]]))
      pal <- if (input$color_pal == "cat") get_palette(n_groups)
             else if (input$color_pal %in% c("Set1", "Dark2")) colorRampPalette(brewer.pal(8, input$color_pal))(n_groups)
             else if (input$color_pal == "viridis") viridis(n_groups)
             else plasma(n_groups)

      DimPlot(sobj, reduction = found[1], group.by = color_by, pt.size = input$pt_size,
              cols = pal, label = TRUE, repel = TRUE) +
        theme_minimal() +
        theme(legend.position = "right") +
        ggtitle(paste0(toupper(reduction), " — ", color_by))
    })
  }

  output$plot_umap <- renderPlot({ make_dim_plot("umap") })
  output$plot_tsne <- renderPlot({ make_dim_plot("tsne") })
  output$plot_pca  <- renderPlot({ make_dim_plot("pca") })

  # --- Spatial Projection ---------------------------------------------------
  observeEvent(input$btn_refresh_annot, {
    sobj <- active_obj(); req(sobj)
    cols <- colnames(sobj@meta.data)
    updateSelectInput(session, "annot_col", choices = cols,
                      selected = if ("cell_type" %in% cols) "cell_type" else cols[1])
    clusters <- sort(unique(sobj$seurat_clusters))
    updateSelectInput(session, "manual_cluster", choices = as.character(clusters))
  })

  output$plot_spatial_proj <- renderPlot({
    input$refresh_plot
    sobj <- active_obj(); req(sobj)
    col <- input$annot_col
    if (is.null(col) || !col %in% colnames(sobj@meta.data)) col <- "seurat_clusters"
    if (!col %in% colnames(sobj@meta.data)) return(NULL)
    isolate({
      tryCatch({
        SpatialDimPlot(sobj, group.by = col, label = TRUE, repel = TRUE,
                       pt.size.factor = input$pt_size * 1.5) +
          theme(legend.position = "right") + ggtitle(paste("Spatial Projection —", col))
      }, error = function(e) {
        # Fallback if no spatial image
        n_groups <- length(unique(sobj@meta.data[[col]]))
        pal <- get_palette(n_groups)
        coords <- tryCatch(GetTissueCoordinates(sobj), error = function(e2) NULL)
        if (is.null(coords)) return(NULL)
        coords[[col]] <- sobj@meta.data[rownames(coords), col]
        cnames <- colnames(coords)
        ggplot(coords, aes_string(x = cnames[1], y = cnames[2], color = col)) +
          geom_point(size = input$pt_size, alpha = input$pt_alpha) +
          scale_color_manual(values = pal) +
          scale_y_reverse() + coord_fixed() + theme_void() +
          ggtitle(paste("Spatial Projection —", col))
      })
    })
  })

  output$plot_spatial_gene <- renderPlot({
    sobj <- active_obj(); req(sobj, input$sp_gene)
    gene <- trimws(input$sp_gene)
    if (gene == "" || !gene %in% rownames(sobj)) return(NULL)
    tryCatch({
      SpatialFeaturePlot(sobj, features = gene, pt.size.factor = input$pt_size * 1.5)
    }, error = function(e) NULL)
  })

  # --- Manual Annotation ----------------------------------------------------
  observeEvent(input$btn_manual_annot, {
    sobj <- active_obj(); req(sobj, input$manual_cluster, input$manual_label)
    label <- trimws(input$manual_label)
    if (label == "") { showNotification("請輸入標註名稱", type = "warning"); return() }
    if (!"cell_type" %in% colnames(sobj@meta.data)) {
      sobj[["cell_type"]] <- as.character(sobj$seurat_clusters)
    }
    cells <- WhichCells(sobj, idents = input$manual_cluster)
    sobj$cell_type[cells] <- label
    rv$datasets[[rv$active_dataset]]$seurat <- sobj
    showNotification(paste0("已標註 Cluster ", input$manual_cluster, " → ", label,
                            " (", length(cells), " cells)"), type = "message")
  })

  # --- Marker-Based Cell Type Annotation ------------------------------------
  observeEvent(input$btn_marker_annot, {
    sobj <- active_obj(); req(sobj)
    gene_text <- trimws(input$marker_genes)
    cell_name <- trimws(input$marker_cell_name)
    if (gene_text == "" || cell_name == "") {
      showNotification("請輸入 marker genes 和細胞群名稱", type = "warning"); return()
    }
    genes <- trimws(unlist(strsplit(gene_text, "[,;\\s]+")))
    genes <- genes[genes != ""]

    # --- Resolve gene symbols to feature names in the assay -----------------
    assay_feats <- rownames(sobj)
    # Direct match first (if rownames are already gene symbols)
    valid_genes <- intersect(genes, assay_feats)

    # If no match, try converting SYMBOL → Ensembl ID via org.*.eg.db
    if (length(valid_genes) < length(genes)) {
      unmatched <- setdiff(genes, valid_genes)
      org_key <- if (rv$organism == "mouse") "Mm" else "Hs"
      map_fn <- paste0("sym2eg", org_key)

      if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
          requireNamespace(org_key, quietly = TRUE)) {
        tryCatch({
          db_pkg <- get(paste0("org.", org_key, ".eg.db"), envir = asNamespace(org_key))
          map_res <- AnnotationDbi::select(db_pkg,
                                           keys = unmatched,
                                           columns = "ENSEMBL",
                                           keytype = "SYMBOL")
          # Keep unique, non-NA Ensembl IDs that exist in the assay
          ensembl_matches <- map_res$ENSEMBL[!is.na(map_res$ENSEMBL) & map_res$ENSEMBL %in% assay_feats]
          ensembl_matches <- unique(ensembl_matches)
          new_hits <- setdiff(ensembl_matches, valid_genes)
          if (length(new_hits) > 0) {
            valid_genes <- c(valid_genes, new_hits)
          }
        }, error = function(e) NULL)
      }
    }

    # If still no match, try reverse: assume input is Ensembl ID
    if (length(valid_genes) == 0) {
      valid_genes <- intersect(genes, assay_feats)
    }

    if (length(valid_genes) == 0) {
      sample_feats <- head(assay_feats, 5)
      showNotification(
        paste0("找不到匹配基因。輸入: ", paste(genes, collapse = ", "),
               "\nAssay 特徵範例: ", paste(sample_feats, collapse = ", "),
               "\n請確認使用正確的基因命名 (SYMBOL 或 Ensembl ID)"),
        type = "error", duration = 15)
      return()
    }
    tryCatch({
      sobj <- AddModuleScore(sobj, features = list(valid_genes), name = "marker_score")
      if (!"cell_type" %in% colnames(sobj@meta.data)) {
        sobj[["cell_type"]] <- as.character(sobj$seurat_clusters)
      }
      positive_cells <- sobj$marker_score1 > input$marker_threshold
      sobj$cell_type[positive_cells] <- cell_name
      rv$datasets[[rv$active_dataset]]$seurat <- sobj
      n_pos <- sum(positive_cells)
      showNotification(paste0("已定義 '", cell_name, "': ", n_pos, " cells (",
                              length(valid_genes), "/", length(genes), " genes matched)"),
                       type = "message")
    }, error = function(e) {
      showNotification(paste("Marker 定義錯誤:", e$message), type = "error")
    })
  })

  # --- RCTD Deconvolution ---------------------------------------------------
  observeEvent(input$btn_run_rctd, {
    sobj <- active_obj(); req(sobj, rv$ref_seurat)
    if (!requireNamespace("spacexr", quietly = TRUE)) {
      showNotification("請安裝 spacexr: devtools::install_github('dmcable/spacexr')", type = "error")
      return()
    }
    withProgress(message = "RCTD 分析中...", value = 0, {
      tryCatch({
        ref <- rv$ref_seurat
        label_col <- trimws(input$ref_label_col)
        if (!label_col %in% colnames(ref@meta.data)) {
          stop(paste0("Reference 中找不到欄位: ", label_col))
        }
        incProgress(0.2, detail = "Preparing reference...")
        counts_ref <- ref[["RNA"]]$counts
        cluster_ref <- as.factor(ref@meta.data[[label_col]])
        levels(cluster_ref) <- gsub("/", "-", levels(cluster_ref))
        cluster_ref <- droplevels(cluster_ref)
        nUMI_ref <- ref$nCount_RNA
        reference <- spacexr::Reference(counts_ref, cluster_ref, nUMI_ref)

        incProgress(0.2, detail = "Preparing query...")
        counts_q <- GetAssayData(sobj, layer = "counts")
        coords_q <- tryCatch(GetTissueCoordinates(sobj)[, 1:2], error = function(e) {
          emb <- Embeddings(sobj, reduction = names(sobj@reductions)[1])[, 1:2]
          data.frame(x = emb[, 1], y = emb[, 2])
        })
        query <- spacexr::SpatialRNA(coords_q, counts_q, colSums(counts_q))

        incProgress(0.4, detail = "Running RCTD...")
        rctd <- spacexr::create.RCTD(query, reference, max_cores = 2)
        rctd <- spacexr::run.RCTD(rctd, doublet_mode = "doublet")

        incProgress(0.1, detail = "Adding results...")
        sobj <- AddMetaData(sobj, metadata = rctd@results$results_df)
        if ("first_type" %in% colnames(sobj@meta.data)) {
          sobj$cell_type <- as.character(sobj$first_type)
        }
        rv$datasets[[rv$active_dataset]]$seurat <- sobj
        showNotification("RCTD 完成!", type = "message")
      }, error = function(e) {
        showNotification(paste("RCTD 錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # --- Differential Expression ----------------------------------------------
  observeEvent(input$btn_run_de, {
    sobj <- active_obj(); req(sobj)
    withProgress(message = "差異分析中...", value = 0, {
      tryCatch({
        if (input$de_mode == "roi") {
          req(rv$roi_a_cells, rv$roi_b_cells)
          sobj$roi_group <- "Other"
          sobj$roi_group[colnames(sobj) %in% rv$roi_a_cells] <- "ROI_A"
          sobj$roi_group[colnames(sobj) %in% rv$roi_b_cells] <- "ROI_B"
          Idents(sobj) <- "roi_group"
          incProgress(0.3, detail = "FindMarkers ROI-A vs ROI-B...")
          markers <- FindMarkers(sobj, ident.1 = "ROI_A", ident.2 = "ROI_B",
                                 logfc.threshold = input$de_logfc)
        } else if (input$de_mode == "roi_celltype") {
          req(rv$roi_a_cells, rv$roi_b_cells, input$de_celltype)
          ct <- input$de_celltype
          if (!"cell_type" %in% colnames(sobj@meta.data)) stop("請先定義細胞型態")
          # Find cells of this type in each ROI
          roi_a_ct <- intersect(rv$roi_a_cells, colnames(sobj)[sobj$cell_type == ct])
          roi_b_ct <- intersect(rv$roi_b_cells, colnames(sobj)[sobj$cell_type == ct])
          if (length(roi_a_ct) < 3 || length(roi_b_ct) < 3)
            stop(paste0("ROI 中 '", ct, "' 數量不足 (A:", length(roi_a_ct), ", B:", length(roi_b_ct), ")"))
          sobj$cross_roi_ct <- "Other"
          sobj$cross_roi_ct[colnames(sobj) %in% roi_a_ct] <- "ROI_A"
          sobj$cross_roi_ct[colnames(sobj) %in% roi_b_ct] <- "ROI_B"
          Idents(sobj) <- "cross_roi_ct"
          incProgress(0.3, detail = paste0("FindMarkers ", ct, " ROI-A vs ROI-B..."))
          markers <- FindMarkers(sobj, ident.1 = "ROI_A", ident.2 = "ROI_B",
                                 logfc.threshold = input$de_logfc)
        } else if (input$de_mode == "group") {
          req(input$de_ident_col, input$de_group1, input$de_group2)
          Idents(sobj) <- input$de_ident_col
          incProgress(0.3, detail = paste("FindMarkers", input$de_group1, "vs", input$de_group2))
          markers <- FindMarkers(sobj, ident.1 = input$de_group1, ident.2 = input$de_group2,
                                 logfc.threshold = input$de_logfc)
        } else {
          incProgress(0.3, detail = "FindAllMarkers...")
          markers <- FindAllMarkers(sobj, only.pos = TRUE, logfc.threshold = input$de_logfc)
        }
        incProgress(0.5, detail = "Done!")
        markers$gene <- rownames(markers)
        rv$de_results <- markers
        rv$datasets[[rv$active_dataset]]$seurat <- sobj
        showNotification(paste0("找到 ", nrow(markers), " 個差異基因"), type = "message")
      }, error = function(e) {
        showNotification(paste("DE 錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  # Update DE group selectors
  observe({
    sobj <- active_obj(); req(sobj)
    cols <- colnames(sobj@meta.data)
    updateSelectInput(session, "de_ident_col", choices = cols)
    # Update cross-ROI celltype choices
    if ("cell_type" %in% cols) {
      ct_choices <- sort(unique(sobj$cell_type))
      updateSelectInput(session, "de_celltype", choices = ct_choices)
    }
  })
  observeEvent(input$de_ident_col, {
    sobj <- active_obj(); req(sobj, input$de_ident_col)
    if (input$de_ident_col %in% colnames(sobj@meta.data)) {
      groups <- sort(unique(as.character(sobj@meta.data[[input$de_ident_col]])))
      updateSelectInput(session, "de_group1", choices = groups)
      updateSelectInput(session, "de_group2", choices = groups)
    }
  })

  output$de_table <- DT::renderDataTable({
    req(rv$de_results)
    rv$de_results
  }, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)

  output$de_heatmap <- renderPlot({
    sobj <- active_obj(); req(sobj, rv$de_results)
    markers <- rv$de_results
    if (input$de_mode == "all" && "cluster" %in% colnames(markers)) {
      top_genes <- markers %>% group_by(cluster) %>%
        dplyr::filter(p_val_adj < 0.05) %>% slice_head(n = 5) %>% ungroup()
      genes <- unique(top_genes$gene)
    } else {
      genes <- head(rownames(markers[order(markers$p_val_adj), ]), input$de_top_n)
    }
    genes <- intersect(genes, rownames(sobj))
    if (length(genes) == 0) return(NULL)
    tryCatch({
      sobj_sub <- subset(sobj, downsample = 500)
      sobj_sub <- ScaleData(sobj_sub, features = genes)
      DoHeatmap(sobj_sub, features = genes, size = 3) +
        theme(axis.text.y = element_text(size = 7)) + NoLegend()
    }, error = function(e) NULL)
  })

  # --- Enrichment Analysis --------------------------------------------------
  observeEvent(input$btn_run_enrich, {
    req(rv$de_results)
    withProgress(message = "富集分析中...", value = 0, {
      tryCatch({
        markers <- rv$de_results
        if ("avg_log2FC" %in% colnames(markers)) {
          genes_up <- markers %>% dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
            arrange(desc(avg_log2FC))
        } else {
          genes_up <- markers %>% dplyr::filter(p_val_adj < 0.05)
        }
        gene_symbols <- if ("gene" %in% colnames(genes_up)) genes_up$gene else rownames(genes_up)

        org_db <- ORG_DB[[rv$organism]]$orgDb
        incProgress(0.3, detail = "Converting gene IDs...")
        entrez <- tryCatch({
          AnnotationDbi::mapIds(get(org_db), keys = gene_symbols,
                                keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
        }, error = function(e) setNames(character(0), character(0)))
        entrez <- entrez[!is.na(entrez)]

        incProgress(0.4, detail = "Running enrichment...")
        if (input$enrich_type == "kegg") {
          kegg_code <- if (rv$organism == "human") "hsa" else "mmu"
          res <- enrichKEGG(gene = entrez, organism = kegg_code, pvalueCutoff = input$enrich_pval)
        } else if (startsWith(input$enrich_type, "go_")) {
          ont <- switch(input$enrich_type, go_bp = "BP", go_mf = "MF", go_cc = "CC")
          res <- enrichGO(gene = entrez, OrgDb = get(org_db), ont = ont,
                          pvalueCutoff = input$enrich_pval, readable = TRUE)
        } else {
          # GSEA
          gene_list <- setNames(genes_up$avg_log2FC, gene_symbols)
          gene_list <- sort(gene_list, decreasing = TRUE)
          # Convert to entrez for GSEA
          names_entrez <- tryCatch({
            AnnotationDbi::mapIds(get(org_db), keys = names(gene_list),
                                  keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
          }, error = function(e) NULL)
          if (!is.null(names_entrez)) {
            valid <- !is.na(names_entrez)
            gl <- gene_list[valid]; names(gl) <- names_entrez[valid]
            gl <- sort(gl, decreasing = TRUE)
            res <- gseGO(geneList = gl, OrgDb = get(org_db), ont = "BP",
                         pvalueCutoff = input$enrich_pval)
          } else {
            stop("基因 ID 轉換失敗")
          }
        }
        incProgress(0.2, detail = "Done!")
        rv$enrich_results <- res
        showNotification("富集分析完成!", type = "message")
      }, error = function(e) {
        showNotification(paste("富集分析錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  output$enrich_plot <- renderPlot({
    req(rv$enrich_results)
    res <- rv$enrich_results
    if (nrow(as.data.frame(res)) == 0) return(NULL)
    if (inherits(res, "gseaResult")) {
      dotplot(res, showCategory = input$enrich_show) + ggtitle("GSEA Results")
    } else {
      dotplot(res, showCategory = input$enrich_show) + ggtitle("Enrichment Results")
    }
  })

  output$enrich_table <- DT::renderDataTable({
    req(rv$enrich_results)
    as.data.frame(rv$enrich_results)
  }, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)

  # --- Trajectory Analysis --------------------------------------------------
  observeEvent(input$btn_refresh_annot, {
    # Also update trajectory selectors
    sobj <- active_obj(); req(sobj)
    cols <- colnames(sobj@meta.data)
    updateSelectInput(session, "traj_cluster_col", choices = cols,
                      selected = if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1])
    clusters <- sort(unique(sobj$seurat_clusters))
    updateSelectInput(session, "traj_start", choices = c("Auto" = "", as.character(clusters)))
  })

  observeEvent(input$btn_run_traj, {
    sobj <- active_obj(); req(sobj)
    withProgress(message = "軌跡分析中...", value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Preparing data...")

        # ROI + cell type subset
        if (input$traj_roi_subset && !is.null(input$traj_celltype) && input$traj_celltype != "") {
          roi_cells <- character(0)
          if (input$traj_roi_choice == "a") roi_cells <- rv$roi_a_cells
          else if (input$traj_roi_choice == "b") roi_cells <- rv$roi_b_cells
          else roi_cells <- union(rv$roi_a_cells, rv$roi_b_cells)
          if (length(roi_cells) == 0) stop("請先設定 ROI")
          ct_cells <- colnames(sobj)[sobj$cell_type == input$traj_celltype]
          subset_cells <- intersect(roi_cells, ct_cells)
          if (length(subset_cells) < 10) stop(paste0("符合條件的 cells 太少: ", length(subset_cells)))
          sobj <- subset(sobj, cells = subset_cells)
        }

        # Get PCA embeddings
        avail_red <- names(sobj@reductions)
        pca_red <- intersect(c("pca", "pca.sketch", "full.pca.sketch"), avail_red)
        if (length(pca_red) == 0) stop("請先執行 PCA")
        pca_emb <- Embeddings(sobj, reduction = pca_red[1])
        dims <- min(input$traj_dims, ncol(pca_emb))
        pca_sub <- pca_emb[, 1:dims]

        cluster_col <- input$traj_cluster_col
        if (is.null(cluster_col) || !cluster_col %in% colnames(sobj@meta.data))
          cluster_col <- "seurat_clusters"
        cl <- sobj@meta.data[[cluster_col]]

        incProgress(0.4, detail = "Running slingshot...")
        start_cl <- if (input$traj_start == "") NULL else input$traj_start
        sds <- slingshot(SlingshotDataSet(
          reducedDim = pca_sub,
          clusterLabels = cl
        ), start.clus = start_cl)

        rv$slingshot_res <- sds
        incProgress(0.2, detail = "Done!")
        showNotification("軌跡分析完成!", type = "message")
      }, error = function(e) {
        showNotification(paste("軌跡分析錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  output$traj_plot <- renderPlot({
    sobj <- active_obj(); req(sobj, rv$slingshot_res)
    sds <- rv$slingshot_res
    # Get UMAP for visualization
    avail_red <- names(sobj@reductions)
    umap_red <- intersect(c("umap", "umap.sketch", "full.umap.sketch"), avail_red)
    if (length(umap_red) == 0) return(NULL)
    emb <- Embeddings(sobj, reduction = umap_red[1])[, 1:2]
    df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2])

    col_name <- input$traj_cluster_col
    if (!is.null(col_name) && col_name %in% colnames(sobj@meta.data)) {
      df$cluster <- sobj@meta.data[[col_name]]
    } else {
      df$cluster <- "1"
    }

    n_cl <- length(unique(df$cluster))
    pal <- get_palette(n_cl)

    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
      geom_point(size = input$pt_size, alpha = input$pt_alpha) +
      scale_color_manual(values = pal) +
      theme_minimal() +
      ggtitle("Trajectory (Slingshot) on UMAP")

    # Try to overlay slingshot curves
    tryCatch({
      curves <- slingCurves(sds)
      for (crv in curves) {
        curve_df <- as.data.frame(crv$s[crv$ord, 1:2])
        # Map PCA curve to UMAP space approximately (just plot PCA curve)
        colnames(curve_df) <- c("x", "y")
      }
      # For proper visualization, show on PCA space instead
      pca_red <- intersect(c("pca", "pca.sketch", "full.pca.sketch"), avail_red)
      pca_emb <- Embeddings(sobj, reduction = pca_red[1])[, 1:2]
      df2 <- data.frame(PC1 = pca_emb[, 1], PC2 = pca_emb[, 2], cluster = df$cluster)
      p2 <- ggplot(df2, aes(x = PC1, y = PC2, color = cluster)) +
        geom_point(size = input$pt_size, alpha = input$pt_alpha) +
        scale_color_manual(values = pal) +
        theme_minimal() + ggtitle("Trajectory on PCA")
      for (crv in curves) {
        curve_df <- as.data.frame(crv$s[crv$ord, 1:2])
        colnames(curve_df) <- c("x", "y")
        p2 <- p2 + geom_path(data = curve_df, aes(x = x, y = y),
                              inherit.aes = FALSE, linewidth = 1.5, color = "black")
      }
      p + p2
    }, error = function(e) p)
  })

  # --- CellChat Analysis ----------------------------------------------------
  observeEvent(input$btn_run_cellchat, {
    sobj <- active_obj(); req(sobj)
    ct_col <- input$cc_celltype_col
    if (is.null(ct_col) || !ct_col %in% colnames(sobj@meta.data)) {
      showNotification("請選擇 cell type 欄位", type = "warning"); return()
    }
    withProgress(message = "CellChat 分析中...", value = 0, {
      tryCatch({
        # Subset by ROI scope
        if (input$cc_scope == "roi_a") {
          req(rv$roi_a_cells)
          sobj <- subset(sobj, cells = intersect(colnames(sobj), rv$roi_a_cells))
        } else if (input$cc_scope == "roi_b") {
          req(rv$roi_b_cells)
          sobj <- subset(sobj, cells = intersect(colnames(sobj), rv$roi_b_cells))
        }

        incProgress(0.2, detail = "Creating CellChat object...")
        data_input <- GetAssayData(sobj, layer = "data")
        meta <- sobj@meta.data
        meta$cell_type_cc <- as.character(meta[[ct_col]])
        cellchat <- createCellChat(object = data_input, meta = meta, group.by = "cell_type_cc")

        incProgress(0.2, detail = "Setting DB...")
        db_species <- if (rv$organism == "human") CellChatDB.human else CellChatDB.mouse
        db_use <- subsetDB(db_species, search = input$cc_db)
        cellchat@DB <- db_use

        incProgress(0.3, detail = "Computing interactions...")
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        cellchat <- computeCommunProb(cellchat, type = "triMean")
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)

        incProgress(0.2, detail = "Done!")
        rv$cellchat_obj <- cellchat
        showNotification("CellChat 分析完成!", type = "message")
      }, error = function(e) {
        showNotification(paste("CellChat 錯誤:", e$message), type = "error", duration = 10)
      })
    })
  })

  output$cc_plot <- renderPlot({
    req(rv$cellchat_obj)
    cc <- rv$cellchat_obj
    tryCatch({
      if (input$cc_plot_type == "circle") {
        groupSize <- as.numeric(table(cc@idents))
        netVisual_circle(cc@net$count, vertex.weight = groupSize,
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = "Number of interactions")
      } else if (input$cc_plot_type == "bubble") {
        netVisual_bubble(cc, remove.isolate = FALSE)
      } else {
        netVisual_heatmap(cc)
      }
    }, error = function(e) NULL)
  })

  output$cc_table <- DT::renderDataTable({
    req(rv$cellchat_obj)
    df <- subsetCommunication(rv$cellchat_obj)
    if (nrow(df) > 0) df <- df[order(df$prob, decreasing = TRUE), ]
    df
  }, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)

  # Update CellChat selectors
  observeEvent(input$btn_refresh_annot, {
    sobj <- active_obj(); req(sobj)
    cols <- colnames(sobj@meta.data)
    updateSelectInput(session, "cc_celltype_col", choices = cols,
                      selected = if ("cell_type" %in% cols) "cell_type" else cols[1])
    # Also update trajectory celltype choices
    if ("cell_type" %in% cols) {
      ct <- sort(unique(sobj$cell_type))
      updateSelectInput(session, "traj_celltype", choices = ct)
    }
  })

  # --- Report Export --------------------------------------------------------
  methods_reactive <- reactiveVal("")

  observeEvent(input$btn_gen_report, {
    sobj <- active_obj()
    # Gather parameters
    ds_name <- if (!is.null(rv$active_dataset)) rv$active_dataset else "N/A"
    n_cells <- if (!is.null(sobj)) ncol(sobj) else "N/A"
    n_genes <- if (!is.null(sobj)) nrow(sobj) else "N/A"
    org <- rv$organism
    norm <- input$norm_method
    sketch <- input$use_sketch
    sketch_n <- input$sketch_n
    res <- input$cluster_res
    pcs <- input$n_pcs
    integ <- input$integ_method
    de_lfc <- input$de_logfc
    qc_nc <- input$qc_ncount
    qc_nf <- input$qc_nfeat
    qc_mt <- input$qc_mt

    norm_text <- if (norm == "sct") {
      "Normalization was performed using SCTransform (Hafemeister & Satija, 2019)."
    } else {
      "Data were log-normalized using NormalizeData with default parameters (scale factor = 10,000)."
    }

    sketch_text <- if (sketch && !is.null(sobj) && ncol(sobj) > 100000) {
      paste0("Due to the large dataset size, the Seurat v5 sketch-based workflow was employed, ",
             "subsampling ", format(sketch_n, big.mark = ","), " cells using LeverageScore method. ",
             "Clustering and dimensional reduction were performed on the sketched subset ",
             "and projected back to the full dataset using ProjectData().")
    } else {
      ""
    }

    integ_text <- ""
    if (length(rv$datasets) > 1) {
      method_name <- switch(integ, harmony = "Harmony (Korsunsky et al., 2019)",
                            cca = "Canonical Correlation Analysis (CCA)",
                            rpca = "Reciprocal PCA (RPCA)")
      integ_text <- paste0("Multi-sample integration was performed using IntegrateLayers() ",
                           "with the ", method_name, " method to correct for batch effects. ")
    }

    methods <- paste0(
      "Spatial transcriptomics data were generated using the 10x Genomics Visium HD platform. ",
      "Raw sequencing data were processed using Space Ranger to generate binned expression matrices. ",
      "Downstream analysis was performed using the Visium HD Analysis System, ",
      "a custom R Shiny application built on Seurat v", packageVersion("Seurat"), 
      " (Hao et al., 2024) in R v", paste(R.version$major, R.version$minor, sep = "."), ". ",
      "\n\n",
      "Quality control filtering was applied to remove low-quality bins ",
      "(nCount_Spatial > ", qc_nc, ", nFeature_Spatial > ", qc_nf,
      ", percent mitochondrial genes < ", qc_mt, "%). ",
      norm_text, " ",
      "Highly variable features were identified using FindVariableFeatures(). ",
      sketch_text,
      "\n\n",
      integ_text,
      "Principal component analysis (PCA) was performed, and the first ", pcs,
      " principal components were used for downstream analysis. ",
      "Cell clustering was performed using the Louvain algorithm via FindClusters() ",
      "at a resolution of ", res, ". ",
      "UMAP dimensional reduction was computed for visualization. ",
      "\n\n",
      "Differential gene expression analysis was performed using FindMarkers() ",
      "with a minimum log2 fold-change threshold of ", de_lfc, " ",
      "(Wilcoxon rank-sum test). ",
      "Gene Ontology and KEGG pathway enrichment analyses were performed using ",
      "clusterProfiler v", packageVersion("clusterProfiler"), " (Wu et al., 2021). ",
      "Trajectory inference was performed using slingshot v", packageVersion("slingshot"),
      " (Street et al., 2018). ",
      "Cell-cell communication analysis was performed using CellChat v",
      packageVersion("CellChat"), " (Jin et al., 2021). ",
      "\n\n",
      "All visualizations were generated using ggplot2 v", packageVersion("ggplot2"),
      " and plotly v", packageVersion("plotly"), "."
    )

    methods_reactive(methods)
    showNotification("分析報告已產生!", type = "message")
  })

  output$methods_text <- renderText({
    methods_reactive()
  })

  output$params_table <- DT::renderDataTable({
    sobj <- active_obj()
    params <- data.frame(
      Parameter = c(
        "Active Dataset", "Cells/Bins", "Genes", "Organism",
        "Normalization", "Sketch Workflow", "Sketch N",
        "Clustering Resolution", "PCA Dimensions",
        "Integration Method", "QC: min nCount", "QC: min nFeature",
        "QC: max percent.mt", "DE: min log2FC",
        "ROI-A cells", "ROI-B cells"
      ),
      Value = c(
        if (!is.null(rv$active_dataset)) rv$active_dataset else "N/A",
        if (!is.null(sobj)) as.character(ncol(sobj)) else "N/A",
        if (!is.null(sobj)) as.character(nrow(sobj)) else "N/A",
        rv$organism,
        if (input$norm_method == "sct") "SCTransform" else "LogNormalize",
        as.character(input$use_sketch),
        as.character(input$sketch_n),
        as.character(input$cluster_res),
        as.character(input$n_pcs),
        input$integ_method,
        as.character(input$qc_ncount),
        as.character(input$qc_nfeat),
        paste0(input$qc_mt, "%"),
        as.character(input$de_logfc),
        if (!is.null(rv$roi_a_cells)) as.character(length(rv$roi_a_cells)) else "N/A",
        if (!is.null(rv$roi_b_cells)) as.character(length(rv$roi_b_cells)) else "N/A"
      ),
      stringsAsFactors = FALSE
    )
    params
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$session_info <- renderPrint({
    cat("=== R Session Info ===", "\n")
    si <- sessionInfo()
    print(si)
    cat("\n=== Key Package Versions ===", "\n")
    pkgs <- c("Seurat", "SeuratObject", "ggplot2", "plotly", "dplyr",
              "clusterProfiler", "enrichplot", "slingshot", "CellChat",
              "harmony", "hdf5r", "shiny", "shinydashboard")
    for (p in pkgs) {
      v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
      cat(sprintf("  %-20s %s\n", p, v))
    }
  })

  output$btn_dl_methods <- downloadHandler(
    filename = function() {
      paste0("visiumHD_methods_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      writeLines(methods_reactive(), file)
    }
  )

} # end server

# =============================================================================
# 啟動應用程式
# =============================================================================
shinyApp(ui = ui, server = server)