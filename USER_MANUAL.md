# Visium HD 分析系統 — 使用說明書

## 目錄

1. [系統需求與安裝](#1-系統需求與安裝)
2. [啟動應用程式](#2-啟動應用程式)
3. [Tab 1：資料載入與管理](#3-tab-1資料載入與管理)
4. [Tab 2：QC 與前處理](#4-tab-2qc-與前處理)
5. [Tab 3：ROI 圈選與降維分析](#5-tab-3roi-圈選與降維分析)
6. [Tab 4：空間投影與細胞型態註解](#6-tab-4空間投影與細胞型態註解)
7. [Tab 5：差異基因表達分析](#7-tab-5差異基因表達分析)
8. [Tab 6：功能富集分析](#8-tab-6功能富集分析)
9. [Tab 7：軌跡推斷分析](#9-tab-7軌跡推斷分析)
10. [Tab 8：細胞通訊分析](#10-tab-8細胞通訊分析)
11. [Tab 9：分析報告匯出](#11-tab-9分析報告匯出)
12. [側邊欄：全域視覺化控制](#12-側邊欄全域視覺化控制)
13. [推薦分析流程](#13-推薦分析流程)
14. [常見問題 FAQ](#14-常見問題-faq)

---

## 1. 系統需求與安裝

### 硬體建議

| 項目 | 最低需求 | 建議配置 |
|---|---|---|
| RAM | 16 GB | 64 GB 以上（Visium HD 8µm 資料集通常 >10GB）|
| CPU | 4 核心 | 8+ 核心 |
| 磁碟 | SSD 建議 | — |

### R 版本

- **R ≥ 4.2.0**（建議 4.3+）

### 必要套件安裝

```r
# --- CRAN 套件 ---
install.packages(c(
  "shiny", "shinydashboard", "shinyFiles", "shinyjs",
  "plotly", "ggplot2", "dplyr", "tidyr", "tibble",
  "RColorBrewer", "viridis", "scales", "pheatmap",
  "DT", "patchwork", "harmony"
))

# --- Bioconductor 套件 ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Seurat", "SeuratObject", "hdf5r",
  "clusterProfiler", "enrichplot", "DOSE",
  "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi",
  "slingshot", "SingleCellExperiment"
))

# --- GitHub 套件 ---
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("jinworks/CellChat")
devtools::install_github("dmcable/spacexr")  # RCTD（選用）
```

---

## 2. 啟動應用程式

```r
# 在 R console 或 RStudio 中執行：
shiny::runApp("visium_HD_analysis.R")

# 或指定 port：
shiny::runApp("visium_HD_analysis.R", port = 3838, host = "0.0.0.0")
```

啟動後瀏覽器會自動開啟，顯示 Dashboard 界面。

---

## 3. Tab 1：資料載入與管理

### 3.1 從 Space Ranger 載入

1. **Dataset 名稱**：輸入一個識別名稱（如 `skin_sample1_8um`）
2. **物種**：選擇 Human 或 Mouse（影響 QC 指標的基因 pattern）
3. **Bin Size**：
   - `8µm`（預設，推薦起始解析度）
   - `16µm`（較低解析度，計算更快）
   - `2µm`（亞細胞解析度，僅用於 cell segmentation）
   - `Custom`：自訂 bin code
4. **選擇目錄**：
   - 點擊 **Browse...** 按鈕使用檔案瀏覽器
   - 或在文字框手動輸入路徑
   - 系統會自動偵測並預覽目錄結構（`spatial/`、`.h5`、`binned_outputs/`）

> **路徑提示**：可以指向 Space Ranger 輸出的任一層級：
> - `/path/to/sample/` — 自動偵測 `outs/`
> - `/path/to/sample/outs/` — 直接使用
> - `/path/to/sample/outs/binned_outputs/square_008um/` — 自動回溯

5. **資料格式**：HDF5（.h5，預設）或 Matrix（.mtx）
6. 點擊 **載入 Dataset** 按鈕

### 3.2 從 RDS 載入

- 上傳已儲存的 Seurat RDS 檔案
- 適用於已完成部分分析的資料集
- 同樣需指定物種

### 3.3 RCTD Reference（選用）

- 上傳 scRNA-seq reference RDS 檔案
- 指定 cell type 欄位名稱（如 `subclass_label`）
- 供 Tab 4 的 RCTD deconvolution 使用

### 3.4 多資料集管理

- 載入的 datasets 顯示在下方表格
- **設為作用中**：選擇一行 → 點擊按鈕
- **移除**：選擇一行 → 點擊移除
- **摘要**：顯示作用中 dataset 的 Seurat 物件資訊

---

## 4. Tab 2：QC 與前處理

### 4.1 QC 過濾

| 參數 | 說明 | 預設值 |
|---|---|---|
| 最低 nCount | 每個 bin 的最低 UMI 總數 | 100 |
| 最低 nFeature | 每個 bin 的最低偵測基因數 | 50 |
| 最高 percent.mt | 粒線體基因比例上限（%）| 25 |

- **QC 視覺化**：Violin plot（nCount, nFeature, percent.mt）+ Scatter plot
- 點擊 **套用 QC 過濾** 後，不符條件的 cells/bins 會被移除

> **注意**：閾值應依組織類型調整。代謝活躍組織（如心臟、肌肉）的 percent.mt 自然較高。

### 4.2 前處理 Pipeline

| 參數 | 說明 |
|---|---|
| Normalization 方式 | `LogNormalize`（預設）或 `SCTransform`（更適合 count depth 差異大的樣品）|
| Sketch Workflow | 當 cells >100,000 時自動啟用，使用 LeverageScore 下採樣 |
| Sketch cells | Sketch 數量（預設 50,000）|
| Clustering Resolution | 越高 = 越多 clusters（預設 0.8）|
| PCA dims | 使用的主成分數（預設 30）|

點擊 **執行 Normalize + Cluster + UMAP** 開始完整分析流程。

### 4.3 多資料集整合

當載入 ≥2 個 datasets 後：

1. 確認各 dataset 已完成 QC 過濾
2. 選擇 Integration 方法：
   - **Harmony**（推薦，Seurat 官方空間資料首選）
   - **CCA**（Canonical Correlation Analysis）
   - **RPCA**（Reciprocal PCA）
3. 點擊 **合併 + 整合所有 Datasets**

> **流程**：各自 Normalize → Merge → PCA → IntegrateLayers → Cluster → UMAP

整合後會產生新的 `integrated_Xds` dataset 並自動設為作用中。

---

## 5. Tab 3：ROI 圈選與降維分析

### 5.1 空間座標圈選

- 空間圖以 Plotly 互動圖呈現
- 使用 **Lasso 工具**（工具列最右側的套索圖示）在空間圖上圈選 ROI
- 圈選後：
  - 點擊 **儲存為 ROI-A** → 儲存第一個區域
  - 重新圈選 → 點擊 **儲存為 ROI-B** → 儲存第二個區域
- 下方顯示各 ROI 的 cell 數量

> **效能提示**：當 cells >50,000 時，空間圖會自動 downsample 以維持互動流暢度。

### 5.2 降維圖

- **UMAP** / **t-SNE** / **PCA** 三個 tab 切換
- 使用 **降維圖上色依據** 下拉選單選擇要上色的 metadata 欄位
- 點擊 **更新欄位選單** 刷新可用欄位

---

## 6. Tab 4：空間投影與細胞型態註解

### 6.1 空間投影

- 顯示 `SpatialDimPlot`（依 cluster 或 cell type 上色）
- 可輸入基因名稱查看單一基因的空間表現分布

### 6.2 細胞型態註解方式

#### A. 手動標註（Cluster-based）

1. 點擊 **更新欄位** 刷新 cluster 清單
2. 選擇一個 cluster 編號
3. 輸入細胞型態名稱（如 `Mast cell`）
4. 點擊 **套用標註**

#### B. 以 Markers 定義細胞群

1. 輸入 marker genes，逗號分隔（如 `TPSAB1, CPA3, KIT`）
2. 輸入細胞群名稱（如 `Mast cell`）
3. 調整 Module Score 閾值（預設 0，越高越嚴格）
4. 點擊 **依 Markers 定義細胞**

> **原理**：使用 `Seurat::AddModuleScore()` 計算每個 cell 的 marker 聯合表現分數，超過閾值的 cells 被標記為指定型態。

#### C. RCTD Deconvolution

1. 確認已在 Tab 1 上傳 scRNA-seq Reference RDS
2. 點擊 **執行 RCTD**
3. 結果自動寫入 `first_type` 和 `cell_type` 欄位

---

## 7. Tab 5：差異基因表達分析

### 比較模式

| 模式 | 說明 |
|---|---|
| **ROI-A vs ROI-B** | 全部 ROI-A cells vs 全部 ROI-B cells |
| **跨 ROI 同一細胞群比較** | ROI-A 中的特定 cell type vs ROI-B 中的同一 cell type |
| **指定群組比較** | 選擇任意 metadata 欄位 + 兩個 group 做比較 |
| **Find All Markers** | 找出每個 cluster 的 marker genes |

### 參數

- **min log2FC**：最低 fold change 閾值（預設 0.25）
- **Heatmap Top N**：Heatmap 上顯示的基因數量

### 輸出

- **結果表格**：p_val, avg_log2FC, pct.1, pct.2, p_val_adj
- **Heatmap**：Top N 差異基因的表現熱圖

---

## 8. Tab 6：功能富集分析

### 分析類型

| 類型 | 說明 |
|---|---|
| GO (BP) | Gene Ontology — Biological Process |
| GO (MF) | Gene Ontology — Molecular Function |
| GO (CC) | Gene Ontology — Cellular Component |
| KEGG | KEGG Pathway |
| GSEA (GO) | Gene Set Enrichment Analysis — GO BP |

### 使用方式

1. 先在 Tab 5 完成差異分析
2. 選擇分析類型
3. 調整 p-value cutoff（預設 0.05）
4. 點擊 **執行富集分析**

### 輸出

- **Dotplot**：顯示富集程度、基因比例、p 值
- **結果表格**：完整富集結果

---

## 9. Tab 7：軌跡推斷分析

### 基本設定

1. 選擇 Cluster 欄位
2. 可選擇起始 Cluster（或留 Auto 自動判斷）
3. 設定使用的 PCA 維度數

### 進階：ROI + 細胞型態篩選

勾選 **「僅分析特定 ROI 的特定細胞群」** 後：
1. 選擇 ROI（A / B / 兩者合併）
2. 選擇要分析的細胞群
3. 執行 → 僅對篩選後的 cells 做 slingshot 軌跡推斷

### 輸出

- UMAP + PCA 雙面板視覺化
- 黑色曲線 = slingshot 推斷的發育軌跡

---

## 10. Tab 8：細胞通訊分析

### 設定

1. 選擇 Cell Type 欄位（需先在 Tab 4 完成細胞註解）
2. 選擇 Interaction Database：
   - **Secreted Signaling**：分泌型訊號分子
   - **ECM-Receptor**：細胞外基質與受體
   - **Cell-Cell Contact**：直接細胞接觸
3. 選擇分析範圍：全部 cells / 僅 ROI-A / 僅 ROI-B
4. 點擊 **執行 CellChat**

### 視覺化

| 類型 | 說明 |
|---|---|
| Circle Plot | 細胞群之間的 interaction 數量/強度網路圖 |
| Bubble Plot | Receptor-Ligand pair 的 interaction 氣泡圖 |
| Heatmap | 細胞群間 signaling 強度熱圖 |

### 輸出

- 互動式圖形
- 詳細結果表格（source, target, ligand, receptor, pathway, probability）

---

## 11. Tab 9：分析報告匯出

### 分析參數記錄

點擊 **產出分析報告** 按鈕，系統會自動生成：

1. **Methods Section（投稿用）**：
   - 完整的分析流程描述（英文）
   - 包含所有使用的參數值
   - 可直接複製到論文 Methods 段落

2. **軟體版本表**：
   - R 版本
   - 所有使用套件的版本號
   - 作業系統資訊

3. **分析參數摘要表**：
   - Normalization 方法
   - QC 閾值
   - Clustering resolution
   - Integration 方法（如適用）
   - DE 統計方法
   - Bin size

點擊 **下載 Methods 文字** 可下載 `.txt` 檔案。

---

## 12. 側邊欄：全域視覺化控制

| 控制項 | 說明 |
|---|---|
| 點大小 | 所有圖形的散點大小（0.1–5）|
| 不透明度 | 散點透明度（0.1–1）|
| 調色盤 | 高對比 / Set1 / Dark2 / Viridis / Plasma |
| 解析度 (Assay) | 當載入多個 bin size 時，切換 Spatial assay |
| 🔄 重繪圖形 | 套用新視覺化設定 |

---

## 13. 推薦分析流程

### 單一樣品分析

```
載入 (Tab 1) → QC 過濾 (Tab 2) → Normalize + Cluster (Tab 2)
    → ROI 圈選 (Tab 3) → 細胞註解 (Tab 4)
    → DE 分析 (Tab 5) → 富集分析 (Tab 6)
    → 軌跡分析 (Tab 7) → CellChat (Tab 8)
    → 匯出報告 (Tab 9)
```

### 多樣品整合分析

```
載入樣品 A (Tab 1) → QC 過濾 A (Tab 2)
載入樣品 B (Tab 1) → QC 過濾 B (Tab 2)
    → 合併 + 整合 (Tab 2, Harmony)
    → ROI 圈選 (Tab 3) → 細胞註解 (Tab 4)
    → 跨 ROI 同一細胞群 DE (Tab 5)
    → 富集分析 (Tab 6) → 軌跡分析 (Tab 7)
    → CellChat (Tab 8) → 匯出報告 (Tab 9)
```

### 皮膚組織特定分析範例

```
1. 載入 8µm 皮膚 Visium HD 資料
2. QC: nCount >200, nFeature >100, percent.mt <20
3. Normalize (LogNormalize) + Cluster (resolution 1.0)
4. ROI 圈選：
   - ROI-A: 毛囊區域
   - ROI-B: 真皮乳頭區域
5. Marker-based 細胞定義：
   - KRT15, CD200, LGR5 → Hair follicle stem cell
   - SOX2, BMP4, WNT5A → Dermal papilla cell
   - TPSAB1, CPA3, KIT → Mast cell
6. 跨 ROI DE：比較兩個區域的同一細胞群基因表現差異
7. GO/KEGG 富集分析
8. Slingshot 軌跡：分析 HFSC 的分化路徑
9. CellChat：分析 DP-HFSC 的 signaling
```

---

## 14. 常見問題 FAQ

### Q: 載入資料時出現 "找不到 spatial/ 目錄" 錯誤？

**A**: 確認路徑指向包含 `spatial/` 子目錄的 `outs/` 資料夾。系統會自動偵測，但若目錄結構非標準，請手動指向正確層級。

### Q: 記憶體不足？

**A**: 
- 使用 **Sketch Workflow**（自動啟用於 >100K cells）
- 選擇 `16µm` 解析度（減少 ~4x 資料量）
- 增加系統 RAM 或使用伺服器

### Q: RCTD 執行很慢？

**A**: RCTD 在大資料集上計算密集。建議：
- 先在 sketch 子集上執行
- 調低 `max_cores` 避免記憶體問題

### Q: 為什麼 CellChat 需要先完成細胞註解？

**A**: CellChat 需要已知的 cell type labels 來建構 interaction 網路。請先在 Tab 4 完成至少一種註解方式。

### Q: 2µm 和 8µm 資料的差異？

**A**: 
- **8µm**：標準分析解析度，每個 bin 包含足夠 UMI 做可靠的基因表現定量
- **2µm**：亞細胞解析度，資料極度稀疏，主要用於結合影像做 cell segmentation
- **建議**：一般分析從 8µm 開始，cell segmentation 才用 2µm
