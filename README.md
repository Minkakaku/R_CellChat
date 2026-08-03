# R-cellchat-workflow-template

一个独立、可移植的 CellChat 分析模板。仓库不采用复杂的 `src/`、`workflow/` 分层，核心分析只有一个文件：

```text
cellchat_workflow.R
```

脚本通过 `analysis_mode` 明确区分三种场景：

- `single`：一个样本的完整 CellChat 分析；
- `two`：两个样本或两个条件的直接比较；
- `multi`：三个及以上样本，并自动生成所有两两比较；
- 配体-受体表、通路表、CellChat 对象及主要图形输出；
- PDF 与 300 dpi PNG 双格式保存。

## 使用

先修改 `cellchat_workflow.R` 开头的参数：

```r
input_rds <- "path/to/seurat_object.rds"
species <- "human"
assay <- "RNA"
celltype_col <- "cell_type"
```

单样本分析：

```r
analysis_mode <- "single"
sample_col <- "orig.ident"
samples <- "Sample1"
```

双样本比较：

```r
analysis_mode <- "two"
sample_col <- "group"
samples <- c("CON", "EXP")
```

多样本比较：

```r
analysis_mode <- "multi"
sample_col <- "orig.ident"
samples <- c("Sample1", "Sample2", "Sample3")
```

运行：

```bash
Rscript cellchat_workflow.R
```

## 依赖

```r
install.packages(c("Seurat", "future", "ggplot2", "remotes"))
BiocManager::install("ComplexHeatmap")
remotes::install_github("jinworks/CellChat")
```

## 输出

运行后生成：

- `results/tables/`：配体-受体及通路通讯表；
- `results/plots/`：网络图、热图、气泡图、信号角色和组间差异图；
- `results/README_results_CN.md`：分析模式、样本和比较方向说明；
- `results/sessionInfo.txt`：软件环境记录；
- `data/processed/cellchat_object.rds`：单样本 CellChat 对象；
- `data/processed/cellchat_objects_by_group.rds` 与 `cellchat_merged.rds`：双样本或多样本对象。

比较目录使用 `EXP_vs_CON` 形式命名。差异图按第二组相对第一组解释：红色表示第二组更强，蓝色表示第一组更强。

多样本分析会自动完成所有两两比较。如果某一对样本的通讯网络完全相同，脚本会跳过没有信息量的空差异热图，并继续输出差异网络、配体-受体气泡图和信息流比较图。
