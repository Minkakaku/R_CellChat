# CellChat 独立分析入口
# single、two、multi 共用一个完整脚本，直接接收主 Seurat 流程输出的 RDS。

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(future)
  library(ggplot2)
  library(Seurat)
})

# 功能：读取 --name value 形式的命令行参数。
arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match(name, args)
  if (is.na(index) || index == length(args)) default else args[[index + 1L]]
}

# 功能：解析逗号分隔的样本或条件名称。
split_values <- function(x) {
  if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

input_rds <- arg_value("--input")
output_dir <- arg_value("--output")
analysis_mode <- arg_value("--mode", "single")
if (is.null(input_rds) || is.null(output_dir)) {
  stop("用法：Rscript cellchat_workflow.R --input annotated.rds --output results/cellchat --mode two --sample-col group --samples CON,EXP")
}
object_dir <- arg_value("--object-dir", output_dir)
species <- arg_value("--species", "human")
assay <- arg_value("--assay", "RNA")
celltype_col <- arg_value("--celltype-col", "cell_type")
sample_col <- arg_value("--sample-col", "group")
samples <- split_values(arg_value("--samples", ""))
min_cells <- as.integer(arg_value("--min-cells", "10"))
workers <- as.integer(arg_value("--workers", "4"))
prob_type <- arg_value("--prob-type", "triMean")
population_size <- as.logical(arg_value("--population-size", "TRUE"))
smooth_ppi <- as.logical(arg_value("--smooth-ppi", "FALSE"))
seed <- as.integer(arg_value("--seed", "123"))

# 使用 CellChatDB 的三类主要配体-受体互作。
database_categories <- c(
  "Secreted Signaling",
  "ECM-Receptor",
  "Cell-Cell Contact"
)

# ------------------------------------------------------------------------------
# 2. 输出函数：所有主要图同时保存为 PDF 和 300 dpi PNG
# ------------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- output_dir
plot_dir <- output_dir

# 功能：把标签转换成可安全用于文件名的字符串（函数：safe_name）。
safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

# 功能：把 ggplot 同时保存为矢量 PDF 和 300 dpi PNG（函数：save_ggplot）。
save_ggplot <- function(plot, stem, width = 8, height = 6) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300)
}

# 功能：把基础绘图或 ComplexHeatmap 同时保存为 PDF 和 PNG（函数：save_draw）。
save_draw <- function(draw_fun, stem, width = 8, height = 7) {
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw_fun()
  dev.off()

  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300)
  draw_fun()
  dev.off()
}

# 两组网络完全相同时不画空差异热图，但保留其余比较结果。
has_network_difference <- function(cellchat, pair, measure) {
  difference <- cellchat@net[[pair[2]]][[measure]] - cellchat@net[[pair[1]]][[measure]]
  any(difference != 0)
}

# ------------------------------------------------------------------------------
# 3. 单个数据集的 CellChat 推断
# ------------------------------------------------------------------------------

# 功能：对一个数据集完成 CellChat 配体-受体通讯推断（函数：run_cellchat）。
run_cellchat <- function(seurat_obj, dataset_name) {
  message("Running CellChat: ", dataset_name)

  # 将真实细胞类型作为通讯网络节点，不使用技术性 cluster 编号。
  seurat_obj@meta.data[[celltype_col]] <- factor(seurat_obj@meta.data[[celltype_col]])
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = assay)
  seurat_obj@meta.data$samples <- factor(dataset_name)
  Idents(seurat_obj) <- seurat_obj@meta.data[[celltype_col]]

  cellchat <- createCellChat(
    object = seurat_obj,
    group.by = celltype_col,
    assay = assay
  )

  database <- switch(
    tolower(species),
    human = CellChatDB.human,
    mouse = CellChatDB.mouse,
    stop("species 只能是 human 或 mouse")
  )

  cellchat@DB <- subsetDB(
    database,
    search = database_categories,
    key = "annotation"
  )

  # 识别高表达配体、受体及其互作。
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # 可选：将表达投射到蛋白互作网络。默认关闭，保留原始表达推断。
  if (smooth_ppi) {
    ppi <- if (tolower(species) == "human") PPI.human else PPI.mouse
    cellchat <- smoothData(cellchat, adj = ppi)
  }

  # 计算细胞群之间的通讯概率，并过滤细胞数过少的互作。
  cellchat <- computeCommunProb(
    cellchat,
    type = prob_type,
    raw.use = !smooth_ppi,
    population.size = population_size
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  label <- safe_name(dataset_name)

  # 导出配体-受体层面和通路层面的通讯表。
  write.csv(
    subsetCommunication(cellchat),
    file.path(table_dir, paste0(label, "_ligand_receptor.csv")),
    row.names = FALSE
  )
  write.csv(
    subsetCommunication(cellchat, slot.name = "netP"),
    file.path(table_dir, paste0(label, "_pathway.csv")),
    row.names = FALSE
  )

  plot_single_cellchat(cellchat, dataset_name)
  cellchat
}

# ------------------------------------------------------------------------------
# 4. 单个 CellChat 对象的主要图
# ------------------------------------------------------------------------------

# 功能：输出单个 CellChat 对象的网络、热图、气泡图和信号角色图（函数：plot_single_cellchat）。
plot_single_cellchat <- function(cellchat, dataset_name) {
  label <- safe_name(dataset_name)
  out <- file.path(plot_dir, label)
  dir.create(out, showWarnings = FALSE)

  celltypes <- levels(cellchat@idents)
  celltype_colors <- setNames(hcl.colors(length(celltypes), "Dark 3"), celltypes)
  group_size <- as.numeric(table(cellchat@idents))

  # 全局通讯数量网络。
  save_draw(
    function() {
      par(mar = c(1, 1, 3, 1), xpd = TRUE)
      netVisual_circle(
        cellchat@net$count,
        vertex.weight = group_size,
        color.use = celltype_colors,
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0(dataset_name, " - interaction count")
      )
    },
    file.path(out, "01_network_count")
  )

  # 全局通讯强度网络。
  save_draw(
    function() {
      par(mar = c(1, 1, 3, 1), xpd = TRUE)
      netVisual_circle(
        cellchat@net$weight,
        vertex.weight = group_size,
        color.use = celltype_colors,
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0(dataset_name, " - interaction strength")
      )
    },
    file.path(out, "02_network_weight")
  )

  # 细胞类型之间通讯数量和强度热图。
  save_draw(
    function() ComplexHeatmap::draw(
      netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")
    ),
    file.path(out, "03_heatmap_count")
  )
  save_draw(
    function() ComplexHeatmap::draw(
      netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
    ),
    file.path(out, "04_heatmap_weight")
  )

  # 配体-受体气泡图。
  bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
  save_ggplot(bubble, file.path(out, "05_ligand_receptor_bubble"), 10, 8)

  # 发送和接收信号的中心性散点图。
  role_scatter <- netAnalysis_signalingRole_scatter(cellchat)
  save_ggplot(role_scatter, file.path(out, "06_signaling_role_scatter"), 7, 6)

  # 通路层面的 outgoing / incoming 信号热图。
  heatmap_height <- min(18, max(6, length(cellchat@netP$pathways) * 0.24))
  for (pattern in c("outgoing", "incoming")) {
    save_draw(
      function() ComplexHeatmap::draw(
        netAnalysis_signalingRole_heatmap(cellchat, pattern = pattern)
      ),
      file.path(out, paste0("07_signaling_role_", pattern)),
      width = 9,
      height = heatmap_height
    )
  }
}

# ------------------------------------------------------------------------------
# 5. 多组对象合并与比较
# ------------------------------------------------------------------------------

# 功能：比较双样本或多样本 CellChat 网络并明确比较方向（函数：plot_cellchat_comparison）。
plot_cellchat_comparison <- function(cellchat_list, merged_cellchat) {
  groups <- names(cellchat_list)
  out <- file.path(plot_dir, "comparison")
  dir.create(out, showWarnings = FALSE)

  dataset_colors <- setNames(hcl.colors(length(groups), "Dark 3"), groups)
  if ("CON" %in% groups) dataset_colors["CON"] <- "#72b1b9"
  if ("EXP" %in% groups) dataset_colors["EXP"] <- "#f6ae00"

  # 总通讯数量和总通讯强度比较。
  p_count <- compareInteractions(
    merged_cellchat,
    show.legend = TRUE,
    group = seq_along(groups),
    measure = "count",
    color.use = unname(dataset_colors)
  )
  p_weight <- compareInteractions(
    merged_cellchat,
    show.legend = TRUE,
    group = seq_along(groups),
    measure = "weight",
    color.use = unname(dataset_colors)
  )
  save_ggplot(p_count, file.path(out, "01_total_interaction_count"), 7, 6)
  save_ggplot(p_weight, file.path(out, "02_total_interaction_weight"), 7, 6)

  # 两组时只比较一次；三组及以上时自动生成所有两两比较。
  pairs <- combn(seq_along(groups), 2, simplify = FALSE)
  for (pair in pairs) {
    reference <- pair[1]
    experiment <- pair[2]
    pair_name <- paste0(groups[experiment], "_vs_", groups[reference])
    pair_out <- file.path(out, safe_name(pair_name))
    dir.create(pair_out, showWarnings = FALSE)

    # 差异网络按“第二组相对第一组”解释；例如 EXP_vs_CON。
    save_draw(
      function() netVisual_diffInteraction(
        merged_cellchat,
        comparison = c(reference, experiment),
        weight.scale = TRUE,
        measure = "count"
      ),
      file.path(pair_out, "03_difference_count")
    )
    save_draw(
      function() netVisual_diffInteraction(
        merged_cellchat,
        comparison = c(reference, experiment),
        weight.scale = TRUE,
        measure = "weight"
      ),
      file.path(pair_out, "04_difference_weight")
    )

    # 差异热图使用蓝-白-红配色：红色为第二组更强，蓝色为第一组更强。
    if (has_network_difference(merged_cellchat, pair, "count")) {
      save_draw(
        function() ComplexHeatmap::draw(
          netVisual_heatmap(
            merged_cellchat,
            comparison = c(reference, experiment),
            measure = "count",
            color.heatmap = c("#2166AC", "#B2182B")
          )
        ),
        file.path(pair_out, "05_difference_heatmap_count")
      )
    }
    if (has_network_difference(merged_cellchat, pair, "weight")) {
      save_draw(
        function() ComplexHeatmap::draw(
          netVisual_heatmap(
            merged_cellchat,
            comparison = c(reference, experiment),
            measure = "weight",
            color.heatmap = c("#2166AC", "#B2182B")
          )
        ),
        file.path(pair_out, "06_difference_heatmap_weight")
      )
    }

    # 比较两组的配体-受体互作和通路信息流。
    bubble <- netVisual_bubble(
      merged_cellchat,
      comparison = c(reference, experiment),
      remove.isolate = FALSE,
      angle.x = 45
    )
    save_ggplot(bubble, file.path(pair_out, "07_ligand_receptor_bubble"), 11, 8)

    rank_count <- rankNet(
      merged_cellchat,
      mode = "comparison",
      comparison = c(reference, experiment),
      measure = "count",
      stacked = TRUE,
      do.stat = TRUE
    )
    rank_weight <- rankNet(
      merged_cellchat,
      mode = "comparison",
      comparison = c(reference, experiment),
      measure = "weight",
      stacked = TRUE,
      do.stat = TRUE
    )
    save_ggplot(rank_count, file.path(pair_out, "08_information_flow_count"), 8, 7)
    save_ggplot(rank_weight, file.path(pair_out, "09_information_flow_weight"), 8, 7)
  }
}

# ------------------------------------------------------------------------------
# 6. 主流程
# ------------------------------------------------------------------------------

set.seed(seed)
future::plan("multisession", workers = workers)

seurat_obj <- readRDS(input_rds)
stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

analysis_mode <- match.arg(tolower(analysis_mode), c("single", "two", "multi"))

if (analysis_mode == "single") {
  # 单样本模式：samples 为空时使用整个 Seurat 对象；填写 1 个名称时只分析该样本。
  dataset_name <- "All"
  if (!is.null(sample_col) && length(samples) == 1) {
    stopifnot(sample_col %in% colnames(seurat_obj@meta.data))
    keep_cells <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[sample_col]] == samples[1]
    ]
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
    dataset_name <- samples[1]
  }

  cellchat_list <- setNames(list(run_cellchat(seurat_obj, dataset_name)), dataset_name)
  saveRDS(cellchat_list[[1]], file.path(object_dir, "cellchat_object.rds"))

  result_note <- c(
    "# CellChat 分析结果",
    "",
    "- 分析模式：单样本",
    paste0("- 样本：", dataset_name),
    paste0("- 细胞类型列：", celltype_col),
    paste0("- 配体-受体表：", safe_name(dataset_name), "_ligand_receptor.csv"),
    paste0("- CellChat 对象：", file.path(object_dir, "cellchat_object.rds")),
    "- CellChat 的概率和置换显著性基于细胞层面；不能代替生物学重复层面的条件检验。"
  )
} else {
  # 双样本和多样本模式共用同一套推断；multi 会自动生成全部两两比较。
  stopifnot(sample_col %in% colnames(seurat_obj@meta.data))
  if (analysis_mode == "two") stopifnot(length(samples) == 2)
  if (analysis_mode == "multi") stopifnot(length(samples) >= 3)

  keep_cells <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data[[sample_col]] %in% samples
  ]
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  seurat_list <- SplitObject(seurat_obj, split.by = sample_col)[samples]

  # 只比较每组中细胞数均达到阈值的共同细胞类型。
  shared_celltypes <- Reduce(
    intersect,
    lapply(seurat_list, function(obj) {
      counts <- table(obj@meta.data[[celltype_col]])
      names(counts[counts >= min_cells])
    })
  )

  seurat_list <- lapply(seurat_list, function(obj) {
    cells <- rownames(obj@meta.data)[obj@meta.data[[celltype_col]] %in% shared_celltypes]
    subset(obj, cells = cells)
  })

  cellchat_list <- Map(run_cellchat, seurat_list, names(seurat_list))
  names(cellchat_list) <- names(seurat_list)

  merged_cellchat <- mergeCellChat(
    cellchat_list,
    add.names = names(cellchat_list),
    cell.prefix = TRUE
  )

  saveRDS(cellchat_list, file.path(object_dir, "cellchat_objects_by_group.rds"))
  saveRDS(merged_cellchat, file.path(object_dir, "cellchat_merged.rds"))

  combined_communication <- do.call(
    rbind,
    lapply(names(cellchat_list), function(group) {
      data <- subsetCommunication(cellchat_list[[group]])
      data$dataset <- group
      data
    })
  )
  write.csv(
    combined_communication,
    file.path(table_dir, "all_groups_ligand_receptor.csv"),
    row.names = FALSE
  )

  plot_cellchat_comparison(cellchat_list, merged_cellchat)

  result_note <- c(
    "# CellChat 分析结果",
    "",
    paste0("- 分析模式：", ifelse(analysis_mode == "two", "双样本比较", "多样本比较")),
    paste0("- 比较顺序：", paste(names(cellchat_list), collapse = " → ")),
    paste0("- 共同细胞类型：", paste(shared_celltypes, collapse = ", ")),
    "- 差异图命名为 第二组_vs_第一组；红色表示第二组增强，蓝色表示第一组增强。",
    paste0("- 分组对象：", file.path(object_dir, "cellchat_objects_by_group.rds")),
    paste0("- 合并对象：", file.path(object_dir, "cellchat_merged.rds")),
    "- 网络差异图是数据集间描述性比较；正式 CON/EXP 推断应在生物学样本层面复核。"
  )
}

writeLines(result_note, file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
future::plan("sequential")

message("CellChat results: ", normalizePath(output_dir))
message("CellChat objects: ", normalizePath(object_dir))
