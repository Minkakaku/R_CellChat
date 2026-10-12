# CellChat 独立分析入口
# 来源：R-cellchat-workflow-template 当前公共脚本；适配工作区 QS2、集中README及可恢复运行。
# single、two、multi 共用一个完整脚本，直接接收主 Seurat 流程输出的 RDS/QS2。

suppressPackageStartupMessages({
  library(CellChat)
  library(ComplexHeatmap)
  library(future)
  library(ggplot2)
  library(Seurat)
})

NPG_COLORS <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
npg_palette <- function(labels) {
  colors <- NPG_COLORS
  if (length(labels) > length(colors)) {
    extra_n <- length(labels) - length(colors)
    extra <- grDevices::colorRampPalette(NPG_COLORS)(2L * extra_n + 1L)[seq.int(2L, 2L * extra_n, by = 2L)]
    colors <- c(colors, extra)
  }
  stats::setNames(colors[seq_along(labels)], labels)
}

# 功能：读取 --name value 形式的命令行参数。
settings <- list()
arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match(name, args)
  key <- sub("^--", "", name)
  if (is.na(index) || index == length(args)) {
    if (!is.null(settings[[key]])) settings[[key]] else default
  } else args[[index + 1L]]
}

# 功能：解析逗号分隔的样本或条件名称。
split_values <- function(x) {
  if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

config_file <- arg_value("--config")
if (!is.null(config_file)) settings <- source(config_file, local = TRUE)$value
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
nboot <- as.integer(arg_value("--nboot", "100"))
normalize_rna <- as.logical(arg_value("--normalize", "FALSE"))
resume <- as.logical(arg_value("--resume", "FALSE"))
top_lr <- as.integer(arg_value("--top-lr", "30"))
celltypes <- split_values(arg_value("--celltypes", ""))
sender_networks <- as.logical(arg_value("--sender-networks", "FALSE"))
plot_stems <- split_values(arg_value("--plot-stems", ""))
write_note <- as.logical(arg_value("--write-note", "TRUE"))
session_path <- arg_value("--session-file", file.path(output_dir, "sessionInfo.txt"))
group_colors <- split_values(arg_value("--group-colors", ""))
stopifnot(min_cells >= 1L, workers >= 1L, nboot >= 1L, top_lr >= 1L)
options(future.globals.maxSize = 8 * 1024^3)

# 使用 CellChatDB 的三类主要配体-受体互作。
database_categories <- c(
  "Secreted Signaling",
  "ECM-Receptor",
  "Cell-Cell Contact"
)

# ------------------------------------------------------------------------------
# 2. 输出函数：所有主要图同时保存为 PDF 和 600 dpi PNG
# ------------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
table_dir <- output_dir
plot_dir <- output_dir

# 功能：把标签转换成可安全用于文件名的字符串（函数：safe_name）。
safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

# 功能：把 ggplot 同时保存为矢量 PDF 和 600 dpi PNG（函数：save_ggplot）。
save_ggplot <- function(plot, stem, width = 8, height = 6) {
  if (length(plot_stems) && !basename(stem) %in% plot_stems) return(invisible(NULL))
  plot <- plot + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, limitsize = FALSE, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600, limitsize = FALSE, bg = "white")
}

# CellChat 的气泡图只接受内置 palette 名；保留其概率尺度/图例，仅替换颜色函数。
style_cellchat_bubble <- function(plot) {
  colour_scale <- plot$scales$get_scales("colour")
  colour_scale$palette <- scales::gradient_n_pal(c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019"))
  plot + guides(color = guide_colorbar(display = "rectangles"))
}

# 功能：把基础绘图或 ComplexHeatmap 同时保存为 PDF 和 PNG（函数：save_draw）。
save_draw <- function(draw_fun, stem, width = 8, height = 7) {
  if (length(plot_stems) && !basename(stem) %in% plot_stems) return(invisible(NULL))
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw_fun()
  dev.off()

  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 600)
  draw_fun()
  dev.off()
}

# 两组网络完全相同时不画空差异热图，但保留其余比较结果。
has_network_difference <- function(cellchat, pair, measure) {
  difference <- cellchat@net[[pair[2]]][[measure]] - cellchat@net[[pair[1]]][[measure]]
  any(difference != 0)
}

# 保留一个发送细胞类型的原始边权，其他行置零，节点顺序及数值不变。
sender_network_matrix <- function(mat, sender) {
  stopifnot(is.matrix(mat), sender %in% rownames(mat), all(is.finite(mat)))
  selected <- mat
  selected[,] <- 0
  selected[sender,] <- mat[sender,]
  selected
}

# 长表按 source/target 精确对应数量和强度，供新增网络图复核与再绘。
communication_network_table <- function(net) {
  stopifnot(identical(dim(net$count), dim(net$weight)),
    identical(dimnames(net$count), dimnames(net$weight)))
  table <- expand.grid(source = rownames(net$count), target = colnames(net$count),
    stringsAsFactors = FALSE)
  table$count <- as.vector(net$count)
  table$weight <- as.vector(net$weight)
  table
}

# 各发送者使用同一组网络最大边权；PNG/PDF 复用现有输出函数。
plot_sender_networks <- function(cellchat, out) {
  labels <- rownames(cellchat@net$count)
  colors <- npg_palette(labels)
  sizes <- as.numeric(table(cellchat@idents)[labels])
  ncols <- min(3L, length(labels))
  nrows <- ceiling(length(labels) / ncols)
  write.csv(communication_network_table(cellchat@net),
    file.path(out, "sender_network_edges.csv"), row.names = FALSE)
  for (measure in c("count", "weight")) {
    mat <- cellchat@net[[measure]]
    if (max(mat) <= 0) next
    save_draw(function() {
      par(mfrow = c(nrows, ncols), mar = c(0.5, 0.5, 2, 0.5), xpd = TRUE)
      for (sender in labels) {
        netVisual_circle(sender_network_matrix(mat, sender), vertex.weight = sizes,
          color.use = colors, weight.scale = TRUE, edge.weight.max = max(mat),
          label.edge = FALSE, title.name = sender)
      }
    }, file.path(out, paste0("08_sender_network_", measure)), width = 5 * ncols, height = 5 * nrows)
  }
}

# ------------------------------------------------------------------------------
# 3. 单个数据集的 CellChat 推断
# ------------------------------------------------------------------------------

# 功能：对一个数据集完成 CellChat 配体-受体通讯推断（函数：run_cellchat）。
run_cellchat <- function(seurat_obj, dataset_name) {
  message("Running CellChat: ", dataset_name)
  label <- safe_name(dataset_name)
  cache_path <- file.path(object_dir, paste0("cellchat_", label, ".rds"))
  cache_key <- digest::digest(list(input = normalizePath(input_rds),
    source = file.info(input_rds)[, c("size", "mtime")], cells = colnames(seurat_obj),
    metadata = seurat_obj[[]], CellChat = as.character(packageVersion("CellChat")),
    parameters = list(species, assay, celltype_col, min_cells, prob_type, population_size,
      smooth_ppi, seed, nboot, normalize_rna, database_categories)))
  if (resume && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (!identical(cached@options$workspace_cache_key, cache_key)) stop("Cached inputs/parameters differ: ", cache_path)
    message("Reusing verified inference cache: ", dataset_name)
    plot_single_cellchat(cached, dataset_name)
    return(cached)
  }

  # 将真实细胞类型作为通讯网络节点，不使用技术性 cluster 编号。
  seurat_obj@meta.data[[celltype_col]] <- factor(seurat_obj@meta.data[[celltype_col]])
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = assay)
  if (!"samples" %in% names(seurat_obj@meta.data)) seurat_obj@meta.data$samples <- dataset_name
  seurat_obj@meta.data$samples <- factor(seurat_obj@meta.data$samples)
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
    population.size = population_size,
    nboot = nboot,
    seed.use = seed
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  cellchat@options$workspace_cache_key <- cache_key
  saveRDS(cellchat, cache_path)
  message("INFERENCE_OK ", dataset_name, ": ", ncol(seurat_obj), " cells; ",
    length(cellchat@netP$pathways), " pathways; ", sum(cellchat@net$count), " significant edges")

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

  if (sender_networks) plot_sender_networks(cellchat, out)

  celltypes <- levels(cellchat@idents)
  celltype_colors <- npg_palette(celltypes)
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
        title.name = NULL
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
        title.name = NULL
      )
    },
    file.path(out, "02_network_weight")
  )

  # 细胞类型之间通讯数量和强度热图。
  save_draw(
    function() ComplexHeatmap::draw(
      netVisual_heatmap(cellchat, measure = "count", color.use = celltype_colors,
        color.heatmap = c("#FFF5F0", "#E64B35"))
    ),
    file.path(out, "03_heatmap_count")
  )
  save_draw(
    function() ComplexHeatmap::draw(
      netVisual_heatmap(cellchat, measure = "weight", color.use = celltype_colors,
        color.heatmap = c("#FFF5F0", "#E64B35"))
    ),
    file.path(out, "04_heatmap_weight")
  )

  # 配体-受体气泡图。
  lr <- subsetCommunication(cellchat)
  ranked <- sort(tapply(lr$prob, lr$interaction_name, sum), decreasing = TRUE)
  selected <- data.frame(interaction_name = head(names(ranked), top_lr))
  bubble <- netVisual_bubble(cellchat, pairLR.use = selected, remove.isolate = TRUE, angle.x = 90,
    grid.on = FALSE, line.on = FALSE)
  bubble <- style_cellchat_bubble(bubble)
  save_ggplot(bubble, file.path(out, paste0("05_top", top_lr, "_ligand_receptor_bubble")), 26, 11)

  # 发送和接收信号的中心性散点图。
  role_scatter <- netAnalysis_signalingRole_scatter(cellchat, color.use = celltype_colors,
    do.label = FALSE) + ggrepel::geom_text_repel(ggplot2::aes(label = labels),
      color = "black", size = 3, seed = seed, max.overlaps = Inf,
      box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
      show.legend = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15)))
  save_ggplot(role_scatter, file.path(out, "06_signaling_role_scatter"), 7, 6)

  # 通路层面的 outgoing / incoming 信号热图。
  # CellChat内部height单位为cm；仅放大外层画布不会展开通路行。
  body_height_cm <- max(8, length(cellchat@netP$pathways) * 0.34)
  heatmap_height <- body_height_cm / 2.54 + 3.2
  for (pattern in c("outgoing", "incoming")) {
    save_draw(
      function() {
        ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = pattern,
          color.use = celltype_colors, color.heatmap = "Reds",
          width = 12, height = body_height_cm, font.size = 8)
        ht@matrix_color_mapping <- ComplexHeatmap::ColorMapping(
          name = "Relative strength", col_fun = circlize::colorRamp2(
            c(0, 0.5, 1), c("#FFF5F0", "#F39B7F", "#E64B35")))
        ComplexHeatmap::draw(ht)
      },
      file.path(out, paste0("07_signaling_role_", pattern)),
      width = 8,
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

  dataset_colors <- npg_palette(groups)
  if (length(group_colors)) {
    stopifnot(length(group_colors) == length(groups))
    dataset_colors <- setNames(group_colors, groups)
  }

  # 总通讯数量和总通讯强度比较。
  p_count <- compareInteractions(
    merged_cellchat,
    show.legend = FALSE,
    group = groups,
    group.levels = groups,
    measure = "count",
    color.use = unname(dataset_colors)
  )
  p_weight <- compareInteractions(
    merged_cellchat,
    show.legend = FALSE,
    group = groups,
    group.levels = groups,
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
        color.use = npg_palette(levels(cellchat_list[[reference]]@idents)),
        color.edge = c("#E64B35", "#3C5488"),
        measure = "count"
      ),
      file.path(pair_out, "03_difference_count")
    )
    save_draw(
      function() netVisual_diffInteraction(
        merged_cellchat,
        comparison = c(reference, experiment),
        weight.scale = TRUE,
        color.use = npg_palette(levels(cellchat_list[[reference]]@idents)),
        color.edge = c("#E64B35", "#3C5488"),
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
            color.use = npg_palette(levels(cellchat_list[[reference]]@idents)),
            color.heatmap = c("#3C5488", "#E64B35")
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
            color.use = npg_palette(levels(cellchat_list[[reference]]@idents)),
            color.heatmap = c("#3C5488", "#E64B35")
          )
        ),
        file.path(pair_out, "06_difference_heatmap_weight")
      )
    }

    # 比较两组的配体-受体互作和通路信息流。
    lr <- do.call(rbind, lapply(cellchat_list[pair], subsetCommunication))
    ranked <- sort(tapply(lr$prob, lr$interaction_name, sum), decreasing = TRUE)
    selected <- data.frame(interaction_name = head(names(ranked), top_lr))
    bubble <- netVisual_bubble(
      merged_cellchat,
      comparison = c(reference, experiment),
      pairLR.use = selected,
      remove.isolate = TRUE,
      angle.x = 90,
      color.text = unname(dataset_colors[pair]),
      grid.on = FALSE,
      line.on = FALSE
    )
    bubble <- style_cellchat_bubble(bubble)
    save_ggplot(bubble, file.path(pair_out, paste0("07_top", top_lr, "_ligand_receptor_bubble")), 40, 11)

    rank_count <- rankNet(
      merged_cellchat,
      mode = "comparison",
      comparison = c(reference, experiment),
      measure = "count",
      stacked = TRUE,
      do.stat = FALSE,
      color.use = unname(dataset_colors)
    )
    rank_weight <- rankNet(
      merged_cellchat,
      mode = "comparison",
      comparison = c(reference, experiment),
      measure = "weight",
      stacked = TRUE,
      do.stat = FALSE,
      color.use = unname(dataset_colors)
    )
    flow_height <- max(7, length(unique(unlist(lapply(cellchat_list[pair], function(x) x@netP$pathways)))) * 0.18 + 1)
    save_ggplot(rank_count, file.path(pair_out, "08_information_flow_count"), 8, flow_height)
    save_ggplot(rank_weight, file.path(pair_out, "09_information_flow_weight"), 8, flow_height)
  }
}

# ------------------------------------------------------------------------------
# 6. 主流程
# ------------------------------------------------------------------------------

# 仅重绘已完成的分组模型，不重读Seurat、不重新推断或改写模型。
if (identical(arg_value("--stage", "all"), "plot")) {
  stopifnot(length(samples) >= 2L)
  cellchat_list <- setNames(lapply(samples, function(group) {
    x <- readRDS(file.path(object_dir, paste0("cellchat_", safe_name(group), ".rds")))
    stopifnot(all(as.character(x@meta[[sample_col]]) == group),
      identical(x@options$mode, "single"), length(x@netP$pathways) > 0L)
    plot_single_cellchat(x, group)
    x
  }), samples)
  merged_cellchat <- readRDS(file.path(object_dir, "cellchat_merged.rds"))
  stopifnot(identical(names(merged_cellchat@net), samples))
  plot_cellchat_comparison(cellchat_list, merged_cellchat)
  message("CELLCHAT_REPLOT_OK: ", normalizePath(output_dir))
  quit(save = "no", status = 0L)
}

set.seed(seed)
# 大型 RNA 标准化保持串行，避免 Seurat 闭包向并行进程重复发送完整表达对象。
future::plan("sequential")

seurat_obj <- if (grepl("[.]qs2$", input_rds, ignore.case = TRUE)) qs2::qs_read(input_rds) else readRDS(input_rds)
stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))
# 可选筛选一个或多个真实细胞类型；不修改原始输入文件。
if (length(celltypes)) {
  missing <- setdiff(celltypes, as.character(seurat_obj@meta.data[[celltype_col]]))
  if (length(missing)) stop("Unknown cell types: ", paste(missing, collapse = ", "))
  keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] %in% celltypes]
  seurat_obj <- subset(seurat_obj, cells = keep)
}
# 只保留通讯所需 RNA，释放全局 SCT/integrated 大矩阵；只改内存副本。
DefaultAssay(seurat_obj) <- assay
seurat_obj@assays <- seurat_obj@assays[assay]
seurat_obj@reductions <- list(); seurat_obj@graphs <- list(); seurat_obj@neighbors <- list()
gc(verbose = FALSE)
seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = assay)
if (normalize_rna) seurat_obj <- NormalizeData(seurat_obj, assay = assay,
  normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
expression <- GetAssayData(seurat_obj, assay = assay, layer = "data")
stopifnot(nrow(expression) > 0, ncol(expression) == ncol(seurat_obj),
  all(is.finite(expression@x)), all(expression@x >= 0), max(expression@x) > 0)
rm(expression)
message("INPUT_OK: ", ncol(seurat_obj), " cells; assay=", assay, "; RNA LogNormalize=", normalize_rna)
print(table(seurat_obj[[celltype_col, drop = TRUE]], seurat_obj[[sample_col, drop = TRUE]]))
if (workers > 1L) future::plan("multisession", workers = workers)

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

  saveRDS(merged_cellchat, file.path(object_dir, "cellchat_merged.rds"))

  plot_cellchat_comparison(cellchat_list, merged_cellchat)

  result_note <- c(
    "# CellChat 分析结果",
    "",
    paste0("- 分析模式：", ifelse(analysis_mode == "two", "双样本比较", "多样本比较")),
    paste0("- 比较顺序：", paste(names(cellchat_list), collapse = " → ")),
    paste0("- 共同细胞类型：", paste(shared_celltypes, collapse = ", ")),
    "- 差异图命名为 第二组_vs_第一组；红色表示第二组增强，蓝色表示第一组增强。",
    paste0("- 分组对象：", paste(file.path(object_dir,
      paste0("cellchat_", safe_name(names(cellchat_list)), ".rds")), collapse = ", ")),
    paste0("- 合并对象：", file.path(object_dir, "cellchat_merged.rds")),
    "- 网络差异图是数据集间描述性比较；正式 CON/EXP 推断应在生物学样本层面复核。"
  )
}

if (write_note) writeLines(result_note, file.path(output_dir, "README_results_CN.md"))
capture.output(sessionInfo(), file = session_path)
future::plan("sequential")

message("CellChat results: ", normalizePath(output_dir))
message("CellChat objects: ", normalizePath(object_dir))
