# R-cellchat-workflow-template

独立 CellChat 分析点。正式代码只有根目录 cellchat_workflow.R，直接接收主 Seurat 流程生成的注释后 RDS/QS2。

## 模式

- single：整个对象或指定一个数据集的通讯网络。
- two：两个数据集的网络推断和直接比较。
- multi：三个及以上数据集，并生成全部两两比较。

## 运行

单数据集：

    Rscript cellchat_workflow.R --input annotated.rds --output results/cellchat --mode single --species human --celltype-col cell_type

双数据集：

    Rscript cellchat_workflow.R --input annotated.rds --output results/cellchat --mode two --sample-col group --samples CON,EXP

多数据集：

    Rscript cellchat_workflow.R --input annotated.rds --output results/cellchat --mode multi --sample-col orig.ident --samples S1,S2,S3

可通过 --object-dir 指定 CellChat RDS；默认和表格、图、README_results_CN.md 一起放在输出目录。最终图同时保存 PDF 和 600 dpi PNG。

## 解释边界

CellChat 根据配体/受体表达和数据库先验推断通讯概率。两组或多组网络图是数据集间描述性比较；细胞层面置换或概率不能替代生物学重复层面的 CON/EXP 推断。若有多个生物学样本，应分别运行或在样本层面复核网络指标。

方法依据：CellChat 官方教程与论文（PMID: 33597522）。

## 配置、缓存与重绘

`--config settings.R` 读取返回 list 的 R 配置，命令行值优先。`--nboot`/`--seed` 控制推断，`--normalize TRUE` 可显式标准化 RNA；`--resume TRUE` 检查输入与参数缓存键后复用模型。QS2 输入需要 qs2，缓存需要 digest，绘图需要 scales/ggrepel/circlize 及已有 CellChat 依赖。

`--stage plot --samples CON,EXP --object-dir <缓存目录>` 只读取 `cellchat_<组>.rds` 和 `cellchat_merged.rds` 重绘，不重做推断；仍需填写 `--input` 与 `--output`。`--top-lr` 默认 30，`--plot-stems` 可选图名，`--group-colors` 按 samples 顺序传逗号分隔颜色。`--write-note FALSE` 适合集中维护项目说明，`--session-file` 指定会话记录路径。

与早期模板的输出差异：保留每组模型、合并模型、每组 LR/通路 CSV；不再额外输出模型列表 RDS 或所有组拼接 LR CSV。默认数据集配色使用 NPG；需要 CON/EXP 固定颜色时显式传 `--group-colors '#E64B35,#4DBBD5'`。rankNet 关闭细胞层面的组间显著性检验，仍只作描述性网络比较。

正式气泡点图同时关闭内部网格与比较组间虚线（`grid.on = FALSE`、`line.on = FALSE`）。原项目已执行脚本、已记录哈希及历史结果保留。
