# R-cellchat-workflow-template

独立 CellChat 分析点。正式代码只有根目录 cellchat_workflow.R，直接接收主 Seurat 流程生成的注释后 RDS。

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
