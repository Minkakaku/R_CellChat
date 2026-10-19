# R_CellChat

跨项目复用的 R 公共分析代码库。仓库管理公共方法、默认配置、环境和项目示例；实际项目统一在 `/workspace/projects/<项目名>/` 使用这些方法，R/Python 共用同一项目，数据和运行产物不提交 Git。

## 集中项目规范

每个项目在 `logs/README.md` 维护唯一的项目总记录，汇总项目概况并按时间追加所有操作，包括检查、依赖与代码变更、参数调整、人工决策、运行、失败、重跑和结果核对；每条注明时间与时区、目的、脚本或命令、输入与参数、代码和环境版本、状态、日志及固定结果路径。根目录 README 作为入口和运行索引，链接到总记录。总记录纳入项目本地 Git，原始运行日志忽略；项目内容不上传 GitHub。

每个实际项目仅在本地拥有独立 Git 仓库，用于本地版本记录；只有公共代码库放在 GitHub，项目代码、配置和数据保留在本地。初始化时按需从 GitHub 以子模块检出公共库到项目 `deps/`，通过 `.gitmodules` 记录来源、子模块指针固定提交。从本地项目或备份恢复时获取已记录的依赖版本，分析启动不自动拉取最新版；升级依赖须单独验证和记录。

同一实际项目只创建一次，统一维护样本表、分组、人工注释、配色、编号 R/Python 步骤和分析结果；方法参数可按分析点拆分。公共源码继续由各代码库维护。完整规则见本库 [AGENTS.md](AGENTS.md) 与工作区 `/workspace/AGENTS.md`。

当前加载器和步骤仍按仓库内目录定位。集中项目接入须先适配显式代码仓库路径与共用配置，再建立统一模板、验证路径与跨语言数据交接，最后迁移实际项目。下列 `projects/project1/` 和复制命令仅用于仓库内兼容示例，不是集中项目的创建方法。

## 目录

```text
src/R/       公共分析方法与配置加载器
config/            公共默认参数与参考配置
envs/              共享环境定义
projects/project1/ 仓库内兼容示例（未运行）
  config/          项目覆盖、样本表、QC/注释等人工决策
  workflow/        编号短脚本与分阶段 runner
  data/processed/  中间对象
  data/reference/  项目参考资源
  results/         按分析点保存图表与中文说明
  logs/            原始运行日志
    README.md      项目总记录（全部操作历史）
tests/             路径和项目配置检查
```

## 仓库内兼容示例

以下命令仅用于验证旧接口，复制出的目录保持在本仓库内。正式项目按上述集中项目规范接入。

在仓库根目录执行：

```bash
conda env create -f envs/environment.yaml
# 按 CellChat 官方说明安装 CellChat 并记录版本。
cp -R projects/project1 projects/compat_check
```

修改 `projects/compat_check/config/config.yaml` 和样本/注释表，准备输入对象与参考文件。示例中的物种、reference、root、分组和文件名必须按真实项目核对，不能直接作为生物学结论。

每个编号脚本可单独执行；脚本从自身位置定位项目，因此复制后不会继续写入 project1。Python 包项目先完成 editable install。

```bash
Rscript projects/compat_check/workflow/10_cellchat.R
```

可用步骤：
- `10_cellchat.R`

`run_all.sh` 对主流程只运行到 QC 检查点；独立分析方法需显式选择步骤，避免把不同候选方法混为同一条流程。已运行脚本不改写；重跑可复用原步骤，确需修改时新增脚本并记录来源。`results/` 的编号固定对应分析点，由项目配置维护；同一分析点重跑、调整参数或补图仍使用原目录并更新对应产物，不按运行次数新增结果目录。仅新增分析点时分配新的结果编号。运行批次、参数与代码版本追加到 `logs/README.md`，原始运行日志另存 `logs/`。

## 配置与路径

以下为仓库内兼容示例的配置规则。集中项目还需加载共用样本、注释、配色与方法配置，代码库根目录显式指定；接入前按集中项目规范完成适配。

公共 `config/default.yaml` → `config/reference.yaml` → 项目 `config/config.yaml`，后者覆盖前者。嵌套映射按键合并，列表和标量整体替换。共享方法接收配置/参数，不包含实际项目路径与样本名。项目相对输入输出路径以项目目录为基准；参考文件路径写在配置中。

## 结果与解释

最终图保存同名 PNG/PDF，配套 CSV/TSV 和 `README_results_CN.md` 放在分析点结果目录。中间对象放在 `data/processed/`。细胞或空间 bins 不作为独立生物学重复；条件比较用样本级汇总。方法特有约束见 [AGENTS.md](AGENTS.md) 和公共代码注释。

## 验证

```bash
Rscript tests/test_project.R
```

此检查验证配置、路径或安全导入，不运行真实组学分析。完整分析需要对应生信依赖、真实输入数据和人工 QC/注释。

检查独立入口参数而不运行分析：`Rscript projects/compat_check/workflow/10_cellchat.R --dry-run`。
