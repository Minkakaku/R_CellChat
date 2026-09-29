# 项目分层

公共方法在仓库 src/R 或 src/python；公共默认在根 config；项目参数与人工表格在 projects/<项目>/config；实际执行脚本在项目 workflow。编号步骤仅加载配置、调用共享方法与写日志。中间对象在项目 data/processed，结果在项目 results，环境在 envs。

新项目复制 projects/project1，再按仓库 README 调整配置和运行。不要修改公共代码中的路径来适配项目。保留已经执行的步骤，重跑新增编号和输出目录。
