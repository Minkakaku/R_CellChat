script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
workspace <- normalizePath(file.path(project, "../.."))
config <- list()
for (path in c(file.path(workspace, "config/default.yaml"), file.path(workspace, "config/reference.yaml"), file.path(project, "config/config.yaml"))) {
  config <- utils::modifyList(config, yaml::read_yaml(path), keep.null = TRUE)
}
for (key in c("input_rds", "output_dir", "object_dir")) {
  if (!grepl("^/", config[[key]])) config[[key]] <- file.path(project, config[[key]])
}
if ("--dry-run" %in% commandArgs(TRUE)) {
  print(config)
} else {
  source(file.path(workspace, "src/R/cellchat.R"))
  run_cellchat_analysis(config)
}
