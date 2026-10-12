# 验证发送者矩阵及绘图支撑表的方向、数值与节点顺序。
args <- commandArgs(trailingOnly = TRUE)
script <- if (length(args)) args[[1]] else "cellchat_workflow.R"
definitions <- parse(script)
selected <- c("sender_network_matrix", "communication_network_table")
for (expr in definitions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% selected) eval(expr)
}
count <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("T", "B"), c("T", "B")))
weight <- count / 10
sender <- sender_network_matrix(count, "B")
stopifnot(identical(sender["B",], count["B",]), all(sender["T",] == 0),
          identical(dimnames(sender), dimnames(count)))
table <- communication_network_table(list(count = count, weight = weight))
stopifnot(table$count[table$source == "B" & table$target == "T"] == 2,
          table$weight[table$source == "T" & table$target == "B"] == 0.3)
bad <- weight
colnames(bad) <- rev(colnames(bad))
stopifnot(inherits(try(communication_network_table(list(count=count, weight=bad)), silent=TRUE), "try-error"))
message("SENDER_NETWORK_TESTS_OK")
