if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,       # 数据处理
  httr,            # API请求
  jsonlite,        # JSON解析
  visNetwork,      # 交互式网络图
  clusterProfiler, # 富集分析
  org.Hs.eg.db,    # 人类基因数据库
  progress,        # 进度条
  DT               # 交互式表格
)

#DGidb下载tsv文件
library(data.table)
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure6/")
data <- fread("gene_interaction_results.tsv", sep = "\t")
drug_data<-data[data$`interaction score`>=1]
table(drug_data$gene)
write.csv(drug_data,"drug-gene_interaction.csv")
# 网络图可视化
nodes <- data.frame(
  id = unique(c(drug_data$gene, drug_data$drug)),
  label = unique(c(drug_data$gene, drug_data$drug)),
  group = ifelse(
    unique(c(drug_data$gene, drug_data$drug)) %in% drug_data$gene, 
    "Gene", 
    "Drug"
  )
)
library(dplyr)
edges <- data.frame(from=drug_data$gene, to = drug_data$drug)
visNetwork(nodes, edges) %>%
  visGroups(groupname = "Gene", color = "#e41a1c", shape = "diamond") %>%
  visGroups(groupname = "Drug", color = "#377eb8", shape = "dot") %>%
  visEdges(arrows = "to") %>%
  visOptions(
    highlightNearest = TRUE,
    nodesIdSelection = TRUE
  ) %>%
  visLayout(randomSeed = 123)  # 固定布局保证可重复性




