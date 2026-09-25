################GO/KEGG富集分析
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
up<-intersect(dat$SLE_up,dat$MDD_up)
up<-data.frame(up)
write.xlsx(up, 'up.xlsx', rowNames = FALSE, sep = '\t')
down<-intersect(dat$SLE_down,dat$MDD_down)
down<-data.frame(down)
write.xlsx(down, 'down.xlsx', rowNames = FALSE, sep = '\t')
# install.packages("gProfileR")
# library("gProfileR")
# gprofiler(up, organism = "hsapiens")
library(openxlsx)
up<-read.xlsx('up.xlsx')
down<-read.xlsx('down.xlsx')
gene<-c(up$up,down$down)
# 安装必要包（如未安装）
# if (!require("BiocManager")) install.packages("BiocManager")
# options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/")) 
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# BiocManager:: install ("org.Hs.eg.db")
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(dplyr)
load("coDEGgene.RData")
# 读取基因列表（示例文件格式：每行一个gene symbol）
# gene_symbols <- data.frame(gene)
# colnames(gene_symbols)<-"SYMBOL"
gene_symbols <- as.character(gene)
# 基因ID转换
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
# 输出转换结果
cat("成功转换基因数:", nrow(gene_entrez), "\n")
cat("未转换成功的基因:", 
    setdiff(gene_symbols, gene_entrez$SYMBOL), "\n")

# GO分析（BP/MF/CC分别进行）
go_analysis <- function(ontology) {
  enrichGO(gene          = gene_entrez$ENTREZID,
           OrgDb         = org.Hs.eg.db,
           ont           = ontology,
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.2,
           readable      = TRUE)
}

ego_bp <- go_analysis("BP") %>% mutate(Ontology = "BP")
ego_mf <- go_analysis("MF") %>% mutate(Ontology = "MF")
ego_cc <- go_analysis("CC") %>% mutate(Ontology = "CC")
write.csv(ego_bp,"ego_bp.csv")
write.csv(ego_mf,"ego_mf.csv")
write.csv(ego_cc,"ego_cc.csv")
# 合并GO结果（各取前5个显著条目）
go_combined <- bind_rows(
  ego_bp@result %>% slice_min(p.adjust, n = 5),
  ego_mf@result %>% slice_min(p.adjust, n = 5),
  ego_cc@result %>% slice_min(p.adjust, n = 5)
) %>% 
  arrange(Ontology, p.adjust) %>%
  mutate(Term = paste0(Ontology, ": ", Description))

# 绘制GO柱状图
go_plot <- ggplot(go_combined, 
                  aes(x = reorder(Term, -log10(p.adjust)), 
                      y = -log10(p.adjust), 
                      fill = Ontology)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(BP = "#1F77B4", MF = "#FF7F0E", CC = "#2CA02C")) +
  labs(title = "GO Enrichment Analysis (Top5 per Category)",
       x = "GO Terms",
       y = "-log10(Adjusted p-value)") +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right"
  )


# KEGG分析
kegg_res <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)
# 结果预处理
kegg_df <- as.data.frame(kegg_res) %>%
  arrange(p.adjust) %>%       # 按校正p值排序
  mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))  # 计算富集因子
library(ggrepel)
# 可视化1：经典气泡图
p1 <- ggplot(kegg_df %>% head(15), 
             aes(x = richFactor, 
                 y = reorder(Description, richFactor),
                 size = Count,
                 color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", 
                       name = "-log10(adj.p)") +
  scale_size_continuous(range = c(3, 8), 
                        name = "Gene Count") +
  labs(title = "KEGG Pathway Enrichment", 
       x = "Rich Factor", 
       y = "") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        legend.position = "right") +
  geom_text_repel(aes(label = Description),  # 自动调整标签位置
                  size = 4, 
                  box.padding = 0.5)
p1
# 可视化2：柱状图（按富集因子排序）
p2 <- ggplot(kegg_df %>% head(15), 
             aes(x = reorder(Description, -p.adjust), 
                 y = -log10(p.adjust),
                 fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  coord_flip() +  # 横向显示
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", 
                      name = "-log10(adj.p)") +
  labs(title = "", 
       x = "", 
       y = "-log10(Adjusted p-value)") +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
        legend.position = "none")
p2
# 可视化3：通路图（以第一个显著通路为例）
library(pathview)
hsa_pathway_id <- sub("hsa", "", kegg_df$ID[1])  # 提取通路编号
p<-pathview(gene.data = gene_entrez$ENTREZID,
         pathway.id = hsa_pathway_id,
         species = "hsa",
         kegg.native = TRUE,  # 显示KEGG原生风格
         key.pos = "topright")
hsa_pathway_id <- sub("hsa", "", kegg_df$ID[2])  # 提取通路编号
p<-pathview(gene.data = gene_entrez$ENTREZID,
            pathway.id = hsa_pathway_id,
            species = "hsa",
            kegg.native = TRUE,  # 显示KEGG原生风格
            key.pos = "topright")
# 保存结果
ggsave("kegg_bubble.pdf", p1, width = 10, height = 8, dpi = 300)
ggsave("kegg_bar.pdf", p2, width = 10, height = 6, dpi = 300)
write.csv(kegg_df, "kegg_results.csv")

# 显示图形
print(go_plot)
print(kegg_plot)

# 保存结果
ggsave("GO_Top5.pdf", go_plot, width = 10, height = 8)
ggsave("KEGG_dotplot.pdf", kegg_plot, width = 10, height = 8)

# 保存分析结果
write.csv(go_combined, "GO_Top5_results.csv")
write.csv(kegg_res@result, "KEGG_results.csv")


###########疾病富集分析
# 安装必要的包（如果尚未安装）
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install(c("DOSE", "clusterProfiler", "org.Hs.eg.db", "ggplot2"))

# # 安装必要的包（如果尚未安装）
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("DOSE", "clusterProfiler", "org.Hs.eg.db", "ggplot2"),force = TRUE)
BiocManager::install( "clusterProfiler",force = TRUE)
# 加载所需的库
# 加载包
# library(DOSE)#不加载才能运行DO分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# DO富集分析
do_result <- enrichDO(gene  = gene_entrez$ENTREZID,
                      ont   = "HDO",           # 固定参数
                      pvalueCutoff  = 0.05,           # p值阈值
                      pAdjustMethod = "BH",
                      readable      = TRUE)           # 显示基因符号

# 结果可视化
if (nrow(do_result) > 0) {
  # 点图可视化
  p <- dotplot(do_result, 
               showCategory = 10,        # 显示前15个显著条目
               font.size = 12,
               title = "Disease Ontology Enrichment") +
    scale_color_gradient(low = "red", high = "blue", 
                         name = "p.adjust") +       # 修改颜色映射
    theme_bw() +
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
  
  # 保存结果
  write.csv(do_result@result, "DO_analysis_results.csv")
  ggsave("DO_dotplot.pdf", p, width = 10, height = 8)
  
} else {
  message("没有显著富集的疾病条目，建议：")
  message("1. 检查基因是否与疾病相关")
  message("2. 放宽pvalueCutoff阈值（当前阈值：0.05）")
}

