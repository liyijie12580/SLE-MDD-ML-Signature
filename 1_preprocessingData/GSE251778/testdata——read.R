
#############GSE251778#####数据下载有问题
# #特征明确的 MDD （n=80） 和健康对照 （n=89） 个体队列的外周血样本中
setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/test/")
library(GEOquery)
library(limma)
library(affy)
# gset <- getGEO('GSE251778', destdir="./",
#                AnnotGPL = T,     ## 注释文件
#                getGPL = T)       ## 平台文件## 平台文件
# ## 获取临床信息
# cli1<-pData(gset[[1]])
# cli2<-pData(gset[[2]])
# cli<-rbind(cli1,cli2)
# meta2<-data.frame(cli$geo_accession,cli$title)
# colnames(meta2)<-c("name","case")
# library(tidyr)
# meta2_split <- meta2 %>%
#   separate(
#     col = case,           # 要拆分的列名
#     into = c("tissue", "case", "gender","other"),  # 新列名
#     sep = "-",            # 分隔符
#     remove = TRUE         # 是否删除原始列
#   )
# print(meta2_split)
# meta2_split$group<-ifelse(meta2_split$case == "Control", 0, 1)
# table(meta2_split$group)
# write.csv(meta2_split, file = "GSE251778_meta.csv", row.names = T,quote = F)
## 获取表达矩阵
library(variancePartition)
library(DESeq2)
library(dplyr)
library(ggplot2)
################
###   Data   ###
################
setwd("./GSE289146/")
data = readRDS('data.mRNA.rds')
raw.data.Female = data$raw.data.Female
raw.data.Male = data$raw.data.Male
cov.Female = data$cov.Female
cov.Male = data$cov.Male

rm(data)

###################
###   Quality   ###
###################

Male_CTRL = (row.names(cov.Male[cov.Male$Group == "Control",]))
Male_MDD = (row.names(cov.Male[cov.Male$Group == "Patient",]))

Female_CTRL = (row.names(cov.Female[cov.Female$Group == "Control",]))
Female_MDD = (row.names(cov.Female[cov.Female$Group == "Patient",]))

Male_filter = raw.data.Male[which(rowMeans(raw.data.Male) > 10), ]
Female_filter = raw.data.Female[which(rowMeans(raw.data.Female) > 10), ]

cov.Male = cov.Male[rownames(cov.Male) %in% colnames(Male_filter),]
cov.Female = cov.Female[rownames(cov.Female) %in% colnames(Female_filter),]

library(dplyr)
exp1<-Male_filter
exp2<-Female_filter
MDD<-merge(exp1,exp2,by="row.names",all = T)
rownames(MDD)<-MDD[,1]
MDD<-MDD[,-1]
MDD[is.na(MDD)]<-0
summary(colSums(MDD))
cov.Male$gender<-"M"
cov.Female$gender<-"F"
rownames(cov.Female)<-paste0("X",rownames(cov.Female))
meta<-rbind(cov.Male,cov.Female)
meta$group<-ifelse(meta$Group == "Control", 0, 1)
colnames(meta)
meta1 <- meta %>%
  dplyr::select("Age", "BMI","group","gender" )
meta1$Cohort2<-"GSE251778"
meta1$Cohort<-paste0("GSE251778_",meta1$gender)
setwd("./new")
#write.csv(MDD_test1, file = "meta1_expr.csv", row.names = T,quote = F)
write.csv(meta1, file = "meta1.csv", row.names = T,quote = F)

library(org.Hs.eg.db)
# 示例 Ensembl ID 列表
ensembl_ids <-rownames(MDD)

# 去除版本号（如果ID形如 ENSG00000141510.5）
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# 转换 Ensembl ID 到 Gene Symbol
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"  # 处理重复ID（可选：first/list/Filter/CharacterList）
)

# 转换为数据框
result <- data.frame(Ensembl_ID = ensembl_ids, Gene_Symbol = gene_symbols)
head(result)

# 将行名转换为 Gene Symbol
# rownames(MDD) <- result$Gene_Symbol[match(rownames(MDD), result$Ensembl_ID)]

# 删除无对应 Symbol 的行（可选）
# expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]

# 添加 Gene Symbol 到 MDD 并处理重复项
library(dplyr)
library(dplyr)
library(tibble)
# 将 Gene Symbol 添加到数据中
MDD_with_symbol <- MDD %>%
  as.data.frame() %>%
  mutate(Gene_Symbol = result$Gene_Symbol[match(rownames(MDD), result$Ensembl_ID)]) %>%
  filter(!is.na(Gene_Symbol))  # 删除无对应 Symbol 的行

# 按 Gene Symbol 合并重复项（取均值）
MDD_merged <- MDD_with_symbol %>%
  group_by(Gene_Symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames(var = "Gene_Symbol")

# 转换回矩阵（如果需要）
MDD_merged <- as.matrix(MDD_merged)
write.csv(MDD_merged, file = "meta1_expr.csv", row.names = T,quote = F)
MDD_merged_log <-log2(MDD_merged/colSums(MDD_merged)*1000000+1)
MDD_merged_log[is.infinite(MDD_merged_log)]<-0
summary(colSums(MDD_merged_log))
MDD_merged_log<-MDD_merged_log[rowSums(MDD_merged_log>0)>=3,]
write.csv(MDD_merged_log, file = "meta1_expr_log.csv", row.names = T,quote = F)