rm(list=ls())  #清空环境内变量
options(stringsAsFactors = F)  #避免自动将字符串转换为R语言因子
library(stringr)
setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
#################################测试集导入#################################################
setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/test/")

#####################################GSE251778#############################数据GEO下载有问题
## 获取表达矩阵
BiocManager::install("biomaRt")
library(variancePartition)
library(DESeq2)
library(dplyr)
library(ggplot2)
################
###   Data   ###
################GSE289146包含GSE251778
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
# MDD_merged_log <-log2(MDD_merged/colSums(MDD_merged)*100000000)
# MDD_merged_log[is.infinite(MDD_merged_log)]<-0
# summary(colSums(MDD_merged_log))
# MDD_merged_log<-MDD_merged_log[rowSums(MDD_merged_log>0)>=3,]
# write.csv(MDD_merged_log, file = "meta1_expr_log.csv", row.names = T,quote = F)



#################################################GSE52790####################################
library(GEOquery)
library(limma)
library(affy)
gset <- getGEO('GSE52790', destdir="./",
               AnnotGPL = T,     ## 注释文件
               getGPL = T)       ## 平台文件## 平台文件
gset[[1]]
exp<-exprs(gset[[1]])
## 获取临床信息
cli<-pData(gset[[1]])
meta2<-data.frame(cli$geo_accession,cli$title)
colnames(meta2)<-c("name","case")
meta2$case<-substr(meta2$case,1,2)
table(meta2$case)
meta2$group<-ifelse(meta2$case == "NC", 0, 1)
table(meta2$group)

#geneid 转换
GPL<-fData(gset[[1]])
gpl<-data.frame(GPL$ID,GPL$gene_symbol)
colnames(gpl)<-c("ID","Gene symbol")
gpl<-gpl[!(gpl$`Gene symbol`==""),]
gpl<-gpl[!(gpl$`Gene symbol`=="---"),]
exp<-as.data.frame(exp)
exp$ID<-rownames(exp)
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
table(duplicated(exp_symbol$`Gene symbol`))
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`Gene symbol`)
MDD_test2<-exp_unique[,meta2$name]
write.csv(MDD_test2, file = "meta2_expr.csv", row.names = T,quote = F)
write.csv(meta2, file = "meta2.csv", row.names = T,quote = F)




##################################	合并test data##################################
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/1_preprocessingData/")
MDD1<-read.csv("./GSE251778//meta1_expr.csv",row.names = 1)
MDD1_meta<-read.csv("GSE251778//meta1.csv",row.names = 1)
MDD1_logTPM <-log2(MDD1/colSums(MDD1)*1000000+1)
MDD1_logTPM[is.na(MDD1_logTPM)]<-0
summary(colSums(MDD1_logTPM))
MDD1_logTPM<-MDD1_logTPM[rowSums(MDD1_logTPM>0)>=3,]
write.csv(MDD1_logTPM, file = "./GSE251778/meta1_expr_log2TPM.csv", row.names = T,quote = F)

MDD2<-read.csv("./GSE52790/meta2_expr.csv",row.names = 1)
MDD2_meta<-read.csv("./GSE52790/meta2.csv",row.names = 1)
rownames(MDD2_meta)<-MDD2_meta$name
table(colnames(MDD1) %in% rownames(MDD1_meta))
table(colnames(MDD2) %in% rownames(MDD2_meta))
#检查是否有批次效应
MDD1_meta$batch=MDD1_meta$Cohort
MDD2_meta$batch="GSE52790"
MDD_meta<-bind_rows(MDD1_meta,MDD2_meta)
write.csv(MDD_meta, file = "MDDtest_meta.csv", row.names = T,quote = F)
MDD<-merge(MDD1_logTPM,MDD2,by="row.names",all=T)
rownames(MDD)<-MDD[,1]
MDD<-MDD[,-1]
MDD[is.na(MDD)]<-0
summary(colSums(MDD1_logTPM))
summary(colSums(MDD2))
summary(colSums(MDD))
write.csv(MDD, file = "MDDtest_expr.csv", row.names = T,quote = F)
table(colnames(MDD) %in% rownames(MDD_meta))
