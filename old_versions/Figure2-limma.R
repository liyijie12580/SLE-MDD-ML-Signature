setwd("C:/Users/sunshiny/Downloads/SLEMDD/")
MDD<-read.csv("MDD_expr_combat.csv",row.names = 1)
SLE<-read.csv("SLE_expr.csv",row.names = 1)
MDD_meta<-read.csv("MDD_meta.csv",row.names = 1)
SLE_meta<-read.csv("SLE_meta.csv",row.names = 1)

## 导入R包
library(limma)
library(dplyr)
######################################MDD
job="MDD"
expr_data<-MDD#输入文件TPM原始值，行名是基因，列名是样本
expr_data <- expr_data[rowSums(expr_data)>3,] #删除表达量为0的基因
expr_data<-expr_data[,colSums(expr_data)>=200]
#expr_data = log2(expr_data) #log化处理
#expr_data[expr_data == -Inf] = 0 #将log化后的负无穷值替换为0
group<-MDD_meta
levels(factor(group$case))
# 构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$case))
colnames(design) <- levels(factor(group$case))
rownames(design) <- colnames(expr_data)
# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(CTRL-MDD,levels = design) #根据实际的样本分组修改，这里对照组CK，处理组HT
# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
colnames(DEG)
# [1] "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        
## 导出所有的差异结果
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 0.25, "up-regulated",
                              ifelse(DEG$logFC < -0.25, "down-regulated", "unchanged")))

# DEG$regulate <- ifelse(DEG$adj.P.Val >0.05, "unchanged",
#                        ifelse(DEG$logFC > 0.25, "up-regulated",
#                               ifelse(DEG$logFC < -0.25, "down-regulated", "unchanged")))
table(DEG$regulate)
# tempOutput<-DEG
# nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
# diffsig <- nrDEG  
# write.csv(diffsig, "MDD.limmaOut.csv") 
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure2")
write.table(table(DEG$regulate),file = paste0(job,"_","DEG_result_025_005.txt"),
            sep = "\t",quote = F,row.names = T,col.names = T)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = paste0(job,"_","DEG_result.csv"),
            sep = "\t",quote = F,row.names = F,col.names = T)

# 火山图的绘制 ------------------------------------------------------------------
head(DEG)
DEG$Genes <- rownames(DEG)
pdf(paste0(job,"_","volcano2.pdf"),width = 10,height = 7)
library(ggplot2)
p<-ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
  geom_point(alpha=0.5,size=2,aes(color=regulate))+ #点的透明度，大小
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
  geom_vline(xintercept = c(-0.25,0.25),lty=4,col ="black",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ #adj.p.val分界线
  theme_bw()  #火山图绘制
DEG$Genesadd<-DEG$Genes
DEG$Genesadd[DEG$regulate=="unchanged"]<-NA
p+
  #geom_text(data=DEG,mapping = aes(label=Genesadd))+
  ggrepel::geom_text_repel(aes(label=Genesadd),DEG)
dev.off()


#################################SLE
job="SLE"
head(SLE)
expr_data<-SLE#输入文件TPM原始值，行名是基因，列名是样本
expr_data <- expr_data[rowSums(expr_data)>3,] #删除表达量为0的基因
expr_data<-expr_data[,colSums(expr_data>0)>=200]
dim(expr_data)
#expr_data = log2(expr_data) #log化处理
#expr_data[expr_data == -Inf] = 0 #将log化后的负无穷值替换为0
group<-SLE_meta[colnames(expr_data),]
levels(factor(group$case))
group$case<-gsub(" Case","",gsub("case/control: ","",group$case))
levels(factor(group$case))
# 构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$case))
colnames(design) <- levels(factor(group$case))
rownames(design) <- colnames(expr_data)
# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(Control-SLE,levels = design) #根据实际的样本分组修改，这里对照组CK，处理组HT
# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
colnames(DEG)
# [1] "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        
## 导出所有的差异结果
DEG$regulate <- ifelse(DEG$P.Value >0.05, "unchanged",
                       ifelse(DEG$logFC > 0.25, "up-regulated",
                              ifelse(DEG$logFC < -0.25, "down-regulated", "unchanged")))
# DEG$regulate <- ifelse(DEG$adj.P.Val >0.05, "unchanged",
#                        ifelse(DEG$logFC > 1, "up-regulated",
#                               ifelse(DEG$logFC < -0.25, "down-regulated", "unchanged")))
table(DEG$regulate)
# tempOutput<-DEG
# nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
# diffsig <- nrDEG  
# write.csv(diffsig, "MDD.limmaOut.csv") 

write.table(table(DEG$regulate),file = paste0(job,"_","DEG_result_025_005.txt"),
            sep = "\t",quote = F,row.names = T,col.names = T)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = paste0(job,"_","DEG_result.csv"),
            sep = "\t",quote = F,row.names = F,col.names = T)

# 火山图的绘制 ------------------------------------------------------------------

DEG$Genes <- rownames(DEG)
pdf(paste0(job,"_","volcano2.pdf"),width = 7,height = 7)
library(ggplot2)
p<-ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
  geom_point(alpha=0.5,size=2,aes(color=regulate))+ #点的透明度，大小
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
  geom_vline(xintercept = c(-0.25,0.25),lty=4,col ="black",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ #adj.p.val分界线
  theme_bw()  #火山图绘制
DEG$Genesadd<-DEG$Genes
DEG$Genesadd[DEG$regulate=="unchanged"]<-NA
p+
  ggrepel::geom_text_repel(aes(label=Genesadd),DEG)
dev.off()

###############韦恩图（VennDiagram 包，适用样本数 2-5）
# install.packages("VennDiagram")
# 
# install.packages("openxlsx")
library(VennDiagram)
library("openxlsx")
#读入作图文件，all.txt即上述提到的记录group1-4的元素名称的文件
dat <- read.table('gene_veen.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
colnames(dat)
#以2个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
###############up
venn_list <- list(SLE_up =dat$SLE_up[1:997],MDD_up = dat$MDD_up[1:51])
venn.diagram(venn_list, filename = 'venn_up2.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 1.5, cat.fontfamily = 'serif', 
             output=TRUE)
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                        col = 'black', cex = 1.5, fontfamily = 'serif', 
                        cat.cex = 1.5, cat.fontfamily = 'serif',
                        output=TRUE)
pdf('venn2_up.pdf' ,width = 7, height = 7)
grid.draw(venn.plot)
dev.off()

inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
gene<-data.frame(inter[-c(5, 6)])
#SLE_up∩MDD_up
# "KLRG1", "KLRB1", "SYTL2"
write.xlsx(gene, 'venn_up2_inter.xlsx', rowNames = FALSE, sep = '\t')


###############down
venn_list <- list(SLE_down = dat$SLE_down[1:925],MDD_down = dat$MDD_down[1:112])
venn.diagram(venn_list, filename = 'venn_down2.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 1.5, cat.fontfamily = 'serif',
             output=TRUE)
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                        col = 'black', cex = 1.5, fontfamily = 'serif', 
                        cat.cex = 1.5, cat.fontfamily = 'serif',
                        output=TRUE)
pdf('venn2_down.pdf' ,width = 7, height = 7)
grid.draw(venn.plot)
dev.off()

inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
gene<-data.frame(inter[-c(5, 6)])
#SLE_down∩MDD_down
# IFI27, RSAD2, OLFM4, TNFAIP6, CAMP, SIPA1L2, DEFA4, MPO, TCN1, SLC26A8, CEACAM8, DDX60L, CEACAM6, HP, S100A12, CACNA1E, SLPI, RNASE3, C9orf66, ELANE, WDR41, SERPING1, CLEC5A, MMP8, NSUN7, TSHZ3
write.xlsx(gene, 'venn_down2_inter.xlsx', rowNames = FALSE, sep = '\t')


#以4个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(SLE_up =dat$SLE_up[1:997],MDD_up = dat$MDD_up[1:52],
                  SLE_down = dat$SLE_down[1:925],MDD_down = dat$MDD_down[1:112])

venn.diagram(venn_list, filename = 'venn4.png', imagetype = 'png', 
             fill = c('red', 'blue', 'green', 'orange'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green', 'orange'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green', 'orange'), cex = 1.5, fontfamily = 'serif')

venn.plot <- venn.diagram(venn_list, filename = NULL, fill = c('red', 'blue', 'green', 'orange'), alpha = 0.50, 
                          cat.col = c('red', 'blue', 'green', 'orange'), cat.cex = 1.5, cat.fontfamily = 'serif',
                          col = c('red', 'blue', 'green', 'orange'), cex = 1.5, fontfamily = 'serif')
pdf('venn4.pdf' ,width = 7, height = 7)
grid.draw(venn.plot)
dev.off()

intersect(dat$SLE_up,dat$MDD_down)
#[1] "ALAS2"
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
gene<-data.frame(inter[-c(5, 6)])
write.xlsx(gene, 'venn4_inter.xlsx', rowNames = FALSE, sep = '\t')

################GO/KEGG富集分析基因提取

up<-intersect(dat$SLE_up,dat$MDD_up)
up<-data.frame(up)
write.xlsx(up, 'up.xlsx', rowNames = FALSE, sep = '\t')
down<-intersect(dat$SLE_down,dat$MDD_down)
down<-data.frame(down)
write.xlsx(down, 'down.xlsx', rowNames = FALSE, sep = '\t')
# install.packages("gProfileR")
# library("gProfileR")
# gprofiler(up, organism = "hsapiens")
gene<-c(up$up,down$down)
gene<-gene[1:29]
save(gene,file = "coDEGgene.RData")

