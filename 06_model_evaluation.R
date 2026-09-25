setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/")
load("AUC_logisticmodel.RData")

library(pROC)

# 锁定最终基因
final_genes <- names(coef(logisticmodel$`glmBoost+Stepglm[both]`))[-1] #去掉截距
final_genes
selected_genes = c("KLRB1","IFI27","CEACAM6","HP","C9orf66","WDR41","TSHZ3")
train_pred = predict(logisticmodel$`glmBoost+Stepglm[both]`, newdata = as.data.frame(Train_set[,selected_genes]), type = "response")
roc_train = roc(Train_class[[classVar]], train_pred)
auc(roc_train)
plot(roc_train, main="Train ROC")
test_pred = predict(logisticmodel$`glmBoost+Stepglm[both]`, newdata = as.data.frame(Test_set[,selected_genes]), type = "response")
roc_test = roc(Test_class[[classVar]], test_pred)
auc(roc_test)
plot(roc_test, col="red", add=TRUE)




library(car)
# 调取你的最终logistic模型
fit <- logisticmodel$`glmBoost+Stepglm[both]`
# 计算7个基因VIF
vif_result <- vif(fit)
vif_result
# KLRB1    IFI27  CEACAM6       HP  C9orf66    WDR41    TSHZ3 
# 1.035281 1.100093 1.125770 1.119358 1.049207 1.047843 1.109534

Train_class <- read.table(file.path(data.path, "Training_class_more_info.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class$Cohort<-"Train-GSE98793"
Train_class$Cohort3<-"Training"
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
table(Test_class$Cohort)
Test_class$Cohort4<-paste0("Test-",Test_class$Cohort)
# Test_class2<-data.frame(Test_class$outcome,row.names = rownames(Test_class))
# colnames(Test_class2)<-"outcome"
library(dplyr)


train_cov <- data.frame(
  y = Train_class[[classVar]],
  sex = Train_class$gender,
  Train_set[,final_genes]
)
# 协变量校正模型
fit_cov <- glm(y ~ sex + ., family = binomial, data = train_cov)
# 对比校正前后基因OR、P值
cbind(原始OR=exp(coef(fit)), 校正后OR=exp(coef(fit_cov))[names(coef(fit))])




###################
# 固定7个核心基因
final_genes <- c("KLRB1","IFI27","CEACAM6","HP","C9orf66","WDR41","TSHZ3")
library(tidyverse)
library(pROC)
library(ggroc)
library(ggplot2)
library(rms)
library(corrplot)
library(reshape2)
library(ggforestplot)
library(RColorBrewer)
library(ggpubr)
# 整理全部需要绘图的队列
# 1 训练集
train_pred <- predict(fit, data.frame(Train_set[,final_genes]), type="response")
roc_train <- roc(Train_class[[classVar]], train_pred)

# 2 拆分3个测试队列
cohort_list <- unique(Test_class$Cohort)
roc_list <- list()
roc_list[["Training cohort"]] <- roc_train

for(co in cohort_list){
  idx <- Test_class$Cohort == co
  sub_x <- Test_set[idx, final_genes]
  sub_y <- Test_class[[classVar]][idx]
  sub_p <- predict(fit, newdata=data.frame(sub_x), type="response")
  roc_list[[co]] <- roc(sub_y, sub_p)
}

# ggplot绘制叠加ROC
ggroc(roc_list, size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
  labs(x="1-Specificity", y="Sensitivity", title="ROC curves across all cohorts") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  annotate("text", x=0.75, y=0.2,
           label=paste0(names(roc_list),": AUC=",
                        round(sapply(roc_list,auc),3),collapse="\n"))
ggsave("All_cohort_ROC.pdf", width=8, height=6)


# 提取原始、校正OR与置信区间
ori_beta <- coef(fit)[final_genes]
ori_or <- exp(ori_beta)
ori_ci <- exp(confint(fit)[final_genes,])

cov_beta <- coef(fit_cov)[final_genes]
cov_or <- exp(cov_beta)
cov_ci <- exp(confint(fit_cov)[final_genes,])

# 整理绘图数据框
forest_df <- data.frame(
  Gene = rep(final_genes,2),
  Model = rep(c("Unadjusted","Gender-adjusted"),each=7),
  OR = c(ori_or, cov_or),
  lower = c(ori_ci[,1], cov_ci[,1]),
  upper = c(ori_ci[,2], cov_ci[,2])
)

# 绘制森林图
ggplot(forest_df,aes(x=OR,y=Gene,color=Model))+
  geom_errorbarh(aes(xmin=lower,xmax=upper),height=0.2,linewidth=1)+
  geom_point(size=2)+
  geom_vline(xintercept=1,linetype="dashed",color="black")+
  scale_x_log10()+
  labs(x="Odds Ratio (log scale)",y="Gene",title="OR before and after gender adjustment")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5))
ggsave("Forestplot_OR.pdf",width=9,height=5)


train_df <- data.frame(y=Train_class[[classVar]], Train_set[,final_genes])
dd <- datadist(train_df)
options(datadist="dd")
lrm_fit <- lrm(y~., data=train_df, x=T, y=T)

pdf("Nomogram.pdf",width=10,height=7)
plot(nomogram(lrm_fit, fun=plogis, funlabel="Predicted disease probability"))
dev.off()
# 
# 
# #模块 4：分队列校准曲线 + 批量 HL 检验
# # 加载全部跑完的环境，一键补齐所有变量
# load("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/AUC_logisticmodel.RData")
# 
# # 锁定最终logistic模型
# fit <- logisticmodel$`glmBoost+Stepglm[both]`
# # 锁定7个目标基因
# final_genes <- names(coef(fit))[-1]
# final_genes
# 
# # 固定必备4组核心数据（你的代码已经完成标准化）
# # Train_set：行=样本，列=7基因；训练集标准化表达
# # Train_class：训练集标签，列名outcome，0/1
# # Test_set：行=样本，列=7基因；测试集标准化表达
# # Test_class：测试集，包含outcome、Cohort（GSE分组）
# 
# library(ggplot2)
# library(ResourceSelection)
# 
# # 通用校准曲线函数，参数固定：模型、表达矩阵、真实分组、标题
# cal_simple <- function(glm_model, expr_mat, true_y, plot_title){
#   df <- data.frame(
#     true_y = true_y,
#     pred_p = predict(glm_model, newdata = data.frame(expr_mat[,final_genes]), type = "response")
#   )
#   # Hosmer-Lemeshow检验
#   hl_p <- hoslem.test(df$true_y, df$pred_p, g=10)$p.value
#   # 概率分箱
#   df$bin <- cut(df$pred_p, breaks = seq(0,1,0.2), include.lowest = TRUE)
#   bin_summary <- aggregate(cbind(pred_p, true_y) ~ bin, data = df, mean)
#   
#   p <- ggplot(bin_summary, aes(x = pred_p, y = true_y)) +
#     geom_line(color = "#E63946", linewidth = 1) +
#     geom_point(size = 2, color = "#E63946") +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
#     labs(x = "Predicted probability", y = "Observed probability",
#          title = paste0(plot_title, " | HL P = ", round(hl_p, 3))) +
#     xlim(0, 1) + ylim(0, 1) +
#     theme_bw() + theme(plot.title = element_text(hjust = 0.5))
#   print(p)
#   return(hl_p)
# }
# 
# cal_simple(
#   glm_model = fit,
#   expr_mat = Train_set,
#   true_y = Train_class$outcome,
#   plot_title = "Training cohort calibration"
# )
# 
# # 批量导出全部校准图到PDF
# pdf("All_cohort_calibration.pdf", width=8, height=6)
# # 训练集
# cal_simple(fit, Train_set, Train_class$outcome, "Training cohort")
# # 遍历所有测试队列
# all_cohort <- unique(Test_class$Cohort)
# for (co in all_cohort) {
#   idx <- Test_class$Cohort == co
#   cal_simple(
#     glm_model = fit,
#     expr_mat = Test_set[idx, ],
#     true_y = Test_class$outcome[idx],
#     plot_title = paste0(co, " calibration")
#   )
# }
# dev.off()

# 
# library(rms)
# # 训练集构建lrm
# train_df <- data.frame(y = Train_class$outcome, Train_set[,final_genes])
# dd <- datadist(train_df)
# options(datadist = "dd")
# lrm_fit <- lrm(y ~ ., data = train_df, x=T, y=T)
# 
# # 绘制GSE52790 bootstrap校准
# pdf("GSE52790_bootstrap_cal.pdf", width=7, height=6)
# plot(calibrate(lrm_fit, newdata = Test_set[Test_class$Cohort=="GSE52790",], B=100),
#      xlab="Predicted probability", ylab="Observed probability",
#      main="GSE52790 Bootstrap Calibration")
# dev.off()



setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model/")

library(rms)
library(ggplot2)

# 锁定模型与7个基因
fit <- logisticmodel$`glmBoost+Stepglm[both]`
final_genes <- names(coef(fit))[-1]

# 仅用训练集构建lrm模型（全局唯一，所有队列共用）
train_df <- data.frame(y = Train_class$outcome, Train_set[,final_genes])
dd <- datadist(train_df)
options(datadist = dd)
lrm_fit <- lrm(y ~ ., data = train_df, x = TRUE, y = TRUE)
setwd("./InputData/test/")
plot_cali_boot <- function(sub_expr, sub_y, cohort_name, B=100){
  plot_df <- data.frame(y = sub_y, sub_expr[,final_genes])
  # Bootstrap校准计算
  cal_obj <- calibrate(lrm_fit, newdata = plot_df, B = B)
  cal_df <- as.data.frame(cal_obj)
  
  # 自动匹配预测概率列、观测概率列，兼容所有rms版本
  pred_col <- grep("pred|mean", colnames(cal_df), value = TRUE)[1]
  obs_col  <- grep("KM|obs", colnames(cal_df), value = TRUE)[1]
  
  # 计算MAE
  pred_vals <- cal_df[[pred_col]]
  obs_vals  <- cal_df[[obs_col]]
  mae <- round(mean(abs(pred_vals - obs_vals)),3)
  
  # 导出PDF
  pdf(paste0("BootCal_",cohort_name,".pdf"), width=7, height=6)
  plot(cal_obj,
       xlab = "Predicted probability",
       ylab = "Observed probability",
       main = paste0(cohort_name, "\nMAE = ", mae, ", n = ", nrow(plot_df)),
       legend = TRUE)
  mtext(paste0("Bootstrap repetitions = ", B), side=1, line=4)
  dev.off()
  
  # 返回汇总信息
  return(data.frame(Cohort=cohort_name, SampleSize=nrow(plot_df), MAE=mae))
}


plot_cali_boot <- function(sub_expr, sub_y, cohort_name, B=100){
  plot_df <- data.frame(y = sub_y, sub_expr[,final_genes])
  # Bootstrap校准计算
  cal_obj <- calibrate(lrm_fit, newdata = plot_df, B = B)
  cal_df <- as.data.frame(cal_obj)
  
  # 自动匹配预测概率列、观测概率列，兼容所有rms版本
  pred_col <- grep("pred|mean", colnames(cal_df), value = TRUE)[1]
  obs_col  <- grep("KM|obs", colnames(cal_df), value = TRUE)[1]
  
  # 计算MAE，NA过滤防止NaN
  pred_vals <- cal_df[[pred_col]]
  obs_vals  <- cal_df[[obs_col]]
  valid_idx <- !is.na(pred_vals) & !is.na(obs_vals)
  mae <- if(sum(valid_idx) > 0) round(mean(abs(pred_vals[valid_idx] - obs_vals[valid_idx])), 3) else NA
  
  # 导出PDF（用绝对路径目录）
  pdf(paste0("BootCal_", cohort_name, ".pdf"), width = 7, height = 6)
  plot(cal_obj,
       xlab = "Model-predicted values",
       ylab = "Observed outcome proportion",
       main = paste0(cohort_name, "\nMAE = ", mae, ", n = ", nrow(plot_df)),
       legend = TRUE)
  mtext(paste0("Bootstrap repetitions = ", B), side = 1, line = 4)
  dev.off()
  
  return(data.frame(Cohort = cohort_name, SampleSize = nrow(plot_df), MAE = mae))
}






mae_all <- list()

# 1. 训练集
mae_all[[1]] <- plot_cali_boot(Train_set, Train_class$outcome, "Training_Cohort")

# 2. GSE251778 男性单独
idx_m <- Test_class$Cohort == "GSE251778_M"
mae_all[[2]] <- plot_cali_boot(Test_set[idx_m,], Test_class$outcome[idx_m], "GSE251778_Male")

# 3. GSE251778 女性单独
idx_f <- Test_class$Cohort == "GSE251778_F"
mae_all[[3]] <- plot_cali_boot(Test_set[idx_f,], Test_class$outcome[idx_f], "GSE251778_Female")

# 4. GSE251778 男女合并
idx_all_251778 <- Test_class$Cohort %in% c("GSE251778_M","GSE251778_F")
mae_all[[4]] <- plot_cali_boot(Test_set[idx_all_251778,], Test_class$outcome[idx_all_251778], "GSE251778_Combined")

# 5. 独立外部队列 GSE52790
idx_52790 <- Test_class$Cohort == "GSE52790"
mae_all[[5]] <- plot_cali_boot(Test_set[idx_52790,], Test_class$outcome[idx_52790], "GSE52790")

# 汇总所有队列样本量、MAE表格
mae_summary_table <- do.call(rbind, mae_all)
write.table(mae_summary_table, "All_Calibration_MAE_Summary.txt", sep="\t", row.names=F, quote=F)
print(mae_summary_table)





cor_matrix <- cor(Train_set[,final_genes], method="pearson")
pdf("Gene_cor_heatmap.pdf", width=7, height=6)
corrplot(cor_matrix, method="color", addCoef.col="black", number.cex=0.7,
         col = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
         tl.col="black", tl.srt=45)
dev.off()





set.seed(666)
train_y <- Train_class$outcome
# 5折拆分
fold_idx <- createFolds(train_y, k=5, returnTrain = TRUE)
auc_vec <- c()

for(i in 1:5){
  tr_idx <- fold_idx[[i]]
  val_idx <- setdiff(1:nrow(Train_set), tr_idx)
  # 训练折拟合logistic
  fold_fit <- glm(train_y[tr_idx] ~ ., family = binomial, data = data.frame(Train_set[tr_idx,final_genes]))
  pred_val <- predict(fold_fit, newdata=data.frame(Train_set[val_idx,final_genes]), type="response")
  auc_vec[i] <- auc(roc(train_y[val_idx], pred_val))
}

# 绘图
cv_df <- data.frame(AUC = auc_vec)
ggplot(cv_df, aes(x="5-fold cross validation", y=AUC)) +
  geom_boxplot(fill="#74A9CF") +
  geom_jitter(width=0.15, size=2) +
  ylim(0.5,1) +
  labs(y="AUC", title="AUC distribution of 5-fold cross validation") +
  theme_bw()
ggsave("Fivefold_AUC_box.pdf", width=4, height=5)







# 加载依赖
library(ggplot2)
library(tidyr)
library(ggpubr)

# 整合训练集表达+分组标签
expr_df <- as.data.frame(Train_set)
expr_df$outcome <- Train_class$outcome

# 正确长格式转换（把expr改成expr_df）
expr_long <- pivot_longer(
  data = expr_df, 
  cols = all_of(final_genes), 
  names_to = "Gene", 
  values_to = "Expression"
)

expr_long$Group <- factor(expr_long$outcome, labels = c("Control","Disease"))

# 绘图
ggplot(expr_long, aes(x=Group, y=Expression, fill=Group)) +
  geom_violin(scale="width", alpha=0.7) +
  geom_boxplot(width=0.2, fill="white") +
  facet_wrap(~Gene, ncol=4, scales="free_y") +
  stat_compare_means(comparisons = list(c("Control","Disease")), method="wilcox.test", label="p.signif") +
  scale_fill_manual(values=c("#89CFF0","#FF6B6B")) +
  labs(x="Group", y="Z-score normalized expression") +
  theme_bw() + theme(legend.position = "none")

ggsave("Gene_violin_case_ctrl.pdf", width=12, height=6)









##########################################
# ============================================================
# 二分类模型性能指标批量计算
# 包含：灵敏度、特异度、PPV、NPV、95% CI、Brier Score
# 加载包
# ============================================================预处理
logisticmodel$`glmBoost+Stepglm[both]`
# 加载包
library(rms)
library(caret)
library(ModelMetrics)
library(writexl)

# 读取你已有的全部文件
setwd("C:/Users/sunshiny/Downloads/SLEMDD/R代码/Figure5/model//")
load("C:/Users/sunshiny/Downloads/SLEMDD/coDEGgene.RData")
logist_mod <- readRDS("./Result/logisticmodel.rds")
logist_mod<-logisticmodel$`glmBoost+Stepglm[both]`


RS_mat <- read.table("./Result/RS_mat.txt", header = T, sep = "\t")
RS_mat <- data.frame(RS_mat$X,RS_mat$glmBoost.Stepglm.both.)
colnames(RS_mat)<-c("samplename","prob")

Train_class <- read.table(file.path(data.path, "Training_class_more_info.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class$Cohort<-"Train-GSE98793"
Train_class$samplename<-rownames(Train_class)
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
table(Test_class$Cohort)
Test_class$samplename<-rownames(Test_class)

train_sub <- Train_class[, c("outcome", "Cohort", "samplename")]
test_sub <- Test_class[, c("outcome", "Cohort", "samplename")]

# 行合并
fea_df_all <- rbind(train_sub, test_sub)


# 左合并，以fea_df_all全部样本为准
df_merge <- merge(fea_df_all, RS_mat, by = "samplename", all.x = TRUE)

# 规范列名：outcome重命名为y，最终四列：samplename, Cohort, y, prob
colnames(df_merge)[colnames(df_merge)=="outcome"] <- "y"
df_final <- df_merge[, c("samplename", "Cohort", "y", "prob")]

# 检查匹配是否正常，有没有概率缺失
table(is.na(df_final$prob))
table(df_final$Cohort)

write.table(df_final, "finalmodel_cohort_prob_label.txt", sep = "\t", row.names = FALSE)




##########################################
# ============================================================
# 二分类模型性能指标批量计算
# 包含：灵敏度、特异度、PPV、NPV、95% CI、Brier Score
# ============================================================


######################################################################################
# ##################Favorable AUC indicates stable risk ranking, while zero sensitivity 
# at cutoff=0.5 was caused by global downward probability shift in external cohorts; 
# cohort-specific Youden cutoffs were used for final classification calculation.


#install.packages(c("pROC","binom","openxlsx"))
# ============================================================
# 二分类模型性能指标批量计算
# 支持：固定阈值(0.5) + 各队列专属最优截断值(ROC约登指数)
# 包含：灵敏度、特异度、PPV、NPV、95% CI、Brier Score
# ============================================================

library(dplyr)

# ---------- 辅助函数：二项比例95% CI (Wald方法) ----------
ci_binomial <- function(p, n, conf.level = 0.95) {
  if (n == 0 || is.na(p) || is.nan(p) || is.infinite(p)) {
    return(c(NA, NA))
  }
  z <- qnorm(1 - (1 - conf.level) / 2)
  se <- sqrt(p * (1 - p) / n)
  low <- max(0, p - z * se)
  high <- min(1, p + z * se)
  return(c(low, high))
}

# ---------- 核心函数：计算单个数据集的所有指标 ----------
calc_metrics <- function(df, threshold = 0.5) {
  y_true <- df$y
  y_prob <- df$prob
  y_pred <- ifelse(y_prob >= threshold, 1, 0)
  
  TP <- sum(y_true == 1 & y_pred == 1)
  TN <- sum(y_true == 0 & y_pred == 0)
  FP <- sum(y_true == 0 & y_pred == 1)
  FN <- sum(y_true == 1 & y_pred == 0)
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
  
  sens_ci <- ci_binomial(sensitivity, TP + FN)
  spec_ci <- ci_binomial(specificity, TN + FP)
  ppv_ci <- ci_binomial(PPV, TP + FP)
  npv_ci <- ci_binomial(NPV, TN + FN)
  
  brier <- mean((y_prob - y_true)^2)
  
  n_total <- nrow(df)
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  
  return(data.frame(
    N_total = n_total,
    N_positive = n_pos,
    N_negative = n_neg,
    Threshold = threshold,
    Sensitivity = round(sensitivity, 4),
    Sensitivity_95CI = paste0("(", round(sens_ci[1], 4), ", ", round(sens_ci[2], 4), ")"),
    Specificity = round(specificity, 4),
    Specificity_95CI = paste0("(", round(spec_ci[1], 4), ", ", round(spec_ci[2], 4), ")"),
    PPV = round(PPV, 4),
    PPV_95CI = paste0("(", round(ppv_ci[1], 4), ", ", round(ppv_ci[2], 4), ")"),
    NPV = round(NPV, 4),
    NPV_95CI = paste0("(", round(npv_ci[1], 4), ", ", round(npv_ci[2], 4), ")"),
    Brier_Score = round(brier, 4),
    TP = TP, TN = TN, FP = FP, FN = FN,
    stringsAsFactors = FALSE
  ))
}

# ---------- 计算约登指数最优阈值 ----------
calc_optimal_threshold <- function(df) {
  y_true <- df$y
  y_prob <- df$prob
  thresholds <- sort(unique(y_prob))
  
  best_youden <- -Inf
  best_threshold <- 0.5
  
  for (th in thresholds) {
    y_pred <- ifelse(y_prob >= th, 1, 0)
    TP <- sum(y_true == 1 & y_pred == 1)
    TN <- sum(y_true == 0 & y_pred == 0)
    FP <- sum(y_true == 0 & y_pred == 1)
    FN <- sum(y_true == 1 & y_pred == 0)
    
    sens <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    spec <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
    youden <- sens + spec - 1
    
    if (youden > best_youden) {
      best_youden <- youden
      best_threshold <- th
    }
  }
  return(best_threshold)
}

# ---------- 固定阈值批量计算 ----------
calc_fixed_threshold <- function(df, threshold = 0.5) {
  overall <- calc_metrics(df, threshold)
  overall$Cohort <- "Overall"
  
  cohorts <- unique(df$Cohort)
  results <- lapply(cohorts, function(c) {
    sub_df <- df %>% filter(Cohort == c)
    res <- calc_metrics(sub_df, threshold)
    res$Cohort <- c
    return(res)
  })
  
  result_df <- do.call(rbind, c(list(overall), results))
  result_df <- result_df %>% 
    select(Cohort, N_total, N_positive, N_negative, Threshold,
           Sensitivity, Sensitivity_95CI, 
           Specificity, Specificity_95CI,
           PPV, PPV_95CI, 
           NPV, NPV_95CI, 
           Brier_Score,
           TP, TN, FP, FN)
  
  return(result_df)
}

# ---------- 各队列最优阈值批量计算 ----------
calc_optimal_threshold_all <- function(df) {
  overall_th <- calc_optimal_threshold(df)
  overall <- calc_metrics(df, overall_th)
  overall$Cohort <- "Overall"
  overall$Optimal_Threshold <- overall_th
  
  cohorts <- unique(df$Cohort)
  results <- lapply(cohorts, function(c) {
    sub_df <- df %>% filter(Cohort == c)
    opt_th <- calc_optimal_threshold(sub_df)
    res <- calc_metrics(sub_df, opt_th)
    res$Cohort <- c
    res$Optimal_Threshold <- opt_th
    return(res)
  })
  
  result_df <- do.call(rbind, c(list(overall), results))
  result_df <- result_df %>% 
    select(Cohort, N_total, N_positive, N_negative, Optimal_Threshold,
           Sensitivity, Sensitivity_95CI, 
           Specificity, Specificity_95CI,
           PPV, PPV_95CI, 
           NPV, NPV_95CI, 
           Brier_Score,
           TP, TN, FP, FN)
  
  return(result_df)
}

# ============================================================
# 直接运行
# ============================================================

# 1. 固定阈值0.5
results_fixed <- calc_fixed_threshold(df_final, threshold = 0.5)
cat("\n========== 固定阈值(0.5)结果 ==========\n")
print(results_fixed)

# 2. 各队列专属最优阈值
results_optimal <- calc_optimal_threshold_all(df_final)
cat("\n========== 各队列最优阈值结果 ==========\n")
print(results_optimal)

# 3. 保存结果
write.csv(results_fixed, "metrics_fixed_threshold.csv", row.names = FALSE)
write.csv(results_optimal, "metrics_optimal_threshold.csv", row.names = FALSE)
cat("\n结果已保存！\n")


