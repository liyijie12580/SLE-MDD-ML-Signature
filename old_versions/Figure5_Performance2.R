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
