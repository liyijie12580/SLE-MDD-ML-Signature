library(pROC)
code.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/Codes"
data.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/InputData/"
res.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/Result/"
setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/Result/")
model <- readRDS(file.path(res.path, "model.rds"))
RS_list<-read.table(file ="RS_mat.txt",header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
fea_df<-read.table(file ="fea_df.txt",header = T, sep = "\t", check.names = F,stringsAsFactors = F)
Class_mat<-read.table(file ="Class_mat.txt",header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
AUC_mat <-read.table(file ="AUC_mat.txt",header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# MDDtrain_class <- read.csv("MDD_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
# MDDtrain_class$group<-ifelse(MDDtrain_class$case == "CTRL", 0, 1)
# MDDtrain_class<-data.frame(MDDtrain_class$group,MDDtrain_class$gender,MDDtrain_class$studygroup,row.names = MDDtrain_class$name)
# colnames(MDDtrain_class)<-c("outcome","gender","studygroup")
# 
# MDDtest_class <- read.csv("MDDtest_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
# MDDtest_class<-data.frame(MDDtest_class$group,MDDtest_class$batch,row.names = MDDtest_class$name)
# colnames(MDDtest_class)<-c("outcome","source")
Train_class <- read.table(file.path(data.path, "Training_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
table(Test_class$Cohort)
Test_class2<-data.frame(Test_class$outcome,row.names = rownames(Test_class))
colnames(Test_class2)<-"outcome"
class<-rbind(Train_class,Test_class2)

# 
setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/ ")
load("C:/Users/sunshiny/Downloads/SLEMDD/coDEGgene.RData")
MDDtrain <- read.csv("MDD_expr_combat.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtrain <-MDDtrain[gene,]
MDDtrain_class <- read.csv("MDD_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
MDDtrain_class$group<-ifelse(MDDtrain_class$case == "CTRL", 0, 1)
# MDDtrain<-data.frame(t(MDDtrain))
# MDDtrain$group<- factor(MDDtrain_class$group)
# MDDtrain$batch <-factor(MDDtrain_class$batch)
# MDDtrain$studygroup<-factor(MDDtrain_class$studygroup)
# MDDtrain$gender<-factor(MDDtrain_class$gender)
# MDDtrain$age<-factor(MDDtrain_class$age)
library(DMwR)
set.seed(123)
data_balanced <- SMOTE(group ~ ., MDDtrain, perc.over = 200)
table(data_balanced$group)
MDDtrain_class<-data.frame(data_balanced$group,row.names = rownames(data_balanced))
colnames(MDDtrain_class)<-c("outcome")
MDDtrain <-t(as.matrix(data_balanced[, -ncol(data_balanced)]))
# 
# RS_list<-read.table(file ="RS_mat.txt",header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# 
# ##############trian
# setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/model/")
# RS_list<-RS_list[rownames(MDDtrain_class),]
# MDDtrain_class <- read.csv("MDD_meta.csv", row.names = 1,check.names = F,stringsAsFactors = F)
# MDDtrain_class$group<-ifelse(MDDtrain_class$case == "CTRL", 0, 1)
# MDDtrain_class<-data.frame(MDDtrain_class$group,MDDtrain_class$studygroup,
#                            MDDtrain_class$age,MDDtrain_class$gender,MDDtrain_class$batch,row.names = MDDtrain_class$name)


trainname<-rownames(MDDtrain_class)
intersect()
class[trainname,]
single_model_dataset<-data.frame(outcomes=class[trainname,]$outcome,predictions=RS_list[,"XGBoost"],
                         Cohort=c(rep("Training GSE98793",448),rep("Testing GSE39653",45),
                         rep("Testing GSE52790",22)))
library(runway)
threshperf_plot(single_model_dataset,
                outcome = 'outcomes',
                positive = '1',
                prediction = 'predictions')
cal_plot(single_model_dataset,
         outcome = 'outcomes',
         positive = '1',
         prediction = 'predictions',
         n_bins = 0,
         show_loess = TRUE)


multi_model_dataset<-data.frame(outcomes=class$outcome[449:515],predictions=RS_list[,"XGBoost"][449:515],
                                 Cohort=c(rep("Testing GSE39653",45),
                                          rep("Testing GSE52790",22)))
threshperf_plot_multi(multi_model_dataset,
                      outcome = 'outcomes',
                      positive = '1',
                      prediction = 'predictions',
                      model = 'Cohort')
cal_plot_multi(multi_model_dataset,
               outcome = 'outcomes',
               positive = '1',
               prediction = 'predictions',
               model = 'Cohort',
               n_bins = 0,
               show_loess = TRUE)

cal_plot_multi(multi_model_dataset,
               outcome = 'outcomes',
               positive = '1',
               prediction = 'predictions',
               model = 'Cohort',
               n_bins = 10)


single_model_dataset<-data.frame(outcomes=class$outcome[449:515],predictions=RS_list[,"XGBoost"][449:515])

                                library(runway)
                                threshperf_plot(single_model_dataset,
                                                outcome = 'outcomes',
                                                positive = '1',
                                                prediction = 'predictions')
                                cal_plot(single_model_dataset,
                                         outcome = 'outcomes',
                                         positive = '1',
                                         prediction = 'predictions',
                                         n_bins = 0,
                                         show_loess = TRUE)





#train
roc_obj <- roc(
  response = class$outcome[1:448],    # 真实标签
  predictor = RS_list[,"XGBoost"][1:448],  # 预测概率
  levels = c("0", "1")       # 指定因子水平顺序（负类在前，正类在后）
)

plot(roc_obj,
     col="red",
     print.auc = TRUE,       # 显示 AUC 值
     auc.polygon = TRUE,     # 填充 AUC 区域
     auc.polygon.col = "lightblue",
     max.auc.polygon = TRUE, # 显示最大 AUC 区域
     grid = TRUE,            # 显示网格线
     legacy.axes = TRUE,     # 横轴为假阳性率（FPR），纵轴为真阳性率（TPR）
main = "Training GSE98793 ROC Curve"
)




#GSE39653
rownames(class)[449:493]

roc_obj <- roc(
  response = class$outcome[449:493],    # 真实标签
  predictor = RS_list[,"XGBoost"][449:493],  # 预测概率
  levels = c("0", "1")       # 指定因子水平顺序（负类在前，正类在后）
)


plot(roc_obj,
     col="red",
     print.auc = TRUE,       # 显示 AUC 值
     auc.polygon = TRUE,     # 填充 AUC 区域
     auc.polygon.col = "lightblue",
     max.auc.polygon = TRUE, # 显示最大 AUC 区域
     grid = TRUE,            # 显示网格线
     legacy.axes = TRUE,     # 横轴为假阳性率（FPR），纵轴为真阳性率（TPR）
     main = "Testing GSE39653 ROC Curve"
)



#GSE52790

roc_obj <- roc(
  response = class$outcome[494:515],    # 真实标签
  predictor = RS_list[,"XGBoost"][494:515],  # 预测概率
  levels = c("0", "1")       # 指定因子水平顺序（负类在前，正类在后）
)


plot(roc_obj,
     col="red",
     print.auc = TRUE,       # 显示 AUC 值
     auc.polygon = TRUE,     # 填充 AUC 区域
     auc.polygon.col = "lightblue",
     max.auc.polygon = TRUE, # 显示最大 AUC 区域
     grid = TRUE,            # 显示网格线
     legacy.axes = TRUE,     # 横轴为假阳性率（FPR），纵轴为真阳性率（TPR）
     main = "Testing GSE52790 ROC Curve"
)


