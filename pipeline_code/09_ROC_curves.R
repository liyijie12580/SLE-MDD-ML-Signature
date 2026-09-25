library(pROC)
code.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model2/Codes"
data.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model2/InputData/"
res.path="C:/Users/sunshiny/Downloads/SLEMDD/MDD/model2/Result/"
setwd("C:/Users/sunshiny/Downloads/SLEMDD/MDD/model2/Result/")
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
#Train_class <- read.table(file.path(data.path, "Training_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.table(file.path(data.path, "Training_class_more_info.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class$age<-as.numeric(gsub("age: ","",Train_class$age))
Train_class$gender<-gsub("gender: ","",Train_class$gender)
Train_class$Cohort<-"GSE98793"
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
table(Test_class$Cohort)
Test_class2<-data.frame(Test_class$outcome,row.names = rownames(Test_class))
colnames(Test_class2)<-"outcome"
library(dplyr)
class<-bind_rows(Train_class,Test_class)
Class_mat$outcome<-class$outcome

aaa<-data.frame(Class_mat$`Enet[alpha=0.1]`,Class_mat$outcome)[1:192,]
aaa<-data.frame(Class_mat$`Enet[alpha=0.1]`,Class_mat$outcome)[193:297,]
aaa<-data.frame(Class_mat$`Enet[alpha=0.1]`,Class_mat$outcome)[298:361,]
aaa<-data.frame(Class_mat$`Enet[alpha=0.1]`,Class_mat$outcome)[362:383,]
a1<-c(1:192)
a2<-c(193:297)
a3<-c(298:361)
a4<-c(362:383)


rownames(AUC_mat)

library(runway)
single_model_dataset<-data.frame(outcomes=class$outcome,predictions=RS_list$`RF+SVM`)
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
cal_plot(single_model_dataset,
         outcome = 'outcomes',
         positive = '1',
         prediction = 'predictions')


multi_model_dataset<-data.frame(outcomes=class$outcome,predictions=RS_list$`RF+SVM`,
                                Cohort=class$Cohort)
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
               show_loess = TRUE,
               n_bins = 0
               )

cal_plot_multi(multi_model_dataset,
               outcome = 'outcomes',
               positive = '1',
               prediction = 'predictions',
               model = 'Cohort',
               n_bins = 5)







a1



#train
roc_obj <- roc(
  response = class$outcome[1:192],    # 真实标签
  predictor = RS_list$`Enet[alpha=0.1]`[a1],  # 预测概率
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


roc_obj <- roc(
  response = class$outcome[193:237],    # 真实标签
  predictor = RS_list$`Enet[alpha=0.1]`[193:237],  # 预测概率
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
  response = class$outcome[238:259],    # 真实标签
  predictor = RS_list$`Enet[alpha=0.1]`[238:259],  # 预测概率
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


