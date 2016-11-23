install.packages("vcd")
require(vcd)
require(MASS)

gene <- read.delim("~/Dropbox1/Dropbox/CMPE 239/Project/Dataset/GENETOX_Bacterial mutagenicity_NTP (1)/GENETOX_Bacterial mutagenicity_NTP.txt")

# take subset of data only for year 2015
gene.subset <- droplevels(subset(gene,  format(as.Date(gene$START_DATE),"%y")== "15"  ))
View(gene.subset)
dim(gene.subset)
#[1] 5260   25
str(gene.subset)

# Convert chemical name, study title, subject name,accession number to charachers
# start date to date
gene.subset$CHEMICAL_NAME <- as.character(gene.subset$CHEMICAL_NAME)
gene.subset$STUDY_TITLE <- as.character(gene.subset$STUDY_TITLE)
gene.subset$START_DATE <- as.Date(gene.subset$START_DATE)
gene.subset$SUBJECT_NAME <- as.character(gene.subset$SUBJECT_NAME)
gene.subset$ACCESSION_NUMBER <- as.character(gene.subset$ACCESSION_NUMBER)
gene.subset$DEPOSITOR_STUDY_NUMBER <- as.character(gene.subset$DEPOSITOR_STUDY_NUMBER)
gene.subset$ORGANIZATION_NAME <- as.character(gene.subset$ORGANIZATION_NAME)

# create two new variables to classify as (non)mutagen 
mutagen <- 1
nonmutagen <- 0


# data looks exponentially distrbuted
hist(gene.subset$DOSE, breaks=500)
hist(gene.subset$COUNT, breaks=50)

str(gene.subset)

# Count has NA's , for all the rows that have NA are replace with COUNT_MEAN value

gene.subset$newCount <- ifelse(is.na(gene.subset$COUNT), gene.subset$COUNT_MEAN, gene.subset$COUNT)

dose <- gene.subset$DOSE
count <-  gene.subset$newCount

cor(dose,count, use="complete")
#-0.2706109, they shuld not be highly corelated which is true in this case

# There are 74 entries with no COUNT, COUNT_SEM, COUNT_MEAN
# Not sure if we can remove these rows as there count is not there
# for now removing these rows assuming they provide no information without COUNT, but need to look into it again
# remove rest of the rows with no count
dim(gene.subset)
gene.subset <-gene.subset[-which(is.na(count)),]
dim(gene.subset)
#[1] 5186   26

mean( gene.subset$DOSE)
#[1] 1217.302
mean(gene.subset$newCount)
#[1] 184.5499
min(gene.subset$newCount) #0
max(gene.subset$newCount) #2242

# Add new column if study conclusion is Poistive , its Mutage
# if its Equivocal or Negative its Non mutagen
gene.subset$STUDY_CONCLUSION <- as.character(gene.subset$STUDY_CONCLUSION)
gene.subset$RESULT <- ifelse(gene.subset$STUDY_CONCLUSION == "Positive", 1,0)

# convert back to factor variable
gene.subset$STUDY_CONCLUSION <- as.factor(gene.subset$STUDY_CONCLUSION)
gene.subset$RESULT <- as.factor(gene.subset$RESULT)


#gene.subset$OUT <- relevel(gene.subset$STUDY_CONCLUSION, ref="Negative")
#str(gene.subset)

# Split the data to train and test
set.seed(1234)

gene.holdout.variable <- rbinom(n = dim(gene.subset)[1], size = 1, prob = .2)

# training data
training.data <- gene.subset[gene.holdout.variable==0, ]
dim(training.data)
# test data
test.data <- gene.subset[gene.holdout.variable==1, ]
dim(test.data)
#View(training.data)
#further split train data into train and validation data
train.holdout.variable <- rbinom(n = dim(training.data)[1], size = 1, prob = .2)
train.holdout.variable

# final train data
gene.train <- training.data[train.holdout.variable==0, ]
dim(gene.train)
str(gene.train)
View(gene.train)
write.csv(gene.train,"gene.train.csv")
# validation data
gene.valid <- training.data[train.holdout.variable==1, ]
dim(gene.valid)
View(gene.valid)

contrasts(gene.train$RESULT)

#################################################

plot(gene.train$STUDY_CONCLUSION, gene.train$DOSE)
  
plot(gene.train$STUDY_CONCLUSION, gene.train$COUNT) 
plot(dose,count)

#################### 1. Logistic Regression#############################
# perform logistic regression and find the best model based on AIC value

#View(gene.train)
model1 <- glm(RESULT ~ DOSE + newCount, family = binomial(link='logit'), data=gene.train)
model2 <- glm(RESULT ~ DOSE + newCount + STRAIN, family = binomial(link='logit'), data=gene.train)
model3 <- glm(RESULT ~ DOSE + newCount + STRAIN + TRIAL_RESULT,family=binomial(link='logit'), data=gene.train)
model4 <- glm(RESULT ~ DOSE + newCount + STRAIN + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model5 <- glm(RESULT ~ DOSE + newCount + STRAIN +  + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model6 <- glm(RESULT ~ DOSE + STRAIN + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model7 <- glm(RESULT ~ newCount + STRAIN +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model8 <- glm(RESULT ~ DOSE +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)
model9 <- glm(RESULT ~ DOSE + newCount + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,family=binomial(link='logit'),data=gene.train)

model10 <- glm(RESULT ~ DOSE +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED + TREATMENT_GROUP_TYPE,family=binomial(link='logit'),data=gene.train)
model11 <- glm(RESULT ~ DOSE + newCount + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED + TREATMENT_GROUP_TYPE,family=binomial(link='logit'),data=gene.train)
model12 <- glm(RESULT ~ newCount + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED + TREATMENT_GROUP_TYPE,family=binomial(link='logit'),data=gene.train)
model13 <- glm(RESULT ~  newCount + STRAIN, family = binomial(link='logit'), data=gene.train)
model14 <- glm(RESULT ~ newCount + STRAIN + TRIAL_RESULT,family=binomial(link='logit'), data=gene.train)

str(gene.train)

summary(model1)  # AIC 1836.5
summary(model2) #AIC: 1838.7

summary(model3) #AIC: 1352.3
summary(model4)# AIC: 1394.5
summary(model5) #AIC: 1394.5
summary(model6) #AIC: 1395.5
summary(model7)
summary(model8)#AIC: 1402.4
summary(model9) #AIC: 1402.6

summary(model10) #AIC: 1402.6
summary(model11) #AIC: 1404.6
summary(model12) #AIC: 1426.4
summary(model13) #AIC: 1875.2
summary(model14) #AAIC: 1431.4


# ANOVA
anova(model1, test="Chisq") # newcount, DOSE significant
anova(model2, test="Chisq") #  STRAIN not significant. DOSE , newcountsignificant
anova(model3, test="Chisq") #   STRAIN not significant. DOSE,newcount, TRIAL_RESULT significant
anova(model4, test="Chisq") #   STRAIN not significant. DOSE, TRIAL_RESULT,newcount,MICROSOMAL_ACTIVATION_USED significant
anova(model5, test="Chisq")
anova(model6, test="Chisq") # STRAIN not significant
anova(model7, test="Chisq") # newCOunt , STRAIN not dignificant
anova(model8, test="Chisq") # DOSE, TRIAL_RESULT, MICROSOMAL_ACTIVATION_USED all significant
anova(model9, test="Chisq") #DOSE, newCount, TRIAL_RESULT, MICROSOMAL_ACTIVATION_USED all significant

anova(model10, test="Chisq") #

anova(model11, test="Chisq")
anova(model12, test="Chisq")
anova(model13, test="Chisq")
anova(model14, test="Chisq")

View(gene.valid)

# based on p values and AIC model9 is significant
# validate model against valid data
fitted.results9 <- predict(model9,newdata=gene.valid,type='response')
final.fitted.results9 <- ifelse(fitted.results9 > 0.5,1,0)

fitted.results4 <- predict(model4,newdata=gene.valid,type='response')
final.fitted.results4 <- ifelse(fitted.results4 > 0.5,1,0)

# for model 9 find the accuracy , model is 94% accurate
misClasificError9 <- mean(final.fitted.results9 != gene.valid$RESULT)
print(paste('Accuracy',1-misClasificError9))
#[1] "Accuracy 0.953028430160692

###on test data###

test.results9 <- predict(model9,newdata=test.data,type='response')
test.fitted.results9 <- ifelse(test.results9 > 0.5,1,0)
# for model 9 find the accuracy , model is 94% accurate
test_misClasificError9 <- mean(test.fitted.results9 != test.data$RESULT)
print(paste('Accuracy',1-test_misClasificError9))
#[1] "Accuracy 0.943702290076336"

test.results4 <- predict(model4,newdata=test.data,type='response')
test.fitted.results4 <- ifelse(test.results4 > 0.5,1,0)
# for model 9 find the accuracy , model is 94% accurate
test_misClasificError4 <- mean(test.fitted.results4 != test.data$RESULT)
print(paste('Accuracy',1-test_misClasificError4))





contr.treatment(2)
contrasts(gene.valid$RESULT) = contr.treatment(2)

#The ROC is a curve generated by plotting the true positive rate (TPR)
#against the false positive rate (FPR) at various threshold settings 
#while the AUC is the area under the ROC curve. As a rule of thumb, 
#a model with good predictive ability should have an AUC closer to 1 (1 is ideal) than to 0.5.


install.packages("ROCR")
library(ROCR)
p <-  predict(model9,newdata=test.data,type='response')
pr <- prediction(p, test.data$RESULT)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
#[1] 0.8100697

# need to perform  cross validation such as k-fold cross validation

###################################


#################Multinomial Regression #################
installed.packages("foreign")
installed.packages("nnet")
installed.packages("ggplot2")
installed.packages("reshape2")

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

library(foreign)
library(nnet)
library(ggplot2)
library(reshape2)

str(gene.train$STUDY_CONCLUSION)

#gene.train$OUT <- relevel(gene.valid)
mmodel1 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount,data=gene.train)
summary(mmodel1)
z1 <- summary(mmodel1)$coefficients/summary(mmodel1)$standard.errors
p1 <- (1 - pnorm(abs(z1), 0, 1))*2
p1  # newCount is not significant


mmodel2 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount + STRAIN,data=gene.train)
summary(mmodel2)
z2 <- summary(mmodel2)$coefficients/summary(mmodel2)$standard.errors
p2 <- (1 - pnorm(abs(z2), 0, 1))*2
p2  #newCount is not Significant

mmodel3 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount + STRAIN + MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel3)
#AIC: 4405.767 
z3 <- summary(mmodel3)$coefficients/summary(mmodel3)$standard.errors
p3 <- (1 - pnorm(abs(z3), 0, 1))*2
p3


mmodel4 <- multinom(STUDY_CONCLUSION ~ DOSE +STRAIN + newCount+TRIAL_RESULT+MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel4)
#AIC: 3489.88
z4 <- summary(mmodel4)$coefficients/summary(mmodel4)$standard.errors
p4 <- (1 - pnorm(abs(z4), 0, 1))*2
p4

mmodel5 <- multinom(STUDY_CONCLUSION ~ DOSE  + STRAIN + MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel5)
#AIC: 4409.717
z5 <- summary(mmodel5)$coefficients/summary(mmodel5)$standard.errors
#2-tailed z test
p5 <- (1 - pnorm(abs(z5), 0, 1))*2
p5 # good

mmodel6 <- multinom(STUDY_CONCLUSION ~ DOSE +STRAIN +TRIAL_RESULT+MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel6)
#Residual Deviance: 3454.682 
#AIC: 3490.682 
z6 <- summary(mmodel6)$coefficients/summary(mmodel6)$standard.errors
#2-tailed z test
p6 <- (1 - pnorm(abs(z6), 0, 1))*2
p6 # good

mmodel7 <- multinom(STUDY_CONCLUSION ~ newCount + STRAIN +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel7)
#Residual Deviance: 3496.214
#AIC: 3532.214 
z7 <- summary(mmodel7)$coefficients/summary(mmodel7)$standard.errors
#2-tailed z test
p7 <- (1 - pnorm(abs(z7), 0, 1))*2
p7

mmodel8 <- multinom(STUDY_CONCLUSION ~ DOSE +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel8)
#Residual Deviance: 3488.153
#AIC: 3516.153
z8 <- summary(mmodel8)$coefficients/summary(mmodel8)$standard.errors
#2-tailed z test
p8 <- (1 - pnorm(abs(z8), 0, 1))*2
p8 # good

mmodel9 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount + TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,data=gene.train)
#Residual Deviance: 3486.138 
#AIC: 3518.138 
summary(mmodel9)
z9 <- summary(mmodel9)$coefficients/summary(mmodel9)$standard.errors
p9 <- (1 - pnorm(abs(z9), 0, 1))*2
p9


mmodel10 <- multinom(STUDY_CONCLUSION ~ DOSE + newCount + STRAIN +TRIAL_RESULT + MICROSOMAL_ACTIVATION_USED,data=gene.train)
summary(mmodel10)
#Residual Deviance: 3486.138 
#AIC: 3518.138 
summary(mmodel10)
z10 <- summary(mmodel10)$coefficients/summary(mmodel10)$standard.errors
p10 <- (1 - pnorm(abs(z10), 0, 1))*2
p10



# models 6 ad model 8 are good based on p value, AIC , Residual

model6.predict <- predict(mmodel6,gene.valid)
model8.predict <- predict(mmodel8,gene.valid)
model10.predict <- predict(mmodel8,gene.valid)

require(ROCR)
x.logit.prob <- predict(mmodel10, type="prob", newdata=gene.valid, probability = T)
x.logit.prob.rocr <- prediction(x.logit.prob, gene.valid$STUDY_CONCLUSION)

cm6 <- table(model6.predict, gene.valid$STUDY_CONCLUSION)
print(cm6)

attr(x.logit.prob, "probabilities")
cm8 <- table(model8.predict, gene.valid$STUDY_CONCLUSION)
print(cm8)




1- sum(diag(cm6))/sum(cm6) # 0.1625113
1- sum(diag(cm8))/sum(cm8)

# with test data

tst.mp6 <-  predict(mmodel6,newdata=test.data)
tmm6 <- table(tst.mp6, test.data$STUDY_CONCLUSION)
print(tmm6)


tst.mp8 <-  predict(mmodel8,newdata=test.data)
tmm8 <- table(tst.mp8, test.data$STUDY_CONCLUSION)
print(tmm8)


#############SVM classification ####################
library("e1071")


plot(gene.train$DOSE, gene.train$newCount, col=gene.train$STUDY_CONCLUSION)
plot(gene.train$DOSE, col=gene.train$STUDY_CONCLUSION)
plot(newCount, col=gene.train$STUDY_CONCLUSION)

?svm
svm.model <- svm(STUDY_CONCLUSION ~ DOSE, data = gene.train, cost = 100, gamma = 1)
summary(svm.model)
plot(svm.model, gene.train)
svm.pred <- predict(svm.model, gene.valid)
table(pred = svm.pred, true = gene.valid$STUDY_CONCLUSION)

View(gene.train)

svm.model1 <- svm(STUDY_CONCLUSION ~ DOSE+newCount, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model1)
plot(svm.model1)
svm.pred1 <- predict(svm.model1, gene.valid)
table(pred = svm.pred1, true = gene.valid$STUDY_CONCLUSION)

svm.model2 <- svm(STUDY_CONCLUSION ~ DOSE+STRAIN, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model2)
plot(svm.model2)
svm.pred2 <- predict(svm.model2, gene.valid)
table(pred = svm.pred2, true = gene.valid$STUDY_CONCLUSION)

svm.model3 <- svm(STUDY_CONCLUSION ~ DOSE+newCount+STRAIN, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model3)
plot(svm.model3)
svm.pred3 <- predict(svm.model3, gene.valid)
table(pred = svm.pred3, true = gene.valid$STUDY_CONCLUSION)

svm.model4 <- svm(STUDY_CONCLUSION ~ DOSE+newCount+STRAIN+MICROSOMAL_ACTIVATION_USED, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model4)
plot(svm.model4)
svm.pred4 <- predict(svm.model4, gene.valid)
table(pred = svm.pred4, true = gene.valid$STUDY_CONCLUSION)

svm.model5 <- svm(STUDY_CONCLUSION ~ DOSE+STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model5)
plot(svm.model55)
svm.pred5 <- predict(svm.model5, gene.valid)
table(pred = svm.pred5, true = gene.valid$STUDY_CONCLUSION)
str(gene.valid)

svm.model6 <- svm(STUDY_CONCLUSION ~ DOSE+STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT+TREATMENT_GROUP_TYPE, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model6)
plot(svm.model6)
svm.pred6 <- predict(svm.model6, gene.valid)
table(pred = svm.pred6, true = gene.valid$STUDY_CONCLUSION)

svm.model7 <- svm(STUDY_CONCLUSION ~ DOSE+newCount+STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, cost = 100, gamma = 1, kernel="radial")
summary(svm.model7)
plot(svm.model57)
svm.pred7 <- predict(svm.model7, gene.valid)
table(pred = svm.pred7, true = gene.valid$STUDY_CONCLUSION)
str(gene.valid)
misClasificError7 <- (82+1+7+4+25+2)/(38+82+1+7+620+4+2+25+30)
#[1] 0.1495674

svm.model8 <- svm(STUDY_CONCLUSION ~ DOSE+newCount+STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT+TREATMENT_GROUP_TYPE, data = gene.train, cost = 100, gamma = 2
                  , kernel="radial")
summary(svm.model8)
plot(svm.model8, data=gene.valid)
svm.pred8 <- predict(svm.model8, gene.valid)
table(pred = svm.pred8, true = gene.valid$STUDY_CONCLUSION)

misClasificError8 <- (71+1+11+7+3+22)/(49+71+1+11+613+7+3+22+32)
#[1] 0.1421508

#final models
svm.pred.final7 <- predict(svm.model7, test.data)
summary(svm.pred.final7)
table(pred = svm.pred.final7, true = test.data$STUDY_CONCLUSION)

svm.pred.final8 <- predict(svm.model8, test.data)
summary(svm.pred.final8)
table(pred = svm.pred.final8, true = test.data$STUDY_CONCLUSION)

misClasificError <- (81+71+6+7+1+48)/(57+81+7+16+793+7+1+48+38)
#[1] 0.2041985
# 20%

##################Random forest ######################
# CART
set.seed(1234)

library(randomForest)
?randomForest
rf1 <- randomForest(STUDY_CONCLUSION ~ DOSE+newCount, data = gene.train, importance=T, ntree=1000)
rf1
varImpPlot(rf1)

rf2 <- randomForest(STUDY_CONCLUSION ~ DOSE, data = gene.train, importance=T, ntree=1000)
rf2
varImpPlot(rf2)

rf3 <- randomForest(STUDY_CONCLUSION ~ DOSE+STRAIN, data = gene.train, importance=T, ntree=1000)
rf3
varImpPlot(rf3)


rf4 <- randomForest(STUDY_CONCLUSION ~ DOSE+newCount+ STRAIN+MICROSOMAL_ACTIVATION_USED, data = gene.train, importance=T, ntree=1000)
rf4
varImpPlot(rf4)

svm.model5 <- svm(STUDY_CONCLUSION ~ DOSE+STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, cost = 100, gamma = 1, kernel="radial")


rf7 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT, data = gene.train, importance=T, ntree=1000)
rf7
varImpPlot(rf7)


rf8 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+TREATMENT_GROUP_TYPE, data = gene.train, importance=T, ntree=1000)
rf8
varImpPlot(rf8)

rf9 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+TREATMENT_GROUP_TYPE+STRAIN, data = gene.train, importance=T, ntree=1000)
rf9
varImpPlot(rf9)

rf10 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+STRAIN, data = gene.train, importance=T, ntree=1000)
rf10
varImpPlot(rf10)

rf11 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+newCount, data = gene.train, importance=T, ntree=1000)
rf11
varImpPlot(rf11)


rf12 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+STRAIN+newCount, data = gene.train, importance=T, ntree=1000)
rf12
varImpPlot(rf12)

rf12 <- randomForest(STUDY_CONCLUSION ~ DOSE++TRIAL_RESULT+STRAIN+newCount+TREATMENT_GROUP_TYPE, data = gene.train, importance=T, ntree=1000)
rf12
varImpPlot(rf12)

rf5 <- randomForest(STUDY_CONCLUSION ~ DOSE+newCount+ STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, importance=T, ntree=1000)
rf5
varImpPlot(rf5)
#best model


####################### rpart ####################
install.packages("caret")
library(rpart)
library(rpart.plot)
library(caret)
library(doSNOW)
# prepare training scheme

set.seed(1234)
rf.label <- gene.train$STUDY_CONCLUSION
# create 10 fold cross validation
cv.10.folds  <- createMultiFolds(rf.label, k =10, times=10)

#setup caret's trainControl object
control <- trainControl(method="repeatedcv", number=10, repeats=3)

#setup doSNOW package for multi-core training , This is helpful as 
# we are going to train a lot of trees

features = c("DOSE","newCount", "STRAIN","MICROSOMAL_ACTIVATION_USED","TRIAL_RESULT")
rpart.train.1 <- gene.train[,features]

set.seed(1234)
fir.rf <- train(STUDY_CONCLUSION ~ DOSE+newCount+ STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, method="rf", tuneLength = 3,
                 trControl=control)

set.seed(1234)
fit.rpart <- train(STUDY_CONCLUSION ~ DOSE+newCount+ STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, 
                 method="rpart", tuneLength = 3,
                 trControl=control)


set.seed(1234)
fit.svm <- train(STUDY_CONCLUSION ~ DOSE+newCount+ STRAIN+MICROSOMAL_ACTIVATION_USED+TRIAL_RESULT, data = gene.train, 
                 method="rf", tuneLength = 3,
                 trControl=control)
results <- resamples(list(CART=fit.rpart, SVM=fit.svm, RF=fir.rf))

# compare the models
# summarize differences between modes
summary(results)


#Box and Whisker Plots
#This is a useful way to look at the spread of the estimated accuracies for different methods and how they relate.
# look at the mean and the max columns.
# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

#Density Plots
#show the distribution of model accuracy as density plots. 
#This is a useful way to evaluate the overlap in the estimated behavior of algorithms.

# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))
densityplot(results, scales=scales, pch = "|")

#Dot Plot
#These are useful plots as the show both the mean estimated accuracy as well as the 95% confidence interval 
#(e.g. the range in which 95% of observed scores fell).


# dot plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))
dotplot(results, scales=scales)

#Scatterplot Matrix

#This create a scatterplot matrix of all fold-trial results for an algorithm compared 
#to the same fold-trial results for all other algorithms. All pairs are compared.
# pair-wise scatterplots of predictions to compare models
splom(results)

#For example, eye-balling the graphs it looks like RF and SVM look strongly correlated, 
# SVM and CART look weekly correlated.

#Pairwise xyPlots
#You can zoom in on one pair-wise comparison of the accuracy of trial-folds for two machine learning algorithms with an xyplot.
# xyplot plots to compare models
xyplot(results, models=c("RF", "SVM"))

#In this case we can see the seemingly correlated accuracy of the RF and SVM models

#Statistical Significance Tests
#calculate the significance of the differences between the metric distributions of different machine
#learning algorithms. We can summarize the results directly by calling the summary() function.
# difference in model predictions
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs)

## calculating the values for ROC curve
#pred <- prediction(target_pred, target_class)
#perf <- performance(pred,"tpr","fpr"

