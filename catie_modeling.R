# Prediction of Drug Response in Schizophrenia 
# Authors: D. Siljak and A.M. Chekroud, July 2015

#### Housekeeping
# Data manipulation and plotting
library("ggplot2"); library("RColorBrewer")
# library("readr"); library("stringr")
# Statistical packages
library("glmnet") # http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
library("caret"); library("gbm");library("rpart"); library("lars"); library("MASS"); library("e1071")
# library("mi"); library("mice")
# library("hydroGOF"); library("pROC")
# library("permute"); library("klaR")
# library("safeBinaryRegression")
# Parallel Computing (default:use n/2 cores during model build)
library("doMC")
registerDoMC(14)
library("dplyr")

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/home/CATIE/")
workDir <- getwd()
dataDir <- paste0(workDir, "/data/") 

##can this (if, else) just be an ifelse statement?
lookup <- function(nm, summary){
  index <- which(summary$nm==nm)
  if(is.na(summary$meaning[index])){
    nm
  }
  else{
    summary$meaning[index]
  }
}
  
## Read in pre-processed, imputed data 

CATIE_ind_vars <- read.csv(paste(workDir, "/CATIE_imputed_data.csv", sep=""),
                           stringsAsFactors = FALSE, check.names=FALSE)
CATIE_ind_vars[,1] <- NULL
CATIE_dep_vars <- read.csv(paste(workDir, "/CATIE_dep_vars.csv", sep=""),
                           stringsAsFactors = FALSE, check.names=FALSE)
CATIE_dep_vars[,1] <- NULL

#   sort by subject id
CATIE_ind_vars <- CATIE_ind_vars %>% dplyr::arrange(src_subject_id)
CATIE_dep_vars <- CATIE_dep_vars %>% dplyr::arrange(src_subject_id)


CATIE_dep_vars$outcome <- CATIE_dep_vars$dcr_eff1+2*CATIE_dep_vars$dcr_tae1+3*CATIE_dep_vars$dcr_pat1
CATIE_dep_vars$total   <- ifelse(CATIE_dep_vars$dcr_eff1==1|
                                   CATIE_dep_vars$dcr_tae1==1|
                                   CATIE_dep_vars$dcr_pat1==1, 1, 0)               
### Make 3 separate dfs for the 3 categories of models 
# predicting: dropouts due to lack of effectiveness; tolerance
# issues; or patient decision
effectiveness <- subset(CATIE_ind_vars, select=-c(src_subject_id))
tolerance     <- effectiveness
pdecision     <- effectiveness
effectiveness$dcr_eff1 <- CATIE_dep_vars$dcr_eff1
tolerance$dcr_tae1 <- CATIE_dep_vars$dcr_tae1
pdecision$dcr_pat1 <- CATIE_dep_vars$dcr_pat1

### Convert to numeric matrices as lars takes only those
predictors <- CATIE_ind_vars
predictors$src_subject_id <- NULL
predictors <- predictors %>% as.matrix() %>% scale() 

 ### Fit logistic regression (breaks)
# glm_eff <- glm(dcr_eff1 ~., data=effectiveness, family="binomial")
# glm_tol <- glm(dcr_tae1 ~., data=tolerance, family="binomial")
# glm_pat <- glm(dcr_pat1 ~., data=pdecision, family="binomial")
# 
### Fit recursive partitioning tree
part_eff <- rpart(dcr_eff1 ~., data=effectiveness)
part_tol <- rpart(dcr_tae1 ~., data=tolerance)
part_dec <- rpart(dcr_pat1 ~., data=pdecision)

### Fit elastic net models

## This is a collection of functions I have written. I cloned the directory from my github.
## You may need to install packages locally
source("/home/CATIE/functions/functions.R")

eff_target <- effectiveness$dcr_eff1 %>% factor(labels = c("neg", "pos"))
tol_target <- tolerance$dcr_tae1 %>% factor(labels = c("neg", "pos"))
dec_target <- pdecision$dcr_pat1 %>% factor(labels = c("neg", "pos"))

### Build df for efficacy best model
eff_results <- data.frame(variable=rep(NA,20), beta=rep(0,20), direction=rep(NA, 20))
tol_results <- data.frame(variable=rep(NA,20), beta=rep(0,20), direction=rep(NA, 20))
dec_results <- data.frame(variable=rep(NA,20), beta=rep(0,20), direction=rep(NA, 20))
## check these targets

fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           savePredictions = TRUE, summaryFunction = custom.summary)
eGrid      <- expand.grid(alpha = seq(0,1,by=0.5), .lambda = seq(0,2,by=0.2))
eGrid2     <- expand.grid(alpha = seq(0,1,by=0.1), .lambda = seq(0,2,by=0.2))

train.index <- createDataPartition(y=CATIE_ind_vars$treat11a, p=0.7, list=FALSE)
train_data <-  effectiveness[train.index, ]
test_data <- effectiveness[-train_index,]
set.seed(1)
eff.enet <- train(x = predictors, y = eff_target,
                  method = "glmnet",
                  metric = "ROC",
                  tuneGrid = eGrid,
                  trControl = fitControl)
getTrainPerf(eff.enet)
coefs <- varImpRanker(eff.enet, 20)
eff_results$variable <- sapply(coefs[[1]], function(x) lookup(x, summarystats))
eff_results$beta     <- coefs[[2]]
eff_results$direction <- coefs[[3]]
write.table(eff_results, file=paste(workDir, "/eff_results.csv", sep=""), sep=",")

set.seed(1)
tol.enet <- train(x = predictors, y = tol_target,
                  method = "glmnet",
                  metric = "ROC",
                  tuneGrid = eGrid,
                  trControl = fitControl)
getTrainPerf(tol.enet)
coefs <- varImpRanker(tol.enet, 20)
tol_results$variable <- sapply(coefs[[1]], function(x) lookup(x, summarystats))
tol_results$beta     <- coefs[[2]]
tol_results$direction <- coefs[[3]]
write.table(tol_results, file=paste(workDir, "/tol_results.csv", sep=""), sep=",", row.names=FALSE)

set.seed(1)
dec.enet <- train(x = predictors, y = dec_target,
                  method = "glmnet",
                  metric = "ROC",
                  tuneGrid = eGrid,
                  trControl = fitControl)
getTrainPerf(dec.enet)
coefs <- varImpRanker(dec.enet, 20)
dec_results$variable <- sapply(coefs[[1]], function(x) lookup(x, summarystats))
dec_results$beta     <- coefs[[2]]
dec_results$direction <- coefs[[3]]
write.table(dec_results, file=paste(workDir, "/dec_results.csv", sep=""), sep=",")

### Fit lasso
lasso_eff  <- lars(predictors, CATIE_dep_vars$dcr_eff1)
lasso_tol  <- lars(predictors, CATIE_dep_vars$dcr_tae1)
lasso_pat  <- lars(predictors, CATIE_dep_vars$dcr_pat1)

### Scale variables to prepare for ridge
scaled_eff <- as.data.frame(scale(effectiveness))
scaled_tol <- as.data.frame(scale(tolerance))
scaled_pat <- as.data.frame(scale(pdecision))

### Fit ridge
ridge_eff <- lm.ridge(dcr_eff1 ~., data=scaled_eff)
ridge_tol <- lm.ridge(dcr_tae1~., data=scaled_tol)
ridge_pat <- lm.ridge(dcr_pat1~., data=scaled_pat)


##New outcome measure: medication 

###10fold cross-validaiton
cv_eff <- cv.lars(predictors, CATIE_dep_vars$dcr_eff1)
cv_tol <- cv.lars(predictors, CATIE_dep_vars$dcr_tae1)
cv_pat <- cv.lars(predictors, CATIE_dep_vars$dcr_pat1)

### Build df for efficacy best model
eff_results <- data.frame(variable=rep(NA,20), beta=rep(0,20))

#find optimal fraction of L1 norm
best <- which(cv_eff$cv==min(cv_eff$cv))/100

#generate predictions using path with optimal s value
#and its coefficients (requires two calls bc lars is dumb)
efficacy_preds <- predict.lars(lasso_eff, s=best, newx=predictors,mode="fraction", type="fit")
pred_coefs <- predict.lars(lasso_eff, s=best, mode="fraction", type="coefficients")
abs_coefficients <- abs(pred_coefs$coefficients)
cs <- abs_coefficients[order(abs_coefficients, decreasing=TRUE)][1:20]
eff_results$beta <- cs
eff_results$variable <- names(cs)

#get accuracy given predicted and real values
accuracy <- postResample(pred=efficacy_preds$fit, CATIE_dep_vars$dcr_eff1)


# ## Set up 5-fold cross validation structure
fitControl <- trainControl(method = "cv", number = 5, allowParallel =TRUE)
# ## Create a grid of parameters that we will try out (combinations of alpha and lambda)
eGrid          <- expand.grid(.alpha = seq(0,1,by=0.1), .lambda = seq(0,1,by=0.1))
# ## set seed again (locally)
# set.seed(1)
predictors <- scale(predictors)
## Build the model (elastic net), doing 5 fold cross validation, and also try out every parameter combination of alpha and lambda (picking the one with the best accuracy)
eff_mod <- train(x = predictors, y = eff_target, 
                                    method = "glmnet",
                                    tuneGrid = eGrid,
                                    trControl = fitControl)
# ## This function will just tell you the performance of the model
# getTrainPerf(enetModel)
# ## This function will give you a plot of the relative variable importance (ranked coefficients, scaled to 100)
# plot(varImp(firstModelYouWantToTrain, scale=TRUE), top=20, main="MyModelsVariableImportanceAveragedDuringCrossValidation")

lm3 <- lm(CATIE_dep_vars$stayed ~ I(CATIE_ind_vars$any_psy==0) + I(wholeset$treat11a==1) + CATIE_ind_vars$b1_panss + CATIE_ind_vars$b1_cgis + CATIE_ind_vars$avgcaps1_1 +
            CATIE_ind_vars$n_nmed + CATIE_ind_vars$postflg2 + CATIE_ind_vars$n_psynew)

summary(lm3)

lm4 <- lm(I(CATIE_dep_vars$outcome==2) ~ I(CATIE_ind_vars$any_psy==0) + I(wholeset$treat11a==1) + CATIE_ind_vars$b1_panss + CATIE_ind_vars$avgcaps1_1 +
            +             log(1+CATIE_ind_vars$n_nmed) + CATIE_ind_vars$postflg2 + CATIE_ind_vars$n_psynew)
summary(lm4)

# test LASSO/EN on half-samples
# feature crafting