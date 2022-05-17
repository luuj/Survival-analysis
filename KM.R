##########################
# BST 244 Lab 2
# 02/15/20
# R code for Kaplan-Meier estimator and log-rank tests
# Yang Wang
# Source:
# https://www.ncbi.nlm.nih.gov/pubmed/21611958
# https://www.rdocumentation.org/packages/survival/versions/2.43-3/topics/print.survfit
# https://rdrr.io/cran/survminer/f/vignettes/Specifiying_weights_in_log-rank_comparisons.Rmd
# https://cran.r-project.org/web/packages/survminer/vignettes/Specifiying_weights_in_log-rank_comparisons.html
##########################
rm(list=ls())
## Prepare the dataset
set.seed(1234) # choose your favorite seed

dataset<-read.csv("/Users/yangwang/Desktop/BST244TA/Besdataset/dataset.csv")

head(dataset)
dim(dataset)

dataset_trt <- dataset[dataset$trt == 1,]
dataset_ctrl <- dataset[dataset$trt == 0,]
dim(dataset_trt)
dim(dataset_ctrl)

# can save this dataset_sub for future use
dataset_sub <- rbind(dataset_trt[sample(1:500,400,replace = FALSE),],
                     dataset_ctrl[sample(1:500,400,replace = FALSE),])
dim(dataset_sub)

dataset_sub_trt <- dataset_sub[dataset_sub$trt == 1,]
dataset_sub_ctrl <- dataset_sub[dataset_sub$trt == 0,]
dim(dataset_sub_trt)
dim(dataset_sub_ctrl)


## Load the package containing special survival analysis functions
library(survival)
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("survminer")
library(survminer)
library(survMisc)


## Create survival object for sub treatment group
## (all the following code can be modified for control group)
y_trt = Surv(dataset_sub_trt$t2dth, dataset_sub_trt$c4dth)
y_trt


## Fit survival data using Kaplan-Meier
y_KM_trt = survfit(y_trt ~ 1, type="kaplan-meier")
summary(y_KM_trt)

print(y_KM_trt)

print(y_KM_trt, rmean = 1496)


## Generate Kaplan-Meier survival curve
plot(y_KM_trt, xlab="Days", ylab="Survival Probability for Treatment Group",
     main = "Kaplan-Meier Survival Curve")
abline(h=0.5,col=2)

## Fit survival data using Fleming-Harrington
y_FH_trt = survfit(y_trt ~ 1, type="fleming-harrington")
summary(y_FH_trt)
print(y_FH_trt, rmean = 1496)

## Generate Fleming-Harrington survival curve
plot(y_FH_trt, xlab="Days", ylab="Survival Probability for Treatment Group",
     main = "Fleming-Harriton Survival Curve")
abline(h=0.5,col=2)



## Log-rank test comparing two groups
survdiff(Surv(t2dth, c4dth) ~ trt, data = dataset_sub)

## Change different weights
fit <- survfit(Surv(t2dth, c4dth) ~ trt, data = dataset_sub)
fit
# Original log-rank test (equal weight)
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE)

ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv")
# Sanity check
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "1")
# Gehan-Breslow (generalized Wilcoxon) 
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "n")
# Peto-Peto’s modified survival estimate
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "S1")
# modified Peto-Peto’s (by Andersen)
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "S2")
# Fleming-Harrington (p=1, q=1)
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "FH_p=1_q=1")
# Tharone-Ware
ggsurvplot(fit, data = dataset_sub, pval = TRUE, pval.method = TRUE, surv.median.line = "hv",
           log.rank.weights = "sqrtN")


### Display all weighted log-rank results
 tenfit <- ten(Surv(t2dth, c4dth) ~ trt, data=dataset_sub)
 tenfit
 null_dev <- capture.output(comp(tenfit))
 null_dev
 tests <- attributes(tenfit)$lrt
 tests


