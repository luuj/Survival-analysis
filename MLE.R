##########################
# BST 244 Lab 1
# 02/04/2021
# R code for Weibull MLE 
# Yang Wang
# Source: https://www.itl.nist.gov/div898/handbook/apr/section4/apr413.htm
##########################
rm(list=ls())

## Prepare the dataset
set.seed(2022) # choose your favorite seed

dataset <- read.csv("dataset.csv")
head(dataset)
dim(dataset)

##group dataset by treatment
dataset_trt <- dataset[dataset$trt == 1,]
dataset_ctrl <- dataset[dataset$trt == 0,]
dim(dataset_trt)
dim(dataset_ctrl)

# can save this dataset_sub for future use
dataset_sub <- rbind(dataset_trt[sample(1:500,400,replace = FALSE),],
                     dataset_ctrl[sample(1:500,400,replace = FALSE),])
dataset_sub_trt <- dataset_sub[dataset_sub$trt == 1,]
dataset_sub_ctrl <- dataset_sub[dataset_sub$trt == 0,]
dim(dataset_sub_trt)
dim(dataset_sub_ctrl)


## Load the package containing special survival analysis functions
# install.packages("survival")
library(survival)

## Create survival object for treatment group 
## (all the following code can be modified for control group)
y_trt = Surv(dataset_sub_trt$t2dth, dataset_sub_trt$c4dth)


## Estimate parameters for Weibull distribution
yw_trt = survreg(y_trt ~ 1, dist="weibull")
summary(yw_trt)


## Maximum likelihood estimates:
## For the Weibull model, survreg fits log(T) = log(1/lambda) + (1/p)*log(X), 
## where X has an exponential distribution with mean 1.
#   survreg's scale  =    1/(rweibull shape)   1/p
#   log scale= log (1/p)    beta
#   survreg's intercept = log(rweibull scale)  log(1/lambda)  alpha
lambda <- exp(-coefficients(yw_trt)[1])
p <- 1/yw_trt$scale
names(lambda ) <- names(p ) <- NULL
signif(c(lambda =lambda , p = p ), 6)

se_lambda=sqrt( exp(-coefficients(yw_trt)[1])^2*diag(yw_trt$var)[1] )
var_scale=yw_trt$scale^2*diag(yw_trt$var)[2]
se_p=sqrt( (1/yw_trt$scale^2)^2*var_scale )

## Survival function
t=unique(sort(dataset_sub_trt$t2dth))
St=exp(-(lambda*t)^p )

alpha=coefficients(yw_trt)[1]; beta=yw_trt$icoef[2]

se_St=sqrt( (St*p*(lambda*t)^(p-1)*t*exp(-alpha))^2*diag(yw_trt$var)[1]
            +(St*(lambda*t)^p*log(lambda*t)*exp(-beta))^2*diag(yw_trt$var)[2]
            +2*(St*p*(lambda*t)^(p-1)*t*exp(-alpha))*(St*(lambda*t)^p*log(lambda*t)*exp(-beta))*yw_trt$var[1,2] )
St_up=St+1.96*se_St
St_low=St-1.96*se_St

# install.packages("ggplot2")
library("ggplot2")
ggplot() +  geom_line(aes(t,y=St)) +
  geom_line(aes(x=t,y=St_up),linetype =2 ) +
  geom_line(aes(x=t,y=St_low),linetype =2 ) +
  xlab("time")+ylab("Survival estimate")


# ## Lifetime: expected value and standard deviation
#lambda_wb_trt<-lambda
#p_wb_trt<-p
#mu_wb_trt <- 1/lambda_wb_trt * gamma(1 + 1/p_wb_trt)
#var_wb_trt <- 1/lambda_wb_trt^2 * (gamma(1+2/p_wb_trt) - (gamma(1+1/p_wb_trt))^2)
#sigma_wb_trt <- sqrt(var_wb_trt)
#names(mu_wb_trt) <- names(var_wb_trt) <- names(sigma_wb_trt) <- NULL
#signif(c(mu_wb_trt=mu_wb_trt, sigma_wb_trt=sigma_wb_trt), 6)
# Probability density of fitted model.
# curve(dweibull(x, shape=p_wb_trt, scale=1/lambda_wb_trt),
#       from=0, to=mu_wb_trt+6*sigma_wb_trt, col="blue",
#       xlab="Days", ylab="Probability Density for Treatment Group")



## Fit a lognormal model
yl_trt = survreg(y_trt ~ 1, dist="lognormal")
summary(yl_trt)


## Fit an exponential model
ye_trt = survreg(y_trt ~ 1, dist="exponential")
summary(ye_trt)
lambda_exp_trt <- exp(-coefficients(ye_trt)[1])
names(lambda_exp_trt) <- NULL
signif(c(lambda_exp_trt=lambda_exp_trt), 6)
## Sanity check with closed form solution, uncensored numbers/total observed time
sum(dataset_sub_trt$c4dth)/sum(dataset_sub_trt$t2dth)


### Weibull vs. lognormal vs. exponential models (When comparing values of AIC, smaller is better)
signif(c(weibullAIC_trt=extractAIC(yw_trt)[2],
          lognormalAIC_trt=extractAIC(yl_trt)[2],
          expAIC_trt=extractAIC(ye_trt)[2]), 5)




