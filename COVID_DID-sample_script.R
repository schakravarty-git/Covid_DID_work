
library(data.table)
library(plyr)
library(plm)
library(tidyverse)
library(dplyr)
library(biostat3)
library(broom)
library(multcomp)
library(lme4)
library(lmtest)
library(sandwich)
library(multiwayvcov)
 

##############################
## Global ##
## SE's clustered for the Global level estimates ##

DID_allsteps_precip_v7 <- fread("..../DID_allsteps_precip_v7.csv")
n_distinct(DID_allsteps_precip_v7$Region_ID)
#133 unique regions, 88 countries

#function for linear combination of coefficients with clustered standard errors
#function obtained from Stackexchange
LinearCombTest <- function (lmObject, vars, .vcov = NULL) {
  ## if `.vcov` missing, use the one returned by `lm`
  if (is.null(.vcov)) .vcov <- vcov(lmObject)
  ## estimated coefficients
  beta <- coef(lmObject)
  ## sum of `vars`
  sumvars <- sum(beta[vars])
  ## get standard errors for sum of `vars`
  se <- sum(.vcov[vars, vars]) ^ 0.5
  ## perform t-test on `sumvars`
  tscore <- sumvars / se
  pvalue <- 2 * pt(abs(tscore), lmObject$df.residual, lower.tail = FALSE)
  ## return a matrix
  matrix(c(sumvars, se, tscore, pvalue), nrow = 1L,
         dimnames = list(paste0(paste0(vars, collapse = " + "), " = 0"),
                         c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
}
##

## Global overall ##
model_global <- lm(defor_rate ~ postXtreat + precip + factor(Country) + factor(group) + factor(period), data = DID_allsteps_precip_v7)
summary(model_global)

coeftest(model_global, vcov=vcovCL(x=model_global, cluster=DID_allsteps_precip_v7$Region_ID))

vcv <- cluster.vcov(model_global, cluster = DID_allsteps_precip_v7$Region_ID)
LinearCombTest(model_global, c("postXtreat"), vcv)
# confirms that using LinearCombTest with vcv gives same coeffs and SE's as using coeftest

## Global by protection status
# postXtreat coefficients are estimates for non-indigenous/non-protected/non-hotspot gridcells
# Indigenous
model <- lm(defor_rate ~ postXtreat + precip + periodXindigen + postXtreatXindigen + factor(Country) + factor(group) + factor(period), data = DID_allsteps_precip_v7)
summary(model)
#### VERY IMPORTANT CLUSTERED (same as the coeftest).
glht_model <- glht(model, linfct = c("postXtreat + postXtreatXindigen = 0"), vcov=vcovCL(x=model, cluster=DID_allsteps_precip_v7$Region_ID))
summary(glht_model)

#### VERY IMP: alternative, works perfectly
vcv <- cluster.vcov(model, cluster = DID_allsteps_precip_v7$Region_ID)
LinearCombTest(model, c("postXtreat", "postXtreatXindigen"), vcv)

coeftest(model, vcov=vcovCL(x=model, cluster=DID_allsteps_precip_v7$Region_ID))

########################################################

## Bootstrapped estimates for Africa

covid_4weeks_25km_defor <- read.csv(file = ".../covid_4weeks_25km_defor.csv")

africa_25km_defor <- covid_4weeks_25km_defor[(covid_4weeks_25km_defor$continent=='Africa'),]
table(africa_25km_defor$country)
str(africa_25km_defor)
summary(africa_25km_defor)

set.seed(12345)

myfunc <- function(n,df) {      #define function
  unique_grids <- unique(n)      #unique firms
  sample_grids <- sample(unique_grids, size=length(unique_grids), replace=T ) 
  new_df <- do.call(rbind, lapply(sample_grids, function(x)  df[df$id_grid==x,] ))
}

B = 1000                             
coeffs_b_africa <- matrix(0, nrow = B, ncol = 1)

for (i in 1:B) {
  t <- myfunc(africa_25km_defor$id_grid, africa_25km_defor) 
  reg_t <- lm(Area_cum ~ postXtreat + group + period + factor(country), data = t) 
  coeffs_b_africa[i,] <- summary(reg_t)$coefficients[2]
}

#Confidence intervals
#CI5 <- c(quantile(coeffs_b_africa,0.025), quantile(coeffs_b_africa,0.975))

 ######################################################################
# Country level

length(unique(DID_allsteps_precip_v7$Country))
#88 countries

mean_precip <- DID_allsteps_precip_v7 %>%
  group_by(Country) %>%
  summarize(Mean = mean(precip, na.rm=TRUE))

DID_allsteps_precip_v7 <- DID_allsteps_precip_v7 %>% group_by(Country) %>% mutate(n_country = n())
table(DID_allsteps_precip_v7$Country)

DID_v7_country <- DID_allsteps_precip_v7[(DID_allsteps_precip_v7$n_country >= 40),]
length(unique(DID_v7_country$Country))
#75 countries with >= 40 gridcells each
table(DID_v7_country$Country, DID_v7_country$time_period)

estimates_DID_v7_country <- DID_v7_country %>% 
  group_by(Country) %>%
  do({model = lm(defor_rate ~ postXtreat + precip + group + period, data=.) 
  data.frame(tidy(model),             
             glance(model))})  

table(estimates_DID_v7_country$term)

