# team assignment 3 
# team c10
#resources 
library(data.table)
library(quantmod)
library(ggplot2)
library(moments)
library(xts)
library(MASS)
library(metRology)
library(rugarch)
library(survival)
library(tidyverse)
library(plyr)

#Step 1
load("hpi_msa.rda")
load("mortgage_rates.rda")
load("ue_msa.rda")

#Step 2
load("TeamAssignment3_cdata_Q1.rda")
Data_C$MSA <- as.integer(Data_C$MSA) 
Data_C <- Data_C[MSA %in% ue_msa$MSA]
Data_C <- Data_C[MSA %in% hpi_msa$MSA]
Data_C <- Data_C[!is.na(Data_C$CSCORE_B), ]
Data_C <- Data_C[!is.na(Data_C$ORIG_VAL), ]
Data_C <- Data_C[!is.na(Data_C$ORIG_AMT), ]
Data_C <- Data_C[!is.na(Data_C$ORIG_RT), ]
nrow(Data_C)

Data_C <- Data_C[, ORIG_DTE := as.integer(year(ORIG_DTE))*100 + as.integer(substr(ORIG_DTE,6,7))]
Data_C <- merge(Data_C, rates[,c(1,2)], by.x = 'ORIG_DTE', by.y = 'yearmon', all.x=TRUE)
Data_C <- merge(Data_C, hpi_msa, by.x = c('ORIG_DTE','MSA'), by.y = c('yearmon','MSA'), all.x=TRUE)
Data_C$hpi
names(Data_C)[names(Data_C)== 'hpi'] <- 'hpi0'
Data_C$spread <- Data_C$ORIG_RT- Data_C$rate
Data_C <- Data_C[, c('LOAN_ID','OLTV','CSCORE_B','spread','ORIG_VAL','hpi0','MSA','ORIG_RT','NUM_BO','PURPOSE','PROP_TYP','OCC_STAT','DTI','FTHB_FLG')]
colnames(Data_C)


#Step 3
load("TeamAssignment3_pdata_Q1.rda")
names(Data_P)[names(Data_P)== 'Monthly.Rpt.Prd'] <- 'yearmon'
Data_P <- Data_P[,yearmon := as.integer(year(yearmon))*100 + as.integer(substr(yearmon,6,7))]
Data_P <- Data_P[LOAN_ID %in% Data_C$LOAN_ID,]
nrow(Data_P)
data1 <- merge(Data_P, Data_C, by.x = 'LOAN_ID', by.y= 'LOAN_ID',  all.x=TRUE)
setorder(data1, LOAN_ID, yearmon)
data1$status <- ifelse(data1$Zero.Bal.Code %in% c("02","03","09","15"),"default",
                       ifelse(data1$Zero.Bal.Code %in% c("01"),"prepaid","censored"))

#Step4
data1<- merge(data1, rates, by.x= 'yearmon', by.y= 'yearmon', all.x = TRUE)
data1$cvr <- data1$ORIG_RT/data1$rate

#Step5
data1<- merge(data1, ue_msa,by=c('yearmon','MSA'))

#Step6
data1<- merge(data1, hpi_msa,by=c('yearmon','MSA'))
data1$val <- data1$ORIG_VAL * data1$hpi / data1$hpi0
data1$pneq <- pnorm(log(data1$LAST_UPB/data1$val)/(100*data1$spi))

#Step7
data1$start = data1$Loan.Age
data1$end = data1$Loan.Age + 1

#Estimation
#coxph to estimate a default model
cox.d<- coxph(formula = Surv(start, end, status == "default") ~ CSCORE_B + pneq, data = data1, ties = "efron")
cox.default<- coxph(formula = Surv(start, end, status == "default") ~ CSCORE_B + pneq + spread + ue + cvr, data = data1, ties = "efron")

AIC(cox.d)
AIC(cox.default)

#coxph to estimate a prepayment model
cox.p<- coxph(formula = Surv(start, end, status == 'prepaid') ~ CSCORE_B + pneq, data = data1, ties = "efron")
cox.prepayment<- coxph(formula = Surv(start, end, status == "prepaid") ~ CSCORE_B + pneq + spread + ue + cvr, data = data1, ties = "efron")

AIC(cox.p)
AIC(cox.prepayment)

#prediction

data2<- data1[1,]
data2$CSCORE_B <- 720
data2$pneq <- 0
data2$spread<- median(data1$spread)
data2$ue<- median(data1$ue)
data2$cvr<- median(data1$cvr)
data2$MSA<- median(data1$MSA)
data2$LAST_RT<- median(data1$LAST_RT)
data2$LAST_UPB<- median(data1$LAST_UPB)
data2$OLTV<- median(data1$OLTV)
data2$ORIG_VAL<- median(data2$ORIG_VAL)
data2$hpi0<- median(data2$hpi0)
data2$ORIG_RT<- median(data1$ORIG_RT)
data2$DTI<- median(data1$DTI)
data2$rate<- median(data1$rate)
data2$hpi<- median(data1$hpi)
data2$spi<- median(data1$spi)
data2$val<- median(data1$val)

data2<- data2[rep(1,60),]
data2$start<- seq(0,59)
data2$end<- data2$start +1

#prediciting default probabilites
pred.default<- predict(cox.default, newdata = data2, type = "expected")
h.d <- 1-exp(-pred.default) # hazard rate each month given survival
s.d <- 1-h.d # survival rate each month given survival
cumsurv.d <- cumprod(s.d) # cumulative survival prob each month
cumdef.d <- 1-cumsurv.d # cumulative default prob each month
plot(cumdef.d)

#predicting prepayment probabilities
pred.prepayment<- predict(cox.prepayment, newdata = data2, type = "expected")
h.p <- 1-exp(-pred.prepayment) # hazard rate each month given survival
s.p <- 1-h.p # survival rate each month given survival
cumsurv.p <- cumprod(s.p) # cumulative survival prob each month
cumdef.p <- 1-cumsurv.p # cumulative default prob each month
plot(cumdef.p)
