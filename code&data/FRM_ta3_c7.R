####################################### step 1 #################################
# read Freddie monthly HPI
rm(tita)
library(data.table)
library(zoo)

hpi_master <- fread("http://www.freddiemac.com/fmac-resources/research/docs/fmhpi_master_file.csv")
hpi_master[, yearmon := Year*100+Month]
hpi_msa   <- hpi_master[GEO_Type=="CBSA",]
setorderv(hpi_msa,c("GEO_Code","yearmon"))
hpi_msa[, c("GEO_Type","GEO_Name","Index_SA","Year","Month") := NULL]
names(hpi_msa) <- c("MSA","hpi","yearmon")
hpi_msa[, MSA := as.integer(MSA)]
hpi_msa[, difl := c(0,diff(log(hpi))),by="MSA"]
hpi_msa[, spi := sqrt(12)*c(rep(0,23),rollapply(difl,24,sd)), by="MSA"]
hpi_msa[, difl := NULL]

save(hpi_msa,file="hpi_msa.rda")

##################################

library(quantmod)

# get national 30Y mortgage rates from FRED
mrates <- getSymbols("MORTGAGE30US",src="FRED",auto.assign=FALSE)
mrates <- na.omit(mrates)
str(mrates)
data <- data.table(Date=index(mrates),RT=mrates)
names(data) <- c("Date","RT")
z  <- rollapply(data$RT,4,mean)
data$ma4 <- c(NA,NA,NA,z)
z <- to.period(data,
               period='months',
               OHLC = FALSE)
rates <- data.table(yearmon = index(z),
                    rate    = z$ma4)
names(rates) <- c("yearmon","rate")
rates$yearmon <- as.integer(floor(as.numeric(as.character(format(rates$yearmon,"%Y%m%d")))/100))
setorderv(rates,c("yearmon"))

save(rates,file="mortgage_rates.rda")

#############################
# read unemployment rates and create MSA unemployment rate 
library(zoo)
url1 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU00-04"
url2 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU05-09"
url3 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU10-14"
url4 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU15-19"
url5 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU90-94"
url6 <- "https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU95-99"

dat1 <- fread(url1)
dat2 <- fread(url2)
dat3 <- fread(url3)
dat4 <- fread(url4)
dat5 <- fread(url5)
dat6 <- fread(url6)

UE <- rbind(dat1,dat2,dat3,dat4,dat5,dat6)
UE[, footnote_codes := NULL]
UE[, yearmon := as.integer(year)*100 + as.integer(substr(period,2,3))]

UE[, c("year","period") := NULL]
names(UE) <- c("series_id","ue","yearmon")

MSA <- fread("UE_MSACode.csv")
ue_msa <- UE[series_id %in% MSA$series_id,]
ue_msa <- ue_msa[MSA, on = "series_id"]
ue_msa[, c("MSA_name") := NULL]
ue_msa[, ue := as.numeric(ue)]
ue_msa[, yearmon := as.integer(yearmon)]
ue_msa[, series_id := NULL]
setorderv(ue_msa,c("MSA","yearmon"))

head(rates)
head(ue_msa)
head(hpi_msa)
####################################### step 2 #################################
load("TeamAssignment3_cdata_Q1.rda")

Data_C <- Data_C[Data_C$MSA %in% ue_msa$MSA,]
Data_C <- Data_C[Data_C$MSA %in% hpi_msa$MSA,]
Data_C <- Data_C[!is.na(Data_C$CSCORE_B),]
Data_C <- Data_C[!is.na(Data_C$ORIG_VAL),]
Data_C <- Data_C[!is.na(Data_C$ORIG_AMT),]
Data_C <- Data_C[!is.na(Data_C$ORIG_RT),]
nrow(Data_C)

Data_C$ORIG_DTE <-as.integer(format(Data_C$ORIG_DTE,"%Y%m")) ### !
hpi_msa$yearmon<-as.integer(hpi_msa$yearmon)
Data_C$MSA <- as.integer(Data_C$MSA)

Data_C<-merge(Data_C,rates[,c(1,2)],by.x = "ORIG_DTE",by.y = "yearmon",all.x = TRUE)
Data_C<-merge(Data_C,hpi_msa, by.x =  c("ORIG_DTE","MSA"), by.y = c("yearmon","MSA"), all.x=TRUE)
Data_C$hpi0<-Data_C$hpi

Data_C$spread<-Data_C$ORIG_RT - Data_C$rate

Data_C <- Data_C[,c("LOAN_ID","OLTV","CSCORE_B","spread","ORIG_VAL","hpi0","MSA",
                 "ORIG_RT","NUM_BO","PURPOSE","PROP_TYP","OCC_STAT","DTI","FTHB_FLG")]
Data_C

#############################step3

load("TeamAssignment3_pdata_Q1.rda")

which(colnames(Data_P) == "Monthly.Rpt.Prd")
colnames(Data_P)[2] <- "yearmon"
colnames(Data_P)
head(Data_P$yearmon)

Data_P$yearmon <-as.integer(format(Data_P$yearmon,"%Y%m"))
Data_P <- Data_P[Data_P$LOAN_ID %in% Data_C$LOAN_ID,]

nrow(Data_P)

data1 <-merge(Data_P,Data_C,by = "LOAN_ID",all.x = TRUE)

data1 <- data1[order(data1$LOAN_ID,data1$yearmon),]

data1$status <- ifelse(data1$Zero.Bal.Code %in% c("02","03","09","15"),"default",
                       ifelse(data1$Zero.Bal.Code %in% c("01"),"prepaid","censored"))


#############################step4

data1<-merge(data1,rates,by = "yearmon")
data1$cvr <- data1$ORIG_RT / data1$rate

###################step5

data1 <- merge(data1,ue_msa, by=c("yearmon", "MSA"))

####################step6
data1 <- merge(data1,hpi_msa, by=c("yearmon", "MSA"))
data1$val <- data1$ORIG_VAL * data1$hpi
data1$pneq = pnorm(log(data1$LAST_UPB/data1$val)/(100*data1$spi))

###################step7

data1$start <- data1$Loan.Age
data1$end <- data1$Loan.Age+1

################################################################################
##################Use coxph to estimate a default model, and a prepayment model
library(survival)

unique(data1$PURPOSE)
unique(data1$PROP_TYP)
unique(data1$OCC_STAT)
unique(data1$FTHB_FLG)

# create dummy variables
data1$PURPOSER <- ifelse(data1$PURPOSE =="R",1,0)
data1$PURPOSEP <- ifelse(data1$PURPOSE =="P",1,0)
data1$PURPOSEC <- ifelse(data1$PURPOSE =="C",1,0)

data1$PROP_TYPSF <- ifelse(data1$PROP_TYP == "SF",1,0)
data1$PROP_TYPCO <- ifelse(data1$PROP_TYP == "CO",1,0)
data1$PROP_TYPPU <- ifelse(data1$PROP_TYP == "PU",1,0)
data1$PROP_TYPMH <- ifelse(data1$PROP_TYP == "MH",1,0)
data1$PROP_TYPCP <- ifelse(data1$PROP_TYP == "CP",1,0)

data1$OCC_STATP <- ifelse(data1$OCC_STAT == "P",1,0)
data1$OCC_STATI <- ifelse(data1$OCC_STAT == "I",1,0)
data1$OCC_STATS <- ifelse(data1$OCC_STAT == "S",1,0)

data1$FTHB_FLGY <- ifelse(data1$FTHB_FLG == "Y",1,0)
data1$FTHB_FLGN <- ifelse(data1$FTHB_FLG == "N",1,0)
data1$FTHB_FLGU <- ifelse(data1$FTHB_FLG == "U",1,0)

data1[,NUM_BO := as.integer(NUM_BO)]
head(data1)
data2 <- copy(data1)
data2[is.na(data2)] <- 0

mue <- median(data2$ue)
mcvr <- median(data2$cvr)
mspread <- median(data2$spread)
mnumbo <- median(data2$NUM_BO)
mdti <- median(data2$DTI)
mpropmh <- median(data2$PROP_TYPMH)
mocci <- median(data2$OCC_STATI)
mpurposep <- median(data2$PURPOSEP)
mpropcp <- median(data2$PROP_TYPCP)
mfthby <- median(data2$FTHB_FLGY)
mproppu <- median(data2$PROP_TYPPU)
moccp <- median(data2$OCC_STATP)
moltv <- median(data2$OLTV)
mpropsf <- median(data2$PROP_TYPSF)

install.packages("My.stepwise")
library(My.stepwise)
install.packages("imputeTS")
library(imputeTS)
data1$CSCORE_B <- na_mean(data1$CSCORE_B)
data1$pneq <- na_mean(data1$pneq)
data1$NUM_BO <- na_mean(data1$NUM_BO)
data1$OLTV <- na_mean(data1$OLTV)
data1$spread <- na_mean(data1$spread)
data1$ue <- na_mean(data1$ue)
data1$cvr <- na_mean(data1$cvr)
data1$DTI <- na_mean(data1$DTI)

# stepwise
my.data <- data1
my.data$status1 <- my.data$status=="prepaid"
my.variable.list <- c("OLTV","spread","ue","cvr","NUM_BO","PURPOSER","PURPOSEP","PURPOSEC","PROP_TYPSF","PROP_TYPCO","PROP_TYPPU","PROP_TYPMH","PROP_TYPCP","OCC_STATP","OCC_STATI","OCC_STATS","DTI", "FTHB_FLGY","FTHB_FLGN","FTHB_FLGU")
My.stepwise.coxph(Time = "yearmon", Status = "status1", variable.list = my.variable.list,
                  in.variable = c("CSCORE_B", "pneq"), data = my.data, sle = 0.25, sls = 0.25)

my.data$status2 <- my.data$status=="default"
my.variable.list <- c("OLTV","spread","ue","cvr","NUM_BO","PURPOSER","PURPOSEP","PURPOSEC","PROP_TYPSF","PROP_TYPCO","PROP_TYPPU","PROP_TYPMH","PROP_TYPCP","OCC_STATP","OCC_STATI","OCC_STATS","DTI", "FTHB_FLGY","FTHB_FLGN","FTHB_FLGU")
My.stepwise.coxph(Time = "yearmon", Status = "status2", variable.list = my.variable.list,
                  in.variable = c("CSCORE_B", "pneq"), data = my.data, sle = 0.25, sls = 0.25)

# prepayment model
prepayment <- coxph(formula = Surv(start, end, status == "prepaid") ~ CSCORE_B + pneq + OLTV + 
                    spread + ue + cvr + NUM_BO + PURPOSEP + PROP_TYPMH + PROP_TYPCP + PROP_TYPPU + PROP_TYPSF + OCC_STATI + OCC_STATP + DTI + FTHB_FLGY, data = data1, ties = "efron")
summary(prepayment)

# default model
default <- coxph(formula = Surv(start, end, status == "default") ~CSCORE_B+ pneq + cvr + 
                   spread + OLTV + PURPOSEP + NUM_BO + ue +PURPOSER + DTI + OCC_STATI + PROP_TYPSF + PROP_TYPMH +FTHB_FLGU+FTHB_FLGN, data = data1, ties = "efron")
summary(default)

#######################cumulative probability graph 
# cumulative prepayment probability
start <- seq(0,59)
end <- seq(1,60)
CSCORE_B <- rep(720, times = 60)
pneq <- rep(0, times = 60)
ue <- as.numeric(rep(mue, times = 60))
cvr <- as.numeric(rep(mcvr, times = 60))
spread <- as.numeric(rep(mspread, times = 60))
NUM_BO <- as.numeric(rep(mnumbo, times = 60))
DTI <- as.numeric(rep(mdti, times = 60))
PROP_TYPMH <- as.numeric(rep(mpropmh, times = 60)) 
OCC_STATI <- as.numeric(rep(mocci, times = 60)) 
PURPOSEP <- as.numeric(rep(mpurposep, times = 60)) 
PROP_TYPCP <- as.numeric(rep(mpropcp, times = 60))
FTHB_FLGY <- as.numeric(rep(mfthby, times = 60))
PROP_TYPPU <- as.numeric(rep(mproppu, times = 60))
OCC_STATP <- as.numeric(rep(moccp, times = 60)) 
OLTV <- as.numeric(rep(moltv, times = 60))
PROP_TYPSF <- as.numeric(rep(mpropsf, times = 60))

datapre <- cbind(start,end,CSCORE_B,pneq,ue,cvr,spread,NUM_BO,DTI,PROP_TYPMH,OCC_STATI,PURPOSEP,
                 PROP_TYPCP,FTHB_FLGY,PROP_TYPPU,OCC_STATP,OLTV,PROP_TYPSF)
datapre <- data.frame(datapre)
datapre$status <- "prepaid"
head(datapre)

pred <- predict(prepayment, newdata = datapre, type = "expected")
p <- 1-exp(-pred) # prepayment fraction
d <- 1-p # default & censoring fraction
cumdef <- cumprod(d) # cum default probability
cumprepay <- 1-cumdef # cum prepayment probability

library(ggplot2)
ggplot(data=datapre, aes(x=start, y=cumprepay, group=1)) +
  geom_line() +
  geom_point() +
  expand_limits(y=0) +
  xlab("Time(month)") + ylab("cumulative prepayment rate") +
  ggtitle("cumulative prepayment probability")

# cumulative default probability
start<-c(0:59)
end <- c(1:60)
CSCORE_B<-rep(720,60)
pneq <- rep(0,60)
cvr <- as.numeric(rep(median(data2$cvr),60))
spread <- as.numeric(rep(median(data2$spread),60))
OLTV <- as.numeric(rep(median(data2$OLTV),60))
PURPOSEP <- as.numeric(rep(median(data2$PURPOSEP),60))
NUM_BO <- as.numeric(rep(median(data2$NUM_BO),60))
ue <- as.numeric(rep(median(data2$ue),60))
PURPOSER <- as.numeric(rep(median(data2$PURPOSER),60))
DTI <- as.numeric(rep(median(data2$DTI),60))
OCC_STATI <- as.numeric(rep(median(data2$OCC_STATI),60))
PROP_TYPSF <- as.numeric(rep(median(data2$PROP_TYPSF),60))
PROP_TYPMH <- as.numeric(rep(median(data2$PROP_TYPMH),60))
FTHB_FLGU <- as.numeric(rep(median(data2$FTHB_FLGU),60))
FTHB_FLGN <- as.numeric(rep(median(data2$FTHB_FLGN),60))

ndata<-data.frame(cbind(start,end,CSCORE_B,pneq,cvr,spread,OLTV,PURPOSEP,NUM_BO,ue,PURPOSER,DTI,OCC_STATI,PROP_TYPSF,PROP_TYPMH,FTHB_FLGU,FTHB_FLGN))
ndata$status <- "default"

pred <- predict(default,newdata=ndata,type="expected")
h <- 1-exp(-pred) # hazard rate each month given survival
s <- 1-h # survival rate each month given survival
cumsurv <- cumprod(s) # cumulative survival prob each month
cumdef <- 1-cumsurv # cumulative default prob each month

datadef <- data.frame(ndata$start,cumdef)

ggplot(data=datadef, aes(x=start, y=cumdef, group=1)) +
  geom_line() +
  geom_point() +
  expand_limits(y=0) +
  xlab("Time(month)") + ylab("cumulative default rate") +
  ggtitle("cumulative default probability")

