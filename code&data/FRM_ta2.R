# team assignment 2 _ part 2
# team c10
#resources 

library(data.table)
library(ggplot2)
library(MASS)
library(metRology)
library(quantmod)
library(xts)
library(dplyr)
library(metRology)
library(rugarch)


##### ------ Part 1 -------

###Get Data ----
fund <- fread('Team_10.csv')
com.nms <- c('American Tower','Mastercard',"Moody's",'SBA Communications',"O'Reilly Automotive",
             'CarMax','Roper Technologies','Markel','Capital Senior Living','Dollar Tree','Danaher')
fund$comp.name <- com.nms
slide1 <- fund[,c(3,2)]
slide1$weight <- round(slide1$weight*100,1)
#View(slide1)

# download data
for (i in 1:nrow(fund)) {
  assign(fund$ticker[i], getSymbols(fund$ticker[i],src="yahoo",
                                    from="2007-12-31",
                                    to="2020-02-29",
                                    auto.assign=FALSE))
}

# merge data
df.return <- AMT[,6]
for (i in 2:nrow(fund)) {
  df.return <- cbind(df.return,eval(as.name(fund$ticker[i]))[,6])
}

# price diff
pdiff <- data.frame(matrix(0, nrow=nrow(df.return), ncol=ncol(df.return)))
for (i in 1:11) {
  pdiff[,i] <- diff(df.return[,i])
}

# logret
ret <- data.frame(matrix(0, nrow=nrow(df.return), ncol=ncol(df.return)))
for (i in 1:11) {
  for (k in 2:nrow(df.return)) {
    ret[k,i] <- log(1+as.numeric(pdiff[k,i])/as.numeric(df.return[k-1,i]))
  }
}
ret <- as.xts(ret, order.by=index(AMT))
ret <- ret[2:nrow(ret),]
names(ret) <- fund$ticker

# Part 1 PCA ----

# calculate correlation matrix
corfun <- function(X,ntick=ncol(X),digit=2) {
  ncol  <- ncol(X)
  nobs  <- nrow(X)
  cor.X <- matrix(rep("",ncol*nobs),nrow=ncol,ncol=ncol)
  corx  <- cor(X)
  se    <- sqrt((1-corx)/(nobs-2))
  avg   <- 0
  for (i in 2:ncol) {
    avg <- avg+sum(abs(corx[i,1:(i-1)]))   #  average of |corr| 
    tst <- abs(corx[i,1:(i-1)]/se[i,1:(i-1)])
    cor.X[i,1:(i-1)] <- ifelse(tst>2,as.character(round(corx[i,1:(i-1)],digit)),"")
  }
  avg <- avg/(ncol*(ncol-1)/2)
  for (i in 1:ncol) {
    cor.X[i,i] <- "1"
  }
  rownames(cor.X) <- colnames(X)
  colnames(cor.X) <- colnames(X)
  corout <- list(cor.mat = noquote(cor.X), 
                 num.sig = (sum(cor.X!="")-ntick),
                 avg = avg)
  return(corout)  
}
cor.ret <- corfun(ret)

# correlation matrix
cor.ret$cor.mat
# 55/55 t > 2
cor.ret$num.sig
# avg cor 0.41
cor.ret$avg

pcfun <- function(X) {
  df <- as.data.frame(X)
  pc <- prcomp(X)
  vname <- seq(1:length(pc$sdev))
  eis <- data.frame( vname=vname,var=pc$sdev^2 / sum(pc$sdev^2) )
  if (nrow(eis)>10) { eis <- eis[1:10,] }
  vlevel <- as.character(seq(1:length(pc$sdev)))
  eis$vname <- factor(eis$vname,levels=vlevel)
  p <- ggplot(data=eis,aes(y=var)) +
    geom_bar(aes(x=vname),stat="identity",col="grey",alpha=0.7) + 
    ylim(0,1) +
    xlab("Principal Component")+
    ylab("Explained Variance") 
  pcx <- as.data.frame(pc$x[,1:5])
  cor.table <- matrix(rep(0,5*ncol(df)),nrow=ncol(df),ncol=5)
  nobs <- nrow(X)
  for (i in 1:ncol(df)) {
    cor.pc <- cor(df[,i],pcx)
    #    se.pc  <- sqrt((1-cor.pc)/(nobs-2))
    #    tst.pc <- abs(cor.pc/se.pc)>2
    cor.table[i,] <- cor.pc
    #    tst.table[i,] <- tst.pc
  }
  se.table <- sqrt((1-cor.table)/(nobs-2))
  tst.table<- cor.table/se.table
  cor.table <- ifelse(abs(tst.table)>2,cor.table,0)
  row.names(cor.table) <- colnames(X)
  colnames(cor.table)  <- c("pc1","pc2","pc3","pc4","pc5")
  pcout <- list(pc,eis,p,cor.table)
  return(pcout)
}

pcout <- pcfun(ret)

pcout[[3]]+ggtitle("PC of Fund 11 stocks") # PCA plot
pcout[[2]]$var[1]                          # % of var explained by PC1
round(pcout[[4]],2)                        # cor mat of PC

pc.1 <- pcout[[1]]$x[,1]                   # PC1

# get S&P500 (SPY ETF) data
spy <- data.table(getSymbols('spy',src="yahoo",
                             from="2007-12-31",
                             to="2020-02-29",
                             auto.assign=FALSE))
spy[, Pdiff := c(NA,diff(SPY.Adjusted))]
spy[, ret := as.vector( (Pdiff)/lag(SPY.Adjusted,k=1) )]
spy <- spy[-1]
spy[, logret := log(1+ret)]
spy <- spy$logret 

# correlation PC1 and SPY
round(cor(pc.1,spy),2)

# bootstrap distribution of eigenvalue 1 / total eigenvalues
nsim <- 1000
set.seed(123789)
X2 <- as.matrix(ret)
row.names(X2) <- NULL
e1vec <- rep(0,nsim)
for (i in 1:nsim) {
  Xsim <- apply(X2,2,sample)
  pc.sim <- prcomp(Xsim)
  e1vec[i] <- pc.sim$sdev[1]^2 / sum( pc.sim$sdev^2 )
}
noquote(paste("Bootstrap p-value = ",mean(e1vec>e1),sep=''))

# Part 1 PCA on residuals ----

# regress fund stocks on SPY and save residuals
ntick <- ncol(ret)
resd <- matrix(0,nrow=length(spy),ncol=ntick)
for (j in 1:ntick) {
  df <- data.frame(spy,ret[,j])
  names(df) <- c("xvar","yvar")
  reg <- lm(yvar~xvar,df)
  resd[,j] <- reg$residuals
}
colnames(resd) <- colnames(ret)
cor.res <- corfun(resd)
cor.res[[1]][1:10,1:10]

noquote(sprintf(fmt="%s %6.2f","avg of |correlation| = ",cor.res[[3]]))
noquote(sprintf(fmt="%s %i","no of correlation = ",ntick*(ntick-1)/2))
noquote(sprintf(fmt="%s %i","no of corr |t|>2 = ",cor.res[[2]]))

# results
pcout.res <- pcfun(resd)

pcout.res[[3]]+ggtitle("PC of residuals of Fund 11 stocks") # residual plot
e1.res<-pcout.res[[2]]$var[1]
noquote(sprintf(fmt="%s %6.2f","Variance explained by PC1 = ",e1.res))
z1 <- as.numeric(cor(pcout.res[[1]]$x[,1],spy))
noquote(sprintf(fmt="%s %6.2f","Corr(PC1,spy) = ",z1))

# Bootstrap
nsim <- 1000
set.seed(123789)
X2 <- as.matrix(resd)
row.names(X2) <- NULL
e1vec <- rep(0,nsim)
for (i in 1:nsim) {
  Xsim <- apply(X2,2,sample)
  pc.sim <- prcomp(Xsim)
  e1vec[i] <- pc.sim$sdev[1]^2 / sum( pc.sim$sdev^2 )
}
noquote(sprintf(fmt="%s %6.2f","Bootstrap p-value = ",mean(e1vec>e1.res)))

# cor mat PC and stocks
res.cormat <- round(pcout.res[[4]],2)
res.cormat

#####------part2--------

###data-------
ptf <- fread("Team_10.csv")
ptf = data.table(ptf)
head(ptf)
data1 <- fread("Zacks_2020-02-28.csv")
data1 = data.table(data1)
head(data1)
summary(data1)


###2a--------
#the weight of the portfolio in the six market cap categories:
# million 1,000,000
# billion 1,000,000,000
a = c("Mega","Large","Mid","Small","Micro","Nano")
b = c(300000,10000,2000,300,50,0)
c = c(rep(0,6))
six <- list(a,b,c)
names(six) <- c("capmkt","mktvalue","weight")

for (i in 1:nrow(ptf)){
  # find the ticker cap 
  j = which(data1$Ticker == ptf$ticker[i])   
  mktv = data1$SharesOut[j]*data1$Price[j]
  # find the cap mkt
  k = 1   
  if (mktv < six$mktvalue[5]) {
    k = 6
  } else if (mktv < six$mktvalue[4]) {
    k = 5
  } else if (mktv < six$mktvalue[3]) {
    k = 4
  } else if (mktv < six$mktvalue[2]) {
    k=3
  } else if (mktv < six$mktvalue[1]){
    k=2
  }
  # cap value & mkt value
  six$weight[k] = six$weight[k] + ptf$weight[i]
}
print(six)   # result 
sum(six$weight)  # check


###2b--------
# the weight of the portfolio in the 16 Zacks sectors.

a = unique(data1$ZacksSector)
b = c(rep(0,16))
sectors <- list(a,b)
names(sectors) <- c("zackssector","weight")

for (i in 1:nrow(ptf)){
  j = which(data1$Ticker == ptf$ticker[i])   # match the ticker
  k = which(data1$ZacksSector[j] == sectors$zackssector)
  sectors$weight[k] = sectors$weight[k] + ptf$weight[i]
}

#result
print(sectors$zackssector[which(sectors$weight != 0)]) 
print(sectors$weight[which(sectors$weight != 0)]) 
sum(sectors$weight)  #check


###2c--------
#the CAPM beta of the portfolio
beta1 = 0

for (i in 1:nrow(ptf)){
  j = which(data1$Ticker == ptf$ticker[i])   # match the ticker
  beta1 = beta1 + (ptf$weight[i] * data1$Beta[j])
}

beta1  # result


#####------ part3 --------
lag.dh <- function(x,k=1) {
  tmp <- as.vector(x)
  n   <- length(tmp)
  xlag <- c(rep(NA,k),tmp[1:(n-k)])
  return(xlag)
}

# read the table of mutual fund assets
rawdata <- fread("Team_10.csv")
rawdata <- as.data.frame(rawdata)
ticker <- rawdata$ticker

# download data from Yahoo!Finance
datalist <-  vector("list",1)  
for (i in 1: nrow(rawdata)) {
  datalist[[i]] <- getSymbols(ticker[i],src="yahoo",
                              from="2014-12-29",
                              to="2019-12-29",
                              auto.assign=FALSE)
  
}

pricetable <- data.frame(matrix(ncol = 11, nrow  =nrow(datalist[[1]]) ))
colnames(pricetable) <- ticker
for (i in 1: nrow(rawdata)){
  pricetable[,i] <- datalist[[i]][,6]
}

# get the date from xts object
library(timetk)
date <-tk_index(datalist[[1]], timetk_idx = FALSE, silent = FALSE)

# create function to convert data frame back to xts object
df_xts <- function (df) {
  xts_stock <- xts(df, order.by= date) 
  return(xts_stock)
}
xts_pricetable <- df_xts(pricetable)

# Now compute returns 
returns <- diff(log(xts_pricetable))
returns <- returns[-1,]
returns <- exp(returns) - 1

# Compute portfolio returns 
z <- matrix(returns,nrow=nrow(returns),ncol=11)
w <- matrix(rawdata$weight,11,1)
rp1 <- z %*% w                  # rp1 is the pro forma return of the portfolio
returns$rp1 <- log(1+rp1)

ret <- returns$rp1
names(ret) <- "logret"

# VaR and ES function
CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}

# test for normality
jarque.test(as.vector(ret))  #reject normality => not use normal distribution

# see whether there's volatility clustering 
acf(abs(as.vector(ret)),main="ACF of |daily log return|")   # there's statistically significant autocorrelation of daily log ret
temp <- sample(as.vector(ret),length(ret),replace=FALSE)
acf(temp,main="ACF of randomly permuted |logret|")    # no significant auto correlation if data is randomly permutated 

# Use ARMA-GARCH model 

### Fit an AR(1)--GARCH(1,1) model with Student t innovations ################
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch.t <- ugarchfit(spec = uspec.t, data = ret)
round(garch.t@fit$coef,6)

# simulate 10-day ahead from AR(1)-GARCH(1,1) ~ t 
RNGkind(sample.kind="Rounding")
set.seed(123789)
nper <- 10
nsim <- 100000
boot.garch <- ugarchboot(garch.t,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp",
                         cluster=NULL)
simmat <- boot.garch@fseries
sim <- apply(simmat,1,sum)
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)   # VaR is -0.0259, while ES is -0.0395
