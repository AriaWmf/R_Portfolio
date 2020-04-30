# Lecture 5
library(data.table)
library(quantmod)
library(ggplot2)
library(moments)
library(xts)

setwd("D:/MQM530/DJI/")

# DJIA tickers
dji <- fread("DJI.csv")
Tickers <- dji$Ticker[c(1:7,9:30)]  # exclude "DOW" which starts in 2019

ntick <- length(Tickers)
start_date <- "2008-03-31"   # "V" started on 2008-03-18
end_date <- "2018-12-31"

# retrieve data from Yahoo!Finance
num <- 0
for (i in 1:ntick) {
  tkr <- Tickers[i]
  dat <- getSymbols(tkr,src="yahoo",from=start_date,to=end_date,auto.assign=FALSE)
  print(dat[1,])
  if (num == 0) {
    All <- dat[,6]
    num <- 1
  } else {
    All <- merge(All,dat[,6],join="outer")
  }      
  print(tkr)
}
names(All) <- as.character(Tickers[1:ntick])


####------------- calculate log returns------------
X <- diff(log(All))
X <- X[-1,]

####------------- calculate correlation matrix------------

# calculate correlation matrix
corfun <- function(X,digit=2) {
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
  corout <- list(noquote(cor.X),(sum(cor.X!="")-ntick),avg)
  return(corout)  
}

corout.X <- corfun(X)

# slide 11
noquote(sprintf(fmt="%s %i","no of correlation = ",ntick*(ntick-1)/2))
noquote(sprintf(fmt="%s %i","no of corr |t|>2  = ",corout.X[[2]]))
noquote(sprintf(fmt="%s %6.2f","avg correlation = ",corout.X[[3]]))

# correlation matrix for first 10 DJI compnents
corout.X[[1]][1:10,1:10]   #crrelation matrix


####-------------PCA------------

# principal component analysis
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

# slide 12

###  pcout <- list(pc,eis,p,cor.table)  
###  pcout contain different results 
pcout <- pcfun(X)
pcout[[3]]+ggtitle("PC of DJI 29 stocks")
e1<-pcout[[2]]$var[1]                      
round(pcout[[4]],2)
pc.1 <- pcout[[1]]$x[,1]          # this is the first principal component
pc.1


####-------------------slide13. bootstrap distribution-----------------------

load("vfinx_2018.rda")
F <- vfinx["2008-04-01/2018-12-28"]
z1 <- as.numeric(cor(pc.1,F))
noquote(sprintf(fmt="%s %6.2f","Corr(PC1,vfinx) = ",z1))


# slide 13
# bootstrap distribution of eigenvalue 1 / total eigenvalues
nsim <- 1000
set.seed(123789)
X2 <- as.matrix(X)
row.names(X2) <- NULL
e1vec <- rep(0,nsim)
system.time(
  for (i in 1:nsim) {
    Xsim <- apply(X2,2,sample)
    pc.sim <- prcomp(Xsim)
    e1vec[i] <- pc.sim$sdev[1]^2 / sum( pc.sim$sdev^2 )
  }
)
noquote(sprintf(fmt="%s %6.2f","Bootstrap p-value = ",mean(e1vec>e1)))
## p-value: Bootstrap p-value is 0 ??? Answer: Yes
## the first PC is very important statistically.


####------------- residuals of 1-factor model -----------
# slide 14
# regress DJI stocks on vfinx and save residuals
resd <- matrix(0,nrow=length(F),ncol=ntick)
for (j in 1:ntick) {
  df <- data.frame(F,X[,j])
  names(df) <- c("xvar","yvar")
  reg <- lm(yvar~xvar,df)
  resd[,j] <- reg$residuals
}
colnames(resd) <- colnames(X)
cor.res <- corfun(resd)
cor.res[[1]][1:10,1:10]

noquote(sprintf(fmt="%s %6.2f","avg of |correlation| = ",cor.res[[3]]))
noquote(sprintf(fmt="%s %i","no of correlation = ",ntick*(ntick-1)/2))
noquote(sprintf(fmt="%s %i","no of corr |t|>2 = ",cor.res[[2]]))

####-----run PCA on residuals--slide 15-------------------

pcout.res <- pcfun(resd)
pcout.res[[3]]+ggtitle("PC of residuals of DJI 29 stocks")
e1.res<-pcout.res[[2]]$var[1]
noquote(sprintf(fmt="%s %6.2f","Variance explained by PC1 = ",e1.res))
z1 <- as.numeric(cor(pcout.res[[1]]$x[,1],F))
noquote(sprintf(fmt="%s %6.2f","Corr(PC1,vfinx) = ",z1))


# bootstrap distribution of eigenvalue 1 / total eigenvalues
nsim <- 1000
set.seed(123789)
X2 <- as.matrix(resd)
row.names(X2) <- NULL
e1vec <- rep(0,nsim)
system.time(
  for (i in 1:nsim) {
    Xsim <- apply(X2,2,sample)
    pc.sim <- prcomp(Xsim)
    e1vec[i] <- pc.sim$sdev[1]^2 / sum( pc.sim$sdev^2 )
  }
)
noquote(sprintf(fmt="%s %6.2f","Bootstrap p-value = ",mean(e1vec>e1.res)))

# slide 16
round(pcout.res[[4]],2)
