# lecture 1

# NOTE: if you are using R version 3.6.0 or higher, 
#       add the following line in front of each set.seed() command
#       RNGkind(sample.kind="Rounding"))
#       otherwise the random numbers generated will not be the same as mine
#       and you will not get the same results

library(data.table)
library(ggplot2)
library(MASS)
library(metRology)
library(quantmod)
library(xts)

lag.dh <- function(x,k=1) {
  tmp <- as.vector(x)
  n   <- length(tmp)
  xlag <- c(rep(NA,k),tmp[1:(n-k)])
  return(xlag)
}

setwd("D:/MQM530/R_code/")

#######-------------------------slide.15 read data & format--------------
# read prices from Vanguard
vfiax <- fread("VFIAX_2018_price.csv")
vfiax[,Date := as.Date(Date,format="%m/%d/%Y")]
names(vfiax) <- c("Date","Price")
# read distributions from Vanguard
distr <- fread("VFIAX_2018_div.csv")
distr[,Date := as.Date(Date,format="%Y-%m-%d")]
names(distr) <- c("Date","distr")

# merge distributions into price file     ##### merge 
data.van <- distr[vfiax,on="Date"]
# set NAs to 0
data.van[is.na(data.van)] <- 0 

# slide 16
data.van[, c(1,3,2)]                          ### change column order & pick certain
  
#######-------------------------slide.17 discrete return & log return--------------
# calculate daily (discrete) returns
data.van[, Pdiff := c(NA,diff(Price))]
data.van[, ret := as.vector( (Pdiff+distr)/lag.xts(Price,k=1) )]
data.van <- data.van[-1]    ## delete first row
data.van[,c(1,3,2,4)]

# slide 18
# calculate daily log returns
data.van[, logret := log(1+ret)]
data.van[,c(1,3,2,5,6)]

#######-------------------------slide.19 data from yahoo finance-------------

# download data from Yahoo!Finance
vfiax <- getSymbols("vfiax",src="yahoo",
                        from="2017-12-29",
                        to="2019-01-02",
                        auto.assign=FALSE)
# convert to data.table
tmp <- data.table(date=index(vfiax))
data.yahoo <- cbind(tmp,data.table(vfiax)[,c(4,6)])
head(data.yahoo,3)
tail(data.yahoo,3)

# slide 20
data.yahoo[, Pdiff := c(NA,diff(VFIAX.Adjusted))]
data.yahoo[, ret.yahoo := as.vector(Pdiff/lag.xts(VFIAX.Adjusted,k=1)) ]
data.yahoo <- data.yahoo[-1]
head(data.yahoo[,c(1,2,3,5)],3)
tail(data.yahoo[,c(1,2,3,5)],3)

#######-------------------------slide.21 compare two price--------------
# compare prices 
df <- data.frame(date=data.van$Date,ret.van=data.van$ret,ret.yahoo=data.yahoo$ret.yahoo)
ggplot(df,aes(x=ret.van,y=ret.yahoo)) +
  geom_point() +
  xlab("Returns from Vanguard") +
  ylab("Returns from Yahoo!Finance")

# maximum difference
noquote(sprintf(fmt="%s %12.6f","Maximum difference = ",max( abs(df$ret.van-df$ret.yahoo) ) ))



#######------------------------------------distribution--------------

load("vfinx_2018.rda")

# run data up to end_date
start_date <- "2000-01-03"
end_date   <- "2017-12-31"
vfinx <- vfinx["2000-01-03/2017-12-31"]

#######-----------------------slide.27 plot log_return--------------
# plot log returns
df <- data.frame(Date=index(vfinx),logret=vfinx$logret)
ggplot(df,aes(x=Date,y=logret)) +
  geom_line(lwd=1) +
  ylab("VFINX daily log return") +
  xlab("")

#######-----------------------slide.29 estimate parameters of normal--------------
# N(mu,sig)
# assume 
mu <- mean(vfinx$logret)
sig<- sd(vfinx$logret)
fmt <- "%s %9.6f"
sprintf (fmt,"Mean: ",mu)
sprintf (fmt,"SD:  ",sig)

#######-----------------------slide.30 VaR & ES--------------
# calculate VaR and ES with equations
alpha <- 0.95
VaR <- qnorm(1-alpha,mu,sig)
ES  <- mu - sig*dnorm(qnorm(1-alpha,0,1),0,1)/(1-alpha)
fmt <- "%s %9.6f"
sprintf (fmt,"VaR: ",VaR)
sprintf (fmt,"ES:  ",ES)


#######-----------------------VaR & ES function--------------

# slide.30 calculate VaR and ES for 
CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}



# slide 32
# simulate from the normal distribution
nsim <- 100000
alpha <- 0.95
set.seed(123789)
sim <- rnorm(n=nsim,mean=mu,sd=sig)
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 33
# plot density of returns with normal distribution
ret_min  <- min(vfinx$logret)
ret_max  <- max(vfinx$logret)
# find limits of x-axis in plot, to nearest 0.05
plot_lim <- max(abs(ret_min),abs(ret_max))
plot_lim <- (round(plot_lim*100/5)+1)*5/100
ggplot(df,aes(logret)) +
  geom_histogram(aes(y=..density..),
                 alpha=0.8) +
  scale_x_continuous(name="Returns",limits=c(-plot_lim,plot_lim)) +
  scale_y_continuous(name="Density") + 
  stat_function(fun = dnorm,
                args = list(mean = mu, sd = sig),
                lwd=1,col = 'red') 


#######----------------------slide.35 skewness & kurtosis--------------
rdat <- as.vector(vfinx$logret)
library(moments)
fmt <- "%s %9.2f"
sprintf(fmt,"Skewness: ",skewness(rdat))
sprintf(fmt,"Kurtosis: ",kurtosis(rdat))

#######----------------------slide.36 Jarque-Bera test for normality--------------
# Jarque-Bera test for normality
jarque.test(rdat)   # p < 0.05, reject normality
 
# slide 37
#------------comparing normal and student t------------------
dt2 <- function(x,t_mu,t_sigma,t_df) {
  dt( (x-t_mu)/t_sigma, df=t_df ) / t_sigma
} 
cols <- c("N(0,1)"="black","t(df=10)"="blue","t(df=5)"="green","t(df=3)"="red")
ggplot(data.frame(x=c(-4,4)),aes(x=x)) +
  stat_function(fun=dnorm,args=list(0,1),aes(colour="N(0,1)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((10-2)/10),10),aes(colour="t(df=10)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((5-2)/5),5),aes(colour="t(df=5)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((3-2)/3),3),aes(colour="t(df=3)"),lwd=1) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Outcomes") +
  ggtitle("Normal and Student t densities") +
  scale_colour_manual("Densities", values = cols,
                      breaks=c("t(df=10)","t(df=5)","t(df=3)","N(0,1)"))

# slide 38
ggplot(data.frame(x=c(-4,4)),aes(x=x)) +
  stat_function(fun=dnorm,args=list(0,1),aes(colour="N(0,1)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((10-2)/10),10),aes(colour="t(df=10)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((5-2)/5),5),aes(colour="t(df=5)"),lwd=1) +
  stat_function(fun=dt2,args=list(0,sqrt((3-2)/3),3),aes(colour="t(df=3)"),lwd=1) +
  scale_y_continuous(name="Density") +
  scale_x_continuous(name="Outcomes") +
  ggtitle("Normal and Student t densities") +
  scale_colour_manual("Densities", values = cols,
                      breaks=c("t(df=10)","t(df=5)","t(df=3)","N(0,1)")) +
  xlim(c(-4,-2))

#-----------------------------------------------------------

#######----------------------distribution fit funtcions & tests--------------

# slide 42
## MLE for normal distribution
library(MASS)
rdat <- as.vector(vfinx$logret)
n.fit <- fitdistr(rdat,"normal")
n.fit$loglik                           # log likelihood of the MLE


# slide 43
t.fit <- fitdistr(rdat,"t")           ## MLE for t-distribution
round(t.fit$estimate,6)
t.fit$loglik
m <- t.fit$estimate[1]
s <- t.fit$estimate[2]
tdf <- t.fit$estimate[3]
t.fit$estimate
 
# slide 44                                   
# likelihood ratio test     , which one is better 
likeratio <- -2 * (n.fit$loglik-t.fit$loglik)
round(likeratio,1)     # ratio = 1127, t is better 


#######----------------------distribution AIC -------------
## abs(aic) bigger -> better
# AIC
AIC.n <- 2*2 - 2*n.fit$loglik	# normal distribution has 2 parameters
AIC.t <- 2*3 - 2*t.fit$loglik   # student-t distribution has 3 parameters
noquote(sprintf(fmt="%s %8.0f","AIC of normal = ",AIC.n))
noquote(sprintf(fmt="%s %8.0f","AIC of student t = ",AIC.t))


# slide 46
# plot density of returns with t-distribution
ret_min  <- min(vfinx$logret)
ret_max  <- max(vfinx$logret)
# find limits of x-axis in plot, to nearest 0.05
plot_lim <- max(abs(ret_min),abs(ret_max))
plot_lim <- (round(plot_lim*100/5)+1)*5/100
ggplot(df,aes(logret)) +
  geom_histogram(aes(y=..density..),
                 alpha=0.8) +
  scale_x_continuous(name="Returns",limits=c(-plot_lim,plot_lim)) +
  scale_y_continuous(name="Density") + 
  stat_function(fun = dt2,
                args = list(m, s, tdf),
                lwd=1,col = 'red') 

#######----------------------slide.47 simulating a t-distribution--------------
# simulate from the t distribution
nsim <- 100000
alpha <- 0.95
set.seed(123789)
#sim <- rt.scaled(n=nsim,df=tdf,mean=m,sd=s)
sim <- rt(nsim,tdf)*s + rep(m,nsim)
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 48
# simulate from empirical distribution
# bootstrapping from the actual distribution
nsim <- 100000
alpha <- 0.95
set.seed(123789)
sim <- sample(rdat,nsim,replace=TRUE)
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 50
# 10 day VaR and ES                     #### n-day VaR
# simulate from normal
nsim <- 100000
nper <- 10
alpha <- 0.95
sim <- rep(0,nsim)
set.seed(123789)
for (i in 1:nper) {
  sim <- sim+rnorm(nsim,mean=mu,sd=sig)
}
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 51
# simulate from student-t
nsim <- 100000
nper <- 10
alpha <- 0.95
sim <- rep(0,nsim)
set.seed(123789)
for (i in 1:nper) {
  sim <- sim+rt.scaled(n=nsim,df=tdf,mean=m,sd=s)
}
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

#######----------------------slide.52 bootstrapping -> 10 day--------------
# bootstrap from empirical 1 day
nsim <- 100000
nper <- 10
alpha <- 0.95
sim <- rep(0,nsim)
set.seed(123789)
for (i in 1:nper) {
  sim <- sim+sample(rdat,nsim,replace=TRUE)
}
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 53
# block bootstrap
nsim <- 100000
nper <- 10
alpha <- 0.95
sim <- rep(0,100000)
set.seed(123789)
posn <- seq(from=1,to=length(rdat)-nper+1,by=1)
rpos <- sample(posn,nsim,replace=TRUE)
for (i in 1:nper) {
  sim <- sim+rdat[rpos]
  rpos <- rpos+1
}
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)
