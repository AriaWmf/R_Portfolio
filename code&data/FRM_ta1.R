# R code to generate pro forma returns for Fund L
library(quantmod)
library(ggplot2)
library(moments)
library(xts)
library(data.table)
library(MASS)
library(metRology)
library(nloptr)
library(ggpubr)

### Part 1

### ----------------------------------Q1 ----------------------
logret_1<- fread("MQM530-TeamAssignment1-logret.csv")
logret <- logret_1$logret
rdat <- as.vector(logret)

# VaR (95%,1d) & ES
# fit to normal dist
n.fit <- fitdistr(rdat,"normal")
n.fit$loglik

# fit to student-t distribution for parameters
t.fit <- fitdistr(rdat,"t")
t.fit$loglik

# comparing dist
AIC.n <- 2*2 - 2*n.fit$loglik	
AIC.t <- 2*3 - 2*t.fit$loglik
AIC.t < AIC.n # t has smaller AIC, so pick student t distribution

# parameters  of t-distribution fit 
round(t.fit$estimate,6)
m <- t.fit$estimate[1]
s <- t.fit$estimate[2]
tdf <- t.fit$estimate[3]

# simulation
nsim <- 100000
alpha <- 0.95
RNGkind(sample.kind="Rounding") 
set.seed(123789)
sim <- rt(nsim,tdf)*s + rep(m,nsim)
# calculate VaR and ES for 
CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

####--------------------------------Q2 ------------------
## implied volatility of futures options 
# F is the price of the underlying asset  # price of fut
# X is the strike price of the option     # strike
# t is the time to maturity (in years)    # ttm
# r is the tbill rate (in decimal form)   #  r
# sigma is the volatility of the underlying asset # ?
# BFC is the price of the call option     # call/put=1
# BFP is the price of the put option      # call/put =2
# IVC is the implied volatility of the call option
# IVP is the implied volatility of the put option
BFC <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  C <- exp(-r*t) * (F * N1 - X * N2 )
  return(C)
}
BFP <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  NM1 <- pnorm(-d1)
  NM2 <- pnorm(-d2)
  P <- exp(-r*t) * (X * NM2 - F * NM1 )
  return(P)
}
IVC <- function(F,X,t,r,Call) {
  eval_f_C <- function(sigma) {
    return ( (Call-BFC(F,X,t,r,sigma))^2 )
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_C,
                opts=opts)
  return(es$solution)
}
IVP <- function(F,X,t,r,Put) {
  eval_f_P <- function(sigma) {
    return ( (Put-BFP(F,X,t,r,sigma))^2 )
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_P,
                opts=opts)
  return(es$solution)
}

# data
data2<- fread("MQM530-TeamAssignment1-fund_holdings.csv")
esz <- fread("ESZ2017.csv")
esh <- fread("ESH2018.csv")

# ttm
data2$Date = as.Date(data2$Date,format="%Y-%m-%d")
data2$Expiration = as.Date(data2$Expiration,format="%Y-%m-%d")
data2$'ttm' = as.numeric(data2$Expiration-data2$Date)/365 
data2$'sigma' = 0

# IVP & IVC
for (i in 1:nrow(data2)){
  F = data2$FutPrice[i]
  X = data2$Strike[i]
  t = data2$ttm[i]
  r = data2$r[i]
  m = 0
  if (data2$`Call/Put`[i] == 2) {
    Put = data2$OptPrice[i]
    m = IVP(F,X,t,r,Put)
  } else {
    Call = data2$OptPrice[i]
    m =  IVC(F,X,t,r,Call) 
  }
  data2$sigma[i]=m
}

# call visualization ----
# number of different expiration dates
table(data2$Expiration[data2$`Call/Put`==1])
table(data2$Expiration[data2$`Call/Put`==2])

df1 = data2[data2$Expiration=='2017-11-20' & data2$`Call/Put`==1]
c1120 <- ggplot(df1,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df1 = data2[data2$Expiration=='2017-11-30' & data2$`Call/Put`==1]
c1130 <- ggplot(df1,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df1 = data2[data2$Expiration=='2017-12-15' & data2$`Call/Put`==1]
c1215 <- ggplot(df1,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df1 = data2[data2$Expiration=='2017-12-29' & data2$`Call/Put`==1]
c1229 <- ggplot(df1,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

cplot <- ggarrange(c1120, c1130, c1215, c1229,
          labels = c("11-20", "11-30", "12-15", '12-29'),
          ncol = 2, nrow = 2, font.label=list(size=10))

# put visualization ----
df2 = data2[data2$Expiration=='2017-11-20' & data2$`Call/Put`==2]
p1120 <- ggplot(df2,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df2 = data2[data2$Expiration=='2017-11-30' & data2$`Call/Put`==2]
p1130 <- ggplot(df2,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df2 = data2[data2$Expiration=='2017-12-15' & data2$`Call/Put`==2]
p1215 <- ggplot(df2,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

df2 = data2[data2$Expiration=='2017-12-29' & data2$`Call/Put`==2]
p1229 <- ggplot(df2,aes(x=Strike,y=sigma)) + 
  geom_line() +
  xlab("strike price") +
  ylab("implied volatilities")

pplot <- ggarrange(p1120, p1130, p1215, p1229,
                   labels = c("11-20", "11-30", "12-15", '12-29'),
                   ncol = 2, nrow = 2, font.label=list(size=10))

# plot of call option group by expiration date
cplot

# plot of call option group by expiration date
pplot

###-------------------------------- Q3 ----
vol.inc <- c(1.5, 2, 3, 4)
new.mv <- c(0,0,0,0)

for (i in 1:4) {
  tmp <- data.table(data2)
  tmp[, sigma := sigma*vol.inc[i]] # increase in imp vol
  
  for (k in 1:nrow(tmp)) {
    if (tmp$`Call/Put`[k]==1) {
      tmp[k, new.price := BFC(FutPrice,Strike,ttm,r,sigma)]
    } else {
      tmp[k, new.price := BFP(FutPrice,Strike,ttm,r,sigma)]
    }
  }
  
  tmp[, mv := Contracts*new.price*250] # new market value
  new.mv[i] <- sum(tmp$mv) # new total value of options
}

# new total fund value
total.value <- new.mv+628226078

imp.chg <- 2.4

decrease.half <- function(imp.chg) {
  tmp <- data.table(data2)
  tmp[, sigma := sigma*imp.chg] # increase in imp vol
  
  for (k in 1:nrow(tmp)) {
    if (tmp$`Call/Put`[k]==1) {
      tmp[k, new.price := BFC(FutPrice,Strike,ttm,r,sigma)]
    } else {
      tmp[k, new.price := BFP(FutPrice,Strike,ttm,r,sigma)]
    }
  }
  
  tmp[, mv := Contracts*new.price*250] # new market value
  diff <- abs(-sum(tmp$mv)/628226078-0.5) # object to minimize
  return(diff)
}

opt.par <- optim(imp.chg, decrease.half, method='Nelder-Mead')

# increase 123% for value of fund to decrease 50%
opt.par$par

# ------------------------------Part 2 ----
# portfolio tickers 
Tickers <- c("AAP","CCK","CPRT","CSCO","DECK","EBAY","GOOGL","JCP",
             "KO", "MSI","PFE", "PG",  "SWKS", "TAP", "WBA")
# portfolio weights
Weights <- c(0.04910,0.06043,0.06169,0.10462,0.02770,0.07957,0.10173,0.00627,
             0.07894,0.07918,0.08347,0.06169,0.05792,0.06081,0.08687)

start_date <- "2012-12-31"
end_date <- "2017-12-31"

# retrieve data from Yahoo!Finance
tkr <- Tickers[1]
dat <- getSymbols(tkr,src="yahoo",from=start_date,to=end_date,auto.assign=FALSE)
All <- dat[,6]
ntick <- length(Tickers)
for (tkr in Tickers[2:ntick]) {
  dat <- getSymbols(tkr,src="yahoo",from=start_date,to=end_date,auto.assign=FALSE)
  All <- merge(All,dat[,6],join="outer")
}
names(All) <- as.character(Tickers[1:ntick])

# create returns for each stock in the portfolio
returns <- diff(log(All))
returns <- returns[-1,]
returns <- exp(returns) - 1

# create the pro forma returns for dates with all stock returns
z <- matrix(returns,nrow=nrow(returns),ncol=ntick)
w <- matrix(Weights,ntick,1)
rp1 <- z %*% w                                  # rp1 is the pro forma return of the portfolio
returns$rp1 <- log(1+rp1)

ret <- returns$rp1
names(ret) <- "logret"

# Q4 ----
ret_vec <- as.vector(ret)
jarque.test(ret_vec)
# we can see that since p-value is very low, we reject the assumption of normal distribution

# Q5 ----
acf(abs(ret_vec),main="ACF of |daily log return|")

# randomized order shows no autocorrelation
RNGkind(sample.kind="Rounding")
set.seed(123789)
rvec <- sample(abs(ret_vec),length(abs(ret_vec)),replace=FALSE)
acf(rvec,main="ACF of randomly permuted |logret|")

# From the acf plot, we can see that in most cases there is strong evidence of autocorrelation in absolute value of daily log returns
# And this indicates volatility clustering and contional volatility 

# Q6 ----
# Specify the function that returns VaR and ES
CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}

#Fit an AR(1)--GARCH(1,1) model with Student t innovations 
library(rugarch)
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch.t <- ugarchfit(spec = uspec.t, data = ret_vec)
round(garch.t@fit$coef,6) 

# simulate 1-day ahead from AR(1)-GARCH(1,1) ~ t
RNGkind(sample.kind="Rounding")
set.seed(123789)
nper <- 1
nsim <- 100000
boot.garch <- ugarchboot(garch.t,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp")
sim <- boot.garch@fseries

alpha <- 0.95
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)
