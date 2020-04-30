# option analysis in lecture 4

library(data.table)
library(ggplot2)
library(moments)
library(quantmod)
library(rugarch)
library(xts)
library(nloptr)

# Option functions
# simulate option hedging outcome
# S is the price of the underlying asset
# X is the strike price of the option
# t is the time to maturity (in years)
# r is the tbill rate (in decimal form)
# q is the dividend yield of the underlying asset, paid continuously
# sigma is the volatility of the underlying asset
# BSC is the price of the call option
# BSP is the price of the put option


#####----------- European call price  function slide.4------------------
####---- price , theta, delta, vega, gamma, ivp--------
####-------------------price----------------
BSC <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  C <- S * exp(-q*t) * pnorm(d1) - X * exp(-r*t) * pnorm(d2)
  return(C)
}
BSP <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  P <- X * exp(-r*t) * pnorm(-d2) - S * exp(-q*t) * pnorm(-d1)
  return(P)
}

###-------------- Theta function---------------
# theta of call
BSC.theta <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  thetaC <- q * S * exp(-q*t) * N1 - 0.5 * S * exp(-q*t) * sigma * dnorm(d1) / sqrt(t) - r * X * exp(-r*t) * N2
  return(thetaC)
}
# Theta of put
BSP.theta <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(-d1)
  N2 <- pnorm(-d2)
  thetaP <- r * X * exp(-r*t) * N2 - 0.5 * S * sigma * dnorm(d1) / sqrt(t) - q * S * exp(-r*t) * N1
  return(thetaP)
}
###---------------- delta function---------------

# Delta of call
BSC.delta <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N2 <- pnorm(d2)
  deltaC <- exp(-q*t) * pnorm(d1) 
  return(deltaC)
}
# Delta of put
BSP.delta <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  deltaP <- - exp(-q*t) * pnorm(-d2)
  return(deltaP)
}
###---------------- gamma function---------------

# Gamma of call = Gamma of put
BSC.gamma <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  gammaC <- exp(-q*t) * dnorm(d1) / ( S*sigma*sqrt(t) ) 
  return(gammaC)
}
BSP.gamma <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  gammaP <- exp(-q*t) * dnorm(d1) / ( S*sigma*sqrt(t) )
  return(gammaP)
}
###---------------- vega function---------------

BSC.vega <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  vegaC <- S * exp(-q*t) * dnorm(d1) * sqrt(t) 
  return(vegaC)
}
BSP.vega <- function(S,X,t,r,q,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  vegaP <- S * exp(-q*t) * dnorm(d1) * sqrt(t) 
  return(vegaP)
}
###---------------- IVC function---------------
# implied volatility 
IVC <- function(S,X,t,r,q,Call) {
  if(!require(nloptr)) install.packages('nloptr')
  library(nloptr)
  eval_f_C <- function(sigma) {
    return ( (Call-BSC(S,X,t,r,q,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_C,
                opts=opts)
  return(es$solution)
}

IVP <- function(S,X,t,r,q,Put) {
  if(!require(nloptr)) install.packages('nloptr')
  library(nloptr)
  eval_f_P <- function(sigma) {
    return ( (Put-BSP(S,X,t,r,q,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_P,
                opts=opts)
  return(es$solution)
}



########-------------slide 5--------------
# Call option on SPY
S <- 266.86
X <- 275

###slide 6 , find t 
DATE <- function(yyyy,mm,dd) {
  dte  <- as.Date(sprintf("%i-%i-%i",yyyy,mm,dd),format="%Y-%m-%d")
  return(dte)
}
t <- as.numeric(DATE(2018,3,16)-DATE(2017,12,29))/365
t         # time to maturity in years


### slide 8   find price P
bill.price <- function(settle,mature,discount.rate,redemption=100) {
  P <- redemption*(1-discount.rate*(as.numeric(mature-settle)/360) )
  return(P)
}
discount.rate  <- 0.01275     
# find t bill price, B
B <- bill.price(DATE(2018,1,2),        
                     DATE(2018,3,15),        
                     discount.rate)        

###--------------calcul ???? = ???????????????(????)/????---------
# find continuously compounded interest rate, r
r <- -log(B/100)/( as.numeric(DATE(2018,3,15)-DATE(2017,12,29))/365)
r

###---------- slide 10 Implied Volatility of Options------
# Call option on SPY
S <- 266.86
X <- 275
t <- 0.2109589
r <- 0.01226235
q <- 0.0183                   # div yield of S&P 500
sigma <- 0.07853767
Call <- BSC(S,X,t,r,q,sigma)
round(Call,2)


# slide 13
IVC <- function(S,X,t,r,q,Call) {
  if(!require(nloptr)) install.packages('nloptr')
  library(nloptr)
  eval_f_C <- function(sigma) {
    return ( (Call-BSC(S,X,t,r,q,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_C,
                opts=opts)
  return(es$solution)
}
IVC(266.86,275,t,r,q,1.04) 


###-------------------------slide.14 + VaR & ES of option-------
# SPY and vfinx should have the same GARCH model
setwd("D:/MQM530/R_code/")
load("vfinx_2018.rda")

# run data up to end_date
vfinx <- vfinx["2000-01-03/2017-12-31"]
rdat <- as.vector(vfinx$logret)

CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}
### Fit an AR(1)--GARCH(1,1) model with Student t innovations ################
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch.t <- ugarchfit(spec = uspec.t, data = rdat)
# simulate 1-day ahead from AR(1)-GARCH(1,1) ~ t
set.seed(123789)
nsim <- 100000
nper <- 1
alpha <- 0.95
boot.garch <- ugarchboot(garch.t,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp",
                         cluster=NULL)
sim <- boot.garch@fseries
vol <- boot.garch@fsigma
vol0<- garch.t@fit$sigma[length(rdat)]

####--------- slide.15, method1 - fixed implied vol, 1-day ahead----delta/gamma-----
# method 1: delta/gamma, fixed implied vol, 1-day ahead
alpha <- 0.95
Delta <- BSC.delta(S,X,t,r,q,sigma)
Gamma <- BSC.gamma(S,X,t,r,q,sigma)
delS  <- S * (exp(sim)-1)
Call.chg1 <- Delta * delS + 0.5 * Gamma * delS^2
CalcVaRES(Call.chg1,alpha) 

####--------- slide.16, method2 - fixed implied vol, 1-day ahead-----BS----
# method 2: BS, fixed implied vol

t1 <- t - 1/252
Sn <- S * exp(sim)
Call.chg2 <- BSC(Sn,X,t1,r,q,sigma)-Call 
CalcVaRES(Call.chg2,alpha)

####--------- slide.17, method3 - BS, changed % impl vol -----BS----

# method 3: BS, impl vol changed same % as vfinx vol
t1 <- t - 1/252
sigma1 <- sigma * vol/vol0
Sn <- S * exp(sim)
Call.chg3 <- BSC(Sn,X,t1,r,q,sigma1)-Call 
CalcVaRES(Call.chg3,alpha)

####--------- slide.18, method4 - scenario analysis-------------

# method 4: scenario analysis
Schg <- c(-0.1,0.1)
sigchg<-c(-0.5,0.5)
Callchg<- matrix(0,nrow=2,ncol=2)
for (i in 1:2) {
  for (j in 1:2) {
    Snew <- S * exp(Schg[i])
    signew <- sigma * exp(sigchg[j])
    Cnew <- BSC(Snew,X,t,r,q,signew)
    Callchg[i,j] <- Cnew-Call
  }
}
Callchg


####----------------selected strategies -----------------------
# Slide 19
S <- 266.86
t <- 0.2109589
r <- 0.01226235
q <- 0.0183

# Call option on SPY
X.c <- 275
Call <- 1.04
IV.c <- IVC(S,X.c,t,r,q,Call)

# Put option on SPY
X.p <- 240
Put <- 0.97
IV.p <- IVP(S,X.p,t,r,q,Put)

##-------------1. Unhedged portfolio-------------
V0 <- 10000000
Q.spy <- 9000000/S
Q.cash<- 1000000

# 10-day ahead 95% confidence level
# VaR/ES for unhedged portfolio

vfinx <- vfinx["2000-01-03/2017-12-31"]
rdat <- as.vector(vfinx$logret)
### Fit an AR(1)--GARCH(1,1) model with Student t innovations ################
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch.t <- ugarchfit(spec = uspec.t, data = rdat)
vol0    <- garch.t@fit$sigma[length(rdat)]

set.seed(123789)
nper <- 10                 #10-day VaR
nsim <- 100000
boot.garch <- ugarchboot(garch.t,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp",
                         cluster=NULL)
simmat <- boot.garch@fseries
sim <- apply(simmat,1,sum)       # this sums the simulated 1-day log returns to get 10-day log returns
vol <- boot.garch@fsigma[,nper]  # this vector contains the last day's volatility of each simulation
alpha <- 0.95
VaR_ES <- CalcVaRES(sim,alpha)
VaR1 <- Q.spy * S * (exp(VaR_ES[1])-1)
ES1  <- Q.spy * S * (exp(VaR_ES[2])-1)
round(VaR1,0)
round(ES1,0)


### loss value 
# scenario SPY down 20%    
loss <- Q.spy * S * (exp(-0.2)-1)
loss

##-------------1. covered call-------------
# 9000,000 in SPY + 1000,000 in cash
# sell spy call option, strike at 275, 2018-3-16 due, @$1.04

###---- VAR/ES for covered call
# assume IV for all changes with SPY's vol
#  changes in t,s,iv
Q.call <- Q.spy/100    # number of call potion
tn   <- t - 10/252     # time 
Sn <- S * exp(sim)     #
IV.cn <- IV.c * vol/vol0
sim2 <- Q.spy * (Sn-S) - Q.call * (BSC(Sn,X.c,tn,r,q,IV.cn)-Call)*100

VaR_ES2 <- CalcVaRES(sim2,alpha)
VaR_ES2[1]
VaR_ES2[2]

# scenario SPY down 20%, implied vol rises 50%
Sn <- S * exp(-0.2)
IVn<- IV.c * 1.5
loss <- Q.spy * S * (exp(-0.2)-1) - Q.call * ( BSC(Sn,X.c,tn,r,q,IVn)-Call )*100
loss


