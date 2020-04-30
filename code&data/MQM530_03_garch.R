# Conditional volatility and GARCH
library(data.table)
library(ggplot2)
library(MASS)
library(metRology)
library(quantmod)
library(xts)
library(rugarch)
library(moments)

CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}


setwd("D:/MQM530/R_code/")

############----------- data  prep -----------------
load("vfinx_2018.rda")
# run data up to end_date
vfinx <- vfinx["2000-01-03/2017-12-31"]
rdat <- as.vector(vfinx$logret)    # list of logret
############--------- volatility clustering -----------------

# acf of daily log returns (from Conditional Mean & ARMA)
acf(rdat,main="ACF of daily log return")
# slide 5
acf(abs(rdat),main="ACF of |daily log return|")
##### ---------- we can see strong evidence of volatility clustering

############---------randomize & acf-----------------

# slide 6
# randomly permute |log return| and run acf()
set.seed(123437)
rvec <- sample(rdat,length(rdat),replace=FALSE)
acf(rvec,main="ACF of randomly permuted |logret|")

############--------- slide.7+ ------------

#model table
model.table <- data.frame(model=c("IID~n","IID~t","AR(1)-GARCH(1,1)~n","AR(1)-GARCH(1,1)~t"),
                          loglik=rep(0,4),
                          aic  = rep(0,4) )
model.table$loglik[1:2] <- c(13554.2,14117.73)			# from Lecture 1 slide 42, 43
model.table$aic[1:2] <- 2*c(2,3) - 2*model.table$loglik[1:2]	# from Lecture 1 slide 45


##########-----------AR(1)-GARCH(1,1) (normality)------------

# slide 8
### Fit an AR(1)--GARCH(1,1) model with normal innovations ###################
uspec.N <- ugarchspec(variance.model = list(model = "sGARCH", 
                                            garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), 		
                                        include.mean = TRUE), 		
                      distribution.model = "norm")

# slide 9 
### fit the GARCH in rdat   
### using Maximum Likelihood
rdat <- as.vector(vfinx$logret)			
garch.N <- ugarchfit(spec = uspec.N, data = rdat)
round(garch.N@fit$coef,6)

# slide 10
### Test for normality of standardized residuals
z <- garch.N@fit$z
round(mean(z),6)
round(sd(z),6)
round(skewness(z),2)
round(kurtosis(z),2)
jarque.test(z)      # TB = 511, not normality -> incorrect model 

# slide 11    
acf(z,main="ACF of standardized residuals of AR(1)-GARCH(1,1)~N ")
acf(abs(z),main="ACF of |standardized residuals| of AR(1)-GARCH(1,1)~N ")
## shows that GRACH remove the volatility clustering 

##########-----------GARCH(1,1)~t (t-distribution)------------

# slide 13
### Fit an AR(1)--GARCH(1,1) model with Student t innovations ################
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch.t <- ugarchfit(spec = uspec.t, data = rdat)
round(garch.t@fit$coef,6)

# slide 14
z <- garch.t@fit$z    # z is the standardized residuals
round(mean(z),6)
round(sd(z),6)
round(skewness(z),2)
round(kurtosis(z),2)
3+6/(garch.t@fit$coef[6]-4)
jarque.test(z)

# slide 15
acf(z,main="ACF of standardized residuals of AR(1)-GARCH(1,1)~t ")
acf(abs(z),main="ACF of |standardized residuals| of AR(1)-GARCH(1,1)~t ")

###############################
####------ Likelihood-ratio test for comparing the two models-------------------

### The likelihood ratio statistic is twice the difference between the log-likelihoods
LL.N <- garch.N@fit$LLH # log-likelihood of the model based on normal innovations
LL.t <- garch.t@fit$LLH # log-likelihood of the model based on t innovations
LL.N
LL.t
LRT <- -2*(LL.N-LL.t) # likelihood-ratio test statistic
LRT > qchisq(0.95, 1) # => H0 is rejected
1-pchisq(LRT, df = 1) # p-value (probability of such an extreme result if normal hypothesis were true)

model.table$loglik[3:4] <- c(LL.N,LL.t)
model.table$aic[3:4]    <- c(2*5-2*LL.N,2*6-2*LL.t)

# slide 16
model.table   ### AR-GARCH-t is the best
#############################

####------mdel fit & forecasting -------------------

# slide 17
### plot the fitted volatility
df <- data.frame(Date=index(vfinx),s=garch.t@fit$sigma*sqrt(260))

ggplot(df,aes(x=Date,y=s)) +
  geom_line(lwd=1) +
  ggtitle("Fitted Standard Deviations") +
  xlab("") +
  ylab("annualized volatility")

df[which.max(df$s),]
df[which.min(df$s),]


####------1-day ahead simulate-------------------

# slide 19
# simulate 1-day ahead from AR(1)-GARCH(1,1) ~ t
set.seed(123789)
nper <- 1
nsim <- 100000
boot.garch <- ugarchboot(garch.t,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp")
sim <- boot.garch@fseries             ## simulated outcomes
sim
#########--------------slide 20  1-day ahead VaR based on simulated outcmes---------
alpha <- 0.95
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

# slide 21
# comparison with VIX
vix <- getSymbols(Symbols = "^VIX", src="yahoo",return.class="xts",
                  from=df$Date[1],to="2017-12-31",
                  auto.assign=FALSE)[,6]
names(vix) <- "vix"
vix <- as.data.frame(vix)
dte <- row.names(vix)
vix$Date <- as.Date(dte,format="%Y-%m-%d")
vix$vix  <- vix$vix/100

df2 <- merge(df,vix,by="Date",all.x=TRUE)

# correlation between log(vix) and log(s)
round(cor(log(df2$s),log(df2$vix)),2)

cols <- c("s"="blue","vix"="red")
ggplot(df2,aes(x=Date)) +
  geom_line(aes(x=Date,y=log(s),colour="s"),lwd=1) +
  geom_line(aes(x=Date,y=log(vix),colour="vix"),lwd=1) +
  ggtitle("Fitted Standard Deviations vs VIX (log scale)") +
  ylab("log(volatility)") +
  xlab("") +
  scale_colour_manual("Volatility", values = cols,
                      breaks=c("vix","s"))
# slide 22
ggplot(df2,aes(x=log(vix),y=log(s))) +
  geom_point() + 
  xlim(c(-3,0)) +
  ylim(c(-3,0)) +
  xlab("log(vix)") +
  ylab("log(s)") +
  geom_abline(lwd=1,color="red")



##########----------- slide 23 Estimate from xxx to certain date -----
# fit until Oct 15, 2008
z <- vfinx["2000-01-03/2008-10-15"]            # find data from XXX -> xxx
rdat <- as.vector(z$logret)
garch.t.2008 <- ugarchfit(spec = uspec.t, data = rdat)    # fit garch model
# simulate 1-day ahead
set.seed(123789)
nsim <- 100000
nper <- 1
alpha <- 0.95
boot.garch <- ugarchboot(garch.t.2008,
                         method=c("Partial","Full")[1],
                         sampling="raw",
                         n.ahead=nper,
                         n.bootpred=nsim,
                         solver="solnp",
                         cluster=NULL)
sim <- boot.garch@fseries
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)

####------------------------slide 25  Rolling Garch-------------
# rolling garch in 2018
load("vfinx_2018.rda")
vfinx <- vfinx["2000-01-03/2018-12-31"]
rdat <- as.vector(vfinx$logret)
n2017 <- length(vfinx["2000-01-03/2017-12-31"])
system.time(
  roll.garch <- ugarchroll(spec=uspec.t,
                           data=rdat,
                           n.ahead=1,
                           n.start= n2017,
                           refit.every=1,
                           refit.window="recursive",
                           calculate.VaR=TRUE,
                           VaR.alpha=0.05,
                           cluster=NULL,
                           keep.coef=TRUE)
)  # user 2682.84   system 57.61  elapsed 360.68
#    with d
str(roll.garch@forecast$VaR)
dfvar <- roll.garch@forecast$VaR
dfvar <- dfvar[1:251,]                             # keep only 2018
names(dfvar) <- c("VaR","actual")
dte <- index(vfinx["2018-01-01/"])
dfvar$date <- dte
head(dfvar)

# slide 26
cols <- c("var"="red","actual"="blue")
ggplot(dfvar) +
  geom_line(aes(x=date,y=VaR,colour="var"),lwd=1) +
  geom_col(aes(x=date,y=actual,colour="actual")) +
  scale_colour_manual(NULL,values=cols) +
  theme(legend.position="none") +
  ggtitle("1-day 95% VaR during 2018") +
  xlab("") +
  ylab("") +
  ylim(c(-0.05,0.05)) 

# how many times was var violated?
mean(dfvar$actual<lag(dfvar$VaR))

#######------------------slide 27  10-day ahead from garch-------------------
# simulate 10-day ahead from AR(1)-GARCH(1,1) ~ t
set.seed(123789)
nper <- 10                       # nper = 10
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
round(VaR_ES,6)


