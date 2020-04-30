# Conditional mean
library(data.table)
library(ggplot2)
library(MASS)
library(metRology)
library(quantmod)
library(xts)

setwd("D:/MQM530/R_code/")

# Conditional Mean and ARIMA

data("AirPassengers")

########-------------- slide 5 / seasonality / plot
# levels graph
ts.plot(AirPassengers, xlab="Year", ylab="Number of Passengers", main="Monthly totals of international airline passengers, 1949-1960")
abline(reg=lm(AirPassengers~time(AirPassengers)))
# variance is growing
# there is a seasonal pattern

# slide 6
# log graph
ts.plot(log(AirPassengers), xlab="Year", ylab="Number of Passengers", main="Monthly totals of international airline passengers, 1949-1960")
abline(reg=lm(log(AirPassengers)~time(AirPassengers)))
# variance of log is not growing
# but series is not stationary (it is growing over time)

# slide 7
# log difference graph
ts.plot(diff(log(AirPassengers)), xlab="Year", ylab="Growth of Passengers", main="Monthly growth rates of international airline passengers, 1949-1960")

#######----------------------------slide10: acf and pacf---------------------------------
# autocorrelation coefficients
# partial acf   (residual)
acf(diff(log(AirPassengers)))
pacf(diff(log(AirPassengers)))   # arma(1,1)

#######----------------------------slide11: ARMA model---------------------------------
(fit <- arima(log(AirPassengers), c(1, 1, 1),seasonal = list(order = c(1, 1, 1), period = 12)))

# slide 12
pred <- predict(fit, n.ahead = 10*12)    # apply ARMA to predict
ts.plot(AirPassengers,2.718^pred$pred, log = "y", lty = c(1,3), main="Monthly totals of international airline passengers")


#######---------------------------try this on VFINX---------------------------------

load("vfinx_2018.rda")

# slide 13
# run data up to end_date
start_date <- "2000-01-03"
end_date   <- "2017-12-31"
vfinx <- vfinx["2000-01-03/2017-12-31"]
# plot log returns
df <- data.frame(Date=index(vfinx),logret=vfinx$logret)
ggplot(df,aes(x=Date,y=logret)) +
  geom_line(lwd=1) +
  ylab("VFINX daily log return") +
  xlab("")
# summary
summary(vfinx$logret)

# slide 14
acf(vfinx,main="ACF of daily log returns")
pacf(vfinx,main="PACF of daily log returns",ylim=c(-0.5,0.5))

# slide 15
#######-------------------------AR(1) / MA(1)---------------------------------

# run AR(1) and MA(1)
vfinx.ar1 <- arima(vfinx,c(1,0,0))   #AR(1)
vfinx.ma1 <- arima(vfinx,c(0,0,1))   #MA(1)
vfinx.ar1
vfinx.ma1


