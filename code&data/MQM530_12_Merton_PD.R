# merton's model of structure default
library(ggplot2)
library(reshape2)

# slide 5

V0 <- 1000
mu <- 0.15
sigma <- 0.35
mean <- mu - 0.5 * sigma^2
meanln <- log(V0) + mu - 0.5 * sigma^2

F <- 600

# slide 6

# simulate 100 paths, with time intervals of 0.01
npath <- 100
Vmat <- matrix(rnorm(npath*101),ncol=npath,nrow=101)
st  <- sigma*sqrt(0.01)
mt  <- mean*0.01
Vmat[1,] <- V0
for (j in 2:101) {
  Vmat[j,] <- Vmat[j-1,]*exp(Vmat[j,]*st+mt)
}

graph.df <- melt(Vmat)
graph.df$Var1 <- graph.df$Var1/100
ggplot(graph.df,aes(x=Var1,y=value,group=Var2)) + 
  geom_line() +
  xlab("time") +
  geom_hline(yintercept=600,color="red",lwd=1)

# slide 7
# graph log normal distribution

x<- seq(0,5000,length = 1000)
y <- dlnorm(x, meanlog = meanln, sdlog = sigma, log = FALSE)
graph.df <- data.frame(x=x, y=y)
ggplot(graph.df) +
  geom_line(aes(x=x,y=y),lwd=1,color="blue") +
  geom_vline(xintercept = 600,lwd=1,color="red") +
  xlab("V1") +
  ylab("frequency") +
  scale_x_continuous(breaks=c(0,600,1000,2000,3000,4000,5000)) +
  scale_y_continuous(breaks=seq(0,0.001,0.00025))

# slide 10

d2s <- ( log(V0) - log(F) + mean) / sigma
d1s <- d2s + sigma

PD  <- 1-pnorm(d2s)
RR  <- (V0/F) * exp(mu) * (1-pnorm(d1s)) / (1-pnorm(d2s))
LGD <- 1 - RR

PD
RR
LGD


#-------------------------------------------
#verify by simulation

nsim <- 100000
V1   <- rep(V0,nsim)
rsim <- rnorm(nsim,0,1)
V1   <- V1 * exp(rsim*sigma+mean)
graph_df <- data.frame(Vt=V1)
ggplot(graph_df) +
  geom_density(aes(x=Vt),lwd=1,color="blue") +
  geom_vline(xintercept=600,color="red",lwd=1) +
  scale_x_continuous(breaks=c(0,600,1000,2000,3000,4000,5000),limits=c(0,5000)) +
  scale_y_continuous(breaks=seq(0,0.001,.00025))

PD.sim   <- mean(V1<F)
RR.sim   <- mean(V1[V1<F]) / F

PD.sim
RR.sim
1-RR.sim


#-------------------------------------------
# S is the price of the underlying asset
# X is the strike price of the option
# t is the time to maturity (in years)
# r is the tbill rate (in decimal form)
# q is the dividend yield of the underlying asset, paid continuously
# sigma is the volatility of the underlying asset
# BSC is the price of the call option
# BSP is the price of the put  option
# BSC.delta is the delta pf the call option
# BSP.delta is the delta of the put  option
# BSC.gamma is the gamma of the call option
# BSP.gamma is the gamma of the put  option

BSC <- function(S,X,t,r,q=0,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  C <- S * exp(-q*t) * N1 - X * exp(-r*t) * N2
  C
}
BSP <- function(S,X,t,r,q=0,sigma) {
  d1 <- log(S/X) + (r-q+0.5*sigma^2)*t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  NM1 <- pnorm(-d1)
  NM2 <- pnorm(-d2)
  P <- X * exp(-r*t) * NM2 - S * exp(-q*t) * NM1
  P
}

# slide 15
q  <- 0
rf <- 0.02
V0 <- 1000
F  <- 600
S0 <- BSC(V0,F,1,rf,q,sigma)
B0 <- V0 - S0
V0
S0
B0

# slide 16
d1 <- ( log(V0/F)+(rf+0.5*sigma^2) ) / sigma
sigS <- sigma * pnorm(d1) * (V0 / S0)
sigS

# slide 17
# solve for V0 and sigV
library(rootSolve)
S0 <- 420
sigS <- 0.75
model <- function(x) {
  F1 <- S0 - BSC(x[1],F,1,rf,q,x[2])
  d1 <- ( log(x[1]/F)+(rf+0.5*x[2]^2) )/x[2]
  F2 <- sigS - x[2]*pnorm(d1)*x[1]/S0 
  c(F1 = F1, F2 = F2)
}
(ss <- multiroot(f=model,start=c(800,0.75)))
# x[1] is V0, x[2] is sigV

