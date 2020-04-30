# Lecture 8

# mortgage amortization table
library(data.table)
library(ggplot2)



# loan(EAD)   ->   1-PD: full repay
#                  PD  : default    ->  1-LGD : recovery 
#                                        0
#                                
            
#EAD : exposure of default
#LGD : loss given default
#PD  : probability of default
# LGD = 1- Recovery rate
# Expect Loss = EAD*LGD*PD


P <- 100000   # loan principal
R <- 0.06     # interest rate, APR
T <- 30       # length of loan in years
m <- 12       # number of payments in a year 


mortgage <- function(P,R,T=30,m=12) {
  n <- T * m  # number of periods of the loan
  n1<- n+1
  r <- R/m    # periodic interest rate
  AF<- ( 1 - (1+r)^(-n) ) / r    # annuity factor
  c <- P / AF # payment per period
  df <- data.frame(period = seq(0,n),
                  begin.bal = rep(0,n1),
                  payment   = c(0,rep(c,n)),
                  interest  = rep(0,n1),
                  prin.repay= rep(0,n1),
                  end.bal   = rep(0,n1)
                  )
  df$end.bal[1] <- P
  df$end.bal[2:n1] <- P * (1 - ( (1+r)^df$period[2:n1]-1) / ( (1+r)^n-1) )
  df$interest[2:n1]<- df$end.bal[1:n]*r
  df$prin.repay[2:n1] <- df$payment[2:n1]-df$interest[2:n1]
  df$begin.bal[2:n1] <- df$end.bal[1:n]
  return(df)
}

df <- mortgage(P,R,T,m)
round(head(df),2)
round(tail(df),2)

# slide 8

# graph the interest and principal repayment each month
df2 <- df
df2[1,] <- NA
ggplot(df2) +
  geom_line(aes(x=period,y=payment),lwd=0.8) +
  geom_area(aes(x=period,y=payment,fill="principal")) +
  geom_line(aes(x=period,y=interest),lwd=0.8) +
  geom_area(aes(x=period,y=interest,fill="interest")) +
  xlim(c(0,360)) +
  ylim(c(0,1000))+
  xlab("Months") +
  ylab("") +
  scale_fill_manual(values=c("chartreuse3","cornflowerblue")) +
  theme(legend.title = element_blank()) +
  ggtitle("Monthly payments of a mortgage loan")

df3 <- df
df3$col <- "pos"
ggplot(df3,aes(x=period,y=end.bal)) +
  geom_line(color="black",lwd=0.8) +
  geom_area(fill='cornflowerblue') +
  xlim(c(0,360)) +
  ylim(c(0,P))+
  xlab("Months") +
  ylab("") +
  ggtitle("Monthly UPB of a mortgage loan")

# slide 35
# example of a constant hazard
df <- data.table(t=c(1:1000)/10)
lambda <- 0.02
df[, S := exp(-lambda*t)]
ggplot(df,aes(x=t,y=S)) +
  geom_line() +
  ylim(c(0,1)) +
  xlab("Time to Event") +
  ylab("Porportion Surviving") +
  ggtitle("Survival function for constant hazard rate")
