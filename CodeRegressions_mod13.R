rm(list = ls())
library(MASS)
library(ks)
library(kedd)

set.seed(8675309)
simulateiv <- function(n=1000, size=1000, rhoxz, rhoxe, eevs= 1, exo =1, instrument =1 ){
##rho = correlation of instrumental variable with x
##  Initialize matrices
rhoxz <- matrix(rhoxz, nrow=length(rhoxz), ncol=eevs)
rhoxe <- matrix(rhoxe, nrow= length(rhoxe), ncol =eevs)
r1 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r2 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r3 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r4 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r5 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r6 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r7 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
r8 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
results <- matrix(0, nrow=nrow(rhoxe), ncol= 8)

for (j in 1:(nrow(rhoxe))){

for (i in 1:size){
##Create model vcov matrix
# Let's add in one continuous exogenous regressor - we'll most surely want this later on
# For now, there's no reason to assume multicollinearity between endog var and exog var
sig <- matrix(0, nrow=eevs+exo+instrument+1, ncol = eevs+exo+instrument+1)
for (k in 1:ncol(rhoxe)){
    sig[1, k+1] = rhoxe[j, k]
    sig[k+1, 1] = rhoxe[j, k]
    sig[k+1, (eevs+exo+2):(eevs+exo+instrument+1)] = rhoxz[, k]
    sig[(eevs+exo+2):(eevs+exo+instrument+1), k+1] = rhoxz[, k]
  }
## Add 1 diagonal, initialize vector
  variable = c()
  for (h in 1:(eevs+exo+instrument+1)){
    sig[h,h]=1
    variable[h] = 0
  }
xez <- mvrnorm(n, variable,Matrix::nearPD(sig, TRUE, TRUE)$mat)
##Obtain error, e
e <- xez[, 1]
##Obtain EEV, x
x <- xez[, 2:(1+eevs)]
##Obtain IV z (excluded exogenous regressor)
z <- xez[, (2+eevs+exo):(1+ eevs+ instrument+ exo)]
##Obtain included exogenous regressor
xo <- xez[, (2+eevs):(1+ eevs+exo)]
xo <- as.data.frame(xo)

##Specify and sample from true model
ypre = 1 + 0.5*x + e

for (g in 1:exo){
ypre = ypre - (.3/exo)*xo[,g]
}

y_values = rep(0,n)
y_values=replace(y_values,which(ypre>0),1)

# On second thought, let's "build up" to what was previously r
# We will probably end up taking many of these out, but it's good to double check the code makes sense
# Let's go in this order: OLS on ypre, 2SLS on ypre, OLS on yval, probit on yval, 2SLPM on yval

olsdata = as.data.frame(cbind(ypre, x, xo))
olspyre <- lm(ypre ~., data = olsdata)
r1[i, j] <- olspyre$coefficients[2]

firstdata = as.data.frame(cbind(x, z,xo ))
xHat <- lm(x~., data=firstdata)$fitted.values
secondata = as.data.frame(cbind(ypre, xHat, xo))
twoslsypre <-lm(ypre ~., secondata)
r2[i, j] <- twoslsypre$coefficients[2]

yvaldata = as.data.frame(cbind(y_values, x, xo))
olsyval <- lm(y_values ~., data=yvaldata)
r3[i, j] <- olsyval$coefficients[2]

dat = as.data.frame(cbind(y_values, x,z,xo))
probyval <- glm(y_values ~., family = binomial(link = "probit"), data = yvaldata)
r4[i, j] <- probyval$coefficients[2]

yvaltsdata = as.data.frame(cbind(y_values, xHat, xo))
twoslsyval <- lm(y_values ~., data=yvaltsdata)
r5[i, j] <- twoslsyval$coefficients[2]

##Control function
## Run first stage regression
probitcf1 <- lm(x~., data=firstdata)
## Collect residuals
v = probitcf1$residuals
cfdata = as.data.frame(cbind(y_values, x, xo, v))
probitcf2 <- glm(y_values ~. , family = binomial(link = "probit"), data = cfdata, control=list(epsilon = 1e-8, maxit = 100))
r6[i, j] <- probitcf2$coefficients[2]

## Special regressor
##This is our step 1

##demean our special regressor, the exogenous variable "xo". This is our V
dat$V <- rt(n, 10, 0)
attach(dat)
o <- cor(V,e)
h <- cor(V, z)
i <- cor(V, xo)
j <- cor(V, x)
while (abs(o) >= .001 | abs(i) >= .001 | abs(h) >= .001 | abs(j) >= .001){
  dat$V <- rt(1000,20, 0)
  o <- cor(V,e)
  h <- cor(V, z)
  i <- cor(V, xo)
  j <- cor(V, x)
}


##Obtain residuals, which are our U_i
xo <- xez[, (2+eevs):(1+ eevs+exo)]
fssr <- lm(V ~ xo + z , )
dat$u <- fssr$residuals

##Step 2
fhat <- kde(dat$u,eval.points=dat$u)$estimate


##Step 3
## Create vector that will collect our I(Vi > 0)
dat$idv <- rep(0,nrow(dat))
## Collect identity function of v
dat$idv <- replace(dat$idv,which(V>=0),1)

##Obtain Ti; y_values=Di, fi and idv is defined above. 
dat$t <- (dat$y_values- dat$idv)/fhat

## Trim extreme values
lower = .025*n + 1
upper = .975*n

sortdat <- dat[order(dat$t),]
trimdat <- sortdat[lower:upper,]
firststagespec <- trimdat[, (2 + eevs):(1+ eevs + instrument + exo)]
firststagespec <- as.data.frame(cbind(firststagespec, trimdat$x))
colnames(firststagespec)[instrument + exo ] = "x"

## Step 4
##Conduct 2SLS
specialreg <- lm(x~., data=firststagespec)
specHat <- specialreg$fitted.values
## Estimate beta

secondstagespec <- trimdat[, (2 + eevs):(1+ eevs + instrument + exo)]
secondstagespec$t <- trimdat$t
secondstagespec$specHat  <- specHat
specialreg2 <- lm(t~specHat, data=trimdat)
r7[i, j] <- specialreg2$coefficients[2]

# JM: Look at the bottom of p. 821. Comparability of estimates between the special regressor and probit could be an issue
# We will want to start looking into marginal effects and the AIF 
xo = as.matrix(xo, ncol = exo, nrow= n)
z = as.matrix(z, ncol = instrument, nrow=n)
  ## Maximum Likelihood
  lik=function(theta){
    
    alpha1<-theta[1]
    beta<-theta[2]
    alpha2<-theta[3]
    logsig<-theta[4]
    iht<-theta[5]
    gamma<-as.matrix(theta[6:(5+exo)], ncol=1, nrow= exo)
    pi1<- as.matrix(theta[(6+exo):(5+2*exo)], ncol = 1, nrow=exo)
    pi2<- as.matrix(theta[(6+2*exo):(5+2*exo+instrument)], ncol =1, nrow=instrument)
    
    y = y_values
    
    rho = (exp(2*iht)-1)/(exp(2*iht)+1)
    #rho = tanh(iht)
    m = (alpha1+x*beta+xo%*%gamma+(rho*(x-alpha2-xo%*%pi1-z%*%pi2)/exp(logsig))/(sqrt(1-(rho^2))))
    
    logl = y*log(pnorm(m,mean=0,sd=1))+(1-y)*log(1-pnorm(m,mean=0,sd=1))+log(dnorm(((x-alpha2-xo%*%pi1-z%*%pi2)/exp(logsig)),mean=0,sd=1))-logsig
    
    return(-sum(logl))
  }
  
  # Vector of starting values
  start<-c(0,0,0,1,.5)
  start[6:(5+exo)] <- 0
  start[(6+exo):(5+2*exo)] <- 0
  start[(6+2*exo):(5+2*exo+instrument)] <- 0
  
  # tryCatch silences error messages to ensure that loop doesn't end prematurely 
  tryCatch({
    out<-nlm(lik,start,iterlim=1000)
    eststore<-out$estimate
    estgrad<-out$gradient
    convcode<-out$code}, error=function(e){cat("Error:",conditionMessage(e), "\n")})

r8[i, j] <- out$estimate[2] 
  
results[j, 1] <- mean(abs(r1[, j]-0.5))
results[j, 2] <- mean(abs(r2[, j]-0.5))  
results[j, 3] <- mean(abs(r3[, j]-0.5))  
results[j, 4] <- mean(abs(r4[, j]-0.5))  
results[j, 5] <- mean(abs(r5[, j]-0.5))  
results[j, 6] <- mean(abs(r6[, j]-0.5))  
results[j, 7] <- mean(abs(r7[, j]-0.5))  
results[j, 8] <- mean(abs(r8[, j]-0.5))    
}
}
return(results)  
}

# mad1 = simulateiv(rhoxz = 0.01, rhoxe = 0.5)
# mad1

##Function use is the same as before, except rhox
mad1 = simulateiv(rhoxz = 0.1, rhoxe = c(.1,.2,0.3,.4,.5))
mad1
mad2 = simulateiv(rhoxz = 0.3, rhoxe = c(.1,.2,0.3,.4,.5))
mad2
mad3 = simulateiv(rhoxz = 0.5, rhoxe =c(.1,.2,0.3,.4,.5))
mad3
mad4 = simulateiv(rhoxz = 0.7, rhoxe = c(.1,.2,0.3,.4,.5))
mad4
mad5 = simulateiv(rhoxz = 0.9, rhoxe = c(.1,.2,0.3,.4,.5))
mad5

##two exogenous variables
#mad6 = simulateiv(rhoxz = 0.9, rhoxe = c(.1,.2,0.3,.4,.5), exo = 2)
#mad6
## two instruments, feed a vector into rhoxz, with length equal to the number of instruments
#mad7 = simulateiv(rhoxz = c(0.5, .5), rhoxe =c(.1,.2,0.3,.4,.5), instrument = 2)
#mad7

## two exogenous, two instrument
#mad8 = simulateiv(rhoxz = c(0.5, .5), rhoxe = c(.1,.2,0.3,.4,.5), instrument = 2, exo=2)
#mad8
