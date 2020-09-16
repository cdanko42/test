rm(list = ls())
library(MASS)
library(ks)
library(lmtest)
library(ivpack)

set.seed(8675309)
simulateiv <- function(n=500, size=1000, rhoxz, rhoxe, eevs= 1, exo =1, instrument =1 ){
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
  
  c1 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c2 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c3 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c4 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c5 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c6 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c7 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  c8 <- matrix(0, nrow=n, ncol= nrow(rhoxe))
  coverage <- matrix(0, nrow=nrow(rhoxe), ncol= 8)
  
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
      xo <- xez[, (2+eevs):(1+ eevs+exo)]
      olsdata = as.data.frame(cbind(ypre, x, xo))
      olspyre <- lm(ypre ~., data = olsdata)
      r1[i, j] <- olspyre$coefficients[2]
      cols <- coeftest(olspyre)[2, 2]
      cover <- function(estimate, se){
        upper <- estimate + 1.96*se
        lower <- estimate - 1.96*se
        if (.5 > lower & .5 < upper){
          return(1)}
        else{
          return(0)}
      }
      
      c1[i, j] <- cover(estimate= r1[i,j], se = cols)
      
      ivpre <- ivreg(ypre~x+ xo, ~z + xo)
      r2[i,j] <- ivpre$coefficients[2]
      invisible(ivse <- robust.se(ivpre)[2,2])
      c2[i, j] <- cover(estimate = r2[i,j], se=ivse)
      
      yvaldata = as.data.frame(cbind(y_values, x, xo))
      olsyval <- lm(y_values ~., data=yvaldata)
      r3[i, j] <- olsyval$coefficients[2]
      invisible(cols3 <- coeftest(olsyval)[2, 2])
      c3[i, j] <- cover(estimate = r3[i,j], se=cols3)
      
      dat = as.data.frame(cbind(y_values, x,z,xo))
      probyval <- glm(y_values ~., family = binomial(link = "probit"), data = yvaldata)
      r4[i, j] <- probyval$coefficients[2]
      invisible(seprobit <- coeftest(probyval)[2,2])
      c4[i, j] <- cover(estimate = r4[i,j], se=seprobit)
      
      ivyval <- ivreg(y_values~x+ xo, ~z + xo)
      r5[i, j] <- ivyval$coefficients[2]
      invisible(iv2se <- robust.se(ivyval)[2,2])
      c5[i,j] <- cover(estimate = r5[i,j], se=iv2se)
      
      ##Control function
      ## Run first stage regression
      probitcf1 <- lm(x~z+xo)
      ## Collect residuals
      v = probitcf1$residuals
      cfdata = as.data.frame(cbind(y_values, x, xo, v))
      probitcf2 <- glm(y_values ~. , family = binomial(link = "probit"), data = cfdata, control=list(epsilon = 1e-8, maxit = 100))
      r6[i, j] <- probitcf2$coefficients[2]
      
      ##Bootstrap standard errors
      probitboots <- function(bootsize=300){
        bootse <- c()
        for (p in 1:bootsize){
          rows = sample(1:n, 1000, replace = TRUE)
          bootdat <- dat[rows, ]
          booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
          bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
          bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
          bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
          bootcf1 <- lm(bootx~bootz+bootxo)
          bootv = bootcf1$residuals
          bootcf2 <- glm(booty ~ bootx + bootxo + bootv , family = binomial(link = "probit"), control=list(epsilon = 1e-8, maxit = 100))
          bootse[p] <- bootcf2$coefficients[2]
        }
        return(sd(bootse))
      }
      
      cprobit <- probitboots()
      c6[i,j] <- cover(estimate = r6[i,j], se=cprobit)
      
      
      
      ## Special regressor
      ##This is our step 1
      
      
      xo <- as.data.frame(xo)
      ##demean our special regressor, the exogenous variable "xo". This is our V
      dat$demeanxo <- xo[,1]-mean(xo[,1])
      demeanxo <- dat$demeanxo
      
      ##Obtain residuals, which are our U_i
      specialregdata = as.data.frame(cbind(demeanxo, x, z))
      if (exo > 1){
        specialregdata = cbind(specialregdata, xo[, 2:exo])
      }
      
      fssr <- lm(demeanxo ~., data=specialregdata)
      dat$u <- fssr$residuals
      
      ##Step 2
      fhat <- kde(dat$u,eval.points=dat$u)$estimate
      
      
      ##Step 3
      ## Create vector that will collect our I(Vi > 0)
      dat$idv <- rep(0,nrow(dat))
      ## Collect identity function of v
      dat$idv <- replace(dat$idv,which(dat$demeanxo>=0),1)
      
      ##Obtain Ti; y_values=Di, fi and idv is defined above. 
      dat$t <- (dat$y_values- dat$idv)/fhat
      
      ## Trim extreme values
      lower = .025*n + 1
      upper = .975*n
      
      sortdat <- dat[order(dat$t),]
      trimdat <- sortdat[lower:upper,]
      firststagespec <- trimdat[, (2 + eevs):(1+ eevs + instrument + exo)]
      firststagespec <- firststagespec[, -(instrument +1)]
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
      
      cfboots <- function(bootsize=300){
        bootse <- c()
        for (p in 1:bootsize){
          rows = sample(1:n, 1000, replace = TRUE)
          bootdat <- as.data.frame(dat[rows, ])
          booty <- as.matrix(bootdat[, 1], ncol=1, nrow=n)
          bootx <- as.matrix(bootdat[, 2], ncol=1, nrow=n)
          bootxo <- as.matrix(bootdat[, (3+instrument):(2+exo+instrument)], ncol=exo, nrow=n)
          bootz <- as.matrix(bootdat[, 3:(2+instrument)], ncol=instrument, nrow=n)
          bootdat$demeanbootxo <- (bootxo[,1]-mean(bootxo[,1]))
          demeanbootxo <- bootdat$demeanbootxo
          
          ##Obtain residuals, which are our U_i
          specialregdata = as.data.frame(cbind(demeanbootxo, bootx, bootz))
          if (exo > 1){
            specialregdata = cbind(specialregdata, bootxo[, 2:exo])
          }
          
          fssr <- lm(demeanbootxo ~., data=specialregdata)
          bootdat$u <- fssr$residuals
          
          ##Step 2
          fhat <- kde(bootdat$u,eval.points=bootdat$u)$estimate
          
          
          ##Step 3
          ## Create vector that will collect our I(Vi > 0)
          bootdat$idv <- rep(0,nrow(bootdat))
          ## Collect identity function of v
          bootdat$idv <- replace(bootdat$idv,which(bootdat$demeanbootxo>=0),1)
          
          ##Obtain Ti; y_values=Di, fi and idv is defined above. 
          bootdat$t <- (bootdat[,1]- bootdat$idv)/fhat
          
          ## Trim extreme values
          lower = .025*n + 1
          upper = .975*n
          
          sortdat <- bootdat[order(bootdat$t),]
          trimdat <- sortdat[lower:upper,]
          firststagespec <- trimdat[, (2 + eevs):(1+ eevs + instrument + exo)]
          firststagespec <- firststagespec[, -(instrument +1)]
          firststagespec <- as.data.frame(cbind(firststagespec, trimdat[2]))
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
          bootse[p] <- specialreg2$coefficients[2]
          
        }
        return(sd(bootse))
      }
      
      cols7 <-  cfboots()
      c7[i,j] = cover(estimate= r7[i,j], se= cols7)
      
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
      
      coverage[j, 1] <- sum(c1[, j])
      coverage[j, 2] <- sum(c2[, j])  
      coverage[j, 3] <- sum(c3[, j])  
      coverage[j, 4] <- sum(c4[, j])
      coverage[j, 5] <- sum(c5[, j])  
      coverage[j, 6] <- sum(c6[, j])
      coverage[j, 7] <- sum(c7[, j])  
      coverage[j, 8] <- sum(c8[, j])    
    }
  }
  return(list(results =results, coverage=coverage ))  
}

# mad1 = simulateiv(rhoxz = 0.01, rhoxe = 0.5)
# mad1

##Function use is the same as before, except rhox
mad1 = simulateiv(rhoxz = 0.1, rhoxe = c(.1,.2,0.3,.4,.5))
sink("NUL")

mad2 = simulateiv(rhoxz = 0.3, rhoxe = c(.1,.2,0.3,.4,.5))
mad3 = simulateiv(rhoxz = 0.5, rhoxe =c(.1,.2,0.3,.4,.5))
mad4 = simulateiv(rhoxz = 0.7, rhoxe = c(.1,.2,0.3,.4,.5))
mad5 = simulateiv(rhoxz = 0.9, rhoxe = c(.1,.2,0.3,.4,.5))
sink()

mad1$results
mad1$coverage

mad2$results
mad2$coverage

mad3$results
mad3$coverage

mad4$results
mad4$coverage

mad5$results
mad5$coverage

##two exogenous variables
#mad6 = simulateiv(rhoxz = 0.9, rhoxe = c(.1,.2,0.3,.4,.5), exo = 2)
#mad6$results
#mad6$coverage
## two instruments, feed a vector into rhoxz, with length equal to the number of instruments
#mad7 = simulateiv(rhoxz = c(0.5, .5), rhoxe =c(.1,.2,0.3,.4,.5), instrument = 2)
#mad7$results
#mad7$coverage

## two exogenous, two instrument
#mad8 = simulateiv(rhoxz = c(0.5, .5), rhoxe = c(.1,.2,0.3,.4,.5), instrument = 2, exo=2)
#mad8
