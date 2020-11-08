#devtools::install_github('selbouhaddani/PPLS/Package/PPLS')
library(mvtnorm); library(PPLS); library(pracma); library(Renvlp); library(pls);library(MASS);
library(PPLS);library(matrixStats);library(fBasics); library(Hmisc); library(ggplot2);
library(reshape2); library(simrel); library(future.apply)
#devtools::install_github("selbouhaddani/PO2PLS@RCpp")
library(OmicsPLS); library(PO2PLS)
#__________________________________________________________________________________________

# peach <- read.csv("peach_spectra+brixvalues.csv")
# dim(peach)
# 
# X <- scale(peach[,2:30])
# Y <- scale(peach[,1])
# p <- dim(X)[2]
# fit <- plsr(Y ~ ., data = data.frame(X,Y), ncomp = p, method = "simpls", validation = "LOO")
# summary(fit)
# plot(RMSEP(fit), legendpos = "topright")
#
# eigen(cor(X))$values


#___________ data simulation __________________________________
#___________________________________________________________________________________
rm(list = ls())
set.seed(100)
runs <- 100 #zhu used 200 replications
n <- 30 # zhu used 50 to 1000
p <- 5
r <- 3
q <- 2
delt <- 0.8
rho <- 0.01

#IRn30p5r3q2delta0_8rho0_8

gramschmidt <- function(x) {
  x <- as.matrix(x)
  # Get the number of rows and columns of the matrix
  n <- ncol(x)
  m <- nrow(x)
  
  # Initialize the Q and R matrices
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)
  
  for (j in 1:n) {
    v = x[,j] # Step 1 of the Gram-Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:(j-1)) {
        r[i,j] <- t(q[,i]) %*% x[,j] # Find the inner product (noted to be q^T a earlier)
        # Subtract the projection from v which causes v to become perpendicular to all columns of Q
        v <- v - r[i,j] * q[,i]
      }
    }
    # Find the L2 norm of the jth diagonal of R
    r[j,j] <- sqrt(sum(v^2))
    # The orthogonalized result is found and stored in the ith column of Q.
    q[,j] <- v / r[j,j]
  }
  
  # Collect the Q and R matrices into a list and return
  qrcomp <- list('Q'=q, 'R'=r)
  return(qrcomp)
}

gamma <- gramschmidt(matrix(rnorm(p*p, 0, 1), p, p))$Q

gamma1 <- gamma[,1:q]
gamma0 <- gamma[,(q+1):p]

sig <- matrix(0, p,p)
for(i in 1:p){
  for(j in 1:p){
    sig[i,j] <- rho^(abs(i - j))
  }
}

# omega <- eigen(sig)$values
# sigx <- gamma1%*%diag(omega[1:q], q) %*%t(gamma1) + gamma0%*%diag(omega[(q+1):p], (p-q)) %*%t(gamma0)

omega <- sort(eigen(sig)$values, decreasing = F)
sigx <- gamma1%*%diag(omega[1:q], q) %*%t(gamma1) + gamma0%*%diag(omega[(q+1):p], (p-q)) %*%t(gamma0)

sigy.x <- matrix(0, r,r)
for(i in 1:r){
  for(j in 1:r){
    sigy.x[i,j] <- delt^(abs(i - j))
  }
}

#sigy.x <- diag(c(2,1.5,0.5), r)

Beta <- gamma1%*%matrix(runif(q*r, 0, 2), q, r)

X <- list(); E <- list(); Y <- list()
for(k in 1:runs){
  X[k] <- list(mvrnorm(n, rep(0, p), sigx))
  E[k] <- list(1*mvrnorm(n, rep(0, r), sigy.x))
  Y[k] <- list(scale(X[[k]]%*%Beta + E[[k]]))
  X[k] <- list(scale(X[[k]]))
}
#cor(X[[1]])


plan(multiprocess)
#______________________Using the Envelope model_____________________________________
fold = k = 10
groups <- sample(rep(seq_len(fold), length.out = n))

myenvCv <- function(X,Y,k){
  #Y <- as.matrix(Y[[1]]);  X <- as.matrix(X[[1]])
  Y <- as.matrix(Y);  X <- as.matrix(X)
  a <- dim(Y);  n <- a[1];  r <- a[2];  p <- ncol(X); p <- p-1
  
  fitt <- plsr(Y ~ X, ncomp = p, method = "simpls")$loadings
  
  M <- list()
  efit <- matrix(M,k,p)
  for(i in 1:k){
    for(j in 1:p){
      efit[i,j] <- list(xenv(X[groups!=i,], Y[groups!=i,], u = j, asy = F, as.matrix(fitt[,1:j])))
    }
  }
  
  prederror <- matrix(M,k,p)
  for(i in 1:k){
    for(j in 1:p){
      prederror[i,j] <- list(as.matrix(Y[groups==i,] - as.matrix(X[groups==i,])%*%efit[[i,j]]$beta))
    }
  }
  
  ecv <- matrix(0,k,p)
  for(i in 1:k){
    for(j in 1:p){
      ecv[i,j] <- norm(prederror[i,j][[1]], type = "F")
    }
  }
  
  return(envpls = colMeans(ecv))
  
}
envrepp1 <- matrix(unlist(future_lapply(1:runs, function(i) myenvCv(X=X[[i]], Y=Y[[i]], k=fold))), ncol=p-1, byrow=T)
ecv_mn1 <- colMeans(envrepp1); ecv_sd1 <- colSds(envrepp1)

sink("EIRn30p5r3q2delta0_8rho0_01.txt")
print(envrepp1)
print(colMeans(envrepp1))
print(colSds(envrepp1))
sink()


#_______________________Using the SIMPLS model_______________________________________
mysimplsCv <- function(X,Y,k){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  

  sfit <- NULL
  for(i in 1:k){
    sfit[i] <- list(plsr(Y[groups!=i,] ~ X[groups!=i,], ncomp = p, method = "simpls"))
  }
  
  
  
  cv <- matrix(0,k,p)
  for(i in 1:k){
    for(j in 1:p){
      cv[i,j] <- norm( as.matrix( Y[groups==i,] - X[groups==i,]%*%as.matrix(sfit[[i]]$coefficients[,,j]) ),type="F")
    } 
  }
  
  return(splsCv <- colMeans(cv))
}
simplsrepp1 <- matrix(unlist(future_lapply(1:runs, function(i) mysimplsCv(X=X[[i]], Y=Y[[i]], k=fold))), ncol=p, byrow=T)
scv_mn1 <- colMeans(simplsrepp1); scv_sd1 <- colSds(simplsrepp1)

sink("SIRn30p5r3q2delta0_8rho0_01.txt")
print(simplsrepp1)
print(colMeans(simplsrepp1))
print(colSds(simplsrepp1))
sink()



#_______________________Using the PPLS model_______________________________________
# myPPLSCv <- function(X, Y, k, emnum){
#   Y <- as.matrix(Y);   X <- as.matrix(X)
#   a <- dim(Y)
#   n <- a[1];   r <- a[2]
#   p <- ncol(X)
#   M=list()
#   pfit <- matrix(M,k,r)
#   for(i in 1:k){
#     for(j in 1:r){
#       pfit[i,j] <- list(PPLS_simult(X=as.matrix(X[groups!=i,]),
#                                     Y=as.matrix(Y[groups!=i,]), a = j, EMsteps = emnum, atol = 1e-09, type = "SVD"))
#     }
#   }
#   pp.pred <- matrix(0, k,r)
#   for(i in 1:k){
#     for(j in 1:r){
#       pp.pred[i,j] <- norm(as.matrix(Y[groups==i,] - 
#                                        X[groups==i,]%*%ginv(pfit[i,j][[1]]$estimates$W%*%t(pfit[i,j][[1]]$Expectations$mu_T)%*%
#                                                               pfit[i,j][[1]]$Expectations$mu_T%*%t(pfit[i,j][[1]]$estimates$W))%*%
#                                        pfit[i,j][[1]]$estimates$W%*%t(pfit[i,j][[1]]$Expectations$mu_T)%*%pfit[i,j][[1]]$Expectations$mu_U%*%
#                                        t(pfit[i,j][[1]]$estimates$C)),type="F")
#     }
#   }
#   return <- colMeans(pp.pred)
# }
# PPLSrepp1 <- matrix(unlist(future_lapply(1:runs, function(i) myPPLSCv(X=X[[i]],
#                                                   Y=Y[[i]], k=fold, emnum = 50))), ncol=r, byrow=T)
# pcv_mn1 <- colMeans(PPLSrepp1); pcv_sd1 <- colSds(PPLSrepp1)



myPPLSCv <- function(X, Y, k, emnum){
  Y <- as.matrix(Y);   X <- as.matrix(X)
  a <- dim(Y);  n <- a[1];   r <- a[2];  p <- ncol(X)
  
  # Y <- as.matrix(Y);   X <- as.matrix(X)
  
  M=list();  pfit <- matrix(M,k,r)
  for(i in 1:k){
    for(j in 1:r){
      pfit[i,j] <- list(PO2PLS(X[groups!=i,], Y[groups!=i,], r = j, 0, 0, steps=emnum, tol=0.000001))
    }
  }
  
  
  pp.pred <- matrix(0, k,r)
  for(i in 1:k){
    for(j in 1:r){
      pp.pred[i,j] <- norm(as.matrix(Y[groups==i,] - X[groups==i,] %*% ginv(pfit[i,j][[1]]$params$W
                                                                            %*% pfit[i,j][[1]]$params$SigT %*% t(pfit[i,j][[1]]$params$W))
                                     %*% pfit[i,j][[1]]$params$W %*% pfit[i,j][[1]]$params$SigT %*% 
                                       pfit[i,j][[1]]$params$B %*% t(pfit[i,j][[1]]$params$C)),  type="F")
    }
  }
  
  return <- colMeans(pp.pred)
}
PPLSrepp1 <- matrix(unlist(future_lapply(1:runs, function(i) myPPLSCv(X=X[[i]],
                                                                      Y=Y[[i]], k=fold, emnum = 150))), ncol=r, byrow=T)
pcv_mn1 <- colMeans(PPLSrepp1); pcv_sd1 <- colSds(PPLSrepp1)


sink("PIRn30p5r3q2delta0_8rho0_01.txt")
print(PPLSrepp1)
print(colMeans(PPLSrepp1))
print(colSds(PPLSrepp1))
sink()



#_____________________Using the OLS model________________________________________
mlmCV <- function(X, Y, k){
  Y <- as.matrix(Y); X <- as.matrix(X);  a <- dim(Y);
  n <- a[1];  r <- a[2];  p <- ncol(X)
  #k=10
  mlm1 <- NULL
  for(i in 1:k){
    mlm1[i] <- list(lm(cbind(Y[groups!=i,]) ~ X[groups!=i,]))
  }
  
  Yhat1 <- NULL
  for(i in 1:k){
    Yhat1[i] <- list( X[groups==i,]%*%coef(mlm1[[i]])[-1,] )
  }
  
  mlmnorm <- rep(0, k)
  for(i in 1:k){
    mlmnorm[i]=norm(as.matrix(Y[groups==i,]-Yhat1[[i]]), type = "F")
  }
  return(rep(mean(mlmnorm),p))
}
mlmrepp <- matrix(unlist(future_lapply(1:runs, function(i) mlmCV(X=X[[i]], Y=Y[[i]], k=fold))), ncol=p, byrow=T)
lcv_mn <- colMeans(mlmrepp); lcv_sd <- colSds(mlmrepp)

sink("LIRn30p5r3q2delta0_8rho0_01.txt")
print(mlmrepp)
print(colMeans(mlmrepp))
print(colSds(mlmrepp))
sink()


#________________________ benchmark______________________________________________________
# min_error <- lapply(1:runs, function(i) sum(eigen(cov(Y[[i]]) - cov(Y[[i]], X[[i]])%*%
#                                         solve(cov(X[[i]]))%*%cov(X[[i]],Y[[i]]))$values))
 

#___________________ Plots ______________________________________________________
x <- c(1:p)
plot(x,c(ecv_mn1, lcv_mn[1]), type="o", col="red", ylab="Cross-validated RMSE",
     xlab=expression("Choice of components"), lty=1, pch=20,
     ylim = c(0, 3), main=expression(list(n==30, p==5, r==3, q==2, rho==0.01)))
#errbar(x, ecv_mn1, ecv_mn1 + ecv_sd1, ecv_mn1 - ecv_sd1, add = T, col = "red", cap=0.01, pch = 20, lwd = 0.5, errbar.col = "red")
lines(x,scv_mn1, type="o",pch=20, col="blue", lty=2)
#errbar(x*1.02, scv_mn1, scv_mn1 + scv_sd1, scv_mn1 - scv_sd1, add = T, cap=0.01, col = "blue", pch = 20, lwd = 0.5, errbar.col = "blue")
lines(c(1:r), pcv_mn1, type="o", col="black",pch=20, lty=3)
#errbar(c(1:r), pcv_mn1, pcv_mn1 + pcv_sd1, pcv_mn1 - pcv_sd1, add = T, cap=0.01, col = "black", pch = 20, lwd = 0.5, errbar.col = "black")
lines(x,lcv_mn, type = "o",pch=20, col = "cyan3", lty=4)
#errbar(x*1.04, lcv_mn, lcv_mn + lcv_sd, lcv_mn - lcv_sd, add = T, col = "cyan3", cap=0.01, pch = 20, lwd = 0.5, errbar.col = "cyan3")
# lines(x, rep(mean(unlist(evar1)), p), col = "purple", lty = 5)
abline(v=q, col = "gray", lty = 6)
legend("bottomright", legend=c("EPLS","SIMPLS","EM-PLS", "OLS"), col=c("red", "blue", "black", "cyan3"), lty=1:4,cex=1)
