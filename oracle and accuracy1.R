library(mvtnorm);library(GrassmannOptim);library(MASS);library(Renvlp);library(matrixStats);library(psych);library(Hmisc)
library(pls);library(spls);library(matrixcalc);library(caret);library(simrel);library(glasso);library(future.apply)
#================================================================================================================
rm(list = ls())
set.seed(10)
runs <- 50 #zhu used 200 replications
r <- 3
q <- 2
p <- 12
pa <- 4
pi <- p - pa
pa.q <- pa - q
n <- 10000 # zhu used 50 to 1000
delt <- 0.95

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

gamma.pa <- gramschmidt(matrix(rnorm(pa*pa, 0, 1), pa, pa))$Q

gamma.a1 <- rbind(gamma.pa[,1:q], matrix(0, pi, q))
gamma.a2 <- rbind(gamma.pa[,(q+1):pa], matrix(0, pi, q))
gamma.i <- rbind(matrix(0, pa, pi), diag(1, pi))
gamma0 <- cbind(gamma.a2, gamma.i)

gamma.both <- cbind(gamma.a1, gamma.a2, gamma.i)

omega1 <- diag(4, q)
omega.01 <- diag(0.8, (pa.q))
omega.02 <- diag(1, pi)
omega0 <- as.matrix(bdiag(omega.01, omega.02))

sigx <- gamma.a1%*%omega1%*%t(gamma.a1) + gamma0%*%omega0%*%t(gamma0)

sigy.x <- matrix(0, r,r)
for(i in 1:r){
  for(j in 1:r){
    sigy.x[i,j] <- delt^(abs(i - j))
  }
}

Beta <- rbind(matrix(c(-2,  4,  3,
                       1,  2,  4,
                       2,  1.5,  5,
                       4, -2, -1), pa, r, byrow = T), matrix(0, pi, r))

X <- list(); E <- list(); Y <- list()
for(k in 1:runs){
  X[k] <- list(mvrnorm(n, rep(0, p), sigx))
  E[k] <- list(mvrnorm(n, rep(0, r), sigy.x))
  Y[k] <- list(scale(X[[k]]%*%Beta + E[[k]]))
  X[k] <- list(scale(X[[k]]))
}

zhu.beta <- lapply(1:runs, function(i) spxenv(X[[i]], Y[[i]], u = q, lambda1 = 1, lambda2=0.02)$beta)
# sink("zhu.beta_100.txt")
# print(zhu.beta)
# sink()

#___________________vectorize beta_________________________________________
new.zhu.beta <- NULL
for(i in 1:runs){
  new.zhu.beta[[i]] <- as.vector(round(zhu.beta[[i]], 7))
}

my.table <- function(actual, predicted){
  actual <- as.vector(actual);  predicted <- as.vector(predicted)
  table.new <- table(ifelse(actual==0 & predicted==0, "TN",
                            ifelse(actual !=0 & predicted !=0, "TP",
                                   ifelse(actual==0 & predicted !=0, "FP", "FN"))))
  return(table.new)
}

#___ table ______________________
tabs <- NULL
for(i in 1:runs){
  tabs[[i]] <- list(my.table(as.vector(Beta), new.zhu.beta[[i]]))[[1]]
}

#___ Accuracy rate ______________
acur <- rep(0, runs)
tpr <- rep(0, runs)
tnr <- rep(0, runs)
for(i in 1:runs){
  acur[i] <- (ifelse(!is.na(tabs[[i]]["TP"]), tabs[[i]]["TP"], 0) + ifelse(!is.na(tabs[[i]]["TN"]), tabs[[i]]["TN"], 0)) /
    (ifelse(!is.na(tabs[[i]]["TP"]), tabs[[i]]["TP"], 0) + ifelse(!is.na(tabs[[i]]["TN"]), tabs[[i]]["TN"], 0) + 
       ifelse(!is.na(tabs[[i]]["FP"]), tabs[[i]]["FP"], 0) + ifelse(!is.na(tabs[[i]]["FN"]), tabs[[i]]["FN"], 0))
  
  tnr[i]  <- ifelse(!is.na(tabs[[i]]["TN"]), tabs[[i]]["TN"], 0)/ (ifelse(!is.na(tabs[[i]]["TN"]), tabs[[i]]["TN"], 0) +
                                                                     ifelse(!is.na(tabs[[i]]["FP"]), tabs[[i]]["FP"], 0))
  tpr[i]  <- ifelse(!is.na(tabs[[i]]["TP"]), tabs[[i]]["TP"], 0)/ (ifelse(!is.na(tabs[[i]]["TP"]), tabs[[i]]["TP"], 0) +
                                                                     ifelse(!is.na(tabs[[i]]["FN"]), tabs[[i]]["FN"], 0))
}
acur
mean(acur)
tpr
mean(tpr)
tnr
mean(tnr)


#________________ bias of estimates ______________________
zhu.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  zhu.oracle.beta[i,] <- as.vector(round(zhu.beta[[i]][1:p.a,], 7))
}

zhu.bias <- colMeans(zhu.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.zhu <- norm_vec(as.vector(zhu.bias))^2
bias.zhu


#________________ standard deviation of estimates ______________________
zhu.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  zhu.oracle.beta[i,] <- as.vector(round(zhu.beta[[i]][1:p.a,], 7))
}
#zhu.oracle.beta 

zhu.total.variance <- sum(eigen(cov(zhu.oracle.beta))$values)
zhu.total.variance
zhu.std <- sqrt(zhu.total.variance)
zhu.std






#________________________________________________________________________
#_____________ SPLS _____________________________________________________
spls.beta <- NULL
for(i in 1:runs){
  spls.beta[i] <- list( spls(X[[i]], Y[[i]], q, eta = 0.8, kappa=0.5, select="simpls", fit="simpls", 
                             scale.x=F, scale.y=F, eps=1e-4, maxstep=5000, trace=F)$betamat )
}

new.spls.beta1 <- NULL
for(i in 1:runs){
  new.spls.beta1[[i]] <- round(spls.beta[[i]][[q]], 7)
}
# sink("new.spls.beta1_100.txt")
# print(new.spls.beta1)
# sink()


#___________________vectorize beta_________________________________________
new.spls.beta <- NULL
for(i in 1:runs){
  new.spls.beta[[i]] <- as.vector(round(spls.beta[[i]][[q]], 7))
}

my.table <- function(actual, predicted){
  actual <- as.vector(actual);  predicted <- as.vector(predicted)
  table.new <- table(ifelse(actual==0 & predicted==0, "TN",
                            ifelse(actual !=0 & predicted !=0, "TP",
                                   ifelse(actual==0 & predicted !=0, "FP", "FN"))))
  return(table.new)
}

#___ table ______________________
tabs1 <- NULL
for(i in 1:runs){
  tabs1[[i]] <- list(my.table(as.vector(Beta), new.spls.beta[[i]]))[[1]]
}

#___ Accuracy rate ______________
acur1 <- rep(0, runs)
tpr1 <- rep(0, runs)
tnr1 <- rep(0, runs)
for(i in 1:runs){
  acur1[i] <- (ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0) + ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0)) /
    (ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0) + ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0) + 
       ifelse(!is.na(tabs1[[i]]["FP"]), tabs1[[i]]["FP"], 0) + ifelse(!is.na(tabs1[[i]]["FN"]), tabs1[[i]]["FN"], 0))
  
  tnr1[i]  <- ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0)/ (ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0) +
                                                                        ifelse(!is.na(tabs1[[i]]["FP"]), tabs1[[i]]["FP"], 0))
  tpr1[i]  <- ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0)/ (ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0) +
                                                                        ifelse(!is.na(tabs1[[i]]["FN"]), tabs1[[i]]["FN"], 0))
}
acur1
mean(acur1)
tpr1
mean(tpr1)
tnr1
mean(tnr1)


#________________ bias of estimates ______________________
spls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  spls.oracle.beta[i,] <- as.vector(round(spls.beta[[i]][[q]][1:p.a,], 7))
}

spls.bias <- colMeans(spls.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.spls <- norm_vec(as.vector(spls.bias))^2
bias.spls



#________________ standard deviation of estimates ______________________
spls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  spls.oracle.beta[i,] <- as.vector(round(spls.beta[[i]][[q]][1:p.a,], 7))
}
#spls.oracle.beta 

spls.total.variance <- sum(eigen(cov(spls.oracle.beta))$values)
spls.total.variance
spls.std <- sqrt(spls.total.variance)
spls.std






#_________________________________________________________________________________
#_________ SIMPLS ________________________________________________________________
mysimpls <- lapply(1:runs, function(i) plsr(Y[[i]] ~ X[[i]], ncomp = q, method = "simpls")$coefficients)

mysimpls.beta1 <- NULL
for(i in 1:runs){
  mysimpls.beta1[[i]] <- round(mysimpls[[i]][,,q], 7)
}
# sink("mysimpls.beta1_100.txt")
# print(mysimpls.beta1)
# sink()



#________________ bias of estimates ______________________
simpls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  simpls.oracle.beta[i,] <- as.vector(round(mysimpls[[i]][,,q][1:p.a,], 7))
}

simpls.bias <- colMeans(simpls.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.simpls <- norm_vec(as.vector(simpls.bias))^2
bias.simpls


#________________ standard deviation of estimates ______________________
simpls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  simpls.oracle.beta[i,] <- as.vector(round(mysimpls[[i]][,,q][1:p.a,], 7))
}
#simpls.oracle.beta 

simpls.total.variance <- sum(eigen(cov(simpls.oracle.beta))$values)
simpls.total.variance
simpls.std <- sqrt(simpls.total.variance)
simpls.std






#_________________________________________________________________________________
#______________________Envelope model_____________________________________________
cook <- lapply(1:runs, function(i) xenv(X[[i]], Y[[i]], u = q, asy = F)$beta)

# sink("cook_100.txt")
# print(cook)
# sink()


#________________ bias of estimates ______________________
cook.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  cook.oracle.beta[i,] <- as.vector(round(cook[[i]][1:p.a,], 7))
}

cook.bias <- colMeans(cook.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.cook <- norm_vec(as.vector(cook.bias))^2
bias.cook


#________________ standard deviation of estimates ______________________
cook.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  cook.oracle.beta[i,] <- as.vector(round(cook[[i]][1:p.a,], 7))
}
#cook.oracle.beta 

cook.total.variance <- sum(eigen(cov(cook.oracle.beta))$values)
cook.total.variance
cook.std <- sqrt(cook.total.variance)
cook.std








#_______________________________________________________________________________
#_____________________OLS_______________________________________________
mmlr <- function(Y, X){
  mvmod <- lm(Y ~ X)
  return(coef(mvmod)[-1,])
}
lm.beta <- lapply(1:runs, function(i) mmlr(Y[[i]], X[[i]]))

# sink("lm.beta_100.txt")
# print(lm.beta)
# sink()

#________________ bias of estimates ______________________
lm.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  lm.oracle.beta[i,] <- as.vector(round(lm.beta[[i]][1:p.a,], 7))
}

lm.bias <- colMeans(lm.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.lm <- norm_vec(as.vector(lm.bias))^2
bias.lm

#________________ standard deviation of estimates ______________________
lm.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  lm.oracle.beta[i,] <- as.vector(round(lm.beta[[i]][1:p.a,], 7))
}
#lm.oracle.beta 

lm.total.variance <- sum(eigen(cov(lm.oracle.beta))$values)
lm.total.variance
lm.std <- sqrt(lm.total.variance)
lm.std









#_____________ Fixed Gamma PLS_____________________
myenvCv <- function(X,Y){
  Y <- as.matrix(Y);  X <- as.matrix(X);  a <- dim(Y);  n <- a[1];  r <- a[2];  p <- ncol(X);  M <- list() 
  
  efit <- matrix(M,p)
  for(j in 1:p){
    efit[j] <- list(xenv(X, Y, u = j, asy = F))
  }
  
  return(efit)
}
em <- lapply(1:runs, function(i) myenvCv(X[[i]], Y[[i]]))

#_________________________________________________________________________________________________
Gamma.fixed1 <- function(par, X, Y, G, low, high, intr){
  n <- nrow(X)
  r <- ncol(Y)
  p <- ncol(X)
  obj <- function(par, X, Y, G, pen){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    
    X=scale(X, center=T, scale=F)
    Y=scale(Y, center=T, scale=F)
    
    G <- as.matrix(G)
    n <- nrow(X)
    r <- ncol(Y)
    p <- ncol(X)
    q <- ncol(G)
    
    B <- matrix(par, p, r)
    
    pen <- pen
    Z <- X%*%G
    sigZY <- cov(Z, Y)*(n - 1)/n
    sigZ <- cov(Z)*(n - 1)/n
    sigY <- cov(Y)*(n - 1)/n 
    sigY.Z <- sigY - t(sigZY)%*%solve(sigZ)%*%sigZY
    
    g.norm1 <- function(x){
      r.norm <- matrix(0,p,r)
      for(i in 1:p){
        for(j in 1:r){   r.norm[i,j] <- abs(x[i,j])    }
      }
      return(sum(r.norm))
    }
    
    g.normq <- function(x){
      r.norm1 <- rep(0,q)
      for(i in 1:q){  r.norm1[i] <- sqrt(sum(x[i,]^2))   }
      return(r.norm1)
    }
    
    g.normp <- function(x){
      r.norm1 <- rep(0,p)
      for(i in 1:p){  r.norm1[i] <- sqrt(sum(x[i,]^2))   }
      return(r.norm1)
    }
    
    J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*sum(g.normp(B))
    
    return(J)
  }
  
  gridd <- seq(low, high, intr) 
  fit4 <- lapply(gridd, function(k) nlm(f=obj, p=par, X=X, Y=Y, G=G, pen=k, iterlim = 2000, gradtol = 1e-16))
  
  mini.fit <- rep(0, length(gridd))
  for(k in 1:length(gridd)){
    mini.fit[k] <- fit4[[k]]$minimum
  }
  
  para1 <- NULL
  for(k in 1:length(gridd)){
    para1[k] <- list(round(matrix(fit4[[k]]$estimate,p,r), 7))
  }
  
  bic.fit <- rep(0, length(gridd))
  for(k in 1:length(gridd)){
    bic.fit[k] <- mini.fit[k] + length(as.vector(para1[[k]])[ which(!as.vector(para1[[k]]) == 0)])*log(n)
    }
  
  bic.fit.min <- which.min(bic.fit)
  opt.Beta <- para1[[bic.fit.min]]
  opt.lambda <- gridd[bic.fit.min]
  
  #____________________________________________________________________________________________________
  
  return(list(para1 = para1, mini.fit = mini.fit, bic.fit = bic.fit, bic.fit.min = bic.fit.min,
              opt.Beta = opt.Beta, opt.lambda = opt.lambda))
}
#plan(multiprocess)
gf <- lapply(1:runs, function(i) Gamma.fixed1(em[[i]][[q]]$beta, X[[i]], Y[[i]], em[[i]][[q]]$Gamma,
                                                     low=0, high=0.008, intr=0.001))

gf.beta1 <- NULL
for(i in 1:runs){
  gf.beta1[[i]] <- round(gf[[i]]$opt.Beta, 5)
}
# sink("gf.beta1_100.txt")
# print(gf.beta1)
# sink()

#___________________vectorize beta_________________________________________
gf.beta <- NULL
for(i in 1:runs){
  gf.beta[[i]] <- as.vector(round(gf[[i]]$opt.Beta, 5))
}

my.table <- function(actual, predicted){
  actual <- as.vector(actual);  predicted <- as.vector(predicted)
  table.new <- table(ifelse(actual==0 & predicted==0, "TN",
                            ifelse(actual !=0 & predicted !=0, "TP",
                                   ifelse(actual==0 & predicted !=0, "FP", "FN"))))
  return(table.new)
}

#___ table ______________________
tabs3 <- NULL
for(i in 1:runs){
  tabs3[[i]] <- list(my.table(as.vector(Beta), gf.beta[[i]]))[[1]]
}

#___ Accuracy rate ______________
acur3 <- rep(0, runs)
tpr3 <- rep(0, runs)
tnr3 <- rep(0, runs)
for(i in 1:runs){
  acur3[i] <- (ifelse(!is.na(tabs3[[i]]["TP"]), tabs3[[i]]["TP"], 0) + ifelse(!is.na(tabs3[[i]]["TN"]), tabs3[[i]]["TN"], 0)) /
    (ifelse(!is.na(tabs3[[i]]["TP"]), tabs3[[i]]["TP"], 0) + ifelse(!is.na(tabs3[[i]]["TN"]), tabs3[[i]]["TN"], 0) + 
       ifelse(!is.na(tabs3[[i]]["FP"]), tabs3[[i]]["FP"], 0) + ifelse(!is.na(tabs3[[i]]["FN"]), tabs3[[i]]["FN"], 0))
  
  tnr3[i]  <- ifelse(!is.na(tabs3[[i]]["TN"]), tabs3[[i]]["TN"], 0)/ (ifelse(!is.na(tabs3[[i]]["TN"]), tabs3[[i]]["TN"], 0) +
                                                                        ifelse(!is.na(tabs3[[i]]["FP"]), tabs3[[i]]["FP"], 0))
  tpr3[i]  <- ifelse(!is.na(tabs3[[i]]["TP"]), tabs3[[i]]["TP"], 0)/ (ifelse(!is.na(tabs3[[i]]["TP"]), tabs3[[i]]["TP"], 0) +
                                                                        ifelse(!is.na(tabs3[[i]]["FN"]), tabs3[[i]]["FN"], 0))
}
acur3
mean(acur3)
tpr3
mean(tpr3)
tnr3
mean(tnr3)


#________________ bias of estimates ______________________
gf.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  gf.oracle.beta[i,] <- as.vector(round(gf[[i]]$opt.Beta[1:p.a,], 5))
}

gf.bias <- colMeans(gf.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.gf <- norm_vec(as.vector(gf.bias))^2
bias.gf


#_____________ standard deviation of estimates ______________________
gf.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  gf.oracle.beta[i,] <- as.vector(round(gf[[i]]$opt.Beta[1:p.a,], 5))
}
#gf.oracle.beta 

gf.total.variance <- sum(eigen(cov(gf.oracle.beta))$values)
gf.total.variance
gf.std <- sqrt(gf.total.variance)
gf.std






#__________________________________________________________________________________
#___________ modified E-SPLS ______________________________________________________
new_espls <- function(X, Y, low.g, high.g, intr, maxiter, d, start.value){
  
  X <- as.matrix(X); Y <- as.matrix(Y);  n <- nrow(X); p <- ncol(X); r <- ncol(Y)
  
  sigxx <- cov(X)*(n - 1)/n; sigxxyy <- cov(X,Y)*(n - 1)/n; sigyy <- cov(Y)*(n - 1)/n
  Sxy <- sigxx - sigxxyy%*%solve(sigyy)%*%t(sigxxyy)
  
  Sigs <- list(Sxy = Sxy, sigxx = sigxx, sigxxyy = sigxxyy)
  
  g.norm <- function(x,b){
    r.norm <- matrix(0,p,d)
    for(i in 1:p){
      for(j in 1:d){
        r.norm[i,j] <- abs(x[i,j])
      }
    }
    return(b*sum(r.norm))
  }
  
  g.norm2 <- function(x,b){
    r.norm1 <- rep(0,p)
    for(i in 1:p){  r.norm1[i] <- sqrt(sum(x[i,]^2))   }
    return(b*sum(r.norm1))
  }
  
  logdet =function(A){
    eigtem=eigen(A,only.values=TRUE,symmetric=TRUE)$values
    y=sum(log(eigtem[eigtem > 1e-10])) #[eigtem > 1e-10]
    return(y)
  }
  
  f <- function(W){
    Qt <- W$Qt
    d <- W$dim[1]
    p <- ncol(Qt)
    lambda <- W$lambda
    U <- matrix(Qt[,1:d], ncol=d)
    V <- matrix(Qt[,(d + 1):p], ncol = (p - d))
    #return( - ( log(det(t(U)%*%Sigs$Sxy%*%U)) + log(det(t(V)%*%Sigs$sigxx%*%V)) + g.norm(U,lambda) ) )
    return( - ( logdet(t(U)%*%Sigs$Sxy%*%U) + logdet(t(V)%*%Sigs$sigxx%*%V) + g.norm(U,lambda)  ) )
    #return( - ( logdet(t(U)%*%Sigs$Sxy%*%U) + logdet(t(V)%*%Sigs$sigxx%*%V) + g.norm(U%*%solve(t(U)%*%sigxx%*%U)%*%t(U)%*%sigxxyy,lambda) ) )
  }
  
  objfun <- function(W) list(value = f(W))
  #low.g <- 100; high.g <- 200; intr <- 100; maxiter <- 400 ; start.value <- eigen(sigxx[[1]])$vectors
  W <- NULL; out <- NULL; grid <- seq(low.g, high.g, intr)
  for(i in grid){
    Qt <- start.value
    W[[i]] <- list(Qt=Qt, dim=c(d,p), lambda = i, Sigs = Sigs)
    out[[i]] <- GrassmannOptim(objfun, W[[i]], eps_conv=1e-5, verbose=F, max_iter = maxiter)
  }
  
  gam <- NULL; gam.full <- NULL; Beta.est <- NULL
  for(i in grid){
    gam[[i]] <- out[[i]]$Qt[,1:d] 
    gam.full[[i]] <- out[[i]]$Qt
    Beta.est[[i]] <- round(out[[i]]$Qt[,1:d]%*%solve(t(out[[i]]$Qt[,1:d])%*%sigxx%*%out[[i]]$Qt[,1:d])%*%t(out[[i]]$Qt[,1:d])%*%sigxxyy, 6)
  }
  
  com.gam <- plyr::compact(gam); r.com.gam <- lapply(com.gam,round,6)
  gam.full <- plyr::compact(gam.full); r.gam.full <- lapply(gam.full,round,6)
  Beta.est <- plyr::compact(Beta.est); r.Beta.est <- lapply(Beta.est,round,6)
  
  MSE <- rep(0, length(grid))
  for(i in 1:length(grid)){
    MSE[i] <- norm(Y - X%*%r.Beta.est[[i]], type = "F") + length(as.vector(r.Beta.est[[i]])[ which(!as.vector(r.Beta.est[[i]]) == 0)])*log(n)
  }
  min.MSE <- which.min(MSE)
  opt.beta <- r.Beta.est[[min.MSE]]
  opt.lambda <- grid[min.MSE] 
  
  #------------------------------------------------------------------
  return(list( r.gam.full = r.gam.full, r.com.gam = r.com.gam, r.Beta.est = r.Beta.est,
               opt.lambda = opt.lambda, MSE = MSE, min.MSE = min.MSE, opt.beta = opt.beta))
}
plan(multiprocess)
new.meslps <- lapply(1:2, function(i) new_espls(X[[i]], Y[[i]], 1, 1000, intr = 10,
                                                  maxiter = 2000, d = q, qr.Q(qr(cov(X[[i]])))))
#new.meslps

#___________________vectorize beta_________________________________________
new.mespls.beta1 <- NULL
for(i in 1:runs){
  new.mespls.beta1[[i]] <- round(new.meslps[[i]]$opt.beta, 5)
}
#new.mespls.beta1

# sink("new.mespls.beta1_10000_2.txt")
# print(new.mespls.beta1)
# sink()


new.mespls.beta <- NULL
for(i in 1:runs){
  new.mespls.beta[[i]] <- as.vector(round(new.meslps[[i]]$opt.beta, 4))
}

my.table <- function(actual, predicted){
  actual <- as.vector(actual);  predicted <- as.vector(predicted)
  table.new <- table(ifelse(actual==0 & predicted==0, "TN",
                            ifelse(actual !=0 & predicted !=0, "TP",
                                   ifelse(actual==0 & predicted !=0, "FP", "FN"))))
  return(table.new)
}

#___ table ______________________
tabs2 <- NULL
for(i in 1:runs){
  tabs2[[i]] <- list(my.table(as.vector(Beta), new.mespls.beta[[i]]))[[1]]
}

#___ Accuracy rate ______________
acur2 <- rep(0, runs)
tpr2 <- rep(0, runs)
tnr2 <- rep(0, runs)
for(i in 1:runs){
  acur2[i] <- (ifelse(!is.na(tabs2[[i]]["TP"]), tabs2[[i]]["TP"], 0) + ifelse(!is.na(tabs2[[i]]["TN"]), tabs2[[i]]["TN"], 0)) /
    (ifelse(!is.na(tabs2[[i]]["TP"]), tabs2[[i]]["TP"], 0) + ifelse(!is.na(tabs2[[i]]["TN"]), tabs2[[i]]["TN"], 0) + 
       ifelse(!is.na(tabs2[[i]]["FP"]), tabs2[[i]]["FP"], 0) + ifelse(!is.na(tabs2[[i]]["FN"]), tabs2[[i]]["FN"], 0))
  
  tnr2[i]  <- ifelse(!is.na(tabs2[[i]]["TN"]), tabs2[[i]]["TN"], 0)/ (ifelse(!is.na(tabs2[[i]]["TN"]), tabs2[[i]]["TN"], 0) +
                                                                        ifelse(!is.na(tabs2[[i]]["FP"]), tabs2[[i]]["FP"], 0))
  tpr2[i]  <- ifelse(!is.na(tabs2[[i]]["TP"]), tabs2[[i]]["TP"], 0)/ (ifelse(!is.na(tabs2[[i]]["TP"]), tabs2[[i]]["TP"], 0) +
                                                                        ifelse(!is.na(tabs2[[i]]["FN"]), tabs2[[i]]["FN"], 0))
}
acur2
mean(acur2)
tpr2
mean(tpr2)
tnr2
mean(tnr2)


#________________ bias of estimates ______________________
mespls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  mespls.oracle.beta[i,] <- as.vector(round(new.meslps[[i]]$opt.beta[1:p.a,], 7))
}

mespls.bias <- colMeans(mespls.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.mespls <- norm_vec(as.vector(mespls.bias))^2
bias.mespls


#_____________ standard deviation of estimates ______________________
mespls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  mespls.oracle.beta[i,] <- as.vector(round(new.meslps[[i]]$opt.beta[1:p.a,], 7))
}
#mespls.oracle.beta 

mespls.total.variance <- sum(eigen(cov(mespls.oracle.beta))$values)
mespls.total.variance
mespls.std <- sqrt(mespls.total.variance)
mespls.std
