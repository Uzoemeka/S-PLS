library(mvtnorm);library(GrassmannOptim);library(MASS);library(Renvlp);library(matrixStats);
library(psych);library(Hmisc);library(pls);library(spls);library(matrixcalc);library(caret);
library(simrel);library(glasso);library(future.apply)
#================================================================================================================
rm(list = ls())
set.seed(10)
runs <- 2 #zhu used 200 replications
r <- 3
q <- 3
p <- 12
pa <- 4
pi <- p - pa
pa.q <- pa - q
n <- 100 # zhu used 50 to 1000
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
gamma.a2 <- rbind(matrix(gamma.pa[,(q+1):pa], pa, pa.q), matrix(0, pi, pa.q))
gamma.i <- rbind(matrix(0, pa, pi), diag(1, pi))
gamma0 <- cbind(gamma.a2, gamma.i)

gamma.both <- cbind(gamma.a1, gamma.a2, gamma.i)

omega1 <- diag(0.8, q)
omega.01 <- diag(30, (pa.q))
omega.02 <- diag(100, pi)
omega0 <- as.matrix(bdiag(omega.01, omega.02))



sigx <- gamma.a1%*%omega1%*%t(gamma.a1) + gamma0%*%omega0%*%t(gamma0)

sigy.x <- matrix(0, r,r)
for(i in 1:r){
  for(j in 1:r){
    sigy.x[i,j] <- delt^(abs(i - j))
  }
}

Beta <- rbind(matrix(c(-2,  4,  3,
                       1,  2.5,  4,
                       2,  1.5,  5,
                       4, -2, -1), pa, r, byrow = T), matrix(0, pi, r))

X <- list(); E <- list(); Y <- list()
for(k in 1:runs){
  X[k] <- list(mvrnorm(n, rep(0, p), sigx))
  E[k] <- list(6*mvrnorm(n, rep(0, r), sigy.x))
  Y[k] <- list(scale(X[[k]]%*%Beta + E[[k]]))
  X[k] <- list(scale(X[[k]]))
}


n_folds <- 5
folds_i <- groups <- sample(rep(1:n_folds, length.out = n))
Zhun_cv <- function(X, Y, n_folds, n_comp){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  
  # folds_i <- sample(rep(1:n_folds, length.out = n))
  
  L<-list()
  prac.Zhun <- matrix(L,n_folds,n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      prac.Zhun[k,j]<-list(spxenv(X[folds_i !=k,], Y[folds_i !=k,], u=j, lambda1=1, lambda2=0.5)$beta)
    }
  }
  
  cv_Zhun <- matrix(NA, nrow = n_folds, ncol = n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      cv_Zhun[k, j] <- norm(as.matrix(Y[folds_i == k,] - X[folds_i == k,]%*%as.matrix(prac.Zhun[k, j][[1]])), type="F")
    }
  }
  
  return(list(prac.Zhun = prac.Zhun, cv_Zhun = cv_Zhun, cmcv_Zhun = colMeans(cv_Zhun)))
}
new_comp <- 11
Zhun.rep <- lapply(1:runs, function(i) Zhun_cv(X[[i]], Y[[i]], n_folds=5, new_comp))
Zhun_mat <- matrix(0, runs, new_comp)
for(i in 1:runs){
  Zhun_mat[i,] <- Zhun.rep[[i]]$cmcv_Zhun
}

# sink("Zhun_matr3q3p12n100del0.95e6ome0.8.30.100.txt")
# print(Zhun_mat)
# print(colMeans(Zhun_mat))
# print(colSds(Zhun_mat))
# sink()

x <- c(1:(new_comp))
plot(colMeans(Zhun_mat), main=expression(list(n==100, r==3, q==3, p==12, p.a==4, delta==0.95)), 
     ylim = c(0,5.1),  type = "l", ylab = "Cross-validated MSE", xlab = "Number of components", lty=1)
errbar(x, colMeans(Zhun_mat), colMeans(Zhun_mat)+colSds(Zhun_mat), colMeans(Zhun_mat)-colSds(Zhun_mat),
       cap=0.004, type = "n", add=TRUE, col="black", lwd=0.5, errbar.col="black")




#__________ SPLS _________________________
spls_cv <- function(X, Y, n_folds, n_comp){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  
  #n_folds <- 2;   n_comp <- 3
 # folds_i <- sample(rep(1:n_folds, length.out = n))
  
  prac.spls <- list()
  for(k in 1:n_folds){
    prac.spls[k] <- list(spls(X[folds_i !=k,], Y[folds_i !=k,], n_comp, eta=0.55, kappa=0.5, select="simpls", fit="simpls",
                              scale.x=F, scale.y=F, eps=1e-4, maxstep=3000, trace=F)$betamat)
  }
  
  cv_spls <- matrix(NA, nrow = n_folds, ncol = n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      cv_spls[k,j] <- norm(as.matrix(Y[folds_i == k,]) - as.matrix(X[folds_i == k,])%*%as.matrix(prac.spls[[k]][[j]]), type="F")
    }
  }
  
  return(list(prac.spls = prac.spls, cv_spls = cv_spls, cmcv_spls = colMeans(cv_spls)))
}
keles.rep <- lapply(1:runs, function(i) spls_cv(X[[i]], Y[[i]], 5, new_comp))
keles_mat <- matrix(0, runs, new_comp)
for(i in 1:runs){
  keles_mat[i,] <- keles.rep[[i]]$cmcv_spls
}

# sink("keles_matr3q3p12n100del0.95e6ome0.8.30.100.txt")
# print(keles_mat)
# print(colMeans(keles_mat))
# print(colSds(keles_mat))
# sink()

lines(x*1.01, colMeans(keles_mat), type="l", col="blue", lty = 2)
errbar(x*1.01, colMeans(keles_mat), colMeans(keles_mat)+colSds(keles_mat), colMeans(keles_mat)-colSds(keles_mat),
       type="n", cap=0.004, add=TRUE, col="blue", pch=19, lwd=0.5, errbar.col="blue")




#_________ SIMPLS ________________________
mysimpls <- function(X, Y, n_folds, n_comp){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  
 # folds_i <- sample(rep(1:n_folds, length.out = n))
  
  prac.mysimpls <- list()
  for(k in 1:n_folds){
    prac.mysimpls[k] <- list(plsr(Y[folds_i !=k,] ~ X[folds_i !=k,], ncomp = n_comp, method = "simpls")$coefficients)
  }
  
  cv_simpls <- matrix(NA, nrow = n_folds, ncol = n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      cv_simpls[k,j] <- norm(as.matrix(Y[folds_i == k,]) - as.matrix(X[folds_i == k,])%*%as.matrix(prac.mysimpls[[k]][,,j]), type="F")
    }
  }
  
  return(list(prac.mysimpls = prac.mysimpls, cv_simpls = cv_simpls, cmcv_simpls = colMeans(cv_simpls)))
}
#cv.simpls = mysimpls(X[[1]], Y[[1]], 5, new_comp)
simpls.rep <- lapply(1:runs, function(i) mysimpls(X[[i]], Y[[i]], 5, new_comp))
simpls_mat <- matrix(0, runs, new_comp)
for(i in 1:runs){
  simpls_mat[i,] <- simpls.rep[[i]]$cmcv_simpls
}

# sink("simpls_matr3q3p12n100del0.95e6ome0.8.30.100.txt")
# print(simpls_mat)
# print(colMeans(simpls_mat))
# print(colSds(simpls_mat))
# sink()

lines(x*1.02, colMeans(simpls_mat), type="l", col="purple", lty = 3)
errbar(x*1.02, colMeans(simpls_mat), colMeans(simpls_mat)+colSds(simpls_mat), colMeans(simpls_mat)-colSds(simpls_mat), type="n", cap=0.004, add=TRUE, col="purple", pch=19, lwd=0.5, errbar.col="blue")




#______________________Using the Envelope model_____________________________________
envCv <- function(X, Y, n_folds, n_comp){
  Y <- as.matrix(Y);  X <- as.matrix(X)
  a <- dim(Y);  n <- a[1]
  r <- a[2];  p <- ncol(X)
  
  #n_folds <- 5; n_comp <- 6 
 # folds_i <- sample(rep(seq_len(n_folds), length.out = n))
  
  L <- list()
  prac.efit <- matrix(L,n_folds,n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      prac.efit[k,j] <- list(xenv(X[folds_i!=k,], Y[folds_i!=k,], u = j, asy = F))
    }
  }
  
  cv_efit <- matrix(NA, nrow = n_folds, ncol = n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      cv_efit[k,j] <- norm(as.matrix(Y[folds_i==k,] - as.matrix(X[folds_i==k,])%*%as.matrix(prac.efit[k, j][[1]]$beta)), "F")
    }
  }
  return(list(prac.efit = prac.efit, cv_efit = cv_efit, cmcv_efit = colMeans(cv_efit)))
}
#cv.envCv <- envCv(X[[1]], Y[[1]], n_folds=5, n_comp=new_comp)
env.rep <- lapply(1:runs, function(i) envCv(X[[i]], Y[[i]], 5, new_comp))
env_mat <- matrix(0, runs, new_comp)
for(i in 1:runs){
  env_mat[i,] <- env.rep[[i]]$cmcv_efit
}

# sink("env_matr3q3p12n100del0.95e6ome0.8.30.100.txt")
# print(env_mat)
# print(colMeans(env_mat))
# print(colSds(env_mat))
# sink()

lines(x*1.01, colMeans(env_mat), type="l", col="forest green", lty = 4)
errbar(x*1.01, colMeans(env_mat), colMeans(env_mat)+colSds(env_mat), colMeans(env_mat)-colSds(env_mat), type="n", cap=0.004, add=TRUE, col="forest green", pch=19, lwd=0.5, errbar.col="forest green")




#___________ New E-SPLS ____________________________________________________________________
plan(multiprocess)

#_________________________________________________________________________________________________
Gamma.fixed1 <- function(par, X, Y, G, benv, low, high, intr){
  
  n <- nrow(X)
  r <- ncol(Y)
  p <- ncol(X)
  
  obj <- function(par, X, Y, G, pen, benv){
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
    
    g.normp <- function(x, w){
      r.norm1 <- rep(0,p)
      for(i in 1:p){  r.norm1[i] <- sqrt(sum(x[i,]^2))*(1/sqrt(sum(w[i,]^2)))^0.5   }
      return(r.norm1)
    }
    
    g.norm2 <- function(x){
      r.norm1 <- rep(0,p)
      for(i in 1:p){  r.norm1[i] <- sqrt(sum(x[i,]^2))   }
      return(r.norm1)
    }
    
    J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*sum(g.normp(B, benv))# + pen*sum(g.norm1(B))
    
    return(J)
  }
  
  gridd <- seq(low, high, intr) 
  fit4 <- lapply(gridd, function(k) nlm(f=obj, p=par, X=X, Y=Y, G=G, benv=benv,  pen=k, iterlim = 2000, gradtol = 1e-16))
  
  mini.fit <- rep(0, length(gridd))
  for(k in 1:length(gridd)){
    mini.fit[k] <- fit4[[k]]$minimum
  }
  
  para1 <- NULL
  for(k in 1:length(gridd)){
    para1[k] <- list(round(matrix(fit4[[k]]$estimate,p,r), 5))
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

gf_cv <- function(X, Y, n_folds, n_comp, low, high, intr){
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  n <- nrow(Y)
  
  #n_folds <- 5;   n_comp <- 6
 # folds_i <- sample(rep(1:n_folds, length.out = n))
  
  L <- list()
  prac.efit <- matrix(L,n_folds,n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      prac.efit[k,j] <- list(xenv(X[folds_i!=k,], Y[folds_i!=k,], u = j, asy = F))
    }
  }
  
  M=list()
  prac.gf <- matrix(M,n_folds,n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      prac.gf[k, j] <- list(Gamma.fixed1(prac.efit[k,j][[1]]$beta, X[folds_i !=k,], Y[folds_i !=k,],
                                         prac.efit[k,j][[1]]$Gamma, prac.efit[k,j][[1]]$beta, low, high, intr))
    }
  }
  
  cv_gf <- matrix(NA, nrow = n_folds, ncol = n_comp)
  for(k in 1:n_folds){
    for(j in 1:n_comp){
      cv_gf[k,j] <- norm(as.matrix(Y[folds_i == k,]) - as.matrix(X[folds_i == k,])%*%as.matrix(prac.gf[k,j][[1]]$opt.Beta), type="F")
    }
  }
  
  return(list(prac.gf = prac.gf, cv_gf = cv_gf, cmcv_gf = colMeans(cv_gf)))
}

gf.rep <- future_lapply(1:runs,function(i) gf_cv(X[[i]], Y[[i]], 5, new_comp,
                                               low=0, high=0.00002, intr=0.00001))

gf_mat <- matrix(0, runs, new_comp)
for(i in 1:runs){
  gf_mat[i,] <- gf.rep[[i]]$cmcv_gf
}

# sink("gf_matr3q3p12n100del0.95e6ome0.8.30.100.txt")
# print(gf_mat)
# print(colMeans(gf_mat))
# print(colSds(gf_mat))
# sink()

lines(x*1.02, colMeans(gf_mat), type="l", col="red", lty = 5)
errbar(x*1.02, colMeans(gf_mat), colMeans(gf_mat)+colSds(gf_mat), colMeans(gf_mat)-colSds(gf_mat), type="n",
       cap=0.004, add=TRUE, col="red", pch=19, lwd=0.5, errbar.col="red")

legend("bottomright", legend=c("E-SPLS","SPLS","SIMPLS","EPLS", "2S-SPLS"),col=c("black", "blue", "purple","forest green", "red"), lty=1:5,cex=0.8)
