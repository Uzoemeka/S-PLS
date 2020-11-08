library(mvtnorm);library(GrassmannOptim);library(MASS);library(Renvlp);library(matrixStats);library(psych);library(Hmisc)
library(pls);library(spls);library(matrixcalc);library(caret);library(glmnet);library(simrel);library(glasso);library(future.apply)
#================================================================================================================
rm(list = ls())
set.seed(10)
runs <- 50 #zhu used 200 replications
r <- 3
q <- 3
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
gamma.a2 <- rbind(matrix(gamma.pa[,(q+1):pa], pa, pa.q), matrix(0, pi, pa.q))
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

Beta <- rbind(matrix(c(-3,  5,  1,
                       3,  0,  0,
                       2,  1, -2,
                       4,  0,  3), pa, r, byrow = T), matrix(0, pi, r)); Beta


Beta1 <- rbind(matrix(c(-3,  5,  1,
                        3,  5e-4,  5e-4,
                        2,  1, -2,
                        4,  5e-4,  3), pa, r, byrow = T), matrix(5e-4, pi, r)); Beta1


X <- list(); E <- list(); Y <- list()
sigxx <- list(); sigxxyy <- list(); sigyy <- list(); Sxy <- list(); Sigs <- list()
for(k in 1:runs){
  X[k] <- list(mvrnorm(n, rep(0, p), sigx))
  E[k] <- list(mvrnorm(n, rep(0, r), sigy.x))
  Y[k] <- list(scale(X[[k]]%*%Beta + E[[k]]))
  X[k] <- list(scale(X[[k]]))
  sigxx[k] <- list(cov(X[[k]])*(n - 1)/n)
  sigxxyy[k] <- list(cov(X[[k]],Y[[k]])*(n - 1)/n)
  sigyy[k] <- list(cov(Y[[k]])*(n - 1)/n)
  Sxy[k] <-list(sigxx[[k]] - sigxxyy[[k]] %*%solve(sigyy[[k]] )%*%t(sigxxyy[[k]]))
  Sigs[k] <- list(list(Sxy = Sxy[[k]], sigxx = sigxx[[k]], sigxxyy = sigxxyy[[k]]))
}
#pairs.panels(cbind(X[[1]][,1:4], Y[[1]][,3]), method = "pearson", hist.col = "#00AFBB", density = F, ellipses = F)



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


#_______________________Fixed Gamma__________________________________________________________________________
Gamma.fixed1 <- function(par, X, Y, G, w, low, high, intr){
  n <- nrow(X)
  r <- ncol(Y)
  p <- ncol(X)
  obj <- function(par, X, Y, G, pen, w){
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
    
    g.norma <- function(x1, x2){
      r.norm <- matrix(0,p,r)
      for(i in 1:p){ for(j in 1:r){   r.norm[i,j] <- abs(x1[i,j])*(1/abs(x2[i,j])^0.5)   }  }
      return(sum(r.norm))   }
    
    g.norm1 <- function(x){
      r.norm <- matrix(0,p,r)
      for(i in 1:p){    for(j in 1:r){   r.norm[i,j] <- abs(x[i,j])    }      }
      return(sum(r.norm))    }
    
    g.norm2 <- function(x){
      r.norm1 <- rep(0,p)
      for(i in 1:p){  r.norm1[i] <- sqrt(sum(x[i,]^2))   }
      return(sum(r.norm1))
    }
    
    g.normp <- function(x1, x2){
      r.norm1 <- rep(0,p)
      for(i in 1:p){  r.norm1[i] <- sqrt(sum(x1[i,]^2))*(1/ sqrt(sum(x2[i,]^2)))^2   }
      return(r.norm1)
    }
    
    #J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + 0.5*g.norma(B, w) + pen*g.normp(B, w)
    J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*g.norma(B, w) + 0.0005*g.norm2(B)
    #J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*g.norma(B, w)
    #J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*g.norm1(B)
    #J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)%*%solve(sigY.Z)) + pen*g.norm2(B) + pen*g.norm1(B)
    #J <- tr((1/n)*t(Y - Z%*%t(G)%*%B)%*%(Y - Z%*%t(G)%*%B)) + pen*g.norm2(B)
    #_______________________________________________________________________________________________________________  
    return(J)
  }
  
  gridd <- seq(low, high, intr)
  fit4 <- lapply(gridd, function(k) nlm(f=obj, p=par, X=X, Y=Y, G=G, pen=k, w=w, iterlim = 2000, gradtol = 1e-16))
  
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
  
  return(list(para1 = para1, mini.fit = mini.fit, bic.fit = bic.fit,
              bic.fit.min = bic.fit.min, opt.Beta = opt.Beta, opt.lambda = opt.lambda))
}
gf <- lapply(1:runs, function(i) Gamma.fixed1(em[[i]][[q]]$beta, X[[i]], Y[[i]], em[[i]][[q]]$Gamma,
                                              w=em[[i]][[q]]$beta, low=0, high=0.007, intr=0.001))
gf.beta1 <- NULL
for(i in 1:runs){
  gf.beta1[[i]] <- round(gf[[i]]$opt.Beta, 5)
}
gf.beta1

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
  gf.oracle.beta[i,] <- as.vector(round(gf[[i]]$opt.Beta[1:pa,], 7))
}

gf.bias <- colMeans(gf.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.gf <- norm_vec(as.vector(gf.bias))^2
bias.gf


#_____________ standard deviation of estimates ______________________
gf.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  gf.oracle.beta[i,] <- as.vector(round(gf[[i]]$opt.Beta[1:pa,], 7))
}
#gf.oracle.beta 

gf.total.variance <- sum(eigen(cov(gf.oracle.beta))$values)
gf.total.variance
gf.std <- sqrt(gf.total.variance)
gf.std



#_______________________________________________________________________________________________________
# Returns the cholesky factor, L, of the matrix A such that A=t(L)%*%L
#_______________________________________________________________________________________________________
sechol <- function(A, tol = .Machine$double.eps, silent= TRUE ){
  if (is.complex(A))  {# if A is a complex matrix Return NULL
    warning("complex matrices not permitted at present")
    return(NULL)
  } else if (!is.numeric(A))  {# if A is not a numerix matrix Return NULL
    warning("non-numeric argument to 'sechol'")
    return(NULL)
  }
  if (is.matrix(A)) {# if A is a matrix but not a square matrix Return NULL
    if (nrow(A) != ncol(A)) {# if A is not a square matrix Return NULL
      warning("non-square matrix in 'sechol'")
      return(NULL)
    }
    if (nrow(A)==1&&A>0) return(as.matrix(sqrt(A))) # if A is a positive scalar (single value) Return its square-root
  } else {
    if (length(A) != 1) { # length of A is not 1 Return NULL
      warning("non-matrix argument to 'sechol'")
      return(NULL)
    }
    if (A>0) {# Returns square-root of non-negative values of A but if the leading minor of order 1 is not positive definite Return NULL
      return(as.matrix(sqrt(A)))
    } 
    warning("the leading minor of order 1 is not positive definite")
    return(NULL)
  }
  n <- nrow(A)
  L <- matrix(rep(0,n*n),ncol=ncol(A))
  tau <- tol ^(1/3)  # made to match gauss
  gamm <- max(A) # largest value in the matrix A
  deltaprev <- 0
  Pprod <- diag(n) # an n by n identity matrix
  if (n > 2)  {
    for (k in 1:(n-2))  {
      if( (min(diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm) # if the minimum of the Schur complement of A[k,k] is < tau*gamm (some small value)
          && (min(svd(A[(k+1):n,(k+1):n])$d)) < 0) {# minimum of the singular values of the matrix A[(k+1):n,(k+1):n]
        dmax <- order(diag(A[k:n,k:n]))[(n-(k-1))] # rearranges diag(A[k:n,k:n]) into ascending 
        if (A[(k+dmax-1),(k+dmax-1)] > A[k,k])  {
          if (!silent) {
            print(paste("iteration:",k,"pivot on:",dmax,"with absolute:",(k+dmax-1)))
          }
          P <- diag(n) # an n by n identity matrix
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        g <- rep(0,length=(n-(k-1)))
        for (i in k:n)  {
          if (i == 1) sum1 <- 0
          else sum1 <- sum(abs(A[i,k:(i-1)]))
          if (i == n) sum2 <- 0
          else sum2 <- sum(abs(A[(i+1):n,i]))
          g[i-(k-1)] <- A[i,i] - sum1 - sum2
        }
        gmax <- order(g)[length(g)]
        if (gmax != k)  {
          if (!silent) {
            print(paste("iteration:",k,"gerschgorin pivot on:",gmax,"with absolute:",(k+gmax-1)))
          }
          P <- diag(ncol(A))
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        normj <- sum(abs(A[(k+1):n,k]))
        delta <- max(0,deltaprev,-A[k,k]+max(normj,tau*gamm))
        if (delta > 0)  {
          A[k,k] <- A[k,k] + delta
          deltaprev <- delta
        }
      }
      L[k,k] <- A[k,k] <- sqrt(A[k,k])
      for (i in (k+1):n)  {
        L[i,k] <- A[i,k] <- A[i,k]/L[k,k]
        A[i,(k+1):i] <- A[i,(k+1):i] - L[i,k]*L[(k+1):i,k]
        if(A[i,i] < 0) A[i,i] <- 0
      }
    }
  }
  A[(n-1),n] <- A[n,(n-1)]
  eigvals <- eigen(A[(n-1):n,(n-1):n])$values
  delta <- max(0,deltaprev,
               -min(eigvals)+tau*max((1/(1-tau))*(max(eigvals)-min(eigvals)),gamm))
  if (delta > 0)  {
    if (!silent) {
      print(paste("delta:",delta))
    }
    A[(n-1),(n-1)] <- A[(n-1),(n-1)] + delta
    A[n,n] <- A[n,n] + delta
    deltaprev <- delta
  }
  L[(n-1),(n-1)] <- A[(n-1),(n-1)] <- sqrt(A[(n-1),(n-1)])
  L[n,(n-1)] <- A[n,(n-1)] <- A[n,(n-1)]/L[(n-1),(n-1)]
  L[n,n] <- A[n,n] <- sqrt(A[n,n] - L[n,(n-1)]^2)
  
  r = t(Pprod)%*%t(L)%*%t(Pprod)
  attr(r,"delta")=delta
  return(r)
}

#_______________________________________________________________________________________________________
# Returns the log determinant of a matrix A using the eigenvalues of the matrix A: log(det(A))
#_______________________________________________________________________________________________________
logdet =function(A){
  eigtem=eigen(A,only.values=TRUE,symmetric=TRUE)$values
  y=sum(log(eigtem[eigtem > 1e-10])) #[eigtem > 1e-10]
  return(y)
}

#_______________________________________________________________________________________________________
#   I have no idea what is happening here: HELP!!!!!!  ; check the papar titled "A note on fast envelope estimation"
#_______________________________________________________________________________________________________
spenvlp <- function(b2, b3, A1, A2, A3, lambda, eps, maxiter, weight, a_vec_init){
  #################################################################################
  #data setup
  lmax_A1 <- max(eigen(A1,only.values = TRUE)$values) # maximum eigenvalue of matrix A1
  lmax_A2 <- max(eigen(A2,only.values = TRUE)$values) # maximum eigenvalue of matrix A2
  lmax_A3 <- max(eigen(A3,only.values = TRUE)$values) # maximum eigenvalue of matrix A3
  gamma <- 4*lmax_A1 + 2*lmax_A2 + 2*lmax_A3 # linear combination of maximum eigenvalues of A1, A2, and A3
  
  #################################################################################
  nobs <- NROW(A1)
  dif <- rep(NA, nobs)
  a_vec = a_vec_init
  vecs = a_vec_init 
  olda_vec = a_vec
  npass <- 0
  while(1){ # while(1) is the same while(TRUE)
    olda_vec <- a_vec
    tmp1=A1%*%olda_vec
    tmp2=A2%*%(olda_vec+b2)
    tmp3=A3%*%(olda_vec+b3)
    U <- drop(4*tmp1/drop(1+olda_vec%*%tmp1) - 2*tmp2/drop(1+(olda_vec+b2)%*%tmp2) - 2*tmp3/drop(1+(olda_vec+b3)%*%tmp3))
    # the drop() function is used to delete the dimensions of an array which have only one level
    
    U_working <- (U + gamma * olda_vec)
    U_norm <- drop(sqrt(crossprod(U_working,U_working)))
    t <- U_norm - weight * lambda	
    if(t > 0){
      a_vec <- U_working * t / (gamma * U_norm)
    }else{
      a_vec <- rep(0,nobs)
    }
    vecs=rbind(vecs,a_vec)
    dif <- a_vec - olda_vec
    if(gamma*sum(dif^2) < eps) break
    
    npass = npass + 1
    if(npass > maxiter) break
  }
  ################################################################################
  # output
  outlist <- list(a_vec = a_vec, vecs=vecs)
  outlist
} 

#_______________________________________________________________________________________________________
#   A code for Gram Schmidt process
#_______________________________________________________________________________________________________
grams <- function(A){
  A=as.matrix(A)
  m=nrow(A)
  n=ncol(A)
  if(n==1) {R=sqrt(sum(A^2)); return(Q=A/R)} #return(list(Q=A/R,R=R)) # scaling the matrix A
  Asave = A
  Q=A
  eps=2.2204e-16
  for (j in 2:n){
    for (k in 1:(j-1)){    
      mult = (t(A[, j]) %*% A[,k]) / (sum(A[,k]^2))
      A[,j] = A[,j]-mult*A[,k];
    }
  }
  for (j in 1:n){
    if (sum(A[,j]^2)<sqrt(eps))  stop('Columns of A are linearly dependent.')
    Q[,j] = A[,j]/ sqrt(sum(A[,j]^2))
  }
  R = t(Q) %*% Asave
  return(Q)
  #return(list(Q=Q,R=R))
}

#_______________________________________________________________________________________________________
#   probably for finding the null basis of a matrix
#_______________________________________________________________________________________________________
nulbasis<-function(A){
  temp= rref(A);
  R=temp$A
  pivcol=temp$jb
  m = nrow(A)
  n = ncol(A)
  r = length(pivcol);
  freecol = 1:n;
  freecol=freecol[-pivcol]
  N = matrix(0,n,n-r)
  N[freecol,  ] = diag(n - r)
  N[pivcol, ] = - R[1 : r, freecol]
  return(N)
}

#_______________________________________________________________________________________________________
#   I think this is some kind of Gaussian elimination; row echelon reduction 
#_______________________________________________________________________________________________________
rref<-function(A){
  m = nrow(A)
  n = ncol(A)
  # Loop over the entire matrix.
  i = 1
  j = 1
  jb = c()
  tol = max(m,n)*2.2204e-016*norm(A,'I')
  while((i<= m)&(j<= n)){
    # Find value and index of largest element in the remainder of column j.
    p = max(abs(A[i:m,j]))
    k = which.max(abs(A[i:m,j]))
    k=k+i-1
    if (p<=tol){
      #The column is negligible, zero it out.
      A[i:m,j] = 0
      j = j + 1
    }
    else{
      # Remember column index
      jb = c(jb,j)
      # Swap i-th and k-th rows.
      A[c(i,k),j:n] = A[c(k,i),j:n]
      # Divide the pivot row by the pivot element.
      A[i,j:n] = A[i,j:n]/A[i,j];
      # Subtract multiples of the pivot row from all the other rows.
      for (k in (1:m)[-i]){
        A[k,j:n] = A[k,j:n] - A[k,j]*A[i,j:n];
      }  
      i = i + 1;
      j = j + 1;
    } #end else
  }#end while  
  result<-list()
  result$A=A
  result$jb=jb
  return(result)
}

#_______________________________________________________________________________________________________
# Returns the initial values for relevant subspace of predictor envelope
#_______________________________________________________________________________________________________
get_Init<-function(M,U,u,G=NULL){ 
  MU = M+U
  invMU<- chol2inv(chol(MU))
  obj<-function(G){ # objective function of interest
    logdet(crossprod(G,M)%*%G)+logdet(crossprod(G,invMU)%*%G)
  }
  
  # get init1  
  tmp.MU <- eigen(MU)
  tmp2.MU <- crossprod(tmp.MU$vectors, U) %*% tmp.MU$vectors
  tmp3.MU <- sort(diag(tmp2.MU), decreasing = TRUE, index.return = TRUE)
  init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])   
  obj1 <- obj(init)
  if(!missing(G)) print(subspace(init,G))
  
  # get init2	
  tmp2 <- diag(1/sqrt(tmp.MU$values)) 
  tmp3 <- sort(diag(tmp2 %*% tmp2.MU %*% tmp2), decreasing = TRUE, index.return = TRUE)
  init.MU <- as.matrix(tmp.MU$vectors[, tmp3$ix[1:u]])
  if(!missing(G))print(subspace(init.MU,G))
  
  obj2 <- obj(init.MU)
  if (obj2 < obj1) {
    init <- init.MU
    obj1 <- obj2
  }
  
  # get init3
  tmp.M <- eigen(M)
  tmp2.M <- crossprod(tmp.M$vectors, U) %*% tmp.M$vectors
  tmp3.M <- sort(diag(tmp2.M), decreasing = TRUE, index.return = TRUE)	
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  obj3 <-  obj(init.M)
  if (obj3 < obj1) {
    init <- init.M
    obj1 <- obj3
  }
  if(!missing(G)) print(subspace(init.M,G))
  
  # get init4				
  tmp2.M <- diag(1/sqrt(tmp.M$values)) 
  tmp3.M <- sort(diag(tmp2 %*% tmp2.M %*% tmp2), decreasing = TRUE, index.return = TRUE)
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  if(!missing(G)) print(subspace(init.M,G))
  obj4 <-  obj(init.M)
  if (obj4 < obj1) {
    init <- init.M
    obj1 <- obj4
  }      
  return(init) 
}

#_______________________________________________________________________________________________________
# Returns the initial values for sparse convariance matrices;  init_method=1 for covariance of Y or 
#init_method=2 for sparse convariance matrices of X and X|Y?
#_______________________________________________________________________________________________________
init_spxenv <- function(X, Y, u,spice=NULL,lambda_spice=0.1,init_method=2){
  X <-as.matrix(X)
  Y <-as.matrix(Y)
  p=ncol(X)
  
  if(init_method==1){init=initial_value(Y,X,u)}
  else if(init_method==2) {
    if(is.null(spice)) {spice=spxenv_spice(X,Y,lambda_spice)}
    init=get_Init(spice$sigXcY,spice$sigX-spice$sigXcY,u)
  }
  #else if(init_method==3)  {init=initial_spls(X,Y,u)}
  #else if(init_method==4) {init = unclass(pls::plsr(Y~X,ncomp=u)$loadings)}
  #else if(init_method==5) {
  #  cv2=glmnet::cv.glmnet(X,Y,nfolds=5)
  #  fit2=glmnet::glmnet(X,Y,lambda = cv2$lambda.min)
  #  where1=which(as.vector(fit2$beta)!=0)
  #  XA=X[,where1]
  #  init = matrix(0,p,u)
  #  tmp = unclass(pls::plsr(Y~XA,ncomp=u)$loadings) #initial_value(Y,XA,u)
  #  init[where1,]=tmp
  #  }	  
  return(init)
}

#_______________________________________________________________________________________________________
# Returns the initial values for relevant subspace of response envelope
#_______________________________________________________________________________________________________
initial_value <- function(X, Y, u,verbose=0,choose=NULL,G=NULL) {
  X <-as.matrix(X)
  Y <-as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  #mX = colMeans(X)
  #mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  sigY <- stats::cov(Yc)*(n-1)/n
  sigX <- stats::cov(Xc)*(n-1)/n 
  sigYX <- stats::cov(Yc, Xc)*(n-1)/n
  if(n>r) invsigY <- chol2inv(chol(sigY))
  
  tmp = sechol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
  tmp3 <- sigYX %*% tmp2
  bsxb <- tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  sigRes <- sigY - bsxb
  
  obj<-function(G){
    tmp1=crossprod(G,sigRes)%*%G
    tmp2=crossprod(G,invsigY)%*%G
    logdet(tmp1)+logdet(tmp2)
  }
  
  #=============
  
  tmp.y <- eigen(sigY)
  tmp2.y <- crossprod(tmp.y$vectors, bsxb) %*% tmp.y$vectors
  tmp3.y <- sort(diag(tmp2.y), decreasing = TRUE, index.return = TRUE)
  ##print(sort(tmp3.y$ix[1:u]))
  init <- as.matrix(tmp.y$vectors[, tmp3.y$ix[1:u]]) 
  if(!missing(G)) print(subspace(G,init))
  #if(verbose) {cat('init1\n');print(init)}
  k <- 1	
  #print(init)	
  if (n > r + 1) {
    obj1 <- obj(init)			
    tmp2 <- diag(1 / sqrt(tmp.y$values))
    tmp3 <- sort(diag(tmp2 %*% tmp2.y %*% tmp2), decreasing = TRUE, index.return = TRUE)
    ##print(sort(tmp3$ix[1:u]))
    init.y <- as.matrix(tmp.y$vectors[, tmp3$ix[1:u]])
    if(!missing(G)) print(subspace(G,init.y))
    #if(verbose) {cat('init2\n');print(init.y)}
    obj2 <- sum(obj(init.y))		
    if (obj2 < obj1) {
      init <- init.y
      obj1 <- obj2
      k <- 2
    }
    if (n > r + p + 1) {
      tmp.res <- eigen(sigRes)
      tmp2.res <- crossprod(tmp.res$vectors, bsxb) %*% tmp.res$vectors
      tmp3.res <- sort(diag(tmp2.res), decreasing = TRUE, index.return = TRUE)	
      init.res <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])				
      #if(verbose) {cat('init3\n');print(init.res)}
      if(!missing(G)) print(subspace(G,init.res))
      obj3 <- obj(init.res)		
      if (obj3 < obj1) {
        init <- init.res
        obj1 <- obj3
        k <- 3
      }
      
      tmp2 <- diag(1 / sqrt(tmp.res$values)) 
      tmp3.res <- sort(diag(tmp2 %*% tmp2.res %*% tmp2), decreasing = TRUE, index.return = TRUE)
      init.res2 <- as.matrix(tmp.res$vectors[, tmp3.res$ix[1:u]])
      #if(verbose) {cat('init4\n');print(init.res)}
      if(!missing(G)) print(subspace(G,init.res2))
      obj4 <- obj(init.res2)
      if (obj4 < obj1) {
        init <- init.res2
        obj1 <- obj4
        k <- 4
      }
    }
  }
  if(verbose)	cat('We choose init',k,'\n')		
  as.matrix(init)
}

#_______________________________________________________________________________________________________
#   The base algorithm for sparse envelope; check the papar titled "A note on fast envelope estimation"
#_______________________________________________________________________________________________________
spxenvbase <- function(X, Y, u,lambda, weight=NULL, ftol=1e-2, maxiter=1e2,eps2=1e-4,  init=NULL,verbose=verbose){
  X =as.matrix(X)
  Y =as.matrix(Y)
  n <- nrow(Y)
  r <- ncol(Y)
  p <- ncol(X)
  if(missing(weight)) weight=rep(1,p-u)
  mX = colMeans(X)
  mY = colMeans(Y)
  Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
  Xc <- as.matrix(scale(X, center = T, scale = FALSE))
  
  sigY <- stats::cov(Yc) * (n-1)/n
  sigX <- stats::cov(Xc) * (n-1)/n
  sigXY <- stats::cov(Xc, Yc) * (n-1)/n
  
  tmp = chol(sigX) # t(tmp)%*%tmp =sigX
  invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
  tmp3 <- crossprod(sigXY, tmp2)
  sigYcX = sigY - tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
  
  tmp = chol(sigY) # t(tmp)%*%tmp =sigY
  invsigY <- chol2inv(tmp) # invsigY=tmp2%*%t(tmp2)
  tmp2 <- backsolve(tmp, diag(r)) # tmp2=inv(tmp)  
  tmp3 <- sigXY %*% tmp2
  U <- tcrossprod(tmp3,tmp3) #sigXY %*% invsigY%*% t(sigXY)
  sigXcY <- sigX - U
  
  betaOLS <- invsigX%*%sigXY
  logDetsigY = logdet(sigY)  
  logDetsigX = logdet(sigX)
  ModelOutput=list()  
  if (u == 0) {
    ModelOutput$mu = mY
    ModelOutput$beta = matrix(0,p, r)
    ModelOutput$Gamma = NULL
    ModelOutput$Gamma0 = diag(p)
    ModelOutput$eta = NULL
    ModelOutput$SigX = sigX
    ModelOutput$Omega = NULL
    ModelOutput$Omega0 = sigX
    ModelOutput$sigYcX = sigYcX 
    ModelOutput$loglik = - n * (r + p) / 2 * (1 + log(2 * pi)) - n / 2 *  (logDetsigX + logDetsigY)
    ModelOutput$paramNum = r + p * (p + 1) / 2 + r * (r + 1) / 2
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$p = p
    ModelOutput$u = u
    ModelOutput$q = 0
    ModelOutput$where1 = NULL
    ModelOutput$where0 = 1:p
    ModelOutput$lambda = lambda
  } else if (u == p) {
    beta = betaOLS;
    ModelOutput$mu = mY - t(beta) %*% mX
    ModelOutput$beta = betaOLS
    ModelOutput$Gamma = diag(p)
    ModelOutput$Gamma0 = NULL
    ModelOutput$eta = beta
    ModelOutput$SigX = sigX;
    ModelOutput$Omega = sigX
    ModelOutput$Omega0 = NULL
    ModelOutput$sigYcX = sigYcX 
    ModelOutput$loglik =  -n * (p + r) * (1 + log(2 * pi)) / 2 - n / 2 * (logDetsigY + logdet(sigXcY));
    ModelOutput$paramNum = r + (p + r) * (p + r + 1) / 2;
    ModelOutput$n = n
    ModelOutput$r = r
    ModelOutput$p = p
    ModelOutput$u = u
    ModelOutput$q = p
    ModelOutput$where1 = 1:p
    ModelOutput$where0 = NULL
    ModelOutput$lambda = lambda
  } 
  else{
    if(missing(init)) {init <- initial_value(Y,X,u)}
    Ginit <- init %*% solve(init[1:u, ]) 
    
    obj1 <- logdet(t(init) %*% sigXcY %*% init) + logdet((t(init) %*% invsigX %*% init))
    widx <- is.finite(weight)      
    Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
    Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
    obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
    
    if (u == (p-1)) { 
      U1c2 <- array(0, dim = c(p-1, p-1))
      V1c2 <- array(0, dim = c(p-1, p-1))
      
      U1c2 <- sigXcY[-p, -p] - as.matrix(sigXcY[-p, p]) %*% sigXcY[p, -p] / sigXcY[p, p]
      V1c2 <- invsigX[-p, -p] - as.matrix(invsigX[-p, p]) %*% invsigX[p, -p] / invsigX[p, p]		
      
      t2 <- sigXcY[-p, p] / sigXcY[p, p]
      t3 <- invsigX[-p, p] / invsigX[p, p]
      invC1 <- chol2inv(sechol(U1c2))
      invC2 <- chol2inv(sechol(V1c2))
      
      i <- 1
      while (i < maxiter) {        
        res <- spenvlp(b2=drop(t2), b3=drop(t3), 
                       A1=diag(p-1), A2=sigXcY[p, p]*invC1, A3=invsigX[p, p]*invC2, 
                       lambda=lambda, eps=eps2, maxiter=1e2, 
                       weight=weight[p-u], 
                       a_vec_init=drop(Ginit[p,]))        
        old_Ginit <- Ginit[p, ]
        Ginit[p, ] <- res$a_vec        
        
        
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigXcY %*% Gamma) + logdet((t(Gamma) %*% invsigX %*% Gamma)) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          break
        }  
        else {
          obj1 <- obj5
          i <- i + 1
        }	
        #if(sum((Ginit[p,]-old_Ginit)^2) < eps) break
        #i <- i + 1		
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
      a <- qr.Q(qr(Ginit), complete = TRUE)
      Gamma <- a[, 1:u]
      Gamma0 <- a[, p]
    } 
    else {
      #obj1 <- logdet(t(init) %*% sigXcY %*% init) + logdet((t(init) %*% invsigX %*% init))
      #widx <- is.finite(weight)      
      #Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
      #Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
      #obj1 <- obj1 + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
      
      GUG <- crossprod(Ginit, (sigXcY %*% Ginit))	
      GVG <- crossprod(Ginit, (invsigX %*% Ginit))		
      
      t4 <- crossprod(Ginit[(u+1):p,], Ginit[(u+1):p, ]) + diag(u)
      i <- 1
      while (i < maxiter) {
        #print(i+1)
        for (j in (u+1):p) {
          g <- as.matrix(Ginit[j, ])
          t2 <- crossprod(Ginit[-j, ], as.matrix(sigXcY[-j, j])) / sigXcY[j, j]
          t3 <- crossprod(Ginit[-j, ], as.matrix(invsigX[-j, j])) / invsigX[j, j]
          
          GUGt2 <- g + t2
          GUG <- GUG - tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG - tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
          
          t4 <- t4 - tcrossprod(as.matrix(Ginit[j, ]), as.matrix(Ginit[j, ]))
          invC1 <- chol2inv(sechol(GUG))
          invC2 <- chol2inv(sechol(GVG))
          invt4 <- chol2inv(chol(t4))
          
          res <- spenvlp(b2=drop(t2), b3=drop(t3), A1=invt4, A2=sigXcY[j, j]*invC1, A3=invsigX[j, j]*invC2, lambda=lambda, eps=eps2, maxiter=1e2, weight=weight[j-u], a_vec_init=drop(Ginit[j,]))
          
          
          Ginit[j, ] <- res$a_vec
          g <- as.matrix(Ginit[j, ])
          # print(g)
          t4 <- t4 + tcrossprod(g, g)
          GUGt2 <- g + t2
          GUG <- GUG + tcrossprod(GUGt2, GUGt2) * sigXcY[j, j]
          
          GVGt2 <- g + t3
          GVG <- GVG + tcrossprod(GVGt2, GVGt2) * invsigX[j, j] 
        }
        a <- qr(Ginit)
        Gamma <- qr.Q(a)
        Ginit.tmp <- as.matrix(Ginit[(u+1):p, ])
        Ginit.tmp <- as.matrix(Ginit.tmp[widx,])
        obj5 <- logdet(t(Gamma) %*% sigXcY %*% Gamma) + logdet((t(Gamma) %*% invsigX %*% Gamma)) + lambda * sum(weight[widx] * sqrt(rowSums(Ginit.tmp^2)))
        if (abs(obj1 - obj5) < ftol * abs(obj1)) {
          break
        }  
        else {
          obj1 <- obj5
          i <- i + 1
        }	
      }
      if(verbose) cat('The number of iterations:',i,'.\n',sep='')
    }
    Gamma <- as.matrix(Gamma)
    if(!is.na(sum(Gamma))) idx <- which(rowSums(abs(Gamma))>0)  
    q = length(idx)
    idx_i = setdiff(1:p, idx);
    if(q>u){
      #---Compute the rest of the parameters based on \Gamma---
      Gamma0 = grams(nulbasis(t(Gamma)))
      Omega = crossprod(Gamma,sigX) %*% Gamma
      invOmega = chol2inv(chol(Omega))
      eta = invOmega %*% t(Gamma) %*% sigXY; 
      Omega0 = crossprod(Gamma0,sigX) %*% Gamma0
      SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
      beta = Gamma %*%  eta
      mu = mY - t(beta) %*% mX      
      # calculate likelihhod
      a = logdet(t(Gamma) %*% sigXcY %*% Gamma);
      b = logdet(t(Gamma) %*% invsigX %*% Gamma);
      l = n * (p + r) * (1 + log(2 * pi)) + n * (a + b + logDetsigX + logDetsigY)
      paramNum = r + u * (r - p + q) + p * (p + 1) / 2 + r * (r + 1) / 2
    }
    else{
      sigX_work=sigX[idx,idx]
      invsigX_work=chol2inv(chol(sigX_work))
      sigXY_work=sigXY[idx,,drop=FALSE]
      sigYcX = sigY - t(sigXY_work) %*% invsigX_work %*% sigXY_work
      Gamma_work = diag(u);
      Gamma = matrix(0,p,u)
      Gamma[idx,] = Gamma_work;        
      Gamma0 = grams(nulbasis(t(Gamma)))
      Omega = t(Gamma) %*% sigX %*% Gamma;
      Omega0 = t(Gamma0) %*% sigX %*% Gamma0
      SigX = Gamma %*% Omega %*% t(Gamma) + Gamma0 %*% Omega0 %*% t(Gamma0)
      beta_work= invsigX_work %*% sigXY_work
      eta=beta_work
      beta=matrix(0,p,r)
      beta[idx,]=beta_work      
      mu = mY - t(beta) %*% mX;
      # calculate BIC and likelihhod
      l = n * (p + r) * (1 + log(2 * pi)) + n * (logdet(sigYcX) + logDetsigX)
      paramNum = r + q * (r - p + q) + p * (p + 1) / 2 + r * (r + 1) / 2
    } #q=u
    ModelOutput$mu = mu;
    ModelOutput$beta = beta
    ModelOutput$Gamma = Gamma
    ModelOutput$Gamma0 = Gamma0
    ModelOutput$eta = eta
    ModelOutput$SigX = SigX
    ModelOutput$Omega = Omega
    ModelOutput$Omega0 = Omega0
    ModelOutput$sigYcX = sigYcX
    ModelOutput$loglik = - 0.5 * l
    ModelOutput$paramNum = paramNum
    ModelOutput$n = n
    ModelOutput$q=q
    ModelOutput$iternum=i
    ModelOutput$lambda=lambda    
    ModelOutput$where1=idx
    ModelOutput$where0=idx_i  
  } # 1<u<r
  ModelOutput
}

#_______________________________________________________________________________________________________
#   algorithm for penalty function
#_______________________________________________________________________________________________________
LassoLambda.spxenv <- function(X, Y, u, ftol=1e-2, maxiter=1e2,eps2=1e-4,lambda=NULL, weight = NULL,init=NULL,verbose=1){
  #t1 = proc.time()
  X = as.matrix(X)
  Y = as.matrix(Y)
  p = ncol(X)
  if(missing(init))     init=initial_value(Y,X,u)
  
  if(missing(lambda)) lambda <- exp(seq(log(1),log(1e-5),len=15)) 
  if(missing(weight)) weight <- rep(1,p-u)
  
  model_vec <- rep(NA, length(lambda))
  minBIC = Inf
  res=NULL
  idmin = NULL
  for(l in 1:length(lambda)){
    if(verbose) cat('------lambda =',lambda[l],'\n')
    tmp <- spxenvbase(X, Y, u, lambda=lambda[l], ftol=ftol,	maxiter=maxiter, eps2=eps2,weight=weight,init=init,verbose=verbose)
    init = tmp$Gamma
    BIC <- -2*tmp$loglik + log(tmp$n) * (tmp$q-u) * u
    model_vec[l] <- BIC
    if(BIC<minBIC) {minBIC=BIC;res=tmp;idmin=l}
  }
  lambda.min <- lambda[idmin]
  out <- res
  out$lambda=lambda.min  
  if(verbose) {
    print(model_vec)
    cat('Based on BIC we choose lambda to be',lambda.min,'.\n')
  }
  out$BIC_seq=model_vec
  #out$fit_time=(proc.time()-t1)[3]
  out
}


#___________algorithm for E-SPLS when n is greater than p__________________
spxenv <- function(X, Y, u, lambda1=NULL,lambda2=NULL, ftol=1e-2, maxiter=1e2,eps2=1e-4,init=NULL,asy=0,verbose=0,init_method=1,lambda_spice=0.1){
  t1 = proc.time()
  X=as.matrix(X)
  Y=as.matrix(Y)
  p = ncol(X)
  r = ncol(Y)
  n = nrow(X)
  if(u==0||u==p) out=xenv(X,Y,u,0)
  else{
    if(missing(lambda1)) lambda1 <- exp(seq(log(1),log(1e-3),len=5))
    if(missing(lambda2)) lambda2 <- exp(seq(log(1),log(1e-3),len=5))
    if(missing(init)){
      init=init_spxenv(X,Y,u=u,lambda_spice=lambda_spice,init_method=init_method)
    }
    GEidx = GE(init)	
    newX = X[, GEidx,drop=FALSE]	
    newinit=init[GEidx,,drop=FALSE]
    
    m1 <- LassoLambda.spxenv(newX, Y, u, ftol=ftol, maxiter=maxiter, eps2=eps2, lambda=lambda1, weight = rep(1,p-u),init=newinit,verbose=verbose)    #calculating weight     
    
    #calculating weight
    Gammahat <- m1$Gamma
    w <- Gammahat %*% solve(Gammahat[1:u, ])
    w_norm <- 1/(rowSums(w^2)^2)[(u+1):p]
    #w_norm2 <- 1/(rowSums(Gammahat^2)^2)[(u+1):p]
    if(verbose) {print(w);cat('weights:',w_norm,'\n')}       
    #adaptive lasso step     
    m2 <- LassoLambda.spxenv(newX, Y, u,  lambda=lambda2,ftol=ftol, maxiter=maxiter,eps2=eps2, weight =
                               w_norm,init=Gammahat,verbose=verbose)
    out=m2
    out$beta=m2$beta[order(GEidx),,drop=FALSE]
    out$Gamma=m2$Gamma[order(GEidx),,drop=FALSE]
    out$Gamma0=m2$Gamma0[order(GEidx),,drop=FALSE]    
    rownames(out$beta)=colnames(X)
    rownames(out$Gamma)=colnames(X) 
    out$where1=     sort(GEidx[m2$where1])     
    out$where0= setdiff(1:p,sort(GEidx[m2$where1]))
    out$BIC_seq = NULL
    out$BIC_seq1=m1$BIC_seq
    out$BIC_seq2=m2$BIC_seq
    lambdas=c(m1$lambda,m2$lambda)
    names(lambdas)=c('first stage','second stage')
    out$lambda=lambdas   
  }
  if(asy){
    Yc <- as.matrix(scale(Y, center = T, scale = FALSE))
    Xc <- as.matrix(scale(X, center = T, scale = FALSE))
    sigY <- stats::cov(Yc) * (n-1)/n
    sigX <- stats::cov(Xc) * (n-1)/n
    sigXY <- stats::cov(Xc, Yc) * (n-1)/n
    tmp = chol(sigX) # t(tmp)%*%tmp =sigX
    invsigX <- chol2inv(tmp) # invsigX=tmp2%*%t(tmp2)
    tmp2 <- backsolve(tmp, diag(p)) # tmp2=inv(tmp)  
    tmp3 <- crossprod(sigXY, tmp2)
    sigYcX = sigY - tcrossprod(tmp3,tmp3) #sigYX %*% invsigX%*% t(sigYX)
    betaOLS <- invsigX%*%sigXY
    if(u==0){
      out$covMatrix = NULL
      out$asySE = NULL
      out$ratio = matrix(1,p, r)
    } # u==0
    else if(u==p){
      covMatrix = kronecker(sigYcX, invsigX);
      asyFm = matrix(sqrt(diag(covMatrix)), p,r)
      out$covMatrix = covMatrix
      out$asySE = asyFm
      out$ratio = matrix(1, p, r)
    } # u==p
    else{
      q = out$q
      idx=out$where1
      idx2 = rep(0,q*r)
      for (k in 1:r)    idx2[((k-1)*q+1) : (k*q)]= idx+(k-1)*p    
      if(q>u){
        eta=out$eta
        Omega=out$Omega
        Gamma=out$Gamma
        invsigYcX=chol2inv(chol(sigYcX))
        invOmega=chol2inv(chol(Omega))
        Gamma_work=Gamma[idx,]
        Gamma0_work=grams(nulbasis(t(Gamma_work)))
        #---compute asymptotic variance and get the ratios---
        asyFm = kronecker(sigYcX, invsigX)
        asyFm = matrix(sqrt(diag(asyFm)), p, r)
        Omega0_1 = t(Gamma0_work) %*% sigX[idx, idx] %*% Gamma0_work
        Omega0_2 = sigX[idx_i, idx_i];  # r-q by r-q
        Omega0_12 = t(G0_work) %*% sigX[idx, idx_i]; # q-u by r-q
        Omega0_1c2 = Omega0_1 - Omega0_12 %*% solve(Omega0_2) %*% t(Omega0_12); # q-u by q-u
        invOmega0_1c2 = chol2inv(chol(Omega0_1c2))
        
        temp = kronecker(eta %*% invsigYcX %*% t(eta), Omega0_1)+
          kronecker(Omega, invOmega0_1c2) + 
          kronecker(invOmega, Omega0_1) - 2 * kronecker(diag(u), diag(q - u))
        invtemp=chol2inv(chol(temp ))
        covMatrix_work = kronecker(sigYcX, Gamma_work %*% invOmega %*% t(Gamma_work))+
          kronecker(t(eta), Gamma0_work) %*% invtemp %*% kronecker(eta, t(Gamma0_work))      
        covMatrix = matrix(0,p*r,p*r);   #pr by pr
        covMatrix[idx2,idx2] = covMatrix_work;
        asySE = matrix(sqrt(diag(covMatrix)), p, r)
        out$covMatrix = covMatrix
        out$asySE = asySE
        out$ratio = asyFm / asySE
      } #q>u
      else{
        asyFm = kronecker(sigYcX,invsigX); # rp by rp
        asyFm = matrix(sqrt(diag(asyFm)), p, r); #r by p
        
        sigX_work=sigX[idx,idx]
        invsigX_work=chol2inv(chol(sigX_work))
        sigXY_work=sigXY[idx,,drop=FALSE]
        sigYcX_work = sigY - t(sigXY_work) %*% invsigX_work %*% sigXY_work
        covMatrix_work = kronecker(sigYcX_work, invsigX_work) # qp by qp
        covMatrix=matrix(0,r*p,r*p)  #rp by rp
        covMatrix[idx2,idx2]=covMatrix_work;
        asySE = matrix(sqrt(diag(covMatrix)), p, r); #r by p
        out$covMatrix = covMatrix
        out$asySE = asySE
        out$ratio = asyFm / asySE
      } #q=u
    }
  } # asy==1
  out$fit_time=as.numeric((proc.time()-t1)[3])
  class(out)='spxenv'
  out
}

zhu.beta <- lapply(1:runs, function(i) spxenv(X[[i]], Y[[i]], u = q, lambda1 = 1, lambda2=2000)$beta)
zhu.beta

#___________________vectorize beta_________________________________________
new.zhu.beta <- NULL
for(i in 1:runs){
  new.zhu.beta[[i]] <- as.vector(round(zhu.beta[[i]], 3))
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
  zhu.oracle.beta[i,] <- as.vector(round(zhu.beta[[i]][1:pa,], 7))
}

zhu.bias <- colMeans(zhu.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.zhu <- norm_vec(as.vector(zhu.bias))^2
bias.zhu


#________________ standard deviation of estimates ______________________
zhu.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  zhu.oracle.beta[i,] <- as.vector(round(zhu.beta[[i]][1:pa,], 7))
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
  spls.beta[i] <- list( spls(X[[i]], Y[[i]], q, eta = 0.6, kappa=0.5, select="simpls", fit="simpls", 
                             scale.x=F, scale.y=F, eps=1e-4, maxstep=5000, trace=F)$betamat )
}
#spls.beta

#___________________vectorize beta_________________________________________
new.spls.beta <- NULL
for(i in 1:runs){
  new.spls.beta[[i]] <- as.vector(round(spls.beta[[i]][[q]], 3))
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
  spls.oracle.beta[i,] <- as.vector(round(spls.beta[[i]][[q]][1:pa,], 7))
}

spls.bias <- colMeans(spls.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.spls <- norm_vec(as.vector(spls.bias))^2
bias.spls



#________________ standard deviation of estimates ______________________
spls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  spls.oracle.beta[i,] <- as.vector(round(spls.beta[[i]][[q]][1:pa,], 7))
}
#spls.oracle.beta 

spls.total.variance <- sum(eigen(cov(spls.oracle.beta))$values)
spls.total.variance
spls.std <- sqrt(spls.total.variance)
spls.std



#_________________________________________________________________________________
#_________ SIMPLS ________________________________________________________________
mysimpls <- lapply(1:runs, function(i) plsr(Y[[i]] ~ X[[i]], ncomp = q, method = "simpls")$coefficients)

#________________ bias of estimates ______________________
simpls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  simpls.oracle.beta[i,] <- as.vector(round(mysimpls[[i]][,,q][1:pa,], 7))
}

simpls.bias <- colMeans(simpls.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.simpls <- norm_vec(as.vector(simpls.bias))^2
bias.simpls



#________________ standard deviation of estimates ______________________
simpls.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  simpls.oracle.beta[i,] <- as.vector(round(mysimpls[[i]][,,q][1:pa,], 3))
}
#simpls.oracle.beta 

simpls.total.variance <- sum(eigen(cov(simpls.oracle.beta))$values)
simpls.total.variance
simpls.std <- sqrt(simpls.total.variance)
simpls.std


#_________________________________________________________________________________
#______________________Envelope model_____________________________________________
cook <- lapply(1:runs, function(i) xenv(X[[i]], Y[[i]], u = q, asy = F)$beta)


#________________ bias of estimates ______________________
cook.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  cook.oracle.beta[i,] <- as.vector(round(cook[[i]][1:pa,], 7))
}

cook.bias <- colMeans(cook.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.cook <- norm_vec(as.vector(cook.bias))^2
bias.cook


#________________ standard deviation of estimates ______________________
cook.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  cook.oracle.beta[i,] <- as.vector(round(cook[[i]][1:pa,], 7))
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
  lm.oracle.beta[i,] <- as.vector(round(lm.beta[[i]][1:pa,], 7))
}

lm.bias <- colMeans(lm.oracle.beta) - as.vector(Beta[1:pa,])
norm_vec <- function(x) sqrt(sum(x^2))
bias.lm <- norm_vec(as.vector(lm.bias))^2
bias.lm

#________________ standard deviation of estimates ______________________
lm.oracle.beta <- matrix(NA, runs, r*pa)
for(i in 1:runs){
  lm.oracle.beta[i,] <- as.vector(round(lm.beta[[i]][1:pa,], 7))
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

