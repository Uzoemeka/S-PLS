#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Data generation
rm( list = ls ())
runs <- 50
r <- 3
q <- 3
p <- 12
pa <- 4
pi <- p - pa
pa.q <- pa - q
n <- 10000
delt <- 0.95

gram_schmidt <- function(x) {
  x <- as.matrix(x)
  m <- nrow(x)
  n <- ncol(x)
  q <- matrix(0, m, n)
  r <- matrix(0, n, n)
  
  process_column <- function(j) {
    v <- x[, j]
    
    if (j > 1) {
      lapply(1:(j - 1), function(i) {
        r[i, j] <<- t(q[, i]) %*% x[, j]
        v <<- v - r[i, j] * q[, i]
      })
    }
    
    r[j, j] <<- sqrt(sum(v^2))
    q[, j] <<- v / r[j, j]
  }
  
  lapply(1:n, process_column)
  
  list(Q = q, R = r)
}


# Generating Gamma matrices using Gram-Schmidt
gamma.pa <- gramschmidt(matrix(rnorm(pa * pa, 0, 1), pa, pa))$Q
gamma.a1 <- rbind(gamma.pa[, 1:q], matrix(0, pi, q))
gamma.a2 <- rbind(matrix(gamma.pa[, (q + 1):pa], pa, pa.q), matrix(0, pi, pa.q))
gamma.i <- rbind(matrix(0, pa, pi), diag(1, pi))
gamma0 <- cbind(gamma.a2, gamma.i)
gamma.both <- cbind(gamma.a1, gamma.a2, gamma.i)

# Define omega matrices
omega1 <- diag(4, q)
omega.01 <- diag(0.8, pa.q)
omega.02 <- diag(1, pi)
omega0 <- as.matrix(bdiag(omega.01, omega.02))

# Compute sigx
sigx <- gamma.a1 %*% omega1 %*% t(gamma.a1) + gamma0 %*% omega0 %*% t(gamma0)

# Create sigy.x matrix using mapply
sigy.x <- matrix(0, r, r)
sigy.x <- mapply(function(i, j) delt^(abs(i - j)), rep(1:r, each = r), rep(1:r, times = r))
dim(sigy.x) <- c(r, r)

# Beta matrices
Beta <- rbind(matrix(c(-3, 5, 1,
                       3, 0, 0,
                       2, 1, -2,
                       4, 0, 3), pa, r, byrow = TRUE),
              matrix(0, pi, r))
Beta1 <- rbind(matrix(c(-3, 5, 1,
                        3, 5e-4, 5e-4,
                        2, 1, -2,
                        4, 5e-4, 3), pa, r, byrow = TRUE),
               matrix(5e-4, pi, r))

# Generating lists using lapply
X <- lapply(1:runs, function(k) mvrnorm(n, rep(0, p), sigx))
E <- lapply(1:runs, function(k) mvrnorm(n, rep(0, r), sigy.x))
Y <- lapply(1:runs, function(k) scale(X[[k]] %*% Beta + E[[k]]))
X <- lapply(X, scale)

# Compute covariance and other statistics using lapply
sigxx <- lapply(X, function(x) cov(x) * (n - 1) / n)
sigxxyy <- lapply(1:runs, function(k) cov(X[[k]], Y[[k]]) * (n - 1) / n)
sigyy <- lapply(Y, function(y) cov(y) * (n - 1) / n)

# Calculate Sxy and Sigs using lapply
Sxy <- lapply(1:runs, function(k) {
  sigxx[[k]] - sigxxyy[[k]] %*% solve(sigyy[[k]]) %*% t(sigxxyy[[k]])
})

Sigs <- lapply(1:runs, function(k) {
  list(Sxy = Sxy[[k]], sigxx = sigxx[[k]], sigxxyy = sigxxyy[[k]])
})


#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Initial value for 2S - SPLS
myenvCv <- function(X, Y) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  # Get dimensions of Y
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  # Generate the list of fitted values using lapply
  efit <- lapply(1:p, function(j) xenv(X, Y, u = j, asy = FALSE))
  
  return(efit)
}

# Apply `myenvCv` function for each run
em <- lapply(1:runs, function(i) myenvCv(X[[i]], Y[[i]]))



#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 2S - SPLS
Gamma.fixed1 <- function(par, X, Y, G, w, low, high, intr) {
  n <- nrow(X)
  r <- ncol(Y)
  p <- ncol(X)
  
  obj <- function(par, X, Y, G, pen, w) {
    X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
    Y <- scale(as.matrix(Y), center = TRUE, scale = FALSE)
    G <- as.matrix(G)
    
    B <- matrix(par, p, r)
    Z <- X %*% G
    sigZY <- cov(Z, Y) * (n - 1) / n
    sigZ <- cov(Z) * (n - 1) / n
    sigY <- cov(Y) * (n - 1) / n
    sigY.Z <- sigY - t(sigZY) %*% solve(sigZ) %*% sigZY
    
    g.norm <- function(x1, x2) sum(abs(x1) * (1 / abs(x2)^0.5))
    g.norm1 <- function(x) sum(abs(x))
    g.norm2 <- function(x) sum(sqrt(rowSums(x^2)))
    g.normp <- function(x1, x2) sqrt(rowSums(x1^2)) * (1 / sqrt(rowSums(x2^2)))^2
    
    J <- tr((1 / n) * t(Y - Z %*% t(G) %*% B) %*% (Y - Z %*% t(G) %*% B) %*% solve(sigY.Z)) +
      pen * g.norm(B, w) +
      0.0005 * g.norm2(B)
    
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norma (B , w )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norm1 ( B )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norm2 ( B ) + pen *g . norm1 ( B )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B )) + pen *g . norm2 ( B )
    
    
    return(J)
  }
  
  gridd <- seq(low, high, intr)
  fit4 <- lapply(gridd, function(k) nlm(f = obj, p = par, X = X, Y = Y, G = G, pen = k, w = w,
                                        iterlim = 2000, gradtol = 1e-16))
  
  mini.fit <- sapply(fit4, function(fit) fit$minimum)
  
  para1 <- lapply(fit4, function(fit) round(matrix(fit$estimate, p, r), 5))
  
  bic.fit <- sapply(1:length(gridd), function(k) mini.fit[k] + 
                      length(as.vector(para1[[k]])[which(as.vector(para1[[k]]) != 0)]) * log(n))
  
  bic.fit.min <- which.min(bic.fit)
  
  opt.Beta <- para1[[bic.fit.min]]
  opt.lambda <- gridd[bic.fit.min]
  
  return(list(para1 = para1, mini.fit = mini.fit, bic.fit = bic.fit, bic.fit.min = bic.fit.min, 
              opt.Beta = opt.Beta, opt.lambda = opt.lambda))
}

# Applying the function in parallel (e.g., using lapply)
gf <- lapply(1:runs, function(i) Gamma.fixed1(
  em[[i]][[q]]$beta, X[[i]], Y[[i]], em[[i]][[q]]$Gamma, w = em[[i]][[q]]$beta, 
  low = 0, high = 0.007, intr = 0.001))

gf.beta1 <- sapply(gf, function(g) round(g$opt.Beta, 5))

gf.beta <- lapply(gf, function(g) as.vector(round(g$opt.Beta, 5)))

# Confusion matrix function for classification performance
my.table <- function(actual, predicted) {
  actual <- as.vector(actual)
  predicted <- as.vector(predicted)
  
  table.new <- table(ifelse(actual == 0 & predicted == 0, "TN", 
                            ifelse(actual != 0 & predicted != 0, "TP", 
                                   ifelse(actual == 0 & predicted != 0, "FP", "FN"))))
  return(table.new)
}


# Create the confusion table (TP, TN, FP, FN) for each run
tabs3 <- lapply(1:runs, function(i) my.table(as.vector(Beta), gf.beta[[i]])[[1]])

# Initialize vectors for metrics
acur3 <- numeric(runs)
tpr3 <- numeric(runs)
tnr3 <- numeric(runs)

# Vectorized calculations for accuracy, true positive rate, and true negative rate
for (i in 1:runs) {
  # Extract confusion matrix counts
  TP <- ifelse(!is.na(tabs3[[i]]["TP"]), tabs3[[i]]["TP"], 0)
  TN <- ifelse(!is.na(tabs3[[i]]["TN"]), tabs3[[i]]["TN"], 0)
  FP <- ifelse(!is.na(tabs3[[i]]["FP"]), tabs3[[i]]["FP"], 0)
  FN <- ifelse(!is.na(tabs3[[i]]["FN"]), tabs3[[i]]["FN"], 0)
  
  # Calculate accuracy
  acur3[i] <- (TP + TN) / (TP + TN + FP + FN)
  
  # Calculate true negative rate (specificity)
  tnr3[i] <- TN / (TN + FP)
  
  # Calculate true positive rate (sensitivity)
  tpr3[i] <- TP / (TP + FN)
}

# Print the results
acur3
mean(acur3)

tpr3
mean(tpr3)

tnr3
mean(tnr3)


# Compute the bias of estimates

# Collect oracle beta values for each run
gf.oracle.beta <- sapply(1:runs, function(i) as.vector(round(gf[[i]]$opt.Beta[1:pa, ], 7)))

# Compute bias
gf.bias <- colMeans(gf.oracle.beta) - as.vector(Beta[1:pa,])

# Function to calculate the norm of a vector
norm.vec <- function(x) sqrt(sum(x^2))

# Compute bias squared
bias.gf <- norm.vec(as.vector(gf.bias))^2
bias.gf

# Compute standard deviation of estimates

# Collect oracle beta values again for consistency (can reuse `gf.oracle.beta` from above)
# gf.oracle.beta is already populated, so no need to repeat the loop

# Compute total variance (sum of eigenvalues of the covariance matrix)
gf.total.variance <- sum(eigen(cov(gf.oracle.beta))$values)

# Compute standard deviation
gf.std <- sqrt(gf.total.variance)
gf.std



#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#_ _ _ _ _ _ _ _ _ _ _ _ _ SPLS method _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# Fit the spls model for each run and store the results
spls.beta <- lapply(1:runs, function(i) spls(X[[i]], Y[[i]], q, eta = 0.8, kappa = 0.5, select = "simpls", 
                                             fit = "simpls", scale.x = FALSE, scale.y = FALSE, 
                                             eps = 1e-4, maxstep = 5000, trace = FALSE)$betamat)

# Round and store the beta values from the spls model
new.spls.beta1 <- lapply(spls.beta, function(beta) round(beta[[q]], 7))

# Vectorize beta values for each run
new.spls.beta <- lapply(spls.beta, function(beta) as.vector(round(beta[[q]], 7)))

# Define confusion table function
my.table <- function(actual, predicted) {
  actual <- as.vector(actual)
  predicted <- as.vector(predicted)
  
  table.new <- table(ifelse(actual == 0 & predicted == 0, "TN", 
                            ifelse(actual != 0 & predicted != 0, "TP", 
                                   ifelse(actual == 0 & predicted != 0, "FP", "FN"))))
  
  return(table.new)
}



# Calculate confusion tables for each run
tabs1 <- lapply(1:runs, function(i) my.table(as.vector(Beta), new.spls.beta[[i]])[[1]])

# Accuracy rate, TPR, TNR calculations
acur1 <- sapply(1:runs, function(i) {
  tp <- ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0)
  tn <- ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0)
  fp <- ifelse(!is.na(tabs1[[i]]["FP"]), tabs1[[i]]["FP"], 0)
  fn <- ifelse(!is.na(tabs1[[i]]["FN"]), tabs1[[i]]["FN"], 0)
  
  # Accuracy
  (tp + tn) / (tp + tn + fp + fn)
})

tpr1 <- sapply(1:runs, function(i) {
  tp <- ifelse(!is.na(tabs1[[i]]["TP"]), tabs1[[i]]["TP"], 0)
  fn <- ifelse(!is.na(tabs1[[i]]["FN"]), tabs1[[i]]["FN"], 0)
  
  # True Positive Rate
  tp / (tp + fn)
})

tnr1 <- sapply(1:runs, function(i) {
  tn <- ifelse(!is.na(tabs1[[i]]["TN"]), tabs1[[i]]["TN"], 0)
  fp <- ifelse(!is.na(tabs1[[i]]["FP"]), tabs1[[i]]["FP"], 0)
  
  # True Negative Rate
  tn / (tn + fp)
})

# Calculate means for accuracy, TPR, and TNR
mean_acur1 <- mean(acur1)
mean_tpr1 <- mean(tpr1)
mean_tnr1 <- mean(tnr1)

# Bias of estimates
spls.oracle.beta <- do.call(rbind, lapply(1:runs, function(i) {
  as.vector(round(spls.beta[[i]][[q]][1:pa,], 7))
}))

spls.bias <- colMeans(spls.oracle.beta) - as.vector(Beta[1:pa,])
bias.spls <- sqrt(sum(spls.bias^2))

# Standard deviation of estimates
spls.total.variance <- sum(eigen(cov(spls.oracle.beta))$values)
spls.std <- sqrt(spls.total.variance)

# Print results
list(
  acur1 = acur1,
  mean_acur1 = mean_acur1,
  tpr1 = tpr1,
  mean_tpr1 = mean_tpr1,
  tnr1 = tnr1,
  mean_tnr1 = mean_tnr1,
  bias_spls = bias.spls,
  spls_std = spls.std
)







