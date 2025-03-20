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

gramschmidt <- function (x) {
  x <- as.matrix (x)
  # Get the number of rows and columns of the matrix
  n <- ncol (x)
  m <- nrow (x)
  # Initialize the Q and R matrices
  q <- matrix (0 , m , n )
  r <- matrix (0 , n , n )
  for (j in 1: n) {
    v = x [,j] # Step 1 of the Gram - Schmidt process v1 = a1
    # Skip the first column
    if (j > 1) {
      for (i in 1:( j -1)) {
        r[i ,j] <- t(q[,i ]) %* % x[, j] # Find the inner product ( noted to be q ˆ T a earlier )
        # Subtract the projection from v which causes v to become pe rp end ic ul ar
        # to all columns of Q
        v <- v - r [i ,j] * q[,i ]
      }
    }
    # Find the L2 norm of the jth diagonal of R
    r[j ,j] <- sqrt ( sum (vˆ2))
    # The o r t h o g o n a l i z e d result is found and stored in the ith column of Q .
    q[,j ] <- v / r[j ,j]
  }
  # Collect the Q and R matrices into a list and return
  qrcomp <- list ( ’Q ’=q, ’R ’= r)
  return ( qrcomp )
}
gamma . pa <- gramschmidt ( matrix ( rnorm ( pa*pa , 0, 1) , pa , pa ))$ Q

gamma . a1 <- rbind ( gamma . pa [ ,1:q], matrix (0 , pi , q))
gamma . a2 <- rbind ( matrix ( gamma . pa [ ,(q +1): pa ] , pa , pa .q), matrix (0 , pi , pa .q))
gamma .i <- rbind ( matrix (0 , pa , pi ), diag (1 , pi ))
gamma0 <- cbind ( gamma .a2 , gamma .i)
gamma . both <- cbind ( gamma .a1 , gamma .a2 , gamma .i)
omega1 <- diag (4 , q)
omega .01 <- diag (0.8 , ( pa .q))
omega .02 <- diag (1 , pi )
omega0 <- as. matrix ( bdiag ( omega .01 , omega .02))
sigx <- gamma . a1 %* % omega1 %* % t( gamma . a1 ) + gamma0 %* % omega0 %* % t( gamma0 )
sigy .x <- matrix (0 , r ,r)
for (i in 1: r ){
  for (j in 1: r ){
    sigy .x[i ,j] <- delt ˆ( abs (i - j ))
  }
}
Beta <- rbind ( matrix (c( -3 , 5 , 1,
                           3, 0 , 0,
                           2, 1 , -2,
                           4, 0 , 3) , pa , r , byrow = T) , matrix (0 , pi , r )); Beta
Beta1 <- rbind ( matrix (c( -3 , 5 , 1,
                            3, 5e -4 , 5e -4 ,
                            2, 1 , -2,
                            4, 5e -4 , 3) , pa , r , byrow = T), matrix (5e -4 , pi , r )); Beta1
X <- list (); E <- list (); Y <- list ()
sigxx <- list (); sigxxyy <- list (); sigyy <- list (); Sxy <- list (); Sigs <- list ()
for (k in 1: runs ){
  X[k ] <- list ( mvrnorm (n , rep (0 , p ), sigx ))
  E[k ] <- list ( mvrnorm (n , rep (0 , r ), sigy . x ))
  Y[k ] <- list ( scale (X [[ k ]] %* % Beta + E [[ k ]]))
  X[k ] <- list ( scale (X [[ k ]]))
  sigxx [k] <- list ( cov (X [[ k ]]) *(n - 1)/n)
  sigxxyy [k] <- list ( cov (X [[ k ]] , Y [[ k ]]) *(n - 1)/n)
  sigyy [k] <- list ( cov (Y [[ k ]]) *(n - 1)/n)
  Sxy [k] <- list ( sigxx [[ k ]] - sigxxyy [[ k ]] %* % solve ( sigyy [[ k ]] )%* % t( sigxxyy [[ k ]]))
  Sigs [k] <- list ( list ( Sxy = Sxy [[ k ]] , sigxx = sigxx [[ k ]] , sigxxyy = sigxxyy [[ k ]]))
}
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Initial value for 2S - SPLS
myenvCv <- function (X , Y ){
  Y <- as. matrix (Y );
  X <- as. matrix (X );
  a <- dim (Y );
  n <- a [1];
  r <- a [2];
  p <- ncol (X );
  M <- list ()
  efit <- matrix (M , p)
  for (j in 1: p ){
    efit [j] <- list ( xenv (X , Y , u = j , asy = F ))
  }
  return ( efit )
}


em <- lapply (1: runs , function (i) myenvCv (X [[ i ]] , Y [[ i ]]))
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 2S - SPLS
Gamma . fixed1 <- function (par , X , Y , G , w , low , high , intr ){
  n <- nrow (X)
  r <- ncol (Y)
  p <- ncol (X)
  obj <- function (par , X , Y , G , pen , w ){
    Y <- as. matrix (Y)
    X <- as. matrix (X)
    X= scale (X , center =T , scale =F)
    Y= scale (Y , center =T , scale =F)
    G <- as. matrix (G)
    n <- nrow (X)
    r <- ncol (Y)
    p <- ncol (X)
    q <- ncol (G)
    B <- matrix (par , p , r)
    pen <- pen
    Z <- X%* %G
    sigZY <- cov (Z , Y)*(n - 1)/n
    sigZ <- cov (Z)*(n - 1)/n
    sigY <- cov (Y)*(n - 1)/n
    sigY .Z <- sigY - t( sigZY )%* % solve ( sigZ )%* % sigZY
    g. norma <- function (x1 , x2 ){
      r. norm <- matrix (0 ,p ,r )
      for (i in 1: p ){
        for (j in 1: r ){
          r. norm [i ,j] <- abs ( x1 [i ,j ])*(1/ abs ( x2 [i ,j ])ˆ0.5)
        }
      }
      return ( sum (r. norm ))
    }
    g. norm1 <- function (x ){
      r. norm <- matrix (0 ,p ,r )
      for (i in 1: p ){
        for (j in 1: r ){
          r. norm [i ,j] <- abs (x[i ,j ])
        }
      }
      return ( sum (r. norm ))
    }
    g. norm2 <- function (x ){
      r. norm1 <- rep (0 , p)
      for (i in 1: p ){
        r. norm1 [ i] <- sqrt ( sum (x[i ,]ˆ2))
      }
      return ( sum (r. norm1 ))
    }
    g. normp <- function (x1 , x2 ){
      r. norm1 <- rep (0 , p)
      for (i in 1: p ){
        r. norm1 [ i] <- sqrt ( sum ( x1 [i ,]ˆ2)) *(1/ sqrt ( sum ( x2 [i ,]ˆ2)))ˆ2
      }
      return (r. norm1 )
    }
    # J <-tr ((1 /n )*t (Y - Z %* %t ( G ) %* %B ) %* %(Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z ))+0.5 *g . norma (B , w )+
    # pen *g . normp (B , w )
    
    
    J <-tr ((1 /n)* t(Y - Z%* % t(G) %* %B)%* %(Y - Z%* % t(G) %* %B)%* % solve ( sigY .Z )) + pen *g. norma (B , w) +
      0.0005 *g. norm2 ( B)
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norma (B , w )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norm1 ( B )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B ) %* % solve ( sigY . Z )) +
    # pen *g . norm2 ( B ) + pen *g . norm1 ( B )
    # J <- tr ((1 /n )*t ( Y - Z %* %t ( G ) %* %B ) %* %( Y - Z %* %t ( G ) %* %B )) + pen *g . norm2 ( B )
    #_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _return (J)
  }
  gridd <- seq ( low , high , intr )
  fit4 <- lapply ( gridd , function (k) nlm (f= obj , p=par , X=X , Y=Y , G =G , pen =k , w =w ,
                                             iterlim = 2000 , gradtol = 1e -16))
  mini . fit <- rep (0 , length ( gridd ))
  for (k in 1: length ( gridd )){
    mini . fit [k] <- fit4 [[ k ]]$ minimum
  }
  para1 <- NULL
  for (k in 1: length ( gridd )){
    para1 [k] <- list ( round ( matrix ( fit4 [[ k ]]$ estimate ,p ,r ), 5))
  }
  bic . fit <- rep (0 , length ( gridd ))
  for (k in 1: length ( gridd )){
    bic . fit [k] <- mini . fit [k] + length (as. vector ( para1 [[ k ]])[ which (! as. vector ( para1 [[ k ]]) == 0)]) log (n)
  }
  bic . fit . min <- which . min ( bic . fit )
  opt . Beta <- para1 [[ bic . fit . min ]]
  opt . lambda <- gridd [ bic . fit . min ]
  #_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _return ( list ( para1 = para1 , mini . fit = mini . fit , bic . fit = bic . fit ,
  bic . fit . min = bic . fit .min , opt . Beta = opt . Beta , opt . lambda = opt . lambda ))
}
gf <- lapply (1: runs , function (i) Gamma . fixed1 ( em [[ i ]][[ q]]$ beta , X [[ i ]] , Y [[ i ]] ,
                                                      em [[ i ]][[ q]]$ Gamma , w= em [[ i ]][[ q]]$ beta ,
                                                      low =0 , high =0.007 , intr =0.001))
gf . beta1 <- NULL
for (i in 1: runs ){
  gf . beta1 [[ i ]] <- round ( gf [[ i ]]$ opt . Beta , 5)
}
gf . beta1
gf . beta <- NULL
for (i in 1: runs ){
  gf . beta [[ i ]] <- as. vector ( round ( gf [[ i ]]$ opt . Beta , 5))
}
my . table <- function ( actual , predicted ){
  actual <- as. vector ( actual ); predicted <- as. vector ( predicted )
  table . new <- table ( ifelse ( actual ==0 & predicted ==0 , " TN " ,
                                  ifelse ( actual ! =0 & predicted ! =0, " TP " ,
                                           ifelse ( actual ==0 & predicted ! =0, " FP " , " FN " ))))
  return ( table . new )
}


#_ _ _ table _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
tabs3 <- NULL
for (i in 1: runs ){
  tabs3 [[ i ]] <- list ( my . table (as. vector ( Beta ), gf . beta [[ i ]]))[[1]]
}
#_ _ _ Accuracy rate _ _ _ _ _ _ _ _ _ _ _ _ _ _
acur3 <- rep (0 , runs )
tpr3 <- rep (0 , runs )
tnr3 <- rep (0 , runs )
for (i in 1: runs ){
  acur3 [i] <- ( ifelse (! is.na( tabs3 [[ i ]][ " TP " ]) , tabs3 [[ i ]][ " TP " ], 0) +
                   ifelse (! is.na( tabs3 [[ i ]][ " TN " ]) , tabs3 [[ i ]][ " TN " ], 0)) /
    ( ifelse (! is.na( tabs3 [[ i ]][ " TP " ]) , tabs3 [[ i ]][ " TP " ], 0) +
        ifelse (! is.na( tabs3 [[ i ]][ " TN " ]) , tabs3 [[ i ]][ " TN " ], 0) +
        ifelse (! is.na( tabs3 [[ i ]][ " FP " ]) , tabs3 [[ i ]][ " FP " ] ,0)+ ifelse (! is.na( tabs3 [[ i ]][ " FN " ]) ,
                                                                                         tabs3 [[ i ]][ " FN " ] ,0))
  tnr3 [i] <- ifelse (! is.na( tabs3 [[ i ]][ " TN " ]) , tabs3 [[ i ]][ " TN " ], 0) /
    ( ifelse (! is.na( tabs3 [[ i ]][ " TN " ]) , tabs3 [[ i ]][ " TN " ], 0) +
        ifelse (! is.na( tabs3 [[ i ]][ " FP " ]) , tabs3 [[ i ]][ " FP " ], 0))
  tpr3 [i] <- ifelse (! is.na( tabs3 [[ i ]][ " TP " ]) , tabs3 [[ i ]][ " TP " ], 0) /
    ( ifelse (! is.na( tabs3 [[ i ]][ " TP " ]) , tabs3 [[ i ]][ " TP " ], 0) +
        ifelse (! is.na( tabs3 [[ i ]][ " FN " ]) , tabs3 [[ i ]][ " FN " ], 0))
}
acur3
mean ( acur3 )
tpr3
mean ( tpr3 )
tnr3
mean ( tnr3 )
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ bias of estimates _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
gf . oracle . beta <- matrix (NA , runs , r*pa )
for (i in 1: runs ){
  gf . oracle . beta [i ,] <- as. vector ( round ( gf [[ i ]]$ opt . Beta [1: pa ,] , 7))
}
gf . bias <- colMeans ( gf . oracle . beta ) - as. vector ( Beta [1: pa ,])
norm _ vec <- function (x) sqrt ( sum (x ˆ2))
bias . gf <- norm _ vec (as. vector ( gf . bias ))ˆ2
bias . gf
#_ _ _ _ _ _ _ _ _ _ _ _ _ standard deviation of estimates _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
gf . oracle . beta <- matrix (NA , runs , r*pa )
for (i in 1: runs ){
  gf . oracle . beta [i ,] <- as. vector ( round ( gf [[ i ]]$ opt . Beta [1: pa ,] , 7))
}
# gf . oracle . beta
gf . total . variance <- sum ( eigen ( cov ( gf . oracle . beta ))$ values )
gf . total . variance
gf . std <- sqrt ( gf . total . variance )
gf . std
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#_ _ _ _ _ _ _ _ _ _ _ _ _ SPLS method _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
spls . beta <- NULL
for (i in 1: runs ){
  spls . beta [i] <- list ( spls (X [[ i ]] , Y [[ i ]] , q, eta = 0.8 , kappa =0.5 , select =" simpls " ,
                                  fit =" simpls " , scale .x=F , scale .y=F , eps =1e -4 , maxstep =5000 ,
                                  trace =F)$ betamat )
}


new . spls . beta1 <- NULL
for (i in 1: runs ){
  new . spls . beta1 [[ i ]] <- round ( spls . beta [[ i ]][[ q]] , 7)
}
# sink (" new . spls . beta1 _ 100. txt ")
# print ( new . spls . beta1 )
# sink ()
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ vectorize beta _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
new . spls . beta <- NULL
for (i in 1: runs ){
  new . spls . beta [[ i ]] <- as. vector ( round ( spls . beta [[ i ]][[ q]] , 7))
}
my . table <- function ( actual , predicted ){
  actual <- as. vector ( actual ); predicted <- as. vector ( predicted )
  table . new <- table ( ifelse ( actual ==0 & predicted ==0 , " TN " ,
                                  ifelse ( actual ! =0 & predicted ! =0, " TP " ,
                                           ifelse ( actual ==0 & predicted ! =0, " FP " , " FN " ))))
  return ( table . new )
}
#_ _ _ table _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
tabs1 <- NULL
for (i in 1: runs ){
  tabs1 [[ i ]] <- list ( my . table (as. vector ( Beta ), new . spls . beta [[ i ]]))[[1]]
}
#_ _ _ Accuracy rate _ _ _ _ _ _ _ _ _ _ _ _ _ _
acur1 <- rep (0 , runs )
tpr1 <- rep (0 , runs )
tnr1 <- rep (0 , runs )
for (i in 1: runs ){
  acur1 [i] <- ( ifelse (! is.na( tabs1 [[ i ]][ " TP " ]) , tabs1 [[ i ]][ " TP " ], 0) +
                   ifelse (! is.na( tabs1 [[ i ]][ " TN " ]) , tabs1 [[ i ]][ " TN " ], 0)) /
    ( ifelse (! is.na( tabs1 [[ i ]][ " TP " ]) , tabs1 [[ i ]][ " TP " ], 0) +
        ifelse (! is.na( tabs1 [[ i ]][ " TN " ]) , tabs1 [[ i ]][ " TN " ], 0) +
        ifelse (! is.na( tabs1 [[ i ]][ " FP " ]) , tabs1 [[ i ]][ " FP " ], 0) +
        ifelse (! is.na( tabs1 [[ i ]][ " FN " ]) , tabs1 [[ i ]][ " FN " ], 0))
  tnr1 [i] <- ifelse (! is.na( tabs1 [[ i ]][ " TN " ]) , tabs1 [[ i ]][ " TN " ], 0) /
    ( ifelse (! is.na( tabs1 [[ i ]][ " TN " ]) , tabs1 [[ i ]][ " TN " ], 0) +
        ifelse (! is.na( tabs1 [[ i ]][ " FP " ]) , tabs1 [[ i ]][ " FP " ], 0))
  tpr1 [i] <- ifelse (! is.na( tabs1 [[ i ]][ " TP " ]) , tabs1 [[ i ]][ " TP " ], 0) /
    ( ifelse (! is.na( tabs1 [[ i ]][ " TP " ]) , tabs1 [[ i ]][ " TP " ], 0) +
        ifelse (! is.na( tabs1 [[ i ]][ " FN " ]) , tabs1 [[ i ]][ " FN " ], 0))
}
acur1
mean ( acur1 )
tpr1
mean ( tpr1 )
tnr1
mean ( tnr1 )
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ bias of estimates _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
spls . oracle . beta <- matrix (NA , runs , r*pa )
for (i in 1: runs ){
  spls . oracle . beta [i ,] <- as. vector ( round ( spls . beta [[ i ]][[ q ]][1: p.a ,] , 7))
}
spls . bias <- colMeans ( spls . oracle . beta ) - as. vector ( Beta [1: pa ,])
norm _ vec <- function (x) sqrt ( sum (x ˆ2))
bias . spls <- norm _ vec (as. vector ( spls . bias ))ˆ2
bias . spls


#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ standard deviation of estimates _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
spls . oracle . beta <- matrix (NA , runs , r*pa )
for (i in 1: runs ){
  spls . oracle . beta [i ,] <- as. vector ( round ( spls . beta [[ i ]][[ q ]][1: p.a ,] , 7))
}
# spls . oracle . beta
spls . total . variance <- sum ( eigen ( cov ( spls . oracle . beta ))$ values )
spls . total . variance
spls . std <- sqrt ( spls . total . variance )
spls . std







