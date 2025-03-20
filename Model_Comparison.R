# devtools :: install _ github ( ’ s elb ou ha dd an i / PPLS / Package /PPLS ’)
library ( mvtnorm ); library ( PPLS ); library ( pracma ); library ( Renvlp ); library ( pls );
library ( PPLS ); library ( matrixStats ); library ( fBasics ); library ( Hmisc ); library ( ggplot2 );
library ( reshape2 ); library ( simrel ); library ( future.apply ) library ( MASS );
# devtools :: install _ github (" s el bou ha dd an i / PO2PLS@RCpp ")
library ( OmicsPLS ); library ( PO2PLS )
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#_ _ _ _ _ _ _ _ _ _ _ data simulation _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
rm( list = ls ())
set.seed (100)
runs <- 100 # zhu used 200 replications
n <- 30 # zhu used 50 to 1000
p <- 5
r <- 3
q <- 2
delt <- 0.8
rho <- 0.01
# I R n 3 0 p 5 r 3 q 2 d e l t a 0 _8 rho0 _8
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
        r[i ,j] <- t(q[,i ]) %* % x[, j] # Find the inner product ( noted to be q ˆ T a
        # earlier ) Subtract the projection from v which causes v to become pe rp end ic ul ar
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
  qrcomp <- list ( 'Q'=q, 'R'= r)
  return ( qrcomp )
}
gamma <- gramschmidt ( matrix ( rnorm (p*p , 0, 1) , p , p ))$ Q
gamma1 <- gamma [ ,1:q]
gamma0 <- gamma [ ,(q +1): p]
sig <- matrix (0 , p ,p)
for (i in 1: p ){
  for (j in 1: p ){
    sig [i , j] <- rhoˆ( abs (i - j ))
  }
}
# omega <- eigen ( sig )$ values
# sigx <- gamma1 %* % diag ( omega [1: q ] , q ) %* %t ( gamma1 ) +
# gamma0 %* % diag ( omega [( q +1): p ] , (p - q )) %* %t ( gamma0 )
omega <- sort ( eigen ( sig )$ values , decreasing = F)
sigx <- gamma1 %* % diag ( omega [1: q], q) %* % t( gamma1 ) +
  gamma0 %* % diag ( omega [(q +1): p], (p -q)) %* % t( gamma0 )
sigy.x <- matrix (0 , r ,r)
for (i in 1: r ){
  for (j in 1: r ){
    sigy.x[i ,j] <- deltˆ( abs (i - j ))
  }
}
# sigy . x <- diag ( c (2 ,1.5 ,0.5) , r )
Beta <- gamma1 %* % matrix ( runif (q *r , 0, 2) , q, r)
X <- list (); E <- list (); Y <- list ()
for (k in 1: runs ){
  X[k ] <- list ( mvrnorm (n , rep (0 , p ), sigx ))
  E[k ] <- list (1* mvrnorm (n , rep (0 , r ), sigy.x ))
  Y[k ] <- list ( scale (X [[ k ]] %* % Beta + E [[ k ]]))
  X[k ] <- list ( scale (X [[ k ]]))
}
# cor ( X [[1]])







plan ( multiprocess )
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Using the Envelope
fold = k = 10
groups <- sample ( rep ( seq_len ( fold ), length.out = n ))
myenvCv <- function (X ,Y ,k ){
  # Y <- as . matrix ( Y [[1]]); X <- as . matrix ( X [[1]])
  
  Y <- as.matrix (Y ); X <- as.matrix (X)
  a <- dim (Y ); n <- a [1]; r <- a [2]; p <- ncol (X ); p <- p -1
  fitt <- plsr (Y ~ X , ncomp = p , method = " simpls ")$loadings
  M <- list ()
  efit <- matrix (M ,k ,p)
  for (i in 1: k ){
    for (j in 1: p ){
      efit [i , j] <- list ( xenv (X[ groups !=i ,] , Y[ groups !=i ,] , u = j , asy = F ,
                                   as.matrix ( fitt [ ,1: j ])))
    }
  }
  prederror <- matrix (M ,k ,p)
  for (i in 1: k ){
    for (j in 1: p ){
      prederror[i ,j] <- list (as.matrix (Y[ groups == i ,] - as.matrix (X[ groups == i ,])
                                            %* % efit [[ i ,j ]]$ beta ))
    }
  }
  ecv <- matrix (0 ,k ,p )
  for (i in 1: k ){
    for (j in 1: p ){
      ecv [i , j] <- norm ( prederror [i , j ][[1]] , type = "F")
    }
  }
  return ( envpls = colMeans ( ecv ))
}
envrepp1 <- matrix ( unlist ( future_lapply (1: runs , function (i) myenvCv (X=X [[ i ]] ,
                                                                               Y=Y [[ i ]] , k = fold ))) ,
                     ncol =p -1 , byrow =T)
ecv_mn1 <- colMeans ( envrepp1 ); ecv_sd1 <- colSds ( envrepp1 )
sink (" EIRn30p5r3q2delta0 _8 rho0 _ 01. txt ")
print ( envrepp1 )
print ( colMeans ( envrepp1 ))
print ( colSds ( envrepp1 ))
sink ()






#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Using the SIMPLS
mysimplsCv <- function (X ,Y ,k ){
  Y <- as.matrix (Y)
  X <- as.matrix (X)
  a <- dim (Y)
  n <- a [1]
  r <- a [2]
  p <- ncol (X)
  sfit <- NULL
  for (i in 1: k ){
    sfit [i] <- list ( plsr(Y[ groups !=i ,] ~ X[ groups !=i ,] , ncomp = p , method = " simpls " ))
  }
  cv <- matrix (0 ,k ,p )
  for (i in 1: k ){
    for (j in 1: p ){
      cv [i ,j ] <- norm( as.matrix ( Y[ groups ==i ,] ~ X[ groups ==i ,] %*% 
                                        as.matrix ( sfit [[ i ]]$coefficients [,, j ]) ))
    }
    }
  return ( splsCv <- colMeans ( cv ))
  
}
simplsrepp1 <- matrix ( unlist ( future_lapply (1: runs , function (i) mysimplsCv (X =X [[ i ]] , Y=Y [[ i ]] ,
                                                                                     k= fold ))) ,
                        ncol =p , byrow =T )
scv_mn1 <- colMeans ( simplsrepp1 ); scv_sd1 <- colSds ( simplsrepp1 )
sink (" SIRn30p5r3q2delta0 _8 rho0 _ 01. txt ")
print ( simplsrepp1 )
print ( colMeans ( simplsrepp1 ))
print ( colSds ( simplsrepp1 ))
sink ()




#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Using the EM - PLS
# myPPLSCv <- function (X , Y , k , emnum ){
# Y <- as . matrix ( Y ); X <- as . matrix ( X )
# a <- dim ( Y )
# n <- a [1]; r <- a [2]
# p <- ncol ( X )
# M = list ()
# pfit <- matrix (M ,k , r )
# for ( i in 1: k ){
# for ( j in 1: r ){
# pfit [i , j ] <- list ( PPLS _ simult ( X = as . matrix ( X [ groups ! =i ,]) ,
# Y = as . matrix ( Y [ groups ! =i ,]) , a = j , EMsteps = emnum , atol = 1e -09 ,
# type = " SVD "))
# }
# }
# pp . pred <- matrix (0 , k , r )
# for ( i in 1: k ){
# for ( j in 1: r ){
# pp . pred [i , j ] <- norm ( as . matrix ( Y [ groups == i ,] -
# X [ groups == i ,] %* % ginv ( pfit [i , j ][[1]] $ estimates $W %* %
# t ( pfit [i , j ][[1]] $ Expectations $mu _T ) %* %
# pfit [i , j ][[1]] $ Expectations $mu _T %* %t ( pfit [i , j ][[1]] $ estimates $W )) %* %
# pfit [i , j ][[1]] $ estimates $W %* %t ( pfit [i , j ][[1]] $ Expectations $mu _T ) %* %
# pfit [i , j ][[1]] $ Expectations $mu _U %* %
# t ( pfit [i , j ][[1]] $ estimates $C )) , type =" F ")
# }
# }
# return <- colMeans ( pp . pred )
# }
# PPLSrepp1 <- matrix ( unlist ( future _ lapply (1: runs , function ( i ) myPPLSCv ( X = X [[ i ]] ,
# Y = Y [[ i ]] , k = fold , emnum = 50))) ,
# ncol =r , byrow = T )
# pcv _ mn1 <- colMeans ( PPLSrepp1 ); pcv _ sd1 <- colSds ( PPLSrepp1 )



myPPLSCv <- function (X , Y , k , emnum ){
  Y <- as.matrix (Y ); X <- as.matrix (X)
  a <- dim (Y ); n <- a [1]; r <- a [2]; p <- ncol (X)
  # Y <- as . matrix ( Y ); X <- as . matrix ( X )
  M= list (); pfit <- matrix (M ,k ,r)
  for (i in 1: k ){
    for (j in 1: r ){
      pfit [i , j] <- list ( PO2PLS ( X[ groups !=i ,] , Y[ groups !=i ,] , r = j , 0 , 0, steps = emnum))
                                      
    }
  }
  pp.pred <- matrix (0 , k ,r)
  for (i in 1: k ){
    for (j in 1: r ){
      pp.pred [i ,j ] <- norm (as.matrix (Y[ groups == i ,] - X[ groups ==i ,] %* % ginv ( pfit [i , j ][[1]] $ params $W
                                                                                              %* % pfit [i , j ][[1]] $ params $ SigT %* % t( pfit [i , j ][[1]] $ params $W ))
                                             %* % pfit [i , j ][[1]] $ params $W %* % pfit [i , j ][[1]] $ params $ SigT %* %
                                               pfit [i , j ][[1]] $ params $B %* % t( pfit [i , j ][[1]] $ params $ C)) ,
                                 type ="F ")
    }
  }
  return <- colMeans ( pp.pred )
}
PPLSrepp1 <- matrix ( unlist ( future_lapply (1: runs , function (i) myPPLSCv (X=X [[ i ]] ,
                                                                                 Y=Y [[ i ]] , k = fold , emnum = 150))) , ncol =r , byrow =T )
pcv_mn1 <- colMeans ( PPLSrepp1 ); pcv_sd1 <- colSds ( PPLSrepp1 )
sink (" PIRn30p5r3q2delta0 _8 rho0 _ 01. txt ")
print ( PPLSrepp1 )
print ( colMeans ( PPLSrepp1 ))
print ( colSds ( PPLSrepp1 ))
sink ()




#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Using the OLS
mlmCV <- function (X , Y , k ){
  Y <- as.matrix (Y ); X <- as.matrix (X ); a <- dim (Y );
  n <- a [1]; r <- a [2]; p <- ncol (X)
  # k =10
  mlm1 <- NULL
  for (i in 1: k ){
    mlm1 [i] <- list (lm( cbind (Y[ groups ! =i ,]) ~ X[ groups ! =i ,]))
  }
  Yhat1 <- NULL
  for (i in 1: k ){
    Yhat1 [i] <- list ( X[ groups ==i ,] %* % coef ( mlm1 [[ i ]])[ -1 ,] )
  }
  mlmnorm <- rep (0 , k )
  for (i in 1: k ){
    mlmnorm [i ]= norm( as.matrix (Y[ groups == i ,] - Yhat1 [[ i ]]) , type = "F ")
  }
  return ( rep ( mean ( mlmnorm ),p ))
}
mlmrepp <- matrix ( unlist ( future_lapply (1: runs , function (i) mlmCV (X =X [[ i ]] , Y=Y [[ i ]] , k= fold ))) ,
                    ncol =p , byrow =T )
lcv_mn <- colMeans( mlmrepp ); lcv_sd <- colSds( mlmrepp )
sink (" LIRn30p5r3q2delta0 _8 rho0 _ 01. txt ")
print ( mlmrepp )
print ( colMeans ( mlmrepp ))
print ( colSds ( mlmrepp ))
sink ()




#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Plots _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
x <- c (1: p)

plot (x ,c( ecv_mn1 , lcv_mn[1]) , type ="o " , col =" red " , ylab =" Cross - validated  RMSE " ,
      xlab = expression (" Number ␣ of ␣ components "), lty =1 , pch =20 ,
      ylim = c(0 , 3) , main = expression ( list (n ==30 , p ==5 , r ==3 , q==2 , rho ==0.01)))
# errbar (x , ecv _mn1 , ecv _ mn1 + ecv _sd1 , ecv _ mn1 - ecv _sd1 , add = T , col = " red " ,
# cap =0.01 , pch = 20 , lwd = 0.5 , errbar . col = " red ")
lines (x , scv_mn1 , type = "o" , pch =20 , col =" blue " , lty =2)
# errbar ( x* 1.02 , scv _mn1 , scv _ mn1 + scv _sd1 , scv _ mn1 - scv _sd1 , add = T ,
# cap =0.01 , col = " blue " , pch = 20 , lwd = 0.5 , errbar . col = " blue ")
lines (c (1: r), pcv_mn1 , type ="o" , col =" black " , pch =20 , lty =3)
# errbar ( c (1: r ) , pcv _mn1 , pcv _ mn1 + pcv _sd1 , pcv _ mn1 - pcv _sd1 , add = T ,
# cap =0.01 , col = " black " , pch = 20 , lwd = 0.5 , errbar . col = " black ")
lines (x , lcv_mn , type = "o " , pch =20 , col = " cyan3 " , lty =4)
# errbar ( x* 1.04 , lcv _mn , lcv _mn + lcv _sd , lcv _mn - lcv _sd , add = T ,
# col = " cyan3 " , cap =0.01 , pch = 20 , lwd = 0.5 , errbar . col = " cyan3 ")
# lines (x , rep ( mean ( unlist ( evar1 )) , p ) , col = " purple " , lty = 5)
abline (v=q, col = " gray " , lty = 6)
legend (" bottomright " , legend =c(" EPLS " ," SIMPLS " ,"EM - PLS " , " OLS "),
        col =c(" red " , " blue " , " black " , " cyan3 " ), lty =1:4 , cex =1)