# packages may be redundant..
library(gasper)
library(GNAR)
library(ggplot2)
library(cccd)
library(igraph)
library(igraphdata)
library(MASS)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)
library(psych)


source("Code/utils.R")

###############################################
## Comparison with classical factor analysis ##
###############################################

#############################
### random sensor network ###
#############################
n <- 10
x <- rep(0:(n-1), n) # define coordinates of vertices
y <- rep(0:(n-1), rep(n,n))

set.seed(1)
x <- runif(n^2, min = 0, max = n)
y <- runif(n^2, min = 0, max = n)

irregular02 <- list()
irregular02$xy <- data.frame(x,y)

N.irregular02 <- nrow(irregular02$xy)
rownames(irregular02$xy) <- 1:N.irregular02

distmat.irregular02 <- as.matrix(dist(irregular02$xy))
A.irregular02 <- c()
for(i in 1:(nrow(distmat.irregular02)-1)){
  for(j in (i+1):ncol(distmat.irregular02)){
    val <- distmat.irregular02[i,j]
    A.irregular02 <- rbind(A.irregular02, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.irregular02, k=5), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.irregular02[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=N.irregular02, ncol=N.irregular02)

colnames(wmat) <- 1:N.irregular02
rownames(wmat) <- 1:N.irregular02

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

## make the same irrgular graphs ##

# sparse weight matrix
# weight matrix
irregular01$A <- wmat

# sparse weight matrix
irregular01$sA <- sp.wmat

irregular01$dist <- distmat.irregular01
irregular01$sdist <- A.irregular01

L.irregular01 <- gasper::laplacian_mat(irregular01$A)

plot_graph(irregular01)
plot_graph_custom3(irregular01, e.size=1.3, v.size=6, vertex_color = "black", 
                   value="value", ratio=0.6, signal=FALSE)


val <- eigensort(L.irregular01)
lmax.irregular01 <- max(val$evalues)
evalues.irregular01 <- val$evalues
evectors.irregular01 <- val$evectors


# sparse weight matrix
# weight matrix
irregular02$A <- wmat

# sparse weight matrix
irregular02$sA <- sp.wmat

irregular02$dist <- distmat.irregular02
irregular02$sdist <- A.irregular02

L.irregular02 <- gasper::laplacian_mat(irregular02$A)

plot_graph(irregular02)

val <- eigensort(L.irregular02)
lmax.irregular02 <- max(val$evalues)
evalues.irregular02 <- val$evalues
evectors.irregular02 <- val$evectors




# filter kernel들 polynomial이어야하지만, chebyshev polynomial로 근사 잘 되므로 okay.
filter_kernel <- function(lambda, lmax, tau1, tau2){
  return(sin(3*tau1*lambda/lmax)*exp(-tau2*lambda/lmax))
}

filter_kernel_highpass <- function(lambda, lmax, tau){
  return(lambda/lmax*exp(-tau*lambda/lmax))
}

filter_kernel_lowpass <- function(lambda, lmax, tau){
  return(exp(-tau*(lambda/lmax)^4))
}

filter_kernel_heat <- function(lambda, lmax, tau){
  return(exp(-tau*lambda/lmax))
}

filter_kernel_Mexican <- function(lambda, lmax, tau1=5, tau2=25){
  return(tau1*lambda/lmax * exp(-tau2*lambda^2/lmax^2))
}
filter_poly <- function(lambda, h){
  res <- c()
  for(i in lambda){
    res <- c(res, sum(h*(i^c(0:(length(h)-1)))))
  }
  return(res)
}


## 3 factors, irregular02 ##
tmp.heat.irregular02 <- polyApprox(function(l) filter_kernel_heat(l, lmax=lmax.irregular02, tau=10), 0, lmax.irregular02, 5)
tmp.Mexican.irregular02 <- polyApprox(function(l) filter_kernel_Mexican(l, lmax=lmax.irregular02,  tau1=5, tau2=25), 0, lmax.irregular02, 7)
tmp.sin.exp.irregular02 <- polyApprox(function(l) filter_kernel(l, lmax=lmax.irregular02, tau1=5, tau2=5), 0, lmax.irregular02, 9)
tmp.highpass.irregular02 <- polyApprox(function(l) filter_kernel_highpass(l, lmax=lmax.irregular02, tau=0.7), 0, lmax.irregular02, 4)
tmp.lowpass.irregular02 <- polyApprox(function(l) filter_kernel_lowpass(l, lmax=lmax.irregular02, tau=3), 0, lmax.irregular02, 4)

lambdaseq <- evalues.irregular02
lambdaseq.normalized <- evalues.irregular02 / lmax.irregular02

# plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.heat.irregular02$p)), type="l")
# plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.Mexican.irregular02$p)), type="l")
# plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.sin.exp.irregular02$p)), type="l")
# plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.highpass.irregular02$p)), type="l")
# plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.lowpass.irregular02$p)), type="l")

H.heat.irregular02 <- evectors.irregular02 %*% diag(filter_poly(evalues.irregular02, h=rev(tmp.heat.irregular02$p))) %*% t(evectors.irregular02)
H.Mexican.irregular02 <- evectors.irregular02 %*% diag(filter_poly(evalues.irregular02, h=rev(tmp.Mexican.irregular02$p))) %*% t(evectors.irregular02)
H.sin.exp.irregular02 <- evectors.irregular02 %*% diag(filter_poly(evalues.irregular02, h=rev(tmp.sin.exp.irregular02$p))) %*% t(evectors.irregular02)
H.highpass.irregular02 <- evectors.irregular02 %*% diag(filter_poly(evalues.irregular02, h=rev(tmp.highpass.irregular02$p))) %*% t(evectors.irregular02)
H.lowpass.irregular02 <- evectors.irregular02 %*% diag(filter_poly(evalues.irregular02, h=rev(tmp.lowpass.irregular02$p))) %*% t(evectors.irregular02)


b_func <- function(i,j,l,lmax,p,q){
  # return(exp(-l/lmax)*sin(2*pi*(i-1/2)/p)*j)
  return(exp(-i*l/lmax)*j^2)
}

bmat_func <- function(p,q,lambdalist){
  lmax <- max(lambdalist)
  blist <- array(0, c(p,q,length(lambdalist)))
  for(l in 1:length(lambdalist)){
    tmp <- matrix(0, nrow=p, ncol=q)
    for(i in 1:p){
      for(j in 1:q){
        tmp[i,j] <- b_func(i,j,lambdalist[l],lmax,p,q)
      }
    }
    blist[,,l] <- tmp
  }
  return(blist)
}

# R realizations case
R.list <- c(50,100,250,500)
p.list <- c(10,20,30,40)


#########################################
## 2 factor simul on irregular01 graph ##
#########################################

rmse.array.irregular012.list <- list()
for(ind in 1:100){
  print(paste("iteration:", ind))
  
  set.seed(ind)
  # tmp.var <- rbind(cbind(diag(1,N.irregular01), diag(1/2,N.irregular01), diag(1/2,N.irregular01)), 
  #                  cbind(diag(1/2,N.irregular01), diag(1,N.irregular01), diag(1/2,N.irregular01)),
  #                  cbind(diag(1/2,N.irregular01), diag(1/2,N.irregular01), diag(1,N.irregular01)))
  tmp.var <- rbind(cbind(diag(1,N.irregular01), diag(1/2,N.irregular01)), 
                   cbind(diag(1/2,N.irregular01), diag(1,N.irregular01)))
  # tilde.white.inputs <- t(mvrnorm(n=max(R.list), mu=rep(0,2*N.irregular01), Sigma=tmp.var)) # (3*N.irregular01) x R 
  tilde.white.inputs <- mvrnorm(n=1, mu=rep(0,2*N.irregular01), Sigma=tmp.var) # (3*N.irregular01) x R 
  white.inputs0 <- evectors.irregular01 %*% tilde.white.inputs[1:N.irregular01] # N.irregular01 x R matrix
  white.inputs1 <- evectors.irregular01 %*% tilde.white.inputs[(N.irregular01+1):(2*N.irregular01)] # N.irregular01 x R matrix
  
  
  f1.irregular01 <- H.Mexican.irregular01 %*% white.inputs0
  f2.irregular01 <- H.highpass.irregular01 %*% white.inputs1
  
  # f1.irregular01 <- H.heat.irregular01 %*% white.inputs0
  # f2.irregular01 <- H.Mexican.irregular01 %*% white.inputs1
  # f3.irregular01 <- H.lowpass.irregular01 %*% white.inputs2
  # f3.irregular01 <- H.highpass.irregular01 %*% white.inputs2
  
  # set.seed(10)
  # white.inputs3 <- t(mvrnorm(n=max(R.list), mu=rep(0,N.irregular01), Sigma=diag(1,N.irregular01))) # N.irregular01 x R matrix
  # # white.inputs3 <- white.inputs2
  # err1.irregular01 <- H2.heat.irregular01 %*% white.inputs3
  # err2.irregular01 <- H2.Mexican.irregular01 %*% white.inputs3
  # err3.irregular01 <- H2.sin.exp.irregular01 %*% white.inputs3
  # err4.irregular01 <- H2.highpass.irregular01 %*% white.inputs3
  # err5.irregular01 <- H2.lowpass.irregular01 %*% white.inputs3
  # err6.irregular01 <- white.inputs3
  
  
  
  # B.irregular01 <- bmat_func(p=6,q=1,lambdalist=evalues.irregular01) # p x q x n
  pmax <- max(p.list) ; q <- 2 ; Rmax <- max(R.list)
  set.seed(ind+1)
  # B.irregular01 <- array(runif(p*q*N.irregular01, 5,10), c(p,q,N.irregular01)) # p x q x n
  B.irregular01 <- array(runif(pmax*q*N.irregular01, 10,20), c(pmax,q,N.irregular01)) # p x q x n
  X.irregular01 <- array(0, c(pmax,N.irregular01,Rmax)) # p x n x R
  f.irregular01 <- array(0,c(q,N.irregular01,Rmax)) # q x n x R
  err.irregular01 <- array(0,c(pmax,N.irregular01,Rmax)) # p x n x R
  
  # 3 factors (q=3)
  f.irregular01[1,,] <- replicate(Rmax,f1.irregular01)
  f.irregular01[2,,] <- replicate(Rmax,f2.irregular01)
  # f.irregular01[3,,] <- replicate(Rmax,f3.irregular01)
  
  # p=50
  set.seed(ind+2)
  for(i in 1:pmax){
    err.irregular01[i,,] <- replicate(Rmax,mvrnorm(n=1, mu=rep(0,N.irregular01), Sigma=diag(1,N.irregular01))) # N.irregular01 x R matrix
  }
  
  
  WB <- windowbank.random(N=length(evalues.irregular01), M=max(R.list), V=evectors.irregular01, sigma=0.1, seed=ind+4)
  for(i in 1:pmax){
    X.irregular01[i,,1] <- evectors.irregular01%*%diag(B.irregular01[i,1,])%*%Conj(t(evectors.irregular01))%*%f.irregular01[1,,1] +
      evectors.irregular01%*%diag(B.irregular01[i,2,])%*%Conj(t(evectors.irregular01))%*%f.irregular01[2,,1] +
      err.irregular01[i,,1]
    
    X.irregular01[i,,2:Rmax] <- (t(WB)*as.vector(X.irregular01[i,,1]))[,2:Rmax]
  }
  
  
  # Fourier coefficients
  X.tilde.irregular01 <- array(0, c(pmax,N.irregular01,Rmax)) # p x n x R
  f.tilde.irregular01 <- array(0,c(q,N.irregular01,Rmax)) # q x n x R
  err.tilde.irregular01 <- array(0,c(pmax,N.irregular01,Rmax)) # p x n x R
  
  for(i in 1:pmax){
    X.tilde.irregular01[i,,] <- Conj(t(evectors.irregular01)) %*% X.irregular01[i,,]
    err.tilde.irregular01[i,,] <- Conj(t(evectors.irregular01)) %*% err.irregular01[i,,]
  }
  
  for(j in 1:q){
    f.tilde.irregular01[j,,] <- Conj(t(evectors.irregular01)) %*% f.irregular01[j,,]
  }
  
  # reordering
  X.ut.irregular01 <- array(0,c(N.irregular01,Rmax,pmax)) # n x Rmax x p
  f.ut.irregular01 <- array(0,c(N.irregular01,Rmax,q)) # n x Rmax x q
  err.ut.irregular01 <- array(0,c(N.irregular01,Rmax,pmax)) # n x Rmax x p
  B.ut.irregular01 <- array(0,c(N.irregular01,pmax,q)) # n x p x q
  
  for(l in 1:N.irregular01){
    X.ut.irregular01[l,,] <- t(X.tilde.irregular01[,l,])
    f.ut.irregular01[l,,] <- t(f.tilde.irregular01[,l,])
    err.ut.irregular01[l,,] <- t(err.tilde.irregular01[,l,])
    B.ut.irregular01[l,,] <- B.irregular01[,,l]
  }
  
  rmse.array.irregular012 <- array(0, c(length(p.list),length(R.list)))
  dimnames(rmse.array.irregular012) <- list(paste("p=",p.list,sep=""), paste("R=",R.list,sep=""))
  k <- 2
  for(ind1 in 1:length(p.list)){
    for(ind2 in 1:length(R.list)){
      print(paste("p =",p.list[ind1], "R =", R.list[ind2]))
      p <- p.list[ind1]
      R <- R.list[ind2]
      
      f.bar.irregular01 <- array(0, c(N.irregular01,R,k)) # n x R x k
      B.bar.irregular01 <- array(0, c(N.irregular01,p,k)) # n x p x k
      # f.check.irregular01 <- array(0, c(N.irregular01,R,k)) # n x R x k
      # Hk <- array(0, c(N.irregular01,q,k)) # n x q x k
      f.tilde.irregular01.re <- array(0, c(k,N.irregular01,R)) # k x n x R
      for(l in 1:N.irregular01){
        eigenres <- eigen(X.ut.irregular01[l,1:R,1:p] %*% Conj(t(X.ut.irregular01[l,1:R,1:p])))
        evectors <- eigenres$vectors
        f.bar.irregular01[l,,] <- sqrt(R)*evectors[,1:k]
        B.bar.irregular01[l,,] <- t(X.ut.irregular01[l,1:R,1:p]) %*% Conj(f.bar.irregular01[l,,]) / R
        
        f.tilde.irregular01.re[,l,] <- t(f.bar.irregular01[l,,]) 
        # f.check.irregular01[l,,] <- X.ut.irregular01[l,1:R,1:p] %*% Conj(B.bar.irregular01[l,,]) / p
        # Hk[l,,] <- t( (t(f.bar.irregular01[l,,])%*%Conj(f.ut.irregular01[l,1:R,]) / R) %*% (Conj(t(B.ut.irregular01[l,,]))%*%B.ut.irregular01[l,,] / p) )
        # rmse.array.irregular01[l,ind1,ind2] <- norm((t(f.check.irregular01[l,,]) - t(Hk[l,,])%*%t(f.ut.irregular01[l,1:R,]))[,1], type="2") / 1
      }
      
      f.irregular01.re <- array(0, c(k,N.irregular01,R)) # k x n x R
      for(j in 1:k){
        f.irregular01.re[j,,] <- evectors.irregular01 %*% f.tilde.irregular01.re[j,,]
      }
      
      
      tmp5 <- c()
      for(pp in 1:p){
        tmp5 <- c(tmp5, X.irregular01[pp,,1] - evectors.irregular01 %*% diag(B.bar.irregular01[,pp,1]) %*% t(Conj(evectors.irregular01)) %*% f.irregular01.re[1,,1] -
                    evectors.irregular01 %*% diag(B.bar.irregular01[,pp,2]) %*% t(Conj(evectors.irregular01)) %*% f.irregular01.re[2,,1])
      }
      rmse.array.irregular012[ind1, ind2] <- norm(tmp5, type="2")
    }
  }
  
  
  fa_result1 <- fa(t(X.irregular01[1:p.list[1],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  # fa.parallel(t(X.irregular01[1:10,,1]), fa = "fa", n.iter = 100, show.legend = TRUE)
  # Print results
  
  fa_result2 <- fa(t(X.irregular01[1:p.list[2],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  # fa.parallel(t(X.irregular01[1:10,,1]), fa = "fa", n.iter = 100, show.legend = TRUE)
  # Print results
  
  fa_result3 <- fa(t(X.irregular01[1:p.list[3],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  
  fa_result4 <- fa(t(X.irregular01[1:p.list[4],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  
  rmse.array.irregular012 <- cbind(rmse.array.irregular012,
                                   c(norm(fa_result1$residual, type="F"),norm(fa_result2$residual, type="F"),norm(fa_result3$residual, type="F"),norm(fa_result4$residual, type="F")))
  rmse.array.irregular012.list[[ind]] <- rmse.array.irregular012
}


res1 <- round(apply(simplify2array(rmse.array.irregular012.list), c(1,2), mean),3)
round(apply(simplify2array(rmse.array.irregular012.list), c(1,2), sd),3)

par(mar=c(5,5,4,2)+0.1, mfrow=c(1,2))
plot(res1[1,1:4], type="l", ylim=c(0, max(res1)+0.1), xaxt="n", xlab="", ylab="reconstruction error", cex.lab=1.6, lwd=1.3, main="Two-factor model", cex.main=1.6)
axis(1, at=c(1,2,3,4), labels=paste("M =", c(50,100,250,500)), cex.axis=1.4)
points(res1[1,1:4], pch=2, cex=2)
lines(rep(res1[1,5],4), lty=2, lwd=1.3)

lines(res1[2,1:4], col="red", lwd=1.3)
points(res1[2,1:4], pch=5, cex=2, col="red")
lines(rep(res1[2,5],4), lty=2, col="red", lwd=1.3)

lines(res1[3,1:4], col="blue", lwd=1.3)
points(res1[3,1:4], pch=6, cex=2, col="blue")
lines(rep(res1[3,5],4), lty=2, col="blue", lwd=1.3)

lines(res1[4,1:4], col="chartreuse4", lwd=1.3)
points(res1[4,1:4], pch=21, cex=2.5, col="chartreuse4")
lines(rep(res1[4,5],4), lty=2, col="chartreuse4", lwd=1.3)

legend("topright", 
       legend=paste("p =", c(10,20,30,40)), 
       col=c("black", "red", "blue", "chartreuse4"), 
       lty=1, 
       lwd=2, 
       cex=1.2)




#########################################
## 3 factor simul on irregular02 graph ##
#########################################

rmse.array.irregular023.list <- list()
for(ind in 1:100){
  print(paste("iteration:", ind))
  
  set.seed(ind)
  # tmp.var <- rbind(cbind(diag(1,N.irregular02), diag(1/2,N.irregular02), diag(1/2,N.irregular02)), 
  #                  cbind(diag(1/2,N.irregular02), diag(1,N.irregular02), diag(1/2,N.irregular02)),
  #                  cbind(diag(1/2,N.irregular02), diag(1/2,N.irregular02), diag(1,N.irregular02)))
  tmp.var <- rbind(cbind(diag(1,N.irregular02), diag(1/2,N.irregular02), diag(1/2,N.irregular02)), 
                   cbind(diag(1/2,N.irregular02), diag(1,N.irregular02), diag(1/2,N.irregular02)),
                   cbind(diag(1/2,N.irregular02), diag(1/2,N.irregular02), diag(1,N.irregular02)))
  # tilde.white.inputs <- t(mvrnorm(n=max(R.list), mu=rep(0,2*N.irregular02), Sigma=tmp.var)) # (3*N.irregular02) x R 
  tilde.white.inputs <- mvrnorm(n=1, mu=rep(0,3*N.irregular02), Sigma=tmp.var) # (3*N.irregular02) x R 
  white.inputs0 <- evectors.irregular02 %*% tilde.white.inputs[1:N.irregular02] # N.irregular02 x R matrix
  white.inputs1 <- evectors.irregular02 %*% tilde.white.inputs[(N.irregular02+1):(2*N.irregular02)] # N.irregular02 x R matrix
  white.inputs2 <- evectors.irregular02 %*% tilde.white.inputs[(2*N.irregular02+1):(3*N.irregular02)] # N.karate x R matrix
  
  
  f1.irregular02 <- H.Mexican.irregular02 %*% white.inputs0
  f2.irregular02 <- H.highpass.irregular02 %*% white.inputs1
  f3.irregular02 <- H.lowpass.irregular02 %*% white.inputs2
  
  pmax <- max(p.list) ; q <- 3 ; Rmax <- max(R.list)
  set.seed(ind+1)
  # B.irregular02 <- array(runif(p*q*N.irregular02, 5,10), c(p,q,N.irregular02)) # p x q x n
  B.irregular02 <- array(runif(pmax*q*N.irregular02, 10,20), c(pmax,q,N.irregular02)) # p x q x n
  X.irregular02 <- array(0, c(pmax,N.irregular02,Rmax)) # p x n x R
  f.irregular02 <- array(0,c(q,N.irregular02,Rmax)) # q x n x R
  err.irregular02 <- array(0,c(pmax,N.irregular02,Rmax)) # p x n x R
  
  # 3 factors (q=3)
  f.irregular02[1,,] <- replicate(Rmax,f1.irregular02)
  f.irregular02[2,,] <- replicate(Rmax,f2.irregular02)
  f.irregular02[3,,] <- replicate(Rmax,f3.irregular02)
  
  # p=50
  set.seed(ind+2)
  for(i in 1:pmax){
    err.irregular02[i,,] <- replicate(Rmax,mvrnorm(n=1, mu=rep(0,N.irregular02), Sigma=diag(1,N.irregular02))) # N.irregular02 x R matrix
  }
  
  
  WB <- windowbank.random(N=length(evalues.irregular02), M=max(R.list), V=evectors.irregular02, sigma=0.1, seed=ind+4)
  for(i in 1:pmax){
    X.irregular02[i,,1] <- evectors.irregular02%*%diag(B.irregular02[i,1,])%*%Conj(t(evectors.irregular02))%*%f.irregular02[1,,1] +
      evectors.irregular02%*%diag(B.irregular02[i,2,])%*%Conj(t(evectors.irregular02))%*%f.irregular02[2,,1] +
      err.irregular02[i,,1]
    
    X.irregular02[i,,2:Rmax] <- (t(WB)*as.vector(X.irregular02[i,,1]))[,2:Rmax]
  }
  
  
  # Fourier coefficients
  X.tilde.irregular02 <- array(0, c(pmax,N.irregular02,Rmax)) # p x n x R
  f.tilde.irregular02 <- array(0,c(q,N.irregular02,Rmax)) # q x n x R
  err.tilde.irregular02 <- array(0,c(pmax,N.irregular02,Rmax)) # p x n x R
  
  for(i in 1:pmax){
    X.tilde.irregular02[i,,] <- Conj(t(evectors.irregular02)) %*% X.irregular02[i,,]
    err.tilde.irregular02[i,,] <- Conj(t(evectors.irregular02)) %*% err.irregular02[i,,]
  }
  
  for(j in 1:q){
    f.tilde.irregular02[j,,] <- Conj(t(evectors.irregular02)) %*% f.irregular02[j,,]
  }
  
  # reordering
  X.ut.irregular02 <- array(0,c(N.irregular02,Rmax,pmax)) # n x Rmax x p
  f.ut.irregular02 <- array(0,c(N.irregular02,Rmax,q)) # n x Rmax x q
  err.ut.irregular02 <- array(0,c(N.irregular02,Rmax,pmax)) # n x Rmax x p
  B.ut.irregular02 <- array(0,c(N.irregular02,pmax,q)) # n x p x q
  
  for(l in 1:N.irregular02){
    X.ut.irregular02[l,,] <- t(X.tilde.irregular02[,l,])
    f.ut.irregular02[l,,] <- t(f.tilde.irregular02[,l,])
    err.ut.irregular02[l,,] <- t(err.tilde.irregular02[,l,])
    B.ut.irregular02[l,,] <- B.irregular02[,,l]
  }
  
  
  rmse.array.irregular023 <- array(0, c(length(p.list),length(R.list)))
  dimnames(rmse.array.irregular023) <- list(paste("p=",p.list,sep=""), paste("R=",R.list,sep=""))
  k <- 3
  for(ind1 in 1:length(p.list)){
    for(ind2 in 1:length(R.list)){
      print(paste("p =",p.list[ind1], "R =", R.list[ind2]))
      p <- p.list[ind1]
      R <- R.list[ind2]
      
      f.bar.irregular02 <- array(0, c(N.irregular02,R,k)) # n x R x k
      B.bar.irregular02 <- array(0, c(N.irregular02,p,k)) # n x p x k
      # f.check.irregular02 <- array(0, c(N.irregular02,R,k)) # n x R x k
      # Hk <- array(0, c(N.irregular02,q,k)) # n x q x k
      f.tilde.irregular02.re <- array(0, c(k,N.irregular02,R)) # k x n x R
      for(l in 1:N.irregular02){
        eigenres <- eigen(X.ut.irregular02[l,1:R,1:p] %*% Conj(t(X.ut.irregular02[l,1:R,1:p])))
        evectors <- eigenres$vectors
        f.bar.irregular02[l,,] <- sqrt(R)*evectors[,1:k]
        B.bar.irregular02[l,,] <- t(X.ut.irregular02[l,1:R,1:p]) %*% Conj(f.bar.irregular02[l,,]) / R
        
        f.tilde.irregular02.re[,l,] <- t(f.bar.irregular02[l,,]) 
        # f.check.irregular02[l,,] <- X.ut.irregular02[l,1:R,1:p] %*% Conj(B.bar.irregular02[l,,]) / p
        # Hk[l,,] <- t( (t(f.bar.irregular02[l,,])%*%Conj(f.ut.irregular02[l,1:R,]) / R) %*% (Conj(t(B.ut.irregular02[l,,]))%*%B.ut.irregular02[l,,] / p) )
        # rmse.array.irregular02[l,ind1,ind2] <- norm((t(f.check.irregular02[l,,]) - t(Hk[l,,])%*%t(f.ut.irregular02[l,1:R,]))[,1], type="2") / 1
      }
      
      f.irregular02.re <- array(0, c(k,N.irregular02,R)) # k x n x R
      for(j in 1:k){
        f.irregular02.re[j,,] <- evectors.irregular02 %*% f.tilde.irregular02.re[j,,]
      }
      
      
      tmp5 <- c()
      for(pp in 1:p){
        tmp5 <- c(tmp5, X.irregular02[pp,,1] - evectors.irregular02 %*% diag(B.bar.irregular02[,pp,1]) %*% t(Conj(evectors.irregular02)) %*% f.irregular02.re[1,,1] -
                    evectors.irregular02 %*% diag(B.bar.irregular02[,pp,2]) %*% t(Conj(evectors.irregular02)) %*% f.irregular02.re[2,,1] - 
                    evectors.irregular02 %*% diag(B.bar.irregular02[,pp,3]) %*% t(Conj(evectors.irregular02)) %*% f.irregular02.re[3,,1])
      }
      rmse.array.irregular023[ind1, ind2] <- norm(tmp5, type="2")
    }
  }
  
  fa_result1 <- fa(t(X.irregular02[1:p.list[1],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  # fa.parallel(t(X.irregular01[1:10,,1]), fa = "fa", n.iter = 100, show.legend = TRUE)
  # Print results
  
  fa_result2 <- fa(t(X.irregular02[1:p.list[2],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  # fa.parallel(t(X.irregular01[1:10,,1]), fa = "fa", n.iter = 100, show.legend = TRUE)
  # Print results
  
  fa_result3 <- fa(t(X.irregular02[1:p.list[3],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  
  fa_result4 <- fa(t(X.irregular02[1:p.list[4],,1]), nfactors = k, rotate = "varimax", fm = "ml")
  
  rmse.array.irregular023 <- cbind(rmse.array.irregular023,
                                   c(norm(fa_result1$residual, type="F"),norm(fa_result2$residual, type="F"),norm(fa_result3$residual, type="F"),norm(fa_result4$residual, type="F")))
  rmse.array.irregular023.list[[ind]] <- rmse.array.irregular023
}


res2 <- round(apply(simplify2array(rmse.array.irregular023.list), c(1,2), mean),3)
round(apply(simplify2array(rmse.array.irregular023.list), c(1,2), sd),3)


plot(res2[1,1:4], type="l", ylim=c(0, max(res2)+0.1), xaxt="n", xlab="", ylab="reconstruction error", cex.lab=1.6, lwd=1.3, main="Three-factor model", cex.main=1.6)
axis(1, at=c(1,2,3,4), labels=paste("M =", c(50,100,250,500)), cex.axis=1.4)
points(res2[1,1:4], pch=2, cex=2)
lines(rep(res2[1,5],4), lty=2, lwd=1.3)

lines(res2[2,1:4], col="red", lwd=1.3)
points(res2[2,1:4], pch=5, cex=2, col="red")
lines(rep(res2[2,5],4), lty=2, col="red", lwd=1.3)

lines(res2[3,1:4], col="blue", lwd=1.3)
points(res2[3,1:4], pch=6, cex=2, col="blue")
lines(rep(res2[3,5],4), lty=2, col="blue", lwd=1.3)

lines(res2[4,1:4], col="chartreuse4", lwd=1.3)
points(res2[4,1:4], pch=21, cex=2.5, col="chartreuse4")
lines(rep(res2[4,5],4), lty=2, col="chartreuse4", lwd=1.3)

legend("topright", 
       legend=paste("p =", c(10,20,30,40)), 
       col=c("black", "red", "blue", "chartreuse4"), 
       lty=1, 
       lwd=2, 
       cex=1.2)
