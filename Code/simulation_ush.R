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

source("Code/utils.R")
######################
## US sensor network##
######################
data_ustemp <- readMat("Data/stationary/US_Hourly_2010_August_1st.mat") # hourly data on 2010.08.01
UShourlytemp <- list()
UShourlytemp$xy <- cbind(data_ustemp$x, data_ustemp$y)
rownames(UShourlytemp$xy) <- 1:nrow(UShourlytemp$xy) 
UShourlytemp$f1 <- data_ustemp$signal[,15]
UShourlytemp$f2 <- data_ustemp$signal[,3]

distmat.htemp <- distm(UShourlytemp$xy, fun = distHaversine) / 1000
A.htemp <- c()
for(i in 1:(nrow(distmat.htemp)-1)){
  for(j in (i+1):ncol(distmat.htemp)){
    val <- distmat.htemp[i,j]
    A.htemp <- rbind(A.htemp, c(i,j,val))
  }
}

# G.knn <- nng(dx=distmat.htemp, k=5, mutual=TRUE)
G.knn <- as.undirected(nng(dx=distmat.htemp, k=7), mode="collapse")
edge.wt <- igraph::as_data_frame(G.knn, what="edges")
edge.wt <- sapply(edge.wt, as.numeric)
edge.wt <- cbind(edge.wt, 0)

for(i in 1:nrow(edge.wt)){
  edge.wt[i,3] <- distmat.htemp[edge.wt[i,1], edge.wt[i,2]]
}  

wmat <- matrix(0, nrow=length(UShourlytemp$f1), ncol=length(UShourlytemp$f1))

colnames(wmat) <- 1:length(UShourlytemp$f1)  
rownames(wmat) <- 1:length(UShourlytemp$f1)

for(i in 1:nrow(edge.wt)){
  wmat[edge.wt[i,1], edge.wt[i,2]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
  wmat[edge.wt[i,2], edge.wt[i,1]] <- exp(-edge.wt[i,3]^2/mean(edge.wt[,3])^2)
}  

sp.wmat <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat <- rbind(sp.wmat, c(edge.wt[i,1], edge.wt[i,2], 
                              wmat[edge.wt[i,1], edge.wt[i,2]]))
}

# weight matrix
UShourlytemp$A <- wmat

# sparse weight matrix
UShourlytemp$sA <- sp.wmat

UShourlytemp$dist <- distmat.htemp
UShourlytemp$sdist <- A.htemp

# visualize
plot_graph(UShourlytemp)

L.ush <- gasper::laplacian_mat(UShourlytemp$A)
N.ush <- nrow(UShourlytemp$xy)
val1 <- eigensort(L.ush)
evalues.ush <- val1$evalues
evectors.ush <- val1$evectors
# largest eigenvalue
lmax.ush <- max(evalues.ush)

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


tmp.heat.ush <- polyApprox(function(l) filter_kernel_heat(l, lmax=lmax.ush, tau=10), 0, lmax.ush, 5)
tmp.Mexican.ush <- polyApprox(function(l) filter_kernel_Mexican(l, lmax=lmax.ush,  tau1=5, tau2=25), 0, lmax.ush, 7)
tmp.sin.exp.ush <- polyApprox(function(l) filter_kernel(l, lmax=lmax.ush, tau1=5, tau2=5), 0, lmax.ush, 9)
tmp.highpass.ush <- polyApprox(function(l) filter_kernel_highpass(l, lmax=lmax.ush, tau=0.7), 0, lmax.ush, 4)
tmp.lowpass.ush <- polyApprox(function(l) filter_kernel_lowpass(l, lmax=lmax.ush, tau=3), 0, lmax.ush, 4)

lambdaseq <- evalues.ush
lambdaseq.normalized <- evalues.ush / lmax.ush

plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.heat.ush$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.Mexican.ush$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.sin.exp.ush$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.highpass.ush$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.lowpass.ush$p)), type="l")

H.heat.ush <- evectors.ush %*% diag(filter_poly(evalues.ush, h=rev(tmp.heat.ush$p))) %*% t(evectors.ush)
H.Mexican.ush <- evectors.ush %*% diag(filter_poly(evalues.ush, h=rev(tmp.Mexican.ush$p))) %*% t(evectors.ush)
H.sin.exp.ush <- evectors.ush %*% diag(filter_poly(evalues.ush, h=rev(tmp.sin.exp.ush$p))) %*% t(evectors.ush)
H.highpass.ush <- evectors.ush %*% diag(filter_poly(evalues.ush, h=rev(tmp.highpass.ush$p))) %*% t(evectors.ush)
H.lowpass.ush <- evectors.ush %*% diag(filter_poly(evalues.ush, h=rev(tmp.lowpass.ush$p))) %*% t(evectors.ush)


# R realizations case
R.list <- c(100,500,1000,2000)
p.list <- c(100,500,1000,2000)
# R.list <- c(2000)
set.seed(10)
tmp.var <- rbind(cbind(diag(1,N.ush), diag(1/2,N.ush)), 
                 cbind(diag(1/2,N.ush), diag(1,N.ush)))
tilde.white.inputs <- t(mvrnorm(n=max(R.list), mu=rep(0,2*N.ush), Sigma=tmp.var)) # (3*N.ush) x R 
white.inputs0 <- evectors.ush %*% tilde.white.inputs[1:N.ush,] # N.ush x R matrix
white.inputs1 <- evectors.ush %*% tilde.white.inputs[(N.ush+1):(2*N.ush),] # N.ush x R matrix


f1.ush <- H.sin.exp.ush %*% white.inputs0
f2.ush <- H.highpass.ush %*% white.inputs1


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

# B.ush <- bmat_func(p=6,q=1,lambdalist=evalues.ush) # p x q x n
pmax <- max(p.list) ; q <- 2 ; Rmax <- max(R.list)
set.seed(20)
# B.ush <- array(runif(p*q*N.ush, 2,4), c(p,q,N.ush)) # p x q x n
B.ush <- array(runif(pmax*q*N.ush, 10,20), c(pmax,q,N.ush)) # p x q x n
X.ush <- array(0, c(pmax,N.ush,Rmax)) # p x n x R
f.ush <- array(0,c(q,N.ush,Rmax)) # q x n x R
err.ush <- array(0,c(pmax,N.ush,Rmax)) # p x n x R

# 2 factors (q=2)
f.ush[1,,] <- f1.ush
f.ush[2,,] <- f2.ush

# p=50
set.seed(30)
for(i in 1:pmax){
  err.ush[i,,] <- t(mvrnorm(n=Rmax, mu=rep(0,N.ush), Sigma=diag(1,N.ush))) # N.ush x R matrix
}

for(i in 1:pmax){
  X.ush[i,,] <- evectors.ush%*%diag(B.ush[i,1,])%*%Conj(t(evectors.ush))%*%f.ush[1,,] +
    evectors.ush%*%diag(B.ush[i,2,])%*%Conj(t(evectors.ush))%*%f.ush[2,,] +
    err.ush[i,,]
}

# Fourier coefficients
X.tilde.ush <- array(0, c(pmax,N.ush,Rmax)) # p x n x R
f.tilde.ush <- array(0,c(q,N.ush,Rmax)) # q x n x R
err.tilde.ush <- array(0,c(pmax,N.ush,Rmax)) # p x n x R

for(i in 1:pmax){
  X.tilde.ush[i,,] <- Conj(t(evectors.ush)) %*% X.ush[i,,]
  err.tilde.ush[i,,] <- Conj(t(evectors.ush)) %*% err.ush[i,,]
}

for(j in 1:q){
  f.tilde.ush[j,,] <- Conj(t(evectors.ush)) %*% f.ush[j,,]
}

# reordering
X.ut.ush <- array(0,c(N.ush,Rmax,pmax)) # n x Rmax x p
f.ut.ush <- array(0,c(N.ush,Rmax,q)) # n x Rmax x q
err.ut.ush <- array(0,c(N.ush,Rmax,pmax)) # n x Rmax x p
B.ut.ush <- array(0,c(N.ush,pmax,q)) # n x p x q

for(l in 1:N.ush){
  X.ut.ush[l,,] <- t(X.tilde.ush[,l,])
  f.ut.ush[l,,] <- t(f.tilde.ush[,l,])
  err.ut.ush[l,,] <- t(err.tilde.ush[,l,])
  B.ut.ush[l,,] <- B.ush[,,l]
}


# construction check
range(t(X.ut.ush[1,,]) - (B.ut.ush[1,,]%*%t(f.ut.ush[1,,]) + t(err.ut.ush[1,,])))


kmax <- 5
kmin.array.ush <- array(0, c(N.ush,length(p.list),length(R.list)))
dimnames(kmin.array.ush) <- list(1:N.ush, paste("p=",p.list,sep=""), paste("R=",R.list,sep=""))
for(ind1 in 1:length(p.list)){
  for(ind2 in 1:length(R.list)){
    print(paste("p =",p.list[ind1], "R =", R.list[ind2]))
    p <- p.list[ind1]
    R <- R.list[ind2]
    tmpres <- matrix(0, nrow=N.ush, ncol=kmax)
    for(k in 1:kmax){
      f.bar.ush <- array(0, c(N.ush,R,k)) # n x R x k
      B.bar.ush <- array(0, c(N.ush,p,k)) # n x p x k
      for(l in 1:N.ush){
        eigenres <- eigen(X.ut.ush[l,1:R,1:p] %*% Conj(t(X.ut.ush[l,1:R,1:p])))
        evectors <- eigenres$vectors
        f.bar.ush[l,,] <- sqrt(R)*evectors[,1:k]
        B.bar.ush[l,,] <- t(X.ut.ush[l,1:R,1:p]) %*% Conj(f.bar.ush[l,,]) / R
        
        tmpres[l,k] <- IC_criterion(X.ut.ush[l,1:R,1:p], as.matrix(B.bar.ush[l,,]), 
                                    as.matrix(f.bar.ush[l,,]), g_penality)
      }
    }
    kmin.array.ush[,ind1,ind2] <- apply(tmpres, 1, which.min)
  }
}

apply(kmin.array.ush, c(2,3), FUN=function(x) table(x)[2] / sum(table(x)))


k <- 2

rmse.array.ush <- array(0, c(N.ush,length(p.list),length(R.list)))
dimnames(rmse.array.ush) <- list(1:N.ush, paste("p=",p.list,sep=""), paste("R=",R.list,sep=""))
for(ind1 in 1:length(p.list)){
  for(ind2 in 1:length(R.list)){
    print(paste("p =",p.list[ind1], "R =", R.list[ind2]))
    p <- p.list[ind1]
    R <- R.list[ind2]
    
    f.bar.ush <- array(0, c(N.ush,R,k)) # n x R x k
    B.bar.ush <- array(0, c(N.ush,p,k)) # n x p x k
    f.check.ush <- array(0, c(N.ush,R,k)) # n x R x k
    Hk <- array(0, c(N.ush,q,k)) # n x q x k
    
    for(l in 1:N.ush){
      eigenres <- eigen(X.ut.ush[l,1:R,1:p] %*% Conj(t(X.ut.ush[l,1:R,1:p])))
      evectors <- eigenres$vectors
      f.bar.ush[l,,] <- sqrt(R)*evectors[,1:k]
      B.bar.ush[l,,] <- t(X.ut.ush[l,1:R,1:p]) %*% Conj(f.bar.ush[l,,]) / R
      
      f.check.ush[l,,] <- X.ut.ush[l,1:R,1:p] %*% Conj(B.bar.ush[l,,]) / p
      Hk[l,,] <- t( (t(f.bar.ush[l,,])%*%Conj(f.ut.ush[l,1:R,]) / R) %*% (Conj(t(B.ut.ush[l,,]))%*%B.ut.ush[l,,] / p) )
      rmse.array.ush[l,ind1,ind2] <- norm(t(f.check.ush[l,,]) - t(Hk[l,,])%*%t(f.ut.ush[l,1:R,]), type="F") / R
    }
  }
}

apply(rmse.array.ush, c(2,3), mean)


