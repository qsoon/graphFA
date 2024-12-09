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

# window filterbank mother function for WGF-based estimation
g <- function(lambda, sigma.sq, m, tau){
  return(exp(-(lambda-m*tau)^2 / sigma.sq))
}

#############################
### 01. graph preparation ###
#############################

###########################
### karate club network ###
###########################
data(karate)
V(karate)$name
E(karate)
edge.attributes(karate)
degree(karate)


N.karate <- length(V(karate)$name)
edge.wt <- igraph::as_data_frame(karate, what="edges")
for(i in 1:nrow(edge.wt)){
  edge.wt[i,1] <- which(V(karate)$name == edge.wt[i,1])
  edge.wt[i,2] <- which(V(karate)$name == edge.wt[i,2])
}
edge.wt <- sapply(edge.wt, as.numeric)
wmat.karate <- matrix(0, nrow=gorder(karate), ncol=gorder(karate))
# colnames(wmat.karate) <- V(karate)$name
# rownames(wmat.karate) <- V(karate)$name

colnames(wmat.karate) <- 1:N.karate
rownames(wmat.karate) <- 1:N.karate

for(i in 1:nrow(edge.wt)){
  wmat.karate[edge.wt[i,1], edge.wt[i,2]] <- edge.wt[i,3]
  wmat.karate[edge.wt[i,2], edge.wt[i,1]] <- edge.wt[i,3]
}

sp.wmat.karate <- c()
for(i in 1:nrow(edge.wt)){
  sp.wmat.karate <- rbind(sp.wmat.karate, c(edge.wt[i,1], edge.wt[i,2], 
                                            wmat.karate[edge.wt[i,1], edge.wt[i,2]]))
}

zachary.karate <- list()
set.seed(1)
layout.karate <- layout.fruchterman.reingold(karate)
layout.karate <- layout.norm(layout.karate, -1, 1, -1, 1)
zachary.karate$xy <- data.frame(x=layout.karate[,1],
                                y=layout.karate[,2])

zachary.karate$A <- wmat.karate
zachary.karate$sA <- sp.wmat.karate

plot_graph(zachary.karate)

L.karate <- gasper::laplacian_mat(wmat.karate)
val <- eigensort(L.karate)
lmax.karate <- max(val$evalues)
evalues.karate <- val$evalues
evectors.karate <- val$evectors

eigenres.karate <- eigen(L.karate)


##############
## gFactor ##
##############

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


## 3 factors, karate ##
tmp.heat.karate <- polyApprox(function(l) filter_kernel_heat(l, lmax=lmax.karate, tau=10), 0, lmax.karate, 5)
tmp.Mexican.karate <- polyApprox(function(l) filter_kernel_Mexican(l, lmax=lmax.karate,  tau1=5, tau2=25), 0, lmax.karate, 7)
tmp.sin.exp.karate <- polyApprox(function(l) filter_kernel(l, lmax=lmax.karate, tau1=5, tau2=5), 0, lmax.karate, 9)
tmp.highpass.karate <- polyApprox(function(l) filter_kernel_highpass(l, lmax=lmax.karate, tau=0.7), 0, lmax.karate, 4)
tmp.lowpass.karate <- polyApprox(function(l) filter_kernel_lowpass(l, lmax=lmax.karate, tau=3), 0, lmax.karate, 4)

lambdaseq <- evalues.karate
lambdaseq.normalized <- evalues.karate / lmax.karate

plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.heat.karate$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.Mexican.karate$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.sin.exp.karate$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.highpass.karate$p)), type="l")
plot(lambdaseq.normalized, filter_poly(lambdaseq, h=rev(tmp.lowpass.karate$p)), type="l")

H.heat.karate <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.heat.karate$p))) %*% t(evectors.karate)
H.Mexican.karate <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.Mexican.karate$p))) %*% t(evectors.karate)
H.sin.exp.karate <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.sin.exp.karate$p))) %*% t(evectors.karate)
H.highpass.karate <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.highpass.karate$p))) %*% t(evectors.karate)
H.lowpass.karate <- evectors.karate %*% diag(filter_poly(evalues.karate, h=rev(tmp.lowpass.karate$p))) %*% t(evectors.karate)


# R realizations case
R.list <- c(10,100,1000)
set.seed(1)
# tmp.var <- rbind(cbind(diag(1,N.karate), diag(1/2,N.karate), diag(1/2,N.karate)), 
#                  cbind(diag(1/2,N.karate), diag(1,N.karate), diag(1/2,N.karate)),
#                  cbind(diag(1/2,N.karate), diag(1/2,N.karate), diag(1,N.karate)))
tmp.var <- rbind(cbind(diag(1,N.karate), diag(1/2,N.karate), diag(0,N.karate)), 
                 cbind(diag(1/2,N.karate), diag(1,N.karate), diag(0,N.karate)),
                 cbind(diag(0,N.karate), diag(0,N.karate), diag(1/4,N.karate)))
tilde.white.inputs <- t(mvrnorm(n=max(R.list), mu=rep(0,3*N.karate), Sigma=tmp.var)) # (3*N.karate) x R 
white.inputs0 <- evectors.karate %*% tilde.white.inputs[1:N.karate,] # N.karate x R matrix
white.inputs1 <- evectors.karate %*% tilde.white.inputs[(N.karate+1):(2*N.karate),] # N.karate x R matrix
white.inputs2 <- evectors.karate %*% tilde.white.inputs[(2*N.karate+1):(3*N.karate),] # N.karate x R matrix


f1.karate <- H.heat.karate %*% white.inputs0
f2.karate <- H.Mexican.karate %*% white.inputs1
f3.karate <- H.lowpass.karate %*% white.inputs2

set.seed(10)
# white.inputs3 <- t(mvrnorm(n=max(R.list), mu=rep(0,N.karate), Sigma=diag(1,N.karate))) # N.karate x R matrix
white.inputs3 <- white.inputs2
err1.karate <- H2.heat.karate %*% white.inputs3
err2.karate <- H2.Mexican.karate %*% white.inputs3
err3.karate <- H2.sin.exp.karate %*% white.inputs3
err4.karate <- H2.highpass.karate %*% white.inputs3
err5.karate <- H2.lowpass.karate %*% white.inputs3
err6.karate <- white.inputs3

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

B.karate <- bmat_func(p=6,q=1,lambdalist=evalues.karate) # p x q x n

p <- 6 ; q <- 3 ; Rmax <- max(R.list)
X.karate <- array(0, c(p,N.karate,Rmax)) # p x n x R
f.karate <- array(0,c(q,N.karate,Rmax)) # q x n x R
err.karate <- array(0,c(p,N.karate,Rmax)) # p x n x R

# 3 factors (q=3)
f.karate[1,,] <- f1.karate
f.karate[2,,] <- f2.karate
f.karate[3,,] <- f3.karate

# p=6
err.karate[1,,] <- err1.karate
err.karate[2,,] <- err2.karate
err.karate[3,,] <- err3.karate
err.karate[4,,] <- err4.karate
err.karate[5,,] <- err5.karate
err.karate[6,,] <- err6.karate

for(i in 1:p){
  X.karate[i,,] <- evectors.karate%*%diag(B.karate[i,1,])%*%Conj(t(evectors.karate))%*%f.karate[1,,] +
    evectors.karate%*%diag(B.karate[i,2,])%*%Conj(t(evectors.karate))%*%f.karate[2,,] +
    evectors.karate%*%diag(B.karate[i,3,])%*%Conj(t(evectors.karate))%*%f.karate[3,,] +
    err.karate[i,,]
  # X.karate[i,,] <- evectors.karate%*%diag(B.karate[i,1,])%*%Conj(t(evectors.karate))%*%f.karate[1,,] +
  #   err.karate[i,,]
}

# Fourier coefficients
X.tilde.karate <- array(0, c(p,N.karate,Rmax)) # p x n x R
f.tilde.karate <- array(0,c(q,N.karate,Rmax)) # q x n x R
err.tilde.karate <- array(0,c(p,N.karate,Rmax)) # p x n x R

for(i in 1:p){
  X.tilde.karate[i,,] <- Conj(t(evectors.karate)) %*% X.karate[i,,]
  err.tilde.karate[i,,] <- Conj(t(evectors.karate)) %*% err.karate[i,,]
}

for(j in 1:q){
  f.tilde.karate[j,,] <- Conj(t(evectors.karate)) %*% f.karate[j,,]
}

# reordering
X.ut.karate <- array(0,c(N.karate,Rmax,p)) # n x Rmax x p
f.ut.karate <- array(0,c(N.karate,Rmax,q)) # n x Rmax x q
err.ut.karate <- array(0,c(N.karate,Rmax,p)) # n x Rmax x p
B.ut.karate <- array(0,c(N.karate,p,q)) # n x p x q

for(l in 1:N.karate){
  X.ut.karate[l,,] <- t(X.tilde.karate[,l,])
  f.ut.karate[l,,] <- t(f.tilde.karate[,l,])
  err.ut.karate[l,,] <- t(err.tilde.karate[,l,])
  B.ut.karate[l,,] <- B.karate[,,l]
}


# construction check
range(t(X.ut.karate[1,,]) - (B.ut.karate[1,,]%*%t(f.ut.karate[1,,]) + t(err.ut.karate[1,,])))
range(X.ut.karate[1,,])
range(B.ut.karate[1,,]%*%t(f.ut.karate[1,,]))
range(err.ut.karate[1,,])


# theta_func(X.ut.karate[l,1:R,], B.bar.karate[l,,], f.bar.karate[l,,])

kmax=6
IC.mat.karate <- array(0, c(N.karate,kmax,length(R.list)))
for(k in 1:kmax){
  print(paste("k =",k))
  
  for(ind in 1:length(R.list)){
    R <- R.list[ind]
    f.bar.karate <- array(0, c(N.karate,R,k)) # n x R x k
    B.bar.karate <- array(0, c(N.karate,p,k)) # n x p x k
    
    for(l in 1:N.karate){
      eigenres <- eigen(X.ut.karate[l,1:R,] %*% Conj(t(X.ut.karate[l,1:R,])))
      evectors <- eigenres$vectors
      f.bar.karate[l,,] <- sqrt(R)*evectors[,1:k]
      B.bar.karate[l,,] <- t(X.ut.karate[l,1:R,]) %*% Conj(f.bar.karate[l,,]) / R
      
      IC.mat.karate[l,k,ind] <- IC_criterion(X.ut.karate[l,1:R,], as.matrix(B.bar.karate[l,,]), 
                                             as.matrix(f.bar.karate[l,,]), g_penality)
    }
  }
}

kmin.mat.karate <- matrix(0,nrow=N.karate, ncol=length(R.list))
colnames(kmin.mat.karate) <- R.list
rownames(kmin.mat.karate) <- 1:N.karate
for(l in 1:N.karate){
  kmin.mat.karate[l,] <- apply(IC.mat.karate[l,,], 2, which.min)
}
kmin.mat.karate
table(kmin.mat.karate[,1])
table(kmin.mat.karate[,2])
table(kmin.mat.karate[,3])

k <- 3

f.bar.karate.list <- list()
B.bar.karate.list <- list()
f.check.karate.list <- list()
Hk.karate.list <- list()
rmse.mat.karate <- matrix(0, nrow=N.karate, ncol=length(R.list))
for(ind in 1:length(R.list)){
  R <- R.list[ind]
  f.bar.karate <- array(0, c(N.karate,R,k)) # n x R x k
  B.bar.karate <- array(0, c(N.karate,p,k)) # n x p x k
  f.check.karate <- array(0, c(N.karate,R,k)) # n x R x k
  Hk <- array(0, c(N.karate,q,k)) # n x q x k
  
  for(l in 1:N.karate){
    # k <- kmin.mat[l,ind]
    # R <- R.list[ind]
    # f.bar.karate <- array(0, c(N.karate,R,k)) # n x R x k
    # B.bar.karate <- array(0, c(N.karate,p,k)) # n x p x k
    # f.check.karate <- array(0, c(N.karate,R,k)) # n x R x k
    # Hk <- array(0, c(N.karate,q,k)) # n x q x k
    
    eigenres <- eigen(X.ut.karate[l,1:R,] %*% Conj(t(X.ut.karate[l,1:R,])))
    evectors <- eigenres$vectors
    f.bar.karate[l,,] <- sqrt(R)*evectors[,1:k]
    B.bar.karate[l,,] <- t(X.ut.karate[l,1:R,]) %*% Conj(f.bar.karate[l,,]) / R
    f.check.karate[l,,] <- X.ut.karate[l,1:R,] %*% Conj(B.bar.karate[l,,]) / p
    Hk[l,,] <- t( (t(f.bar.karate[l,,])%*%Conj(f.ut.karate[l,1:R,]) / R) %*% (Conj(t(B.ut.karate[l,,]))%*%B.ut.karate[l,,] / p) )
    rmse.mat.karate[l,ind] <- norm(t(f.check.karate[l,,]) - t(Hk[l,,])%*%t(f.ut.karate[l,1:R,]), type="F") / R
  }
  
  f.bar.karate.list[[ind]] <- f.bar.karate
  B.bar.karate.list[[ind]] <- B.bar.karate
  f.check.karate.list[[ind]] <- f.check.karate
  Hk.karate.list[[ind]] <- Hk
}

rmse.mat.karate
plot(rmse.mat.karate[20,], type="l")


